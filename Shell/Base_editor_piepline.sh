#!/bin/bash
while getopts "o:1:2:R:" opt
do
  case $opt in
    o)
      out_dir=$OPTARG;;
    1)
      fq1=$OPTARG;;
    2)
      fq2=$OPTARG;;
	R)
      ref=$OPTARG;;	
    ?)
      echo "-o outdir -1 read1 -2 read2"
      exit 1;;
  esac
done

if [ $ref = "NET" ]
then
	bowtie_index="/home/qzh/Analysis/Shengwang/ref/NET/bowtie_index/NET"
	fa_ref="/home/qzh/Analysis/Shengwang/ref/NET/NET.fa"
elif [ $ref = "Lag6" ]
then
	bowtie_index="/home/qzh/Analysis/Shengwang/ref/Lag6/bowtie_index/Lag6"
	fa_ref="/home/qzh/Analysis/Shengwang/ref/Lag6/Lag6.fa"
else
	echo "species not recognized"
	exit
fi

function trim_seq {
	local fq1=$1
	local fq2=$2	
    local out_dir=$3
    local fname=$(basename ${fq1%_1*})

	local clean2_fq1="${out_dir}/${fname}_1.fq.gz"
	local clean2_fq2="${out_dir}/${fname}_2.fq.gz"
	local clean1_fq1="${out_dir}/${fname}_A_1.fq.gz"
	local clean1_fq2="${out_dir}/${fname}_A_2.fq.gz"


    [[ -f ${clean2_fq1} ]] || \
        trim_galore -q 20 \
			--clip_R1 2 --clip_R2 2 \
            --length 30 -j 8 \
			-a A{10} -a T{10} \
			-a2 A{10} -a2 T{10} \
            --paired \
            ${fq1} ${fq2} \
            -o ${out_dir}

	[[ -f ${clean2_fq1} ]] || \
		mv ${out_dir}/${fname}_1_val_1.fq.gz ${clean2_fq1}

	[[ -f ${clean2_fq2} ]] || \
		mv ${out_dir}/${fname}_2_val_2.fq.gz ${clean2_fq2}

	[[ -f ${clean1_fq1} ]] && rm ${clean1_fq1} ${clean1_fq2}
}
function align_bowtie_fq1 {
	local fq1=$1
    local out_dir=$2
    local fname=$(basename ${fq1%_1*})

	local map_Q35_file="${out_dir}/${fname}_Q35_1.map"
	local map_file="${out_dir}/${fname}_1.map"
	[[ -f ${map_file} ]] || \
		bowtie \
			-p 10 -n 3 -e 450 -l 16 \
			-x ${bowtie_index} \
			$fq1 > ${map_file} \
			2> ${out_dir}/${fname}_1_bowtie.log
	
}
function align_bowtie_fq2 {
	local fq2=$1
    local out_dir=$2
    local fname=$(basename ${fq2%_2*})

	local map_Q35_file="${out_dir}/${fname}_Q35_2.map"
	local map_file="${out_dir}/${fname}_2.map"
	[[ -f ${map_file} ]] || \
		bowtie \
			-p 10 -n 3 -e 450 -l 16 \
			-x ${bowtie_index} \
			$fq2 > ${map_file} \
			2> ${out_dir}/${fname}_2_bowtie.log

}

function align_bowtie {
	local fq1=$1	
	local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1.fq*})

	local Merge_bam="${out_dir}/${fname}_merge.bam"
	local bam_file1="${out_dir}/${fname}_1.bam"
	local sam_file1="${out_dir}/${fname}_1.sam"

	[[ -f ${bam_file1} ]] || \
		bowtie \
			-p 10 -n 3 -S \
			-e 450 -l 16 \
			-x ${bowtie_index} \
			$fq1 \
			1> ${sam_file1} \
			2> ${out_dir}/${fname}_bowtie1.log

	[[ -f ${bam_file1} ]] || \
		samtools view -Sub <(samtools view -H ${sam_file1}; \
			samtools view -f 0 ${sam_file1}) | \
			samtools sort -@ 15 -o ${bam_file1} -

	local bam_file2="${out_dir}/${fname}_2.bam"
	local sam_file2="${out_dir}/${fname}_2.sam"
	[[ -f ${bam_file2} ]] || \
		bowtie \
			-p 10 -n 3 -S \
			-e 450 -l 16 \
			-x ${bowtie_index} \
			$fq2 \
			1> ${sam_file2} \
			2> ${out_dir}/${fname}_bowtie2.log

	[[ -f ${bam_file2} ]] || \
		samtools view -Sub <(samtools view -H ${sam_file2}; \
			samtools view -f 16 ${sam_file2}) | \
			samtools sort -@ 15 -o ${bam_file2} -

	[[ -f ${Merge_bam} ]] || \
		samtools merge -@ 12 -f ${Merge_bam} \
			${bam_file1} ${bam_file2} 

	local vcf_file1="${out_dir}/${fname}_r1.vcf"
	local vcf_file="${out_dir}/${fname}.vcf"	
	[[ -f ${vcf_file} ]] || \
		samtools mpileup \
			-A --no-output-ends --no-output-ins \
			--no-output-ins-mods --no-output-del -d 1000000 \
			-f ${fa_ref}  \
			${bam_file1} ${bam_file2} > ${vcf_file1}

	[[ -f ${vcf_file} ]] || \
		awk -v OFS='\t' '($4 > $7) {print $1,$2,$3,$4,$5,$6} ($4 <= $7) {print $1,$2,$3,$7,$8,$9}' ${vcf_file1} \
			> ${vcf_file}

	# ~/Analysis/Shengwang/ref/Lag6/Lag6.fa
	[[ -f ${sam_file1} ]] && rm ${sam_file1}
	[[ -f ${sam_file2} ]] && rm ${sam_file2}	
	[[ -f ${vcf_file1} ]] && rm ${vcf_file1}
}
function align_bwa {
	local fq1=$1	
	local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_rmdup*})

	local Merge_bam="${out_dir}/${fname}_merge.bam"
	local bam_file1="${out_dir}/${fname}_1.bam"
	local sam_file1="${out_dir}/${fname}_1.sam"
	[[ -f ${Merge_bam} ]] || \
		bwa mem -t 10 -T 50 -k 15 -Y \
			/home/qzh/Analysis/Shengwang/ref/Lag6/bwa_index/Lag6 \
			${fq1} \
			1> ${sam_file1} \
			2> ${out_dir}/${fname}_bwa1.log
			
	[[ -f ${Merge_bam} ]] || \
		samtools view -Sub <(samtools view -H ${sam_file1}; \
			samtools view ${sam_file1}) | \
			samtools sort -@ 15 -o ${bam_file1} -

	local bam_file2="${out_dir}/${fname}_2.bam"
	local sam_file2="${out_dir}/${fname}_2.sam"
	[[ -f ${Merge_bam} ]] || \
		bwa mem -t 10 -T 50 -k 15 -Y \
			/home/qzh/Analysis/Shengwang/ref/Lag6/bwa_index/Lag6 \
			${fq2} \
			1> ${sam_file2} \
			2> ${out_dir}/${fname}_bwa2.log

	[[ -f ${Merge_bam} ]] || \
		samtools view -Sub <(samtools view -H ${sam_file2}; \
			samtools view ${sam_file2}) | \
			samtools sort -@ 15 -o ${bam_file2} -

	[[ -f ${Merge_bam} ]] || \
		samtools merge -@ 12 -f ${Merge_bam} \
			${bam_file1} ${bam_file2} 

	local vcf_file="${out_dir}/${fname}.vcf"
	[[ -f ${vcf_file} ]] || \
		samtools mpileup \
			--no-output-ends --no-output-ins \
			--no-output-ins-mods --no-output-del -d 1000000 \
			-f ~/Analysis/Shengwang/ref/Lag6/Lag6.fa \
			${Merge_bam} > ${vcf_file}

	# [[ -f ${sam_file1} ]] && rm ${sam_file1}
	# [[ -f ${sam_file2} ]] && rm ${sam_file2}
	[[ -f ${bam_file1} ]] && rm ${bam_file1}
	[[ -f ${bam_file2} ]] && rm ${bam_file2}		
}
function edit_rate_cal {
	local vcf=$1	
    local out_dir=$2
    local fname=$(basename ${vcf%.vcf*})

	local outfile_file1="${out_dir}/${fname}_plot.txt"
	local outfile_file2="${out_dir}/${fname}.txt"
	python ~/Analysis/Shengwang/src/vcf_base_content.py \
		--infile ${vcf} \
		--outfile ${outfile_file1} \
		--outfile_2 ${outfile_file2} \
		--scale 1 \
		--fname ${fname}

	local mutate_file="${out_dir}/${fname}_position_mutate.txt"	
	cat ${outfile_file1} | \
		grep -v 'mutate\|T->T\|G->G\|A->A\|C->C' \
		> ${mutate_file}

	local mutate_per_file="${out_dir}/${fname}_position_mutate.results"	
	cat ${outfile_file1} | \
		grep -v 'mutate\|T->T\|G->G\|A->A\|C->C' | \
		awk '{x[$1"\t"$2"\t"$3]+=$4;} END{for(i in x) print(i"\t"x[i])}' | \
		sort -k2 -n | \
		awk 'BEGIN{print "Chr\tPos\tref\trate\tSample"}{print $0"\t""'"${fname}"'"}' \
		> ${mutate_per_file}

	local all_base=`cat ${outfile_file1} | \
		grep -v 'mutate' | \
		awk '{s[$6]+=$5;}END{for(i in s)print(i"\t"s[i])}' | \
		awk '{sum+=$2}END{print sum}'`

	local mutate_file="${out_dir}/${fname}_mutate.results"
	cat ${outfile_file1} | \
		grep -v 'mutate\|T->T\|G->G\|A->A\|C->C' | \
		awk '{s[$6]+=$4;}END{for(i in s)print(i"\t"s[i])}' | \
		awk 'BEGIN{print "sample\ttype\trate"}{print "'"${fname}"'""\t"$1"\t"$2}' \
		> ${mutate_file}
}
function SlamDunk_align {
	local fq1=$1	
	local fq2=$2	
    local out_dir=$3
    local fname=$(basename ${fq1%_1.fq*})

	local SlamDunk_bam_file1="${out_dir}/$(basename ${fq1%.gz*})_slamdunk_mapped.bam"
	local SlamDunk_bam_file2="${out_dir}/$(basename ${fq2%.gz*})_slamdunk_mapped.bam"
	local Merge_bam="${out_dir}/${fname}_merge.bam"	
	[[ -f ${Merge_bam} ]] || \
		slamdunk map \
			-r ${fa_ref} \
			-t 10 -5 0 \
			-o ${out_dir} \
			${fq1} ${fq2}

	local bam_file1="${out_dir}/${fname}_1_sort.bam"
	[[ -f ${Merge_bam} ]] || \
		samtools sort -@ 10 -o ${bam_file1} \
			${SlamDunk_bam_file1}

	local bam_file2="${out_dir}/${fname}_2_sort.bam"
	[[ -f ${Merge_bam} ]] || \
		samtools sort -@ 10 -o ${bam_file2} \
			${SlamDunk_bam_file2}

	[[ -f ${Merge_bam} ]] || \
		samtools merge -@ 12 -f ${Merge_bam} \
			${bam_file1} ${bam_file2} 

	local vcf_file="${out_dir}/${fname}.vcf"
	[[ -f ${vcf_file} ]] || \
		samtools mpileup \
			--no-output-ends --no-output-ins \
			--no-output-ins-mods --no-output-del -d 1000000 \
			-f ${fa_ref} \
			${Merge_bam} > ${vcf_file}

	local outfile_file1="${out_dir}/${fname}_plot.txt"
	local outfile_file2="${out_dir}/${fname}.txt"
	[[ -f ${outfile_file2} ]] || \	
		python ~/Analysis/Shengwang/src/vcf_base_content.py \
			--infile ${vcf_file} \
			--outfile ${outfile_file1} \
			--outfile_2 ${outfile_file2} \
			--scale 1 \
			--fname ${fname}

	local mutate_file="${out_dir}/${fname}_position_mutate.txt"	
	cat ${mutate_file} | \
		grep -v 'mutate\|T->T\|G->G\|A->A\|C->C' \
		> ${mutate_file}

	local mutate_per_file="${out_dir}/${fname}_position_mutate.results"	
	cat ${mutate_per_file} | \
		grep -v 'mutate\|T->T\|G->G\|A->A\|C->C' | \
		awk '{x[$1"\t"$2"\t"$3]+=$4;} END{for(i in x) print(i"\t"x[i])}' | \
		sort -k2 -n | \
		awk 'BEGIN{print "Chr\tPos\tref\trate\tSample"}{print $0"\t""'"${fname}"'"}' \
		> ${mutate_per_file}

	local all_base=`cat ${outfile_file1} | \
		grep -v 'mutate' | \
		awk '{s[$6]+=$5;}END{for(i in s)print(i"\t"s[i])}' | \
		awk '{sum+=$2}END{print sum}'`

	local mutate_file="${out_dir}/${fname}_mutate.results"
	cat ${mutate_file} | \
		grep -v 'mutate\|T->T\|G->G\|A->A\|C->C' | \
		awk '{s[$6]+=$4;}END{for(i in s)print(i"\t"s[i])}' | \
		awk 'BEGIN{print "sample\ttype\trate"}{print "'"${fname}"'""\t"$1"\t"$2}' \
		> ${mutate_file}

	[[ -f ${bam_file1} ]] && rm ${bam_file1} ${bam_file2}
}

fname=$(basename $fq1 "_1.fq.gz")

mkdir -p $out_dir/${fname}
# trim
mkdir -p $out_dir/${fname}/01-clean
trim_seq $fq1 $fq2 \
	$out_dir/${fname}/01-clean

# align 
mkdir -p $out_dir/${fname}/02-align
align_bowtie_fq1 $out_dir/${fname}/01-clean/${fname}_1.fq.gz \
	$out_dir/${fname}/02-align

align_bowtie_fq2 $out_dir/${fname}/01-clean/${fname}_2.fq.gz \
	$out_dir/${fname}/02-align

[[ -f $out_dir/${fname}/01-clean/${fname}_1.fq ]] || \
	zcat $out_dir/${fname}/01-clean/${fname}_1.fq.gz \
		> $out_dir/${fname}/01-clean/${fname}_1.fq
# align 
mkdir -p $out_dir/${fname}/03-align_paired
align_bowtie $out_dir/${fname}/01-clean/${fname}_1.fq.gz \
	$out_dir/${fname}/01-clean/${fname}_2.fq.gz \
	$out_dir/${fname}/03-align_paired
 
