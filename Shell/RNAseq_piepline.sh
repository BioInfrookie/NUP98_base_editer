getopt_cmd=$(getopt -o g:2:1:o:h --long help -n $(basename $0) -- "$@")
[ $? -ne 0 ] && exit 1
eval set -- "$getopt_cmd"
# 解析选项
while [ -n "$1" ]
do
    case "$1" in
        -1)
            fq1="$2"
            shift ;;
        -2)
            fq2="$2"
            shift ;;
        -g)
            species="$2"
            shift ;;
        -o)
            out_dir="$2"
            shift ;;
		-h|--help)
            echo -e "$help_str"
            exit ;;
        --) shift
            break ;;
         *) echo "$1 is not an option"
            exit 1 ;;
    esac
    shift
done

if [ "$species" ==  Sac ]
then
    genome="/home/qzh/ref/Saccharomyces/STAR_index"
    rRNA_index="/home/qzh/ref/Saccharomyces/STAR_index/rRNA"
    gtf_file="/home/qzh/ref/Saccharomyces/Saccharomyces_cerevisiae.R64-1-1.52.gtf"
    reporter_index="/home/qzh/Analysis/Shengwang/ref/Nup98/STAR_index"
	# Exon_length="/home/qzh/ref/ss11/ss11_exon_length.txt"
else
    exit 1;
fi

function trim_reads() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1.fq*})

    local fq1_clean2="${out_dir}/${fname}_1_val_1.fq.gz"
    [[ -f ${fq1_clean2} ]] || \
        trim_galore -q 20 \
            --length 30 -j 8 \
			-a A{30} -a T{30} \
			-a2 A{30} -a2 T{30} \
            --paired \
            ${fq1} ${fq2} \
            -o ${out_dir} 
}
function align_genome() {
    local fq1=$1
    local fq2=$2
    local out_dir=$3
    local fname=$(basename ${fq1%_1_val*})

    local genome_align_file="${out_dir}/${fname}_genome"    
    local rRNA_align_file="${out_dir}/${fname}_rRNA"
    [[ -f ${genome_align_file}_Aligned.sortedByCoord.out.bam ]] || \
        STAR --runThreadN 15 --genomeDir ${rRNA_index} \
            --runMode alignReads \
            --genomeLoad NoSharedMemory \
            --readFilesIn ${fq1} ${fq2} \
            --outFileNamePrefix ${rRNA_align_file}_ \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outTmpDir ${rRNA_align_file}_rRNA.tmp \
            --outFilterMultimapNmax 20 \
            --outBAMsortingThreadN 10 \
            --limitBAMsortRAM 2019378293 \
            --outReadsUnmapped Fastx

    [[ -f ${genome_align_file}_Aligned.sortedByCoord.out.bam ]] || \
        STAR --runThreadN 15 --genomeDir ${genome} \
            --runMode alignReads \
            --genomeLoad NoSharedMemory \
            --readFilesIn ${rRNA_align_file}_Unmapped.out.mate1 ${rRNA_align_file}_Unmapped.out.mate2 \
            --outFileNamePrefix ${genome_align_file}_ \
            --readFilesCommand cat \
            --outSAMtype BAM SortedByCoordinate \
            --outTmpDir ${genome_align_file}_genome.tmp \
			--outFilterMultimapNmax 1 \
			--outSAMmultNmax 1 \
            --outBAMsortingThreadN 10 \
            --limitBAMsortRAM 2019378293 \
            --outSAMstrandField intronMotif \
            --outReadsUnmapped Fastx   

    [[ -f ${genome_align_file}_Aligned.sortedByCoord.out.bam.bai ]] || \
        samtools index ${genome_align_file}_Aligned.sortedByCoord.out.bam

    local reporter_file="${out_dir}/${fname}_reporter"
    [[ -f ${reporter_file}_Aligned.sortedByCoord.out.bam ]] || \
        STAR --runThreadN 15 --genomeDir ${reporter_index} \
            --runMode alignReads \
            --genomeLoad NoSharedMemory \
            --readFilesIn ${genome_align_file}_Unmapped.out.mate1 ${genome_align_file}_Unmapped.out.mate2 \
            --outFileNamePrefix ${reporter_file}_ \
            --readFilesCommand cat \
            --outSAMtype BAM SortedByCoordinate \
            --outTmpDir ${reporter_file}_reporter.tmp \
			--outFilterMultimapNmax 1 \
			--outSAMmultNmax 1 \
            --outBAMsortingThreadN 10 \
            --limitBAMsortRAM 2019378293 

    [[ -f ${reporter_file}_Aligned.sortedByCoord.out.bam.bai ]] || \
        samtools index ${reporter_file}_Aligned.sortedByCoord.out.bam

    # Remove
    [[ -f ${rRNA_align_file}_Aligned.sortedByCoord.out.bam ]] && rm ${rRNA_align_file}_Aligned.sortedByCoord.out.bam
    [[ -f ${rRNA_align_file}_Unmapped.out.mate1 ]] && rm ${rRNA_align_file}_Unmapped.out.mate1 ${rRNA_align_file}_Unmapped.out.mate2
}
function quantify_counts() {
    local bam=$1
    local reporter_bam=$2
    local out_dir=$3
    local fname=$(basename ${bam%_Aligned*})

    local counts_file="${out_dir}/${fname}_counts.txt"
    [[ -f ${counts_file} ]] || \
        featureCounts -s 0 \
            -o ${counts_file} \
            -T 6 -M -O --fraction \
            -p -C -B \
            -a ${gtf_file} \
            -F GTF -t exon -g gene_id \
            ${bam}

    cat ${counts_file} | tail -n +3 | \
        awk '{print $1"\t"$NF}' | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "id\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_gene.txt

	local Map_reads=`cat ${counts_file} | tail -n +3 | \
		awk '{sum+=$NF}END{print sum}' | \
		awk '{printf "%.f\n",$0}'`
	local RPM_factor=`echo "scale=8;1000000/${Map_reads}" | bc`
	echo ${RPM_factor}

    cat ${counts_file}  | tail -n +3 | \
        awk '{print $NF}' | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_only_count.txt

    cat ${counts_file}  | tail -n +3 | \
        awk '{print $1"\t"$NF*"'"${RPM_factor}"'"}' | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "Gene_id\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_RPM.txt

    cat ${counts_file}  | tail -n +3 | \
        awk '{print $NF*"'"${RPM_factor}"'"}' | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print fname1}{print $0}' \
        > ${out_dir}/${fname}_gene_only_RPM.txt

    samtools idxstats ${reporter_bam} | sed '$d' | cut -f1,3 | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "id\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_reporter_count.txt
    samtools idxstats ${reporter_bam} | sed '$d' | cut -f3 | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print fname1}{print $0}' \
        > ${out_dir}/${fname}_reporter_only_count.txt
    samtools idxstats ${reporter_bam} | sed '$d' | \
        awk '{print $3*"'"${RPM_factor}"'"}' | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print fname1}{print $0}' \
        > ${out_dir}/${fname}_reporter_only_RPM.txt
    
    samtools idxstats ${reporter_bam} | sed '$d' | cut -f1,3 \
        >  ${out_dir}/${fname}_reporter_tmp.txt
    cat ${counts_file} | tail -n +3 | \
        awk '{print $1"\t"$NF}' \
        >  ${out_dir}/${fname}_genome_tmp.txt
    
    cat ${out_dir}/${fname}_reporter_tmp.txt ${out_dir}/${fname}_genome_tmp.txt | \
        awk -v fname1="'"${fname}"'" 'BEGIN{print "id\t"fname1}{print $0}' \
        > ${out_dir}/${fname}_merge_count.txt

    [[ -f ${out_dir}/${fname}_reporter_tmp.txt ]] && rm ${out_dir}/${fname}_reporter_tmp.txt
    [[ -f ${out_dir}/${fname}_genome_tmp.txt ]] && rm ${out_dir}/${fname}_genome_tmp.txt
}
function kraken2_run {
	local fq1=$1
	local fq2=$2
	local out_dir=$3
	local fname=$(basename ${fq1%_genome*})

	local report_file="${out_dir}/${fname}.report"
	[[ -f ${report_file} ]] || \
		kraken2 --db /data/biodata/kraken2/standard \
			--threads 6  \
			--report ${report_file} \
			--paired $fq1 $fq2

	local unmap_pct="${out_dir}/${fname}_pct.out"
	[[ -f ${unmap_pct} ]] || \
		cat ${report_file} | \
			grep $'\t'D$'\t' | \
			sort -k1 -r | \
			awk 'BEGIN{print "Pct\tspecies\ttype"}{print $1"\t"$6"\t""'"${fname}"'"}' \
			> ${unmap_pct}

	cp /home/qzh/src/pollution_monitoring.Rmd $out_dir/
	cp /home/qzh/src/pollution_monitoring_report.R $out_dir/

	work_dir=$(readlink -fs $out_dir)
	local html_file="${out_dir}/${fname}-report.html"	
	[[ -f ${html_file} ]] || \
		Rscript $out_dir/pollution_monitoring_report.R $work_dir/ ${fname}
}
# piepline
fname=$(basename ${fq1%_1.fq*})
# Clean reads
mkdir -p $out_dir/${fname}/01-clean-data
trim_reads $fq1 $fq2 \
    $out_dir/${fname}/01-clean-data

# Star mapping
mkdir -p $out_dir/${fname}/02-align
align_genome $out_dir/${fname}/01-clean-data/${fname}_1_val_1.fq.gz \
    $out_dir/${fname}/01-clean-data/${fname}_2_val_2.fq.gz \
    $out_dir/${fname}/02-align

clean_reads=$(cat ${out_dir}/${fname}/01-clean-data/${fname}_C.log | \
    grep "Pairs written" | awk '{print $(NF-1)}' | sed 's/,//g')
Unmap_reads=$(seqkit stat ${out_dir}/${fname}/02-align/${fname}_genome_Unmapped.out.mate1 | \
    awk '/DNA/{print $4}' | sed 's/,//g')
mapping_reads=$(echo ${clean_reads} ${Unmap_reads} | \
    awk '{print $1-$2}')

# Feacture counts
mkdir -p $out_dir/${fname}/03-counts
quantify_counts ${out_dir}/${fname}/02-align/${fname}_genome_Aligned.sortedByCoord.out.bam \
    ${out_dir}/${fname}/02-align/${fname}_reporter_Aligned.sortedByCoord.out.bam \
    ${out_dir}/${fname}/03-counts 

# echo "Step4 Unmap reads content"
# mkdir -p ${out_dir}/${fname}/04-Unmap	
# kraken2_run \
#     ${out_dir}/${fname}/02-align/${fname}_genome_Unmapped.out.mate1 \
#     ${out_dir}/${fname}/02-align/${fname}_genome_Unmapped.out.mate2 \
#     ${out_dir}/${fname}/04-Unmap

# Step Summary
reporte_summary="$out_dir/${fname}/${fname}.results"
if [[ ! -f ${reporte_summary} ]]
then
    total_reads=$(seqkit stat $fq1 | awk '{print $4}' | tail -n 1 | sed "s/,//g")
    clean_reads=$(cat ${out_dir}/${fname}/01-clean-data/${fname}_C.log | \
        grep "Pairs written" | awk '{print $(NF-1)}' | sed 's/,//g')
    genome_input=$(cat ${out_dir}/${fname}/02-align/${fname}_genome_Log.final.out | \
        grep "Number of input reads" | awk '{print $NF}')
    rRNA_reads=$(echo ${clean_reads} ${genome_input} | \
        awk '{print $1-$2}')
    Unmap_reads=$(seqkit stat ${out_dir}/${fname}/02-align/${fname}_genome_Unmapped.out.mate1 | \
        awk '/DNA/{print $4}' | sed 's/,//g')
    genome_reads=$(echo ${genome_input} ${Unmap_reads} | \
        awk '{print $1-$2}')
    
    # add header
    echo "name,total_reads,clean_reads,rRNA_reads,genome_reads,Unmap_reads,rRNA_pct,genome_pct,unmap_pct,gene_numbers" \
        > ${out_dir}/${fname}/${fname}.stat
    echo "${fname} ${total_reads} ${clean_reads} ${rRNA_reads} ${genome_reads} ${Unmap_reads}" | \
        awk '{printf "%s,%d,%d,%d,%d,%d,%.2f,%.2f,%.2f\n",$1,$2,$3,$4,$5,$6,$4/$3,$5/$3,$6/$3}' \
        >> ${out_dir}/${fname}/${fname}.stat
    # trans to tab
    cat ${out_dir}/${fname}/${fname}.stat | tr "," "\t" > ${reporte_summary}
    [[ -f ${out_dir}/${fname}/${fname}.stat ]] && rm ${out_dir}/${fname}/${fname}.stat
fi


