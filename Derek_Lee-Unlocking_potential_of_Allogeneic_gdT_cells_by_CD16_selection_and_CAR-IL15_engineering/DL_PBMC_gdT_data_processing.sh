# Creator: Wenbin Guo <wbguo[AT]ucla[DOT]edu>
# Date: 2023-06-20
# Description:
# this script is used to trim the raw reads, align them to the reference genome,        \
# and quantify the gene expression for each sample.                                     \
# Each tool's version can be found in the README file or from the original paper titled:\
#       Unlocking the Potential of Allogeneic Vδ2 T Cells for Ovarian Cancer Therapy    \
#       through CD16 Biomarker Selection and CAR/IL-15 Engineering                      \
# 
# the input: 
# - the reference genome and annotation file
# - the fastq files 
# - a job array file that contains sample ID and file names (format described below)
# the output:
# - trimmed fastq files and aligned bam files 
# - the raw gene expression count matrix


# set path for quality control and trimming
working_path=$HOME/project-liliyang/bulk/
raw_reads_path=$working_path/raw_data/
trm_reads_path=$working_path/reads/

log_preQC_path=$working_path/docs/preQC/
log_postQC_path=$working_path/docs/postQC/
log_trimm_path=$working_path/docs/fastp/

# set reference and annotation file, set the job array file
gtf_file=${working_path}/ref/gencode.v32.primary_assembly.annotation.gtf
fasta_file=${working_path}/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
ja=$working_path/jobs_list
# ja is the job array file, each row is a line which has a format as `ID,file_name`
# for example: 
#Number5,Number5_S16_R1_001.fastq.gz
#Number6,Number6_S17_R1_001.fastq.gz


### step 1: build STAR index
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir ${working_path}/STAR_idx \
    --genomeFastaFiles $fasta_file \
    --sjdbGTFfile $gtf_file \
    --sjdbOverhang 49

# the following lines is used when submitting jobs to the hoffman2 computing cluster@UCLA
# one can also use a for loop to replace it if they don't have access to computing clusters
PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
IFS=',' read -r -a PARMS_arr <<< "$PARMS"
prefix=${PARMS_arr[0]}
fastq_file=${PARMS_arr[1]}
fastq_base=$(echo $fastq_file | cut -d. -f1)


### step 2: quality control and trim adpaters (run for each sample)
raw_file=$raw_reads_path/$fastq_file
new_file=$trm_reads_path/${fastq_base}_trimmed.fastq.gz

echo "raw file $raw_file"
date '+%Y-%m-%d %H:%M:%S'

fastqc "$raw_file" -o "$log_preQC_path" -t 16
echo "preQC Finished"
date '+%Y-%m-%d %H:%M:%S'

fastp -i $raw_file -o $new_file -h $log_trimm_path/${prefix}.html -j $log_trimm_path/${prefix}.json
echo "fastp Finished"
date '+%Y-%m-%d %H:%M:%S'

fastqc $new_file -o $log_postQC_path -t 16
echo "postQC Finished"
date '+%Y-%m-%d %H:%M:%S'


### step 3: align (run for each sample)
echo "This is subject $prefix"
date '+%Y-%m-%d %H:%M:%S'

STAR --runThreadN 16 --runMode alignReads --genomeDir $working_path/STAR_idx \
--readFilesIn $working_path/reads/${fastq_base}_trimmed.fastq.gz \
--outFileNamePrefix $working_path/bam/$prefix/ \
--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMultimapNmax 20 \
--readFilesCommand zcat --outReadsUnmapped Fastx_failed --outSAMtype BAM SortedByCoordinate

mv $working_path/bam/$prefix/Aligned.sortedByCoord.out.bam $working_path/bam/${prefix}.sorted.bam
samtools index $working_path/bam/${prefix}.sorted.bam

echo "Alignment and bam index Finished!"s
date '+%Y-%m-%d %H:%M:%S'


### step4: gene expression quantification
out_file=$working_path/data/expr/featureCounts.txt

cd $working_path/bam
job_list=$(cat $ja |cut -d, -f1 | xargs printf "%s.sorted.bam ")
featureCounts -T 16 -a $gtf_file -o $out_file ${job_list[@]}
