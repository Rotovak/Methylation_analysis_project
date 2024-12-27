# Define variables
REFERENCE_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr17.fa.gz"
ACCESSIONS=("SRR1272802" "SRR1272799" "SRR1272795" "SRR1272791")
REFERENCE_FILE="/chr17.fa"

echo "Workflow started on $(date)"

# Step 1: Download SRA files
echo "Downloading SRA files..."
for ACCESSION in "${ACCESSIONS[@]}"; do
    prefetch $ACCESSION
done

# Step 2: Convert SRA to FASTQ
echo "Converting SRA to FASTQ..."
for ACCESSION in "${ACCESSIONS[@]}"; do
    fastq-dump $ACCESSION.sra
done

echo "Indexing reference genome..."
bwa index $REFERENCE_FILE

# Step 3: Align FASTQ files to reference genome
echo "Aligning reads to reference genome..."
for ACCESSION in "${ACCESSIONS[@]}"; do
    bwa mem -t 6 $REFERENCE_FILE ${ACCESSION}_1.fastq > $ACCESSION.sam
done

# Step 4: Convert SAM to BAM, sort, and index
echo "Processing SAM files..."
for ACCESSION in "${ACCESSIONS[@]}"; do
    samtools view -bS $ACCESSION.sam > $ACCESSION.bam
    samtools sort $ACCESSION.bam -o ${ACCESSION}_sorted.bam
    samtools index ${ACCESSION}_sorted.bam
done

# Step 5: Extract methylation counts using mpileup
echo "Extracting methylation counts with samtools mpileup..."
for ACCESSION in "${ACCESSIONS[@]}"; do
    samtools mpileup -f $REFERENCE_FILE ${ACCESSION}_sorted.bam | cut -f1-3,5 > ${ACCESSION}_methylation_counts.txt
done

echo "Workflow completed on $(date)"
