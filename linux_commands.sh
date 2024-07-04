### Linux commands 
# fastqc version=0.12.1
# featureCounts version=2.0.6
# STAR version=2.7.11b
# samtools version=1.18
# bioawk version=1.0

ID='string_for_analysis_identifier'
GENOME=heligmosomoides_polygyrus.PRJEB15396.WBPS18.genomic.fa
ANNOT=heligmosomoides_polygyrus.PRJEB15396.WBPS18.annotations.gff3
TRANSC=heligmosomoides_polygyrus.PRJEB15396.WBPS18.transcripts.fasta
RPATH=path/to/rnaFiles/
GPATH=path/to/genomeFiles/

# Create index to map reads
STAR --runThreadN 14 --runMode genomeGenerate --genomeDir ${GPATH} --genomeFastaFiles ${GENOME} \
    --sjdbGTFfile ${ANNOT} --genomeSAindexNbases 13 --sjdbGTFtagExonParentTranscript Parent

# List of prefixes for each sample
RPREFS=$( ls ${RPATH}/*.fastq.gz | awk -F '/' '{print $NF}' | cut -f 1 -d '_' | sort -u ) 

# Map reads to genome
for RPREF in ${RPREFS}; do
    STAR --runThreadN 14 --genomeDir ${GPATH} \
        --readFilesIn ${RPATH}${RPREF}_L001_R1.fastq.gz,${RPATH}${RPREF}_L002_R1.fastq.gz \
        ${RPATH}${RPREF}_L001_R2.fastq.gz,${RPATH}${RPREF}_L002_R2.fastq.gz \
        --readFilesCommand zcat --outFileNamePrefix ${RPREF}_ \
        --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted
done

# Check bam files 
samtools quickcheck -v *.bam && echo 'bam files ok' || echo 'above file(s) failed check'

# Count reads 
S=$1 #strandedness, 0
D=$2 #min fragment length, 35
featureCounts -p --countReadPairs -T 14 -g Parent -s ${S} -d ${D} -a ${ANNOT} -o ${ID}featureCounts_d${D}s${S}.txt \
    -O --fraction -B -C */*Aligned.out.bam 2> ${ID}featureCounts_d${D}s${S}.screenoutput.txt

# Get gene lengths
bioawk -c fastx '{ print $name, length($seq) }' ${TRANSC}

## Create counts matrix for DESeq2 input
cut -f 1,7,8,9,10,11,12,13,14,15 ${ID}featureCounts_d${D}s${S}.txt | sed 's,path/to/rnaFiles/,,g' | sed 's/_Aligned.out.bam//g' | sed -n '1!p' > hbak_PRJEB15396.WBPS18_featurecounts.Rmatrix.tsv
