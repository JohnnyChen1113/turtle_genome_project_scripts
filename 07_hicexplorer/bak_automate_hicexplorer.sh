#!/usr/bin/env bash
# To run the hicexplorer automatically
# Author: Junhao Chen
# Date: 2022-01-30
# Show the uasge
if [ -z ${1} ]
then
  	echo "Usage: nohup bash ${0} <genome_file> > out.log &"
	exit $E_BADARGS
fi

# Sort the genome

if [ ! -f ${1%.*}/${1%.*}.sorted.fasta ]
then
seqkit sort -r -l ${1} -o ${1%.*}/${1%.*}.sorted.fasta
else
  echo "The ${1} already sorted!"
fi
# check and generate the folder
if [ ! -d ${1%.*} ]
then
  mkdir ${1%.*}
else
  echo "The ${1%.*} folder already exist!"
fi
# Generate the enzyme files
if [ ! -f ${1%.*}/TTAA.bed ]
then
for i in `cat /mnt/e8/01turtle_genome_assembly/07_hicexplorer/enzyme.txt` 
do 
        hicFindRestSite --fasta ${1%.*}/${1%.*}.sorted.fasta --searchPattern ${i} -o ${1%.*}/${i}.bed
done
else
  echo "The enzyme files already done!"
fi

# Build the index
if [ ! -f ${1%.*}/${1%.*}.bwt ]
then
bwa index -p ${1%.*}/${1%.*} ${1%.*}/${1%.*}.sorted.fasta
else
  echo "The bwa index alreadt exist!"
fi

# Mapping with bwa
if [ ! -f ${1%.*}/1.bam ]
then
  bwa mem -A 1 -B 4 -E 50 -L 0 -t 16 ${1%.*}/${1%.*} 1_R1.fastq.gz | samtools view -Shb - > ${1%.*}/1.bam 
else
  echo "The 1.bam already mapped!"
fi
if [ ! -f ${1%.*}/2.bam ]
then
  bwa mem -A 1 -B 4 -E 50 -L 0 -t 16 ${1%.*}/${1%.*} 1_R2.fastq.gz | samtools view -Shb - > ${1%.*}/2.bam
else
  echo "The 2.bam already mapped!"
fi
# Build the hicMatrix
if [ ! -f ${1%.*}/${1%.*}.h5 ]
then
mkdir ${1%.*}/hicQC 

hicBuildMatrix \
--samFiles ${1%.*}/1.bam ${1%.*}/2.bam \
--binSize 10000 \
--restrictionCutFile ${1%.*}/CT.AG.bed  ${1%.*}/GA.TC.bed  ${1%.*}/GATC.bed ${1%.*}/TTAA.bed \
--restrictionSequence GATC GA.TC CT.AG TTAA \
--danglingSequence GATC A.T T.A TA \
--outFileName ${1%.*}/${1%.*}.h5 \
--outBam ${1%.*}/${1%.*}.bam \
--threads 8 \
--QCfolder ${1%.*}/hicQC \
--inputBufferSize 400000
fi

# Merge matrix bins for plotting
if [ ! -f ${1%.*}/${1%.*}.100bins.h5 ]
then
hicMergeMatrixBins \
--matrix ${1%.*}/${1%.*}.h5 \
--numBins 100 \
--outFileName ${1%.*}/${1%.*}.100bins.h5
fi

# Plot the corrected Hi-C matrix
if [ ! -f ${1%.*}/${1%.*}_1Mb_matrix.png ]
then
hicPlotMatrix \
--matrix ${1%.*}/${1%.*}.100bins.h5 \
--log1p \
--dpi 300 \
--clearMaskedBins \
--colorMap jet \
--title "Hi-C matrix for ${1}" \
--outFileName ${1%.*}/${1%.*}_1Mb_matrix.png
fi

# Plot the logest 24 contigs
# Get the ids of this contigs
cat ${1%.*}/${1%.*}.sorted.fasta | grep ">" | head -25 | tr -d '>' | paste -s -d ' ' > ${1%.*}/logest25.id && \
if [ ! -f ${1%.*}/${1%.*}_logest25_1Mb_matrix.png ]
then
hicPlotMatrix \
--matrix ${1%.*}/${1%.*}.100bins.h5 \
--log1p \
--dpi 300 \
--chromosomeOrder `cat ${1%.*}/logest25.id` \
--clearMaskedBins \
--colorMap jet \
--title "Logest25 scaffolds Hi-C matrix for ${1}" \
--outFileName ${1%.*}/${1%.*}_logest25_1Mb_matrix.png
fi


