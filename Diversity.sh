#!/bin/bash
## Hepacivirus diversity analysis


# path to folder with trimmed fastq data
# Folder per individual including all time points
path= # path to your data

# define a output directory
out= # output folder

## Location of meta data (md) folder and tools
md= # reference and ORF data is deposited in this folder
TVx= # TVx = individuum id - will be updated in the loop
ref= # $md/${TVx}/{yourFastaConsensus}.fa

## path to tools & parameters
SAM2CONSENSUS= # path to samtoconsensus tool


# diversi tools
# path to
diversi=# bin/diversiutils.pl
filter=# bin/diversifilter.pl

# path to ORF information
div_ORF=$md/${TVx}/TV_CodingRegion.txt

# vnvs
dnds=# path to DNDS.jar
dnds_ORF=$md/${TVx}/TV_CodingRegiondnds.txt

# CliqueSNV
clique= # path to CliqueSNV .jar
cliqueStart= # ORF min size
cliqueEnd= # ORF max size


## start analysis
for folder in $path/*
	do
		for file in $folder/*
			do
			# select paired end reads per sample
			if [[ $file == *"fwd.paired.fq"* ]]; then
				filename=$(basename "$file" g | cut -d'.' -f1)
				echo $folder
				echo $file
				echo $filename
				echo "==== mapping to Baseline ==="
				## mapping via tanoti
				echo tanoti -r $ref -i $folder $path/${filename}.rev.paired.fq -p 1 -o $filename.mapped.sam;
				echo
				# echo "process data via samtools"
				echo "==="
				# # data processing
				samtools view -bS $filename.mapped.sam > $filename.mapped1.bam;
				samtools view -h -F 4 -b $filename.mapped1.bam > $filename.mapped.bam;
				samtools sort $filename.mapped.bam -o $filename.mapped.sorted.bam;
				samtools index $filename.mapped.sorted.bam

				echo
				# #dedup data
				echo "Remove duplicates via MarkDuplicates"
				echo "==="
				gatk MarkDuplicates -I $filename.mapped.sorted.bam -M $filename.dedup_1.txt -O $filename.mapped.sorted_dedup.bam --REMOVE_DUPLICATES true -MAX_SEQS 50000 -MAX_FILE_HANDLES 8000 --SORTING_COLLECTION_SIZE_RATIO 0.25 --REMOVE_SEQUENCING_DUPLICATES false --TAGGING_POLICY DontTag -AS false --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --USE_JDK_DEFLATER true --USE_JDK_INFLATER true
				samtools index $filename.mapped.sorted_dedup.bam
				echo
				echo "generate consensus sequence via SAM2CONSENSUS"
				echo "==="
				samtools view -h $filename.mapped.sorted_dedup.bam > $filename.mapped.sorted_dedup.sam
				samtools flagstats $filename.mapped.sorted_dedup.bam > $filename.mapped.sorted_dedup_flagstat.txt
				echo
				# rename .fa header
				$SAM2CONSENSUS -i $filename.mapped.sorted_dedup.sam -o $filename.samtocon.fa
				sed -i '' "1s/.*/>ref/" $filename.samtocon.fa
				head -n 1 $filename.samtocon.fa
				echo
				echo "===="
				echo "Gerating haplotypes via CliqueSNV"
				java -Xmx7g -jar $clique -m snv-illumina -os $cliqueStart -in $filename.mapped.sorted_dedup.sam -oe $cliqueEnd -outDir $filename.toref.cliqueSNV
				echo
				echo "==="
				echo "Diversi tool"
				echo $div_ORF
				perl $diversi -bam $filename.mapped.sorted_dedup.bam -ref $ref -orfs $div_ORF -stub $filename.toref_dedup_E1E2;
				perl $filter -in $filename.toref_dedup_ORF1_RdRp -pQ 0.05 -pS 0.05 -stub $filename.toref_dedup_E1E2_fil;
				echo "==="
				echo "dnds tool"
				echo
				echo $dnds_ORF
				java -jar $dnds $ref $dnds_ORF $filename.mapped.sorted_dedup.sam 3 > $filename.dnds_toref_dedup.txt
fi
done
done
