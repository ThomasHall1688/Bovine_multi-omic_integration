########################################################################################
##########################           FastQC            #################################
########################################################################################

#Check md5sums to see if files were dowloaded properly 
md5sum -c md5sum.txt > md5sum_check.txt #should return file name and OK 

#Create fastqc folder and run fastqc on everything
mkdir fastqc
fastqc *.gz

########################################################################################
##########################       STAR alignment        #################################
########################################################################################


#Work was intially done on Stampede due to the Rodeo failure

#genome = /home/thall/scratch/RNAseq/Genomes/UCD_1.2
#reads = /home/thall/scratch/RNAseq/Alv_mac-data/raw_reads/deconv_raw_reads
#Genome file = GCF_002263795.1_ARS-UCD1.2_genomic.fna
#Annotation file = GCF_002263795.1_ARS-UCD1.2_genomic.gff


########## Make index ##########

#--runThreadN: number of threads
#--runMode: genomeGenerate mode
#--genomeDir: /path/to/store/genome_indices
#--genomeFastaFiles: /path/to/FASTA_file
#--sjdbGTFfile: /path/to/GTF_file
#--sjdbOverhang: readlength -1

STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir $HOME/scratch/RNAseq/Genomes/UCD_1.2 \
--genomeFastaFiles \
$HOME/scratch/RNAseq/Genomes/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
--sjdbGTFfile $HOME/scratch/RNAseq/Genomes/GCF_002263795.1_ARS-UCD1.2_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 \
--outFileNamePrefix \
$HOME/scratch/RNAseq/Genomes/ARS-UCD1.2 &

#or in the actual directory. I think there is a weird directory format on Stampede. 
STAR --runThreadN 20 --runMode genomeGenerate \
--genomeDir . \
--genomeFastaFiles \
./GCF_002263795.1_ARS-UCD1.2_genomic.fna \
--sjdbGTFfile ./GCF_002263795.1_ARS-UCD1.2_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 99 \
--outFileNamePrefix \
/home/thall/scratch/RNAseq/Genomes/ARS-UCD1.2 &


########## Mapping ##########

#--runThreadN NumberOfThreads
#--genomeDir /path/to/genomeDir
#--readFilesIn /path/to/read1 [/path/to/read2 ]


#Mapping all raw files 
for i in *_pe1.fastq.gz; do
STAR --runMode alignReads  --genomeLoad  LoadAndKeep --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
--outReadsUnmapped Fastx  --genomeDir $HOME/scratch/RNAseq/Genomes/UCD_1.2 --readFilesIn $HOME/scratch/RNAseq/Alv_mac-data/raw_reads/deconv_raw_reads/$i $HOME/scratch/RNAseq/Alv_mac-data/raw_reads/deconv_raw_reads/${i%_pe1.fastq.gz}_pe2.fastq.gz \
--runThreadN 10 --outFileNamePrefix $HOME/scratch/RNAseq/Alv_mac-data/new_alignment/${i%_pe1.fastq.gz}_
done

for i in *.bam; do
featureCounts -a \
$HOME/scratch/RNAseq/Genomes/UCD_1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff \
-B -p -C -R -s 1 -T 15 -t gene -g Dbxref -o ./counts.txt \
$HOME/scratch/RNAseq/Alv_mac-data/new_alignment/$i
done

#mapping one

STAR --runMode alignReads  --genomeLoad  LoadAndKeep --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 30000000000 \
--outReadsUnmapped Fastx  --genomeDir $HOME/scratch/RNAseq/Genomes/UCD_1.2 --readFilesIn $HOME/scratch/RNAseq/Alv_mac-data/raw_reads/deconv_raw_reads/N1178_CN_48H_pe1.fastq.gz $HOME/scratch/RNAseq/Alv_mac-data/raw_reads/deconv_raw_reads/N1178_CN_48H_pe2.fastq.gz \
--runThreadN 10 --outFileNamePrefix $HOME/scratch/RNAseq/Alv_mac-data/new_alignment/unmapped_reads/N1178_CN_48H
done

#mapping one (unsorted for test with feature counts. Sorted = 79.0% unsorted = 79.0%  for N1178_CN_48)

STAR --runMode alignReads  --genomeLoad  LoadAndKeep --readFilesCommand zcat --outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx  --genomeDir $HOME/scratch/RNAseq/Genomes/UCD_1.2 --readFilesIn $HOME/scratch/RNAseq/Alv_mac-data/raw_reads/deconv_raw_reads/N1178_CN_48H_pe1.fastq.gz $HOME/scratch/RNAseq/Alv_mac-data/raw_reads/deconv_raw_reads/N1178_CN_48H_pe2.fastq.gz \
--runThreadN 10 --outFileNamePrefix $HOME/scratch/RNAseq/Alv_mac-data/new_alignment/unmapped_reads/N1178_CN_48H
done


########################################################################################
##########################        FeatureCounts        #################################
########################################################################################

#input files:									Give the names of input read files that include the read mapping results. The program automatically detects the file format (SAM or BAM). Multiple files can be provided at the same time.
#-a < input > (annot.ext, annot.inbuilt):		Give the name of an annotation file.
#-A (chrAliases) :								Give the name of a file that contains aliases of chromosome names. The file should be a comma delimited text file that includes two columns. The first column gives the chromosome
#												names used in the annotation and the second column gives the chromosome names used by reads. This file should not contain header lines. Names included in this file are case sensitive.
#-B (requireBothEndsMapped)						If specified, only fragments that have both ends successfully aligned will be considered for summarization. This option should be used together with -p (or is PairedEnd in Rsubread featureCounts).
#-C (countChimericFragments)					If specified, the chimeric fragments (those fragments that have their two ends aligned to different chromosomes) will NOT be counted. This option should be used together with -p (or isPairedEnd in Rsubread featureCounts).
# -d < int > (minFragLength)					Minimum fragment/template length, 50 by default. This option must be used together with -p and -P.
#-D < int > (maxFragLength)						Maximum fragment/template length, 600 by default. This option must be used together with -p and -P.
#-f (useMetaFeatures)							If specified, read summarization will be performed at feature level (eg. exon level). Otherwise, it is performed at metafeature level (eg. gene level).
#-F (isGTFAnnotationFile)						Specify the format of the annotation file. Acceptable formats include ‘GTF’ and ‘SAF’ (see Section 6.2.2 for details). The C version of featureCounts program uses a 
#												GTF format annotation by default, but the R version uses a SAF format annotation by default. The R version also includes in-built annotations.
#-g < input > (GTF.attrType)					Specify the attribute type used to group features (eg. exons) into meta-features (eg. genes) when GTF annotation is provided. ‘gene id’ by default. This attribute type is usually the
#												gene identifier. This argument is useful for the meta-feature level summarization.
#-M (countMultiMappingReads)					If specified, multi-mapping reads/fragments will be counted. A multi-mapping read will be counted up to N times if it has N reported mapping locations. The program uses the ‘NH’
#												tag to find multi-mapping reads.
#-o < input > 									Give the name of the output file. The output file contains the number of reads assigned to each meta-feature (or each feature if -f is specified). Note that the featureCounts function
#												in Rsubread does not use this parameter. It returns a list object including read summarization results and other data.
#-O (allowMultiOverlap) 						If specified, reads (or fragments if -p is specified) will be allowed to be assigned to more than one matched meta-feature (or feature if -f is specified). Reads/fragments overlapping
#												with more than one meta-feature/feature will be counted more than once. Note that when performing meta-feature level summarization, a read (or fragment) will still be counted once
#												if it overlaps with multiple features belonging to the same meta-feature but does not overlap with other meta-features.
#-p (isPairedEnd)								If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads.
#-P (checkFragLength)							If specified, the fragment length will be checked when assigning fragments to meta-features or features. This option must be used together with -p. The fragment length thresholds should be specified using -d and -D options.
#-Q < int > (minMQS)							The minimum mapping quality score a read must satisfy in order to be counted. For paired-end reads, at least one end should satisfy this criteria. 0 by default.
#-R 											Output read assignment results for each read (or fragment if paired end). They are saved to a tab-delimited file that contains four columns including read name, status(assigned or the reason if not assigned), name of target feature/metafeature and number of hits if the read/fragment is counted
#												multiple times. Name of the file is the input file name added with a suffix ‘.featureCounts’.
#-s < int > (isStrandSpecific)					Indicate if strand-specific read counting should be performed. It has three possible values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). 0 by default. For paired-end reads, strand of the first read is taken as the strand of the whole
#												fragment and FLAG field of the current read is used to tell if it is the first read in the fragment.
#-t < input > (GTF.featureType) 				Specify the feature type. Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting. ‘exon’ by default.
#-T < int > (nthreads)							Number of the threads. The value should be between 1 and 32. 1 by default.
#-v 											Output version of the program.
#−−countSplitAlignmentsOnly (countSplitAlignmentsOnly)	If specified, only split alignments (CIGAR strings containing letter ‘N’) will be counted. All the other alignments will be ignored. An example of split alignments is the exon-spanning reads in RNA-seq data. If exon-spanning reads need to be
#												assigned to all their overlapping exons, ‘-f’ and ‘-O’ options should be provided as well.
#−−donotsort 	If specified, paired end reads will not be re-ordered even if reads from the same pair were found not to be next to each other in the input.
#−−ignoreDup (ignoreDup)						If specified, reads that were marked as duplicates will be ignored. Bit Ox400 in FLAG field of SAM/BAM file is used for identifying duplicate reads. In paired end data, the entire
# 												read pair will be ignored if at least one end is found to be a duplicate read.
#−−primary (countPrimaryAlignmentsOnly)			If specified, only primary alignments will be counted. Primary and secondary alignments are identified using bit 0x100 in the Flag field of SAM/BAM files. All primary alignments
# 												in a dataset will be counted no matter they are from multimapping reads or not (ie. ‘-M’ is ignored).
#−−minReadOverlap < int > (minReadOverlap)		Specify the minimum number of overlapped bases required to	assign a read to a feature. 1 by default. Negative values are
#												permitted, indicating a gap being allowed between a read and a feature.
#−−read2pos < int > (read2pos)					The read is reduced to its 5’ most base or 3’ most base. Read summarization is then performed based on the single base position to which the read is reduced. By default, no read
#												reduction will be performed.
#−−readExtension5 < int > (readExtension5)		Reads are extended upstream by < int > bases from their 5’ end. 0 by default.
#−−readExtension3 < int > (readExtension3) 		Reads are extended downstream by < int > bases from their 3’ end. 0 by default.

#This needs to be applied to all .bam files. 
featureCounts -a \
$HOME/scratch/RNAseq/Genomes/UCD_1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff \
-B -p -C -R -s 1 -T 15 -t gene -g Dbxref -o ./counts.txt \
$HOME/scratch/RNAseq/Alv_mac-data/new_alignment/N1178_CN_48H_Aligned.sortedByCoord.out.bam
