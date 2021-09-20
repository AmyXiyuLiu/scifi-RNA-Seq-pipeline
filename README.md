# scifi-RNA-Seq-pipeline
Lab in-build pipeline
   
**Remember to use the root folder of the illumina run for BaseCalls (BCL files), containing the runParameters.xml file.**
  
***
## Description  
This pipeline is modified on the basis of the existing scifiRNA-seq pipeline (https://github.com/epigen/scifiRNA-seq) to be suitable for processing our in-lab sequencing data. The input data of this pipeline is the illumina run BCL files, and the output of this pipeline is a gene expression reads counting tsv file which can directly use for the downstream analysis.    
## Datasets  
**`BCL files`**: /data/sci-fi/210914_Zhang-scifiRNA-QC/MiSeqOutput/210914_M00314_0072_000000000-G992H    
**`Genome fasta`**: /afs/crc.nd.edu/user/x/xliu32/refdata-gex-mm10-2020-A/fasta/genome.fa    
**`Genome gtf`**: /afs/crc.nd.edu/user/x/xliu32/refdata-gex-mm10-2020-A/genes/genes.gtf    
## Required tools
**`STAR`**: /afs/crc.nd.edu/user/x/xliu32/STAR-2.7.4a/bin/Linux_x86_64/STAR    
**`featureCounts`**: /afs/crc.nd.edu/user/x/xliu32/subread-2.0.3-Linux-x86_64/bin/featureCounts    
**`Picard`**: /afs/crc.nd.edu/user/x/xliu32/picard-2.19.2-CeMM-all.jar (https://github.com/DanieleBarreca/picard/releases/2.19.2-CeMM)    
**`scifi`**: /afs/crc.nd.edu/user/x/xliu32/anaconda3/bin/scifi     
     
***
## Running the pipeline         
### Demultiplex the scifi-RNAseq data  
The first steps reads the Illumina raw data and creates an unaligned, undemultiplexed bam file:  
```
java \
-Xmx20G \
-Djava.util.concurrent.ForkJoinPool.common.parallelism=2 \
-Djava.io.tmpdir=./tmp \
-jar ${picard_jar} \
IlluminaBasecallsToMultiplexSam \
RUN_DIR= ${illumina_run_folder} \
LANE=1 \
OUTPUT= ${undemultiplexed_file} \
SEQUENCING_CENTER=BSF \
NUM_PROCESSORS=2 \
APPLY_EAMSS_FILTER=false \
INCLUDE_NON_PF_READS=false \
TMP_DIR=tmp \
CREATE_MD5_FILE=false \
FORCE_GC=false \
MAX_READS_IN_RAM_PER_TILE=9000000 \
MINIMUM_QUALITY=2 \
VERBOSITY=INFO \
QUIET=false \
VALIDATION_STRINGENCY=STRICT \
CREATE_INDEX=false \
GA4GH_CLIENT_SECRETS=client_secrets.json  
```
Then, you should prepare a sample_annotation csv file for the next step. You can see that for each sample and well there is one line and barcode 1 is the round 1 (well-specifc) barcode, while barcode 2 is the 701/702 sample index.    

The second step performs the actual demultiplexing:  
```
java \
-Xmx20G \
-Djava.io.tmpdir=./tmp \
-jar ${picard_jar} \
IlluminaSamDemux \
INPUT= ${undemultiplexed_file} \
OUTPUT_DIR= ${output_dir} \
OUTPUT_PREFIX= ${output_prefix} \
LIBRARY_PARAMS= ${sample_annotation} \
METRICS_FILE= ${output_metrics_file} \
TMP_DIR=./tmp \
COMPRESSION_LEVEL=9 \
CREATE_MD5_FILE=true \
OUTPUT_FORMAT=bam \
BARCODE_TAG_NAME=BC \
BARCODE_QUALITY_TAG_NAME=QT \
MAX_MISMATCHES=1 \
MIN_MISMATCH_DELTA=1 \
MAX_NO_CALLS=2 \
MINIMUM_BASE_QUALITY=0 \
VERBOSITY=INFO \
QUIET=false \
VALIDATION_STRINGENCY=STRICT \
MAX_RECORDS_IN_RAM=500000 \
CREATE_INDEX=false \
GA4GH_CLIENT_SECRETS=client_secrets.json \
USE_JDK_DEFLATER=false \
USE_JDK_INFLATER=false \
DEFLATER_THREADS=4 \
MATCHING_THREADS=4 \
READ_STRUCTURE= 8M13B8S8B16M70T \
TAG_PER_MOLECULAR_INDEX=RX \
TAG_PER_MOLECULAR_INDEX=r2
```
The demultiplexer will split this into the structure 8M13B5S8B16M81T:    
 * The first 8 bases are the UMI (tag RX)
 * The next 13 bases are used as a barcode (round 1 * well-specific barcode)
 * The next 5 are skipped
 * The next 8 bases are also used as a barcode (sample specific barcode)
 * The next 16 bases are marked as UMI and are put into the r2 tag (round 2 – 10x barcode)
 * The last 81 bases are the transcriptome read

The outputs of this step are bam files for each well. eg: ***SCIFI_nuclei_BB2182_A01.bam***    

### Alignment  
Here, we should prepare one **annotation.csv** file and one **Round1_SCIFI_nuclei_BB2182.csv** file.     
```
scifi map -c ~/.scifiRNA-seq.config.yaml --input-bam-glob /afs/crc.nd.edu/user/x/xliu32/scifi_map/SCIFI_nuclei_BB2182_A01.bam /afs/crc.nd.edu/user/x/xliu32/scifi_map/annotation.csv
```   
The outputs of this step are aligned bam files for each well. eg: ***SCIFI_nuclei_BB2182_A01_ALL_STAR_Aligned_out.bam***    

### Separation
#### Separate a well's bam file into the bam of each single cell   
* Extract round 1 and round 2 barcodes from bam        
```
samtools view -h SCIFI_nuclei_BB2182_A01_ALL_STAR_Aligned_out.bam > SCIFI_nuclei_BB2182_A01.txt    
sed -n 's/.*\b\(r2:Z:[[:alnum:]]*\).*/\1/'p SCIFI_nuclei_BB2182_A01.txt > r2.txt    
sed -n 's/.*\b\(BC:Z:[[:alnum:]]*\).*/\1/'p SCIFI_nuclei_BB2182_A01.txt > r1.txt    
paste r2.txt r1.txt > SCIFI_nuclei_BB2182_A01_ALL_STAR_Aligned_out.txt    
perl single_cell.pl SCIFI_nuclei_BB2182_A01_ALL_STAR_Aligned_out.txt    
mkdir single_cell
bash res.sh
```
The **`res.sh`** is an executable file containing command lines to separate bam file per well to bam file per cell using round 1 and round 2 barcodes.   
```
samtools view SCIFI_nuclei_BB2182_A01_ALL_STAR_Aligned_out.bam -h | sed -n "/r2:Z:GTCGTAATCCAGGTAT\tBC:Z:AAGTGATTAGCAATAAGGCGA/p" | samtools view -Sb -T /afs/crc.nd.edu/user/x/xliu32/refdata-gex-mm10-2020-A/fasta/genome.fa > ./single_cell/BC1.bam   
```
The outputs of this step are bam files per cell located in the single_cell folder named **BC${counts}.bam**   

### featureCounts
Generate one gene expression read counts tsv file   
```
/afs/crc.nd.edu/user/x/xliu32/subread-2.0.3-Linux-x86_64/bin/featureCounts -T 4 -F GTF -t exon -g gene_id --extraAttributes gene_name -Q 30 -s 0 -R BAM -a /afs/crc.nd.edu/user/x/xliu32/refdata-gex-mm10-2020-A/genes/genes.gtf -o /afs/crc.nd.edu/user/x/xliu32/featureCounts_output/ALL.STAR.featureCounts.quant_gene.exon.tsv *.bam   
```
The output is one tsv file called **`ALL.STAR.featureCounts.quant_gene.exon.tsv`**. Each column represents one single cell.   
You can quickly check any read counts of any single cell.    
```
awk '{if(${column number of your interested single cell}>0) print $0}' ALL.STAR.featureCounts.quant_gene.exon.tsv    
```

All Done!  
