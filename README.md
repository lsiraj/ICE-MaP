# ICE-MaP
## About
Here we introduce ICE-MaP, Inosine CyanoEthylation and Mutational Profiling, building upon previous experimental methods for detecting inosines. ICE-MaP can accurately quantify inosines at single-base resolution, down to 10% frequency. Moreover, ICE-MaP can detect up to 5 inosines at a 13-base editing hotspot, and provides phase data for disentangling multiple editing patterns. ICE-MaP is compatible with a wider range of library construction and sequencing methods, including long read technology, making it possible to phase multiple  hotspots even hundreds of bases apart. 
## Data
Raw data is submitted to SRA (submission: SUB15494827)
## Processing
Processing of the data is done with the align.sh script, which takes .bam files and aligns them to a user-provided index:

```align.sh [-h] [-o output dir -f full file name for input bam -b base name -i index] -- program to align ICEMaP data ```
Where you specify output directory, file name, the base name of the file you are aligning (will be carried through for naming of the intermediate and final output files), and the index. We will walk through an example here:
### Step 0 (optional)
Create a list of file names (useful to align multiple .bam files):

```ls -l ../MiSeq/*.bam | grep -v mapped | awk '{print $9}' > filelist.txt```
### Step 1
Create fastq of the .bam files:

```bedtools bamtofastq -i ${file} -fq ${dir}/${base}.1.fastq -fq2 ${dir}/${base}.2.fastq```
### Step 2
Merge reads using [flash2](https://github.com/dstreett/FLASH2):

```flash2 -O -m 40 -M 260 -o ${dir}/${base} ${dir}/${base}.1.fastq ${dir}/${base}.2.fastq```
### Step 3
Build a [bowtie2](https://github.com/BenLangmead/bowtie2) index. Here, we are using the 5HT2cR gene as our reference, with N's at the editing sites:

```bowtie2-build ~/ICEMaP/data/references/5HT2cR-ligated-5N.fa ~/ICEMaP/data/references/5HT2cR-ligated-5N ```
### Step 4 
Align:

```bowtie2 --very-sensitive-local --no-unal -x ${index} -U ${dir}/${base}.extendedFrags.fastq -S ${dir}/${base}.mapped.raw.sam 2> ${dir}/${base}.MapStats.txt```

### Step 5
Convert from .bam to .sam:

```samtools view -h -b -S -o ${dir}/${base}.mapped.raw.bam ${dir}/${base}.mapped.raw.sam```
### Step 6
Sort:

```samtools sort -o ${dir}/${base}.mapped.sorted.bam ${dir}/${base}.mapped.raw.bam```

### Post-processing: piledriver
After mapping and sorting the reads, count mutations using [piledriver](https://github.com/arq5x/piledriver):

```bamtools piledriver -fasta ~/ICEMaP/data/references//5HT2cR-ligated-5N.fa -in ${base}.sorted.bam -out ${base}.piledriver.txt```

## Mutation proportions
```coverage.R``` estimates the coverage along the oligo.

``` mutprops.R``` estimates the proportion of each base at A-to-I edited bases using SSII and TGIRT enzymes. 


## Mutation rate
``` mut_rate.R``` estimates the mutational rate of A/G/N called bases to T for de novo A-to-I editing detection


## Detecting admixture
```admixture.R``` parses mutational rates at A-to-I edited sites in pool s with different admixtures of inosine vs adenosine at the A site. 

## ICE-PIK
ICE-PIK phases k-mers to detect multiple A to I edits in the same oligo. 

## Contact
For specific questions about the code or the work, please contact [Layla Siraj](layla.siraj@gmail.com) or [Aaron Lin](alin@broadinstitute.org).
