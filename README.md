# ICE-MaP
## About
Here we introduce ICE-MaP, Inosine CyanoEthylation and Mutational Profiling, building upon previous experimental methods for detecting inosines. ICE-MaP can accurately quantify inosines at single-base resolution, down to 10% frequency. Moreover, ICE-MaP can detect up to 5 inosines at a 13-base editing hotspot, and provides phase data for disentangling multiple editing patterns. ICE-MaP is compatible with a wider range of library construction and sequencing methods, including long read technology, making it possible to phase multiple  hotspots even hundreds of bases apart. 
## Data
To come - GEO!
## Processing
### Step 1
Create a list of file names

```ls -l ../MiSeq/*.bam | grep -v mapped | awk '{print $9}' > filelist.txt```
## Contact
For specific questions about the code or the work, please contact [Layla Siraj](layla.siraj@gmail.com) or [Aaron Lin](alin@broadinstitute.org).
