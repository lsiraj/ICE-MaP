setwd('~/ICEMaP/processing/piledriver_output/')
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(vroom)
library(scales)
library(Rsamtools)
library(Biostrings)

# essential sequences
fasta <- readDNAStringSet("~/ICEMaP/data/references/5HT2cR-ligated-5N.fa")
fasta.seq <- paste(fasta)
# get the positions of the kmers we're interested in
unlist(gregexpr('N', fasta.seq))

# function
#### get cigar sublists
cparse <- function(bam, i){
 pos <- bam[[1]]$pos[i]
 cigar <- bam[[1]]$cigar[i]
 cigar <- gsub('M', 'M:', cigar)
 cigar <- gsub('D', 'D:', cigar)
 cigar <- gsub('I', 'I:', cigar)
 cigar <- gsub('S', 'S:', cigar)
 ciglist <- strsplit(cigar,':')
 return(ciglist)
}


#### clean up the read according to the cigar
seqclean <- function(bam, i){
  
  # get sequence
  seq <-paste(bam[[1]]$seq[i])

  
  
  # get cigar list split up by letter
  ciglist <- cparse(bam,i)
  
  # keep track of insertions
  # I don't care if there's an insertion other than at these sites
  asite = 0
  bsite = 0
  esite = 0
  csite = 0
  dsite = 0
  

  print(seq)
  print(ciglist)
  # split sequence up by cigar strings and modify
  seqlist <- c()
  pos = 0
  for (j in 1:length(ciglist[[1]])){
    # get cigar info
    cigar <- ciglist[[1]][j]
    c <- substr(cigar, nchar(cigar),nchar(cigar))
    print(paste0('cigar is ', cigar, ' and c is ', c))

    # if a match, append sequence to seq-list and add sub-sequence length to position marker
    if (c == 'M'){
      seqlist <- append(seqlist, substr(seq, 1,as.numeric(gsub(c,"",cigar))))
      pos <- pos + as.numeric(gsub(c,"",cigar))
      print('it is a match!')
      #get rid of the portion of sequence that cigar string dealt with
      seq <- substr(seq,as.numeric(gsub(c,"",cigar))+1, nchar(seq))

    }
    
    # if insertion, keep track of what position in the final cleaned read the insertion is at, and how long the insertion i
    else if (c == 'I'){
      if ((pos+1) == 21){asite = 1}
      else if ((pos+1) == 23){bsite = 1}
      else if ((pos+1) == 27){esite = 1}
      else if ((pos + 1) == 28){csite = 1}
      else if ((pos + 1) == 33){dsite = 1}
      print('it is an insertion!')
      # get rid of the portion of sequence that cigar string dealt with
      seq <- substr(seq,as.numeric(gsub(c,"",cigar))+1, nchar(seq))
    }
    
    # if deletion, add dashes of the appropriate length to seq-list, and increment position marker appropriately
    else if (c=='D'){
      seqlist <- append(seqlist, paste(rep('-',as.numeric(gsub(c,"",cigar))), collapse = ''))
      pos <- pos + as.numeric(gsub(c,"",cigar))
      print('it is a deletion!')
    }
    
    else if (c == 'S'){
    # if s, just ignore
      print('soft clipping!')
      # get rid of the portion of sequence that cigar string dealt with
      seq <- substr(seq,as.numeric(gsub(c,"",cigar))+1, nchar(seq))
    }
    print(paste0('fastass is ', fasta.seq))
    print(paste0('seqlist is ', paste(seqlist, collapse = '')))
    print(paste0('seq remaining is ', seq))
  }
  
  # collapse the out sequence
  seq_out <- paste(seqlist, collapse = '')
  
  # update progress bar
 #setTxtProgressBar(pb, i)
  
  # return the start position of the read, the cleaned read, and insertions at ABECD
  return(c(as.numeric(bam[[1]]$pos[i]),
           seq_out,
           as.numeric(asite), 
           as.numeric(bsite), 
           as.numeric(esite), 
           as.numeric(csite), 
           as.numeric(dsite)))
}

# check number of unique characters in cigar
uniqchars <- function(x) unique(strsplit(x, "")[[1]])
# cigar = ''
# for(i in 1:128983){
#   cigar <- paste0(cigar, bam[[1]]$cigar[i])
#   if (i%%1000 == 0){
#     print(i)
#   }
# }
# uniqchars(cigar)

files <- c('170812-ICE-MaP_MultipleInosines.5HT2cR.0.CE0.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.0.CE1.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.A.CE0.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.A.CE1.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.AB.CE0.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.AB.CE1.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.ABC.CE0.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.ABC.CE1.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.ABCDE.CE0.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.ABCDE.CE1.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.B.CE0.mapped.sorted.bam',
              '170812-ICE-MaP_MultipleInosines.5HT2cR.B.CE1.mapped.sorted.bam'
)


## to generate DFs: 
# for (file in files){
#   print(file)
#   outname <- str_split(str_split(file, '170812-ICE-MaP_MultipleInosines.5HT2cR.')[[1]][2], '.mapped.sorted.bam')[[1]][1]
#   print(outname)
#   # read in the .bam file
#   bam <- scanBam(paste0('~/ICEMaP/processing/',file))
#   # run over all reads
#   tot <- length(bam[[1]]$seq)
#   # set progression bar
#   pb <- txtProgressBar(min = 0, max = tot, style = 3)
#   # get cleaned reads and make df
#   a <- lapply(1:tot, seqclean, bam=bam, pb = pb)
#   df <- as.data.frame(do.call(rbind,a))
#   colnames(df) <- c('startpos', 'seq', 'asite_i','bsite_i','esite_i','csite_i','dsite_i')
#   df %>% vroom::vroom_write(paste0('~/ICEMaP/processing/MultipleInosines/170812-ICE-MaP_MultipleInosines.5HT2cR.',outname,'.txt'))
# }

