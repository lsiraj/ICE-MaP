setwd('~/ICEMaP/processing/piledriver_output/')
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(vroom)
library(scales)
library(Rsamtools)
library(Biostrings)

# run after ICEPIK
# get DFs from ICEPIK_dfs.R

### functions part 2
# read in dfs
read_in_df <- function(file){
  df <- vroom::vroom(file)
  df <- df %>% mutate_at(c('startpos', 'asite_i', 'bsite_i','csite_i','esite_i','dsite_i'), as.numeric) %>% dplyr::mutate(rownum = 1:dim(df)[1])
  df <- df %>% dplyr::filter(startpos < 21)
}

seq_split <- function(row, df){
  pos <- as.numeric(df[row,'startpos'])
  # print(substr(seq,21-pos+1,33-pos+1))
  # print(substr(fasta.seq, 21,33))
  seq <- df[row,'seq']
  asite <- substr(seq, 21-pos+1, 21-pos+1)
  bsite <- substr(seq, 23-pos+1, 23-pos+1)
  esite <- substr(seq, 27-pos+1, 27-pos+1)
  csite <- substr(seq, 28-pos+1, 28-pos+1)
  dsite <- substr(seq, 33-pos+1, 33-pos+1)
  return(c(asite, bsite, esite, csite, dsite))
}

seq_tog <- function(row, df){
  pos <- as.numeric(df[row,'startpos'])
  # print(substr(seq,21-pos+1,33-pos+1))
  # print(substr(fasta.seq, 21,33))
  seq <- df[row,'seq']
  asite <- substr(seq, 21-pos+1, 21-pos+1)
  bsite <- substr(seq, 23-pos+1, 23-pos+1)
  esite <- substr(seq, 27-pos+1, 27-pos+1)
  csite <- substr(seq, 28-pos+1, 28-pos+1)
  dsite <- substr(seq, 33-pos+1, 33-pos+1)
  return(paste0(asite, bsite, esite, csite, dsite))
}

get_seqdf <- function(df){
  b <- lapply(1:dim(df)[1], seq_split, df =df)
  seqdf <- as.data.frame(do.call(rbind,b))
  colnames(seqdf) <- c('asite', 'bsite', 'esite','csite','dsite')
  return(seqdf)
}

get_seqlist <- function(df){
  b <- lapply(1:dim(df)[1], seq_tog, df =df)
  seqlist <- as.data.frame(do.call(rbind,b))
  colnames(seqlist) <- c('sites')
  return(seqlist)
}

muts_per_site <- function(seqdf, df){
  # mutations per site
  a <- as.data.frame(seqdf %>% count(asite) %>% rbind(c('insertion', sum(df$asite_i)))) 
  colnames(a) = c('nucleotide','asite')
  b <- as.data.frame(seqdf %>% count(bsite) %>% rbind(c('insertion', sum(df$bsite_i)))) 
  colnames(b) = c('nucleotide','bsite')
  e <- as.data.frame(seqdf %>% count(esite) %>% rbind(c('insertion', sum(df$esite_i)))) 
  colnames(e) = c('nucleotide','esite')
  c <- as.data.frame(seqdf %>% count(csite) %>% rbind(c('insertion', sum(df$csite_i)))) 
  colnames(c) = c('nucleotide','csite')
  d <- as.data.frame(seqdf %>% count(dsite) %>% rbind(c('insertion', sum(df$dsite_i)))) 
  colnames(d) = c('nucleotide','dsite')
  overall <- a %>%
    merge(b, by = 'nucleotide', all = T) %>%
    merge(e, by = 'nucleotide', all = T) %>%
    merge(c, by = 'nucleotide', all = T) %>% 
    merge(d, by = 'nucleotide', all = T) 
  overall %>% 
    column_to_rownames("nucleotide") %>% 
    t %>% 
    as.data.frame %>% 
    rownames_to_column("site")
}
  

icepik <- function(file1, file2, df, tbase, title){ 
  df1 <- read_in_df(paste0('~/ICEMaP/processing/MultipleInosines/', file1))
  df2 <- read_in_df(paste0('~/ICEMaP/processing/MultipleInosines/', file2))
  
  seqdf1 <- get_seqdf(df1)
  seqdf2 <- get_seqdf(df2)
  
  seqlist1 <- get_seqlist(df1)
  seqlist2 <- get_seqlist(df2)
  
  # overall1 <- muts_per_site(seqdf1, df1)
  # overall2 <- muts_per_site(seqdf2, df2)
  
  df <- vroom::vroom(paste0('~/ICEMaP/code/',df,'.txt'))
  
  seqlist1 %>% 
    count(sites) %>%
    dplyr::mutate(prop = n/sum(n)) %>% 
    left_join(.,df %>% dplyr::select(seq,muts), by = c("sites" = "seq")) %>%
    vroom::vroom_write(paste0('~/ICEMaP/tables/', Sys.Date(), '.phasing.CE0.',title, '.txt'))
  
  s1 <- seqlist1 %>% 
    count(sites) %>%
    dplyr::mutate(prop = n/sum(n)) %>% 
    left_join(.,df %>% dplyr::select(seq,muts), by = c("sites" = "seq")) %>%
    group_by(muts) %>%
    summarise(tot_prop = sum(prop)) %>%
    arrange(-tot_prop)
  
  seqlist2 %>% 
    count(sites) %>%
    dplyr::mutate(prop = n/sum(n)) %>% 
    left_join(.,df %>% dplyr::select(seq,muts), by = c("sites" = "seq")) %>%
    vroom::vroom_write(paste0('~/ICEMaP/tables/', Sys.Date(), '.phasing.CE1.',title, '.txt'))
  
  s2 <- seqlist2 %>% 
    count(sites) %>%
    dplyr::mutate(prop = n/sum(n)) %>% 
    left_join(.,df %>% dplyr::select(seq,muts), by = c("sites" = "seq")) %>%
    group_by(muts) %>%
    summarise(tot_prop = sum(prop)) %>%
    arrange(-tot_prop)
  
  
  x <- c(5.5, 3.5,4.5,5.5,6.5,7.5,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,3.5,4.5,5.5,6.5,7.5,5.5)
  y <- c(11,9,9,9,9,9,7,7,7,7,7,7,7,7,7,7,5,5,5,5,5,5,5,5,5,5,3,3,3,3,3,1)
  muts <- c('No mutations',
            'A site','B site','E site','C site','D site',
            'AB site','AE site', 'AC site','AD site', 'BE site','BC site','BD site','EC site','ED site','EC site',
            'ABE site','ABC site','ABD site','AEC site','AED site','ACD site','BEC site','BED site','BCD site','ECD site',
            'ABEC site','ABED site','ABCD site','AECD site','BECD site',
            'All mutated')
  overall <- c(11,
               9,9,9,9,9,
               7,7,7,7,7,7,7,7,7,7,
               5,5,5,5,5,5,5,5,5,5,
               3,3,3,3,3,
               1)
  
  df1 <- data.frame(x,y, muts, overall) %>% 
    left_join(.,s1, by = 'muts') %>% 
    dplyr::mutate(tot_prop = ifelse(is.na(tot_prop), 0, tot_prop)) %>%
    arrange(-tot_prop)
  
  p1 <- ggplot(df1, aes(x=x,y=y,fill=tot_prop))+
    geom_point(size=10, shape=22, aes(fill=tot_prop))+
    scale_fill_gradient(low="white", high="maroon",
                        limits = c(0, 1), oob = scales::squish)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    ylim(0,12)+
    ggtitle(paste0(tbase, ', -ACN'))
  
  p1
  
  
  df1 <- df1 %>% group_by(overall) %>%
    summarise(tot_prop = sum(tot_prop)) %>%
    dplyr::mutate(x = c(-1,-1,-1,-1,-1,-1))
  
  p1a <- ggplot(df1, aes(x=x,y=overall,fill=tot_prop))+
    geom_point(size=10, shape=22, aes(fill=tot_prop))+
    scale_fill_gradient(low="white", high="maroon",
                        limits = c(0, 1), oob = scales::squish)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(), legend.position="none")+
    ylim(0,12)+
    xlim(-1.5, -0.5)+
    ggtitle(paste0(tbase, ', -ACN'))
  
  p1a
  
  df2 <- data.frame(x,y, muts, overall) %>% 
    left_join(.,s2, by = 'muts') %>% 
    dplyr::mutate(tot_prop = ifelse(is.na(tot_prop), 0, tot_prop)) %>%
    arrange(-tot_prop)
  
  p2 <- ggplot(df2, aes(x=x,y=y,fill=tot_prop))+
    geom_point(size=10, shape=22, aes(fill=tot_prop))+
    scale_fill_gradient(low="white", high="maroon",
                        limits = c(0, 1), oob = scales::squish)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    ylim(0,12)+
    ggtitle(paste0(tbase, ', +ACN'))
  
  p2
  
  df2 <- df2 %>% group_by(overall) %>%
    summarise(tot_prop = sum(tot_prop)) %>%
    dplyr::mutate(x = c(-1,-1,-1,-1,-1,-1))
  
  p2a <- ggplot(df2, aes(x=x,y=overall,fill=tot_prop))+
    geom_point(size=10, shape=22, aes(fill=tot_prop))+
    scale_fill_gradient(low="white", high="maroon",
                        limits = c(0, 1), oob = scales::squish)+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none")+
    ylim(0,12)+
    xlim(-1.5, -0.5)+
    ggtitle(paste0(tbase, ', +ACN'))
  
  p2a
  
  
  plt_combined <- p1a + p1 + p2a + p2 + plot_layout(ncol = 2, nrow = 2, heights = c(12))
  
  plt_combined
  
  # cowplot::save_plot(paste0('~/ICEMaP/figures/', Sys.Date(), '.phasing.',title, '.pdf'),
  #                    plt_combined,
  #                    base_height = 6,
  #                    base_width = 10,
  #                    device = cairo_pdf)
  return(plt_combined)
}



file1 <- c('170812-ICE-MaP_MultipleInosines.5HT2cR.0.CE0.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.A.CE0.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.AB.CE0.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.ABC.CE0.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.ABCDE.CE0.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.B.CE0.txt')
file2 <- c('170812-ICE-MaP_MultipleInosines.5HT2cR.0.CE1.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.A.CE1.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.AB.CE1.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.ABC.CE1.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.ABCDE.CE1.txt',
           '170812-ICE-MaP_MultipleInosines.5HT2cR.B.CE1.txt')
dfs <- c('df0', 'dfA', 'dfAB','dfABC','dfABECD', 'dfB')
tbase <- c('No mutations', 'A site',' A and B sites',' A,B,C sites',' All sites',' B site')
title <- c('0','A','AB','ABC','ABCDE','B')

input <- data.frame(file1, file2, dfs, tbase, title)

for (row in 1:dim(input)[1]){
  p <- icepik(input[row,'file1'],
              input[row,'file2'],
              input[row,'dfs'],
              input[row,'tbase'],
              input[row,'title'])
  p
}


# I think we are missing things with insertions at the E and C sites
