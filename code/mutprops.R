setwd('~/ICEMaP/processing/piledriver_output/')
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(vroom)
library(scales)

# functions
props <- function(df){
  df<- df %>%
    # position of interest is start==20
    dplyr::filter(start == 20) %>% 
    dplyr::select(depth, num_A, num_C, num_G, num_T, num_D, num_I) %>%
    dplyr::mutate(A = num_A/depth,
                  C = num_C/depth,
                  G = num_G/depth,
                  T = num_T/depth,
                  indel = (num_D+num_I)/depth)
  return(df)
}
plot_props <- function(f1, f2, cond1, cond2, title, outname, outdir_p, outdir_t){
  
  #file read-in
  df1 <- vroom::vroom(f1)
  df2 <- vroom::vroom(f2)
  
  df <- rbind(props(df1) %>% dplyr::mutate(source = cond1),
              props(df2) %>% dplyr::mutate(source = cond2))
  
  df %>% vroom_write(paste0(outdir_t, Sys.Date(), '.mutprops.', outname, '.txt'))
  print('saved!')
  
  df <- df %>% 
    pivot_longer(cols=c(A, C, G, T, indel), names_to = 'Nucleotide') %>%
    dplyr::select(source, Nucleotide, value) %>%
    dplyr::rename(read.prop = value)
  
  p <- ggplot(df, aes(x=source, y=read.prop, group = Nucleotide, color = Nucleotide)) +
    geom_line()+
    geom_point()+
    scale_color_manual(values=c(jdb_palette("china_novice")[4],
                                jdb_palette("china_basics")[5],
                                jdb_palette("brewer_celsius")[1],
                                jdb_palette("corona")[8],
                                jdb_palette("brewer_red")[7]))+
    theme(aspect.ratio=1.5,
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())+ 
    geom_text(data=df, aes(x=source, y = read.prop, label=paste0(Nucleotide, ": ", as.character(round(read.prop*100, 1)), '%'), vjust = -0.5))+
    ggtitle(title)
  
  plt_combined <- p + plot_layout(ncol = 1, nrow = 1, heights = c(6))
  cowplot::save_plot(paste0(outdir_p, Sys.Date(), '.mutprops.', outname, '.pdf'),
                     plt_combined,
                     base_height = 6,
                     base_width = 12,
                     device = cairo_pdf
  )
  return(p)

}

# set variables
df = data.frame(f1 = c('170609-Layla.A_unt_RP1v2.piledriver.txt', 
                       '170609-Layla.A_unt_RP1v2.piledriver.txt', 
                       '170609-Layla.A_unt_RP1v3.piledriver.txt',
                       '170609-Layla.A_unt_RP1v3.piledriver.txt',
                       '170609-Layla.I_unt_RPIv2.piledriver.txt',
                       '170609-Layla.I_unt_RPIv2.piledriver.txt',
                       '170609-Layla.I_unt_RPIv3.piledriver.txt',
                       '170609-Layla.I_unt_RPIv3.piledriver.txt',
                       '170713-Layla.TGIRT-I-oligo-unt.piledriver.txt',
                       '170713-Layla.TGIRT-I-oligo-unt.piledriver.txt'),
                f2 = c('170609-Layla.A_CE1_RPIv2.piledriver.txt',
                       '170609-Layla.A_CE2_RPIv2.piledriver.txt',
                       '170609-Layla.A_CE1_RPIv3.piledriver.txt',
                       '170609-Layla.A_CE2_RPIv3.piledriver.txt',
                       '170609-Layla.I_CE1_RPIv2.piledriver.txt',
                       '170609-Layla.I_CE2_RPIv2.piledriver.txt',
                       '170609-Layla.I_CE1_RPIv3.piledriver.txt',
                       '170609-Layla.I_CE2_RPIv3.piledriver.txt',
                       '170713-Layla.TGIRT-I-oligo-CE1.piledriver.txt',
                       '170713-Layla.TGIRT-I-oligo-CE2.piledriver.txt'),
                cond1 = c('- ACN',
                          '- ACN',
                          '- ACN',
                          '- ACN',
                          '- ACN',
                          '- ACN',
                          '- ACN',
                          '- ACN',
                          '- ACN',
                          '- ACN'),
                cond2 = c('+ ACN',
                          '+ ACN',
                          '+ ACN',
                          '+ ACN',
                          '+ ACN',
                          '+ ACN',
                          '+ ACN',
                          '+ ACN',
                          '+ ACN',
                          '+ ACN'),
                title=c('A at A site',
                      'A at A site',
                      'A at A site',
                      'A at A site',
                      'I at A site',
                      'I at A site',
                      'I at A site',
                      'I at A site',
                      'TGIRT: I at A site',
                      'TGIRT: I at A site'),
                outname = c('adenosine.ce1.rpiv2',
                            'adenosine.ce2.rpiv2',
                            'adenosine.ce1.rpiv3',
                            'adenosine.ce2.rpiv3',
                            'inosine.ce1.rpiv2',
                            'inosine.ce2.rpiv2',
                            'inosine.ce1.rpiv3',
                            'inosine.ce2.rpiv3',
                            'TGIRT.inosine.ce1',
                            'TGIRT.inosine.ce2'))


for (row in 1:nrow(df)){
  print(row)
  f1 <- df$f1[row]
  f2 <- df$f2[row]
  cond1 <- df$cond1[row]
  cond2 <- df$cond2[row]
  title <- df$title[row]
  outname <- df$outname[row]
  outdir_p <- "~/ICEMaP/figures/"
  outdir_t <- "~/ICEMaP/tables/"
  
  p <- plot_props(f1,f2, cond1, cond2, title, outname, outdir_p, outdir_t)
  p
}

# consider admixtures




170622-Layla.0percent-CE1.piledriver.txt
170622-Layla.0percent-CE2.piledriver.txt
170622-Layla.0percent-untreated.piledriver.txt
170622-Layla.100percent-CE1.piledriver.txt
170622-Layla.100percent-CE2.piledriver.txt
170622-Layla.100percent-untreated.piledriver.txt
170622-Layla.10percent-CE1.piledriver.txt
170622-Layla.10percent-CE2.piledriver.txt
170622-Layla.10percent-untreated.piledriver.txt
170622-Layla.50percent-CE1.piledriver.txt
170622-Layla.50percent-CE2.piledriver.txt
170622-Layla.50percent-untreated.piledriver.txt
170629-Layla.v2_000-CE2.piledriver.txt
170629-Layla.v2_000-unt.piledriver.txt
170629-Layla.v2_005-CE2.piledriver.txt
170629-Layla.v2_005-unt.piledriver.txt
170629-Layla.v2_010-CE2.piledriver.txt
170629-Layla.v2_010-unt.piledriver.txt
170629-Layla.v2_100-CE2.piledriver.txt
170629-Layla.v2_100-unt.piledriver.txt
170629-Layla.v4_000-CE2.piledriver.txt
170629-Layla.v4_000-unt.piledriver.txt
170629-Layla.v4_005-CE2.piledriver.txt
170629-Layla.v4_005-unt.piledriver.txt
170629-Layla.v4_010-unt.piledriver.txt
170629-Layla.v4_100-CE2.piledriver.txt
170629-Layla.v4_100-unt.piledriver.txt


