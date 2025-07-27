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

# filenames
z1 <- '170622-Layla.0percent-untreated.piledriver.txt'
z2 <- '170622-Layla.0percent-CE2.piledriver.txt'
h1 <- '170622-Layla.100percent-untreated.piledriver.txt'
h2 <- '170622-Layla.100percent-CE2.piledriver.txt'
t1 <- '170622-Layla.10percent-untreated.piledriver.txt'
t2 <- '170622-Layla.10percent-CE2.piledriver.txt'
f1 <- '170622-Layla.50percent-untreated.piledriver.txt'
f2 <- '170622-Layla.50percent-CE2.piledriver.txt'

#file read-in
df.z1 <- vroom::vroom(z1)
df.z2 <- vroom::vroom(z2)
df.h1 <- vroom::vroom(h1)
df.h2 <- vroom::vroom(h2)
df.t1 <- vroom::vroom(t1)
df.t2 <- vroom::vroom(t2)
df.f1 <- vroom::vroom(f1)
df.f2 <- vroom::vroom(f2)

df <- rbind(props(df.z1) %>% dplyr::mutate(source = '0% untreated'),
            props(df.z2) %>% dplyr::mutate(source = '0% CE2'),
            props(df.h1) %>% dplyr::mutate(source = '100% untreated'),
            props(df.h2) %>% dplyr::mutate(source = '100% CE2'),
            props(df.t1) %>% dplyr::mutate(source = '10% untreated'),
            props(df.t2) %>% dplyr::mutate(source = '10% CE2'),
            props(df.f1) %>% dplyr::mutate(source = '50% untreated'),
            props(df.f2) %>% dplyr::mutate(source = '50% CE2'))

df %>% vroom_write(paste0('~/ICEMaP/tables/', Sys.Date(), '.mutprops.', '170622.admixtures.ce2', '.txt'))
print('saved!')

# get background for expectation
df.ctrl <- df %>% 
  dplyr::mutate(case = case_when(startsWith(source, '0%') ~ 'all A',
                                 startsWith(source, '100%') ~ 'all I',
                                 startsWith(source, '10%') ~ '10% I',
                                 startsWith(source, '50%') ~ '50% I')) %>% 
  dplyr::mutate(treat = ifelse(endsWith(source, "untreated"), 'untreated','CE2')) %>%  
  dplyr::select(-c(depth, num_A, num_C, num_G, num_T, num_D, num_I, source)) %>%
  pivot_wider(names_from = treat, values_from = c(A,C,G,T,indel)) %>% t(.)
df.ctrl <- data.frame(df.ctrl)[-1,] %>% dplyr::rename(all.A = X1, all.I = X2, ten.I = X3, fifty.I = X4)
df.ctrl$all.A <- as.numeric(df.ctrl$all.A)
df.ctrl$all.I <- as.numeric(df.ctrl$all.I)
df.ctrl$ten.I <- as.numeric(df.ctrl$ten.I)
df.ctrl$fifty.I <- as.numeric(df.ctrl$fifty.I)
df.ctrl  <- df.ctrl %>% 
  dplyr::mutate(ten.I.expect = all.A*.9 + all.I*.1,
                fifty.I.expect = all.A*0.5 + all.I * 0.5) %>% t(.)
  
df.ctrl <- data.frame(df.ctrl) %>% 
  tibble::rownames_to_column('source') %>% 
  pivot_longer(colnames(data.frame(df.ctrl))) %>%
  dplyr::mutate(type = ifelse(endsWith(source, "expect"), "expected", "observed")) %>% 
  dplyr::mutate(source = case_when(source == "all.A" ~ '0% I',
                                   source == 'all.I' ~ '100% I',
                                   source == 'fifty.I' ~ '50% I',
                                   source == 'ten.I' ~ '10% I',
                                   source == 'ten.I.expect' ~ '10% I',
                                   source == 'fifty.I.expect' ~ '50% I')) %>%
  dplyr::mutate(cond = ifelse(endsWith(name, 'untreated'), '-ACN', '+ACN')) %>% 
  dplyr::mutate(name = gsub('_untreated',"", name)) %>%
  dplyr::mutate(name = gsub('_CE2', '', name)) %>%
  dplyr::rename(Nucleotide = name)

p1 <- ggplot(df.ctrl %>% dplyr::filter(source == '0% I'), aes(x=cond, y=value, group = Nucleotide, color = Nucleotide)) +
  geom_line()+
  geom_text(aes(label=Nucleotide))+
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
        axis.line = element_line(),
        legend.position = "none",
        plot.title= element_text(size=12))+
  geom_text(data=df.ctrl %>% dplyr::filter(source == '0% I'), aes(x=cond, y = value, label=paste0(as.character(round(value*100, 1)), '%'), vjust = -0.5))+
  xlab('') + ylab('') + 
  ggtitle('0% I')+
  ylim(0,1)

p1

p2 <- ggplot(df.ctrl %>% dplyr::filter(source == '100% I'), aes(x=cond, y=value, group = Nucleotide, color = Nucleotide)) +
  geom_line()+
  geom_text(aes(label=Nucleotide))+
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
        axis.line = element_line(),
        legend.position = "none",
        plot.title = element_text(size=12))+ 
  geom_text(data=df.ctrl %>% dplyr::filter(source == '100% I'), aes(x=cond, y = value, label=paste0(as.character(round(value*100, 1)), '%'), vjust = -0.5))+
  xlab('') + ylab('') + 
  ggtitle('100% I')+
  ylim(0,1)

p2


p3 <- ggplot(df.ctrl %>% 
               dplyr::filter(source == '50% I') %>%
               dplyr::mutate(interact = paste(sep= '.', Nucleotide, type)), aes(x=cond, y=value, group=interact, color = Nucleotide)) +
  geom_line(aes(linetype=type))+
  geom_text(aes(label=Nucleotide))+
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
        axis.line = element_line(),
        legend.position = "none",
        plot.title = element_text(size=12))+ 
  geom_text(data=df.ctrl %>% dplyr::filter(source == '50% I', type == 'observed') %>%
              dplyr::mutate(interact = paste(sep= '.', Nucleotide, type)), aes(x=cond, y = value, label=paste0(as.character(round(value*100, 1)), '%'), hjust =-0.4))+
  scale_linetype_manual(breaks=c("observed","expected"), values = c(1,2))+
  xlab('') + ylab('') + 
  ggtitle('50% I')+
  ylim(0,1)

p3

p4 <- ggplot(df.ctrl %>% 
               dplyr::filter(source == '10% I') %>%
               dplyr::mutate(interact = paste(sep= '.', Nucleotide, type)), aes(x=cond, y=value, group=interact, color = Nucleotide)) +
  geom_line(aes(linetype=type))+
  geom_text(aes(label=Nucleotide))+
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
        axis.line = element_line(),
        legend.position = "none",
        plot.title = element_text(size=12))+ 
  geom_text(data=df.ctrl %>% dplyr::filter(source == '10% I', type == 'observed') %>%
              dplyr::mutate(interact = paste(sep= '.', Nucleotide, type)), aes(x=cond, y = value, label=paste0(as.character(round(value*100, 1)), '%'), hjust =-0.4))+
  scale_linetype_manual(breaks=c("observed","expected"), values = c(1,2))+
  xlab('') + ylab('') + 
  ggtitle('10% I')+
  ylim(0,1)

p4

p5 <- ggplot(df.ctrl %>% 
               dplyr::filter(source == '10% I') %>%
               dplyr::mutate(interact = paste(sep= '.', Nucleotide, type)) %>%
               dplyr::filter(Nucleotide != 'A'), aes(x=cond, y=value, group=interact, color = Nucleotide)) +
  geom_line(aes(linetype=type))+
  geom_text(aes(label=Nucleotide))+
  scale_color_manual(values=c(jdb_palette("china_basics")[5],
                              jdb_palette("brewer_celsius")[1],
                              jdb_palette("corona")[8],
                              jdb_palette("brewer_red")[7]))+
  theme(aspect.ratio=1.5,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none",
        plot.title= element_text(size=12))+ 
  geom_text(data=df.ctrl %>% dplyr::filter(source == '10% I', type == 'observed') %>%
              dplyr::mutate(interact = paste(sep= '.', Nucleotide, type)) %>%
              dplyr::filter(Nucleotide != 'A'), aes(x=cond, y = value, label=paste0(as.character(round(value*100, 1)), '%'), hjust =-0.4))+
  scale_linetype_manual(breaks=c("observed","expected"), values = c(1,2))+
  ylim(0,0.1)+
  xlab('') + ylab('') +
  ggtitle('10% I - zoom')

p5

p6 <- ggplot(df.ctrl %>% dplyr::filter(source == '0% I', Nucleotide != 'A'), aes(x=cond, y=value, group = Nucleotide, color = Nucleotide)) +
  geom_line()+
  geom_text(aes(label=Nucleotide))+
  scale_color_manual(values=c(jdb_palette("china_basics")[5],
                              jdb_palette("brewer_celsius")[1],
                              jdb_palette("corona")[8],
                              jdb_palette("brewer_red")[7]))+
  theme(aspect.ratio=1.5,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.position = "none",
        plot.title= element_text(size=12))+ 
  geom_text(data=df.ctrl %>% dplyr::filter(source == '0% I', Nucleotide != 'A'), aes(x=cond, y = value, label=paste0(as.character(round(value*100, 1)), '%'), vjust = -0.5))+
  ylim(0,0.1)+
  xlab('') + ylab('')+
  ggtitle('0% I - zoom')
p6

plt_combined <- p2 + p3+ p4+ p1+ p5+ p6 + plot_layout(ncol = 2, nrow = 3)
plt_combined

cowplot::save_plot(paste0('~/ICEMaP/figures/', Sys.Date(), '.admixture.', 'exp.170622', '.pdf'),
                   plt_combined,
                   base_height = 14,
                   base_width = 7,
                   device = cairo_pdf
)





