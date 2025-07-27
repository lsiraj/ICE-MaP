setwd('~/ICEMaP/processing/piledriver_output/')
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(vroom)
library(patchwork)
library(BuenColors)
library(sys)
library(Biostrings)
library(stats)
library(binom)

mut_rate <- function(file){
  df <- vroom::vroom(paste0('~/ICEMaP/processing/piledriver_output/',file)) %>%
        dplyr::select(end, ref, depth, num_A, num_C, num_G, num_T, num_D, num_I) %>%
        dplyr::mutate(A = num_A/depth,
                      C = num_C/depth,
                      G = num_G/depth,
                      T = num_T/depth,
                      indel = (num_D+num_I)/depth)
    
      return(df)
}

# read in all experiments
file1 <- '170609-Layla.A_unt_RP1v2.piledriver.txt'
file2 <- '170609-Layla.A_CE2_RPIv2.piledriver.txt'
file3 <- '170609-Layla.I_unt_RPIv3.piledriver.txt'
file4 <- '170609-Layla.I_CE2_RPIv3.piledriver.txt'
file6 <- '170622-Layla.0percent-CE2.piledriver.txt'
file5 <- '170622-Layla.0percent-untreated.piledriver.txt'
file8 <- '170622-Layla.10percent-CE2.piledriver.txt'
file7 <- '170622-Layla.10percent-untreated.piledriver.txt'
file10 <- '170622-Layla.50percent-CE2.piledriver.txt'
file9 <- '170622-Layla.50percent-untreated.piledriver.txt'
file12 <- '170622-Layla.100percent-CE2.piledriver.txt'
file11 <- '170622-Layla.100percent-untreated.piledriver.txt'
file13 <- '170713-Layla.TGIRT-I-oligo-unt.piledriver.txt'
file14 <- '170713-Layla.TGIRT-I-oligo-CE2.piledriver.txt'
file15 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.0.CE0.piledriver.txt'
file16 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.0.CE1.piledriver.txt'
file17 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.A.CE0.piledriver.txt'
file18 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.A.CE1.piledriver.txt'
file19 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.B.CE0.piledriver.txt'
file20 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.B.CE1.piledriver.txt'
file21 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.AB.CE0.piledriver.txt'
file22 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.AB.CE1.piledriver.txt'
file23 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.ABC.CE0.piledriver.txt'
file24 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.ABC.CE1.piledriver.txt'
file25 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.ABCDE.CE0.piledriver.txt'
file26 <- '170812-ICE-MaP_MultipleInosines.5HT2cR.ABCDE.CE1.piledriver.txt'

# get proportions of nucleotides for each position
df1 <- mut_rate(file1)
df2 <- mut_rate(file2)
df3 <- mut_rate(file3)
df4 <- mut_rate(file4)
df5 <- mut_rate(file5)
df6 <- mut_rate(file6)
df7 <- mut_rate(file7)
df8 <- mut_rate(file8)
df9 <- mut_rate(file9)
df10 <- mut_rate(file10)
df11 <- mut_rate(file11)
df12 <- mut_rate(file12)
df13 <- mut_rate(file13)
df14 <- mut_rate(file14)
df15 <- mut_rate(file15)
df16 <- mut_rate(file16)
df17 <- mut_rate(file17)
df18 <- mut_rate(file18)
df19 <- mut_rate(file19)
df20 <- mut_rate(file20)
df21 <- mut_rate(file21)
df22 <- mut_rate(file22)
df23 <- mut_rate(file23)
df24 <- mut_rate(file24)
df25 <- mut_rate(file25)
df26 <- mut_rate(file26)

# bind together in one experimental file
df <- rbind(df1 %>% dplyr::mutate(cond = 'A at A site',
                                  treat = 'untreated'),
            df2 %>% dplyr::mutate(cond = 'A at A site',
                                  treat = 'CE2'),
            df3 %>% dplyr::mutate(cond = 'I at A site',
                                  treat = 'untreated'),
            df4 %>% dplyr::mutate(cond = 'I at A site',
                                  treat = 'CE2'),
            df5 %>% dplyr::mutate(cond = '0% I at A site',
                                  treat = 'untreated'),
            df6 %>% dplyr::mutate(cond = '0% I at A site',
                                  treat = 'CE2'),
            df7 %>% dplyr::mutate(cond = '10% I at A site',
                                  treat = 'untreated'),
            df8 %>% dplyr::mutate(cond = '10% I at A site',
                                  treat = 'CE2'),
            df9 %>% dplyr::mutate(cond = '50% I at A site',
                                  treat = 'untreated'),
            df10 %>% dplyr::mutate(cond = '50% I at A site',
                                   treat = 'CE2'),
            df11 %>% dplyr::mutate(cond = '100% I at A site',
                                   treat = 'untreated'),
            df12 %>% dplyr::mutate(cond = '100% I at A site',
                                   treat = 'CE2'),
            df13 %>% dplyr::mutate(cond = 'I at A site, TGIRT',
                                   treat = 'untreated'),
            df14 %>% dplyr::mutate(cond = 'I at A site, TGIRT',
                                   treat = 'CE2'),
            df15 %>% dplyr::mutate(cond = 'Multiple Inosines: A at A site',
                                   treat = 'untreated'),
            df16 %>% dplyr::mutate(cond = 'Multiple Inosines: A at A site',
                                   treat = 'CE2'),
            df17 %>% dplyr::mutate(cond = 'Multiple Inosines: I at A site',
                                   treat = 'untreated'),
            df18 %>% dplyr::mutate(cond = 'Multiple Inosines: I at A site',
                                   treat = 'CE2'),
            df19 %>% dplyr::mutate(cond = 'Multiple Inosines: I at B site',
                                   treat = 'untreated'),
            df20 %>% dplyr::mutate(cond = 'Multiple Inosines: I at B site',
                                   treat = 'CE2'),
            df21 %>% dplyr::mutate(cond = 'Multiple Inosines: I at A,B sites',
                                   treat = 'untreated'),
            df22 %>% dplyr::mutate(cond = 'Multiple Inosines: I at A,B sites',
                                   treat = 'CE2'),
            df23 %>% dplyr::mutate(cond = 'Multiple Inosines: I at A,B,C sites',
                                   treat = 'untreated'),
            df24 %>% dplyr::mutate(cond = 'Multiple Inosines: I at A,B,C sites',
                                   treat = 'CE2'),
            df25 %>% dplyr::mutate(cond = 'Multiple Inosines: I at all sites',
                                   treat = 'untreated'),
            df26 %>% dplyr::mutate(cond = 'Multiple Inosines: I at all sites',
                                   treat = 'CE2'))

# get one line per mutation type
df <- df %>% pivot_longer(cols = c('A','C','T','G','indel')) %>% 
  dplyr::rename(nuc = name, prop = value) %>%
  dplyr::mutate(lower = case_when(nuc == 'G' ~ binom.confint(num_G, depth, methods="wilson")$lower,
                                  nuc == 'A' ~ binom.confint(num_A, depth, methods="wilson")$lower,
                                  nuc == 'C' ~ binom.confint(num_C, depth, methods="wilson")$lower,
                                  nuc == 'T' ~ binom.confint(num_T, depth, methods="wilson")$lower,
                                  nuc == 'indel' ~ binom.confint(num_I+num_D, depth, methods="wilson")$lower),
                upper = case_when(nuc == 'G' ~ binom.confint(num_G, depth, methods="wilson")$upper,
                                  nuc == 'A' ~ binom.confint(num_A, depth, methods="wilson")$upper,
                                  nuc == 'C' ~ binom.confint(num_C, depth, methods="wilson")$upper,
                                  nuc == 'T' ~ binom.confint(num_T, depth, methods="wilson")$upper,
                                  nuc == 'indel' ~ binom.confint(num_I+num_D, depth, methods="wilson")$upper))

# # get differences between pre and post
# df_diff <- df %>%
#   dplyr::select(end, ref, cond, nuc, treat, prop) %>%
#   pivot_wider(names_from = treat, values_from = prop) %>%
#   dplyr::mutate(diff = CE2-untreated)
# condition = '10% I at A site'
# ggplot(df_diff %>% dplyr::filter(cond == condition, ref %in% c('A','N')), aes(x=end, y = diff, group = nuc, color = nuc)) +
#   geom_text(aes(label = nuc))+
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line())+
#   scale_color_manual(values=c(jdb_palette("china_novice")[4],
#                               jdb_palette("china_basics")[5],
#                               jdb_palette("brewer_celsius")[1],
#                               jdb_palette("corona")[8],
#                               jdb_palette("brewer_red")[7]))+
#   ylab('Difference between proportion of base between treated and untreated')+
#   xlab('Position along 5HT2CR')+
#   geom_vline(xintercept = 21, alpha = 0.4)+
#   geom_vline(xintercept = 23, alpha = 0.4)+
#   geom_vline(xintercept = 27, alpha = 0.4)+
#   geom_vline(xintercept = 28, alpha = 0.4)+
#   geom_vline(xintercept = 33, alpha = 0.4)+
#   ggtitle(condition)
  

# 
# t <- '10% I at A site'
# p <- ggplot(df %>% dplyr::filter(cond == t,treat == 'CE2',( end >= 18 & end <= 38)), aes(x=end, y=prop, group=nuc, color = nuc, shape = cond)) +
#   geom_text(aes(label = nuc))+
#   geom_errorbar(aes(ymin = lower, ymax=upper))+
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line())+
#   scale_color_manual(values=c(jdb_palette("china_novice")[4],
#                               jdb_palette("china_basics")[5],
#                               jdb_palette("brewer_celsius")[1],
#                               jdb_palette("corona")[8],
#                               jdb_palette("brewer_red")[7]))+
#   ylim(0,1)+
#   geom_vline(xintercept=20.5)+
#   geom_vline(xintercept =21.5)+
#   ggtitle(t)
# 
# p

# hists <- df %>% 
#   dplyr::filter(ref %in% c('A','N', 'G'), 
#                 nuc %in% c('T')) %>% 
#   dplyr::mutate(color = case_when(end == 21 ~ 'A-site',
#                                   end == 23 ~ 'B-site',
#                                   end == 27 ~ 'E-site',
#                                   end == 28 ~ 'C-site',
#                                   end == 33 ~ 'D-site'))
# ggplot(hists, aes(x=log10(prop), color = treat, fill = treat))+
#          geom_histogram(alpha = 0.5, aes(y=..density..))+
#   facet_wrap(~cond)+
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line())


# get a plot of the differences
to_plot <- df %>% 
  dplyr::filter(ref %in% c('A','N', 'G'), 
                nuc %in% c("T")) %>% 
  dplyr::mutate(color = case_when(end == 21 ~ 'A-site',
                                  end == 23 ~ 'B-site',
                                  end == 27 ~ 'E-site',
                                  end == 28 ~ 'C-site',
                                  end == 33 ~ 'D-site')) %>%
  # for error propagation
  dplyr::mutate(se = sqrt((prop*(1-prop))/depth),
                sd = sqrt(prop * (1-prop))) %>%
  dplyr::select(end, ref, cond, nuc, color, prop, treat, se) %>%
  pivot_wider(names_from = treat, values_from = c(prop, se)) %>%
  dplyr::mutate(diff = prop_CE2 - prop_untreated,
                se = sqrt(se_CE2**2 + se_untreated**2)) %>%
  dplyr::mutate(lower = diff - 1.96*se,
                upper = diff + 1.96*se)
  
# plot
p <- ggplot(to_plot, aes(x=end, y=diff)) +
  geom_text(aes(label = nuc, color = nuc))+
  geom_errorbar(aes(ymin = lower, ymax=upper))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  scale_color_manual(values=c(jdb_palette("brewer_red")[7]))+
  #scale_color_manual(values = c(jdb_palette("china_basics")[5]))+
  #scale_color_manual(values = c(jdb_palette("china_basics")[5], jdb_palette("brewer_red")[7]))+
  ylim(-0.1,.1)+
  facet_wrap(~cond)

p
# cowplot::save_plot(paste0('~/ICEMaP/figures/', Sys.Date(), '.T.errorrates.zoom.pdf'),
#                    p,
#                    base_height = 7.5,
#                    base_width = 10,
#                    device = cairo_pdf)

p

# # version 1: 
# to_plot <- df %>% 
#   dplyr::filter(ref %in% c('A','N', 'G'), 
#                 nuc %in% c('T','C')) %>% 
#   dplyr::mutate(color = case_when(end == 21 ~ 'A-site',
#                                   end == 23 ~ 'B-site',
#                                   end == 27 ~ 'E-site',
#                                   end == 28 ~ 'C-site',
#                                   end == 33 ~ 'D-site')) %>%
#   group_by(nuc, cond) %>%
#   dplyr::mutate(m = mean(prop), sd = sd(prop)) %>%
#   ungroup() %>%
#   dplyr::mutate(z = (prop - m) / sd) %>%
#   arrange(-abs(z)) %>%
#   dplyr::select(-c(m,sd, prop, lower, upper)) %>%
#   pivot_wider(names_from = nuc, values_from = z) %>%
#   dplyr::mutate(comp = C + `T`)  %>%
#   dplyr::select(end, ref, cond, color, comp, `T`, C, treat) %>%
#   pivot_wider(names_from = treat, values_from = c(comp, `T`, C)) %>% 
#   dplyr::mutate(C = C_CE2 - C_untreated) %>%
#   dplyr::mutate(`T` = T_CE2 - C_untreated) %>%
#   dplyr::select(end, ref, cond, color, `T`,C) %>%
#   pivot_longer(names_to = 'nuc', cols = c(`T`, C)) %>%
#   dplyr::rename('diff' = 'value') %>%
#   # dplyr::filter((color == 'A-site' & cond %in% c('I at A site',
#   #                                              '10% I at A site',
#   #                                              '50% I at A site',
#   #                                              '100% I at A site',
#   #                                              'I at A site, TGIRT',
#   #                                              'Multiple Inosines: I at A site',
#   #                                              'Multiple Inosines: I at A,B sites',
#   #                                              'Multiple Inosines: I at A,B,C sites',
#   #                                              'Multiple Inosines: I at all sites') )| 
#   #                 (color == 'B-site' & cond %in% c('Multiple Inosines: I at B site',
#   #                                                 'Multiple Inosines: I at A,B sites',
#   #                                                 'Multiple Inosines: I at A,B,C sites',
#   #                                                 'Multiple Inosines: I at all sites')) | 
#   #                 (color == 'C-site' & cond %in% c('Multiple Inosines: I at A,B,C sites',
#   #                                                  'Multiple Inosines: I at all sites'))|
#   #                 (color %in% c('D-site','E-site') & cond == 'Multiple Inosines: I at all sites')) %>%
#   arrange(cond, -diff)
#   
# c = '50% I at A site'
# p <- ggplot(to_plot %>% dplyr::filter(nuc=='T'), aes(x=end, y=diff, fill = nuc))+
#   geom_col()+
#   scale_fill_manual(values=c(jdb_palette('corona')[2]))+
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line())+
#   xlab('Position')+
#   ylab('Z score difference between treated and untreated')+geom_hline(yintercept = 5) + geom_hline(yintercept = 10)
# p + facet_wrap(~cond)
# 
# #jdb_palette('corona')[19], 

# # version 2:
# df %>% 
#   dplyr::filter(ref %in% c('A','N', 'G'), 
#                 nuc %in% c('T','C')) %>% 
#   dplyr::mutate(color = case_when(end == 21 ~ 'A-site',
#                                   end == 23 ~ 'B-site',
#                                   end == 27 ~ 'E-site',
#                                   end == 28 ~ 'C-site',
#                                   end == 22 ~ 'D-site')) %>%
#   group_by(nuc, cond) %>%
#   dplyr::mutate(m = mean(prop), sd = sd(prop)) %>%
#   ungroup() %>%
#   dplyr::mutate(z = (prop - m) / sd) %>%
#   arrange(-abs(z)) %>%
#   dplyr::select(-c(m,sd, prop, lower, upper)) %>%
#   arrange(cond, -z) %>%
#   dplyr::filter(cond == '0% I at A site', treat =='CE2') %>%
#   ggplot(., aes(x=end, y=z, fill = nuc))+
#   geom_col(position='stack')+
#   scale_fill_manual(values=jdb_palette('corona'))+
#   geom_vline(xintercept = 21, aes(linetype='dash'))+
#   geom_vline(xintercept = 23, aes(linetype='dash'))+
#   geom_vline(xintercept = 27, aes(linetype='dash'))+
#   geom_vline(xintercept = 28, aes(linetype='dash'))+
#   geom_vline(xintercept = 33, aes(linetype='dash'))
  


 
nucleotide = 'T'
p <- ggplot(df %>% 
              dplyr::filter(ref %in% c('A','N', 'G'), 
                            nuc == nucleotide) %>% 
              dplyr::mutate(color = case_when(end == 21 ~ 'A-site',
                                              end == 23 ~ 'B-site',
                                              end == 27 ~ 'E-site',
                                              end == 28 ~ 'C-site',
                                              end == 33 ~ 'D-site')), 
            aes(x=cond, y=prop))+
  geom_violin()+
  geom_point(data = df %>%
               dplyr::filter(ref %in% c('A','N','G'),
                             nuc == nucleotide) %>% 
             dplyr::mutate(color = case_when((end == 21) ~ 'A-site',
                                             (end == 23) ~ 'B-site',
                                             (end == 27) ~ 'E-site',
                                             (end == 28 )~ 'C-site',
                                             (end == 33) ~ 'D-site')),
             aes(x=cond, y=prop, colour = color, shape = treat), size = 2.5)+
  geom_point(data = df %>%
               dplyr::filter(ref %in% c('A','N','G'),
                             nuc == nucleotide,
                             end %in% c(21,23,27,28,33)) %>% 
               dplyr::mutate(color = case_when((end == 21) ~ 'A-site',
                                               (end == 23) ~ 'B-site',
                                               (end == 27) ~ 'E-site',
                                               (end == 28 )~ 'C-site',
                                               (end == 33) ~ 'D-site')),
             aes(x=cond, y=prop, colour = color, shape = treat), size = 2.5)+
  scale_color_manual(values = jdb_palette('corona'))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  ggtitle(paste0(nucleotide, ' Mutation rates at A/G/N reference nucleotides'))


p
# 
# cowplot::save_plot(paste0('~/ICEMaP/figures/', Sys.Date(), '.C.errorrates.violin.pdf'),
#                    p,
#                    base_height = 5,
#                    base_width = 7,
#                    device = cairo_pdf)

lin <- to_plot %>% 
  dplyr::filter(cond %in% c('10% I at A site', '50% I at A site','100% I at A site'), color == 'A-site') %>% 
  dplyr::mutate(x = case_when(cond == '10% I at A site' ~ 10, 
                              cond == '50% I at A site'~ 50, 
                              cond == '100% I at A site' ~ 100))

lm_eqn <- function(lin){
  m <- lm(diff ~ x, lin);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

ggplot(lin, aes(x=x, y = diff)) + 
  geom_point() + 
  geom_line() + 
  geom_smooth(method='lm', formula= y~x)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  geom_text(x = 50, y = 0.6, label = lm_eqn(lin), parse = TRUE)

to_plot %>% group_by(cond) %>% dplyr::summarise(mean(diff), mean(upper))

to_plot %>% 
  dplyr::filter((color == 'A-site' & cond %in% c('I at A site',
                                                 '10% I at A site',
                                                 '50% I at A site',
                                                 '100% I at A site',
                                                 'I at A site, TGIRT',
                                                 'Multiple Inosines: I at A site',
                                                 'Multiple Inosines: I at A,B sites',
                                                 'Multiple Inosines: I at A,B,C sites',
                                                 'Multiple Inosines: I at all sites') )|
                    (color == 'B-site' & cond %in% c('Multiple Inosines: I at B site',
                                                    'Multiple Inosines: I at A,B sites',
                                                    'Multiple Inosines: I at A,B,C sites',
                                                    'Multiple Inosines: I at all sites')) |
                    (color == 'C-site' & cond %in% c('Multiple Inosines: I at A,B,C sites',
                                                     'Multiple Inosines: I at all sites'))|
                    (color %in% c('D-site','E-site') & cond == 'Multiple Inosines: I at all sites'))



df %>% 
  dplyr::filter(ref %in% c('A','N', 'G'), 
                nuc %in% c("T"))  %>% 
  dplyr::select(end, ref,num_T, depth, cond, treat) %>%
  group_by(cond, treat) %>% dplyr::summarise(Ttot = sum(num_T), tot = sum(depth) ) %>%
  dplyr::mutate(prop = Ttot/tot,
                se = sqrt((prop*(1-prop))/tot)) %>%
  dplyr::select(-c(Ttot, tot)) %>%
  pivot_wider(names_from = treat, values_from = c('prop','se'))


### get baseline T mutation rates at A sites or G sites, across each condition
p<- df %>% 
  dplyr::mutate(treat = ifelse(treat=='CE2', 'Treated','Untreated')) %>%
  dplyr::mutate(name = paste(ref,treat)) %>% 
  dplyr::filter(ref %in% c('G','A','N'), nuc == 'T',!grepl('TGIRT', cond)) %>%
  dplyr::filter(!(cond == 'I at A site' & end == 21 ),
                !(cond %in% c('10% I at A site', '50% I at A site', '100% I at A site') & end == 21),
                !(cond == 'Multiple Inosines: I at A site' & end == 21),
                !(cond == 'Multiple Inosines: I at B site' & end == 23),
                !(cond == 'Multiple Inosines: I at A,B sites' & end %in% c(21,23)),
                !(cond == 'Multiple Inosines: I at A,B,C sites' & end %in% c(21,23,28)),
                !(cond == 'Multiple Inosines: I at all sites' & end %in% c(21,23,27,28,33))) %>% 
  dplyr::mutate(se = sqrt((prop * (1-prop))/depth),
                lc = prop -1.96*se, uc = prop +1.96*se) %>% 
dplyr::select(end, ref,num_T, depth, cond, treat, name, prop, se) %>%
  ggplot(., aes(x = name, y = prop*100)) + 
  geom_violin(trim=, draw_quantiles = c(0.25, 0.5, 0.75))  +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
  ggtitle("Percentage of T mutations at A,G, or N positions without inosines")+
  xlab('') +
  ylab( '% T mutation')+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())

cowplot::save_plot(paste0('~/ICEMaP/figures/', Sys.Date(), '.T.violinplots.treated.untreated.pdf'),
                   p,
                   base_height = 4,
                   base_width = 6,
                   device = cairo_pdf)


to_plot <- df %>% 
  dplyr::mutate(treat = ifelse(treat=='CE2', 'Treated','Untreated')) %>%
  dplyr::mutate(name = paste(ref,treat)) %>% 
  dplyr::filter(ref %in% c('G','A','N'), nuc == 'T',!grepl('TGIRT', cond)) %>%
dplyr::mutate(se = sqrt((prop*(1-prop))/depth),
              sd = sqrt(prop * (1-prop))) %>%
  dplyr::select(end, ref, cond, nuc, prop, treat, se) %>%
  pivot_wider(names_from = treat, values_from = c(prop, se)) %>%
  dplyr::mutate(diff = prop_Treated - prop_Untreated,
                se = sqrt(se_Treated**2 + se_Untreated**2)) %>%
  dplyr::mutate(lower = diff - 1.96*se,
                upper = diff + 1.96*se) %>%
  dplyr::mutate(color = case_when(cond == 'I at A site' & end == 21  ~ 'N: I',
                                  cond %in% c('10% I at A site', '50% I at A site', '100% I at A site') & end == 21 ~ 'N: I',
                                  cond == 'Multiple Inosines: I at A site' & end == 21 ~ 'N: I',
                                  cond == 'Multiple Inosines: I at B site' & end == 23 ~ 'N: I',
                                  cond == 'Multiple Inosines: I at A,B sites' & end %in% c(21,23) ~ 'N: I',
                                  cond == 'Multiple Inosines: I at A,B,C sites' & end %in% c(21,23,28) ~ 'N: I',
                                  cond == 'Multiple Inosines: I at all sites' & end %in% c(21,23,27,28,33) ~ 'N: I',
                                  ref == 'N'~ 'N: no I',
                                  ref == 'G' ~ 'G',
                                  ref == 'A' ~ 'A')) %>%
  dplyr::filter(grepl('%', cond))

to_plot$cond <- factor(to_plot$cond,levels=c("0% I at A site","10% I at A site","50% I at A site", '100% I at A site'))

p <-  ggplot(to_plot, aes(x = ref, y = diff*100)) + 
  geom_violin(trim=T, draw_quantiles = c(0.25, 0.5, 0.75))  +
  geom_jitter(height = 0, width = 0.1, alpha = 0.5)+
  ggtitle("Difference in % of T mutations at A,G, or N positions without inosines")+
  xlab('') +
  ylab( 'Difference in % T mutation')+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  facet_wrap(~cond, nrow = 1)

cowplot::save_plot(paste0('~/ICEMaP/figures/', Sys.Date(), '.T.violinplots.diff.pdf'),
                   p,
                   base_height = 4,
                   base_width = 6,
                   device = cairo_pdf)



 df %>% 
  dplyr::mutate(treat = ifelse(treat=='CE2', 'Treated','Untreated')) %>%
  dplyr::mutate(name = paste(ref,treat)) %>% 
  dplyr::filter(ref %in% c('G','A','N'), nuc == 'T',!grepl('TGIRT', cond)) %>%
  dplyr::mutate(se = sqrt((prop*(1-prop))/depth),
                sd = sqrt(prop * (1-prop))) %>%
  dplyr::select(end, ref, cond, nuc, prop, treat, se) %>%
  pivot_wider(names_from = treat, values_from = c(prop, se)) %>%
  dplyr::mutate(diff = prop_Treated - prop_Untreated,
                se = sqrt(se_Treated**2 + se_Untreated**2)) %>%
  dplyr::mutate(lower = diff - 1.96*se,
                upper = diff + 1.96*se) %>%
  dplyr::mutate(color = case_when(cond == 'I at A site' & end == 21  ~ 'N: I',
                                  cond %in% c('10% I at A site', '50% I at A site', '100% I at A site') & end == 21 ~ 'N: I',
                                  cond == 'Multiple Inosines: I at A site' & end == 21 ~ 'N: I',
                                  cond == 'Multiple Inosines: I at B site' & end == 23 ~ 'N: I',
                                  cond == 'Multiple Inosines: I at A,B sites' & end %in% c(21,23) ~ 'N: I',
                                  cond == 'Multiple Inosines: I at A,B,C sites' & end %in% c(21,23,28) ~ 'N: I',
                                  cond == 'Multiple Inosines: I at all sites' & end %in% c(21,23,27,28,33) ~ 'N: I',
                                  ref == 'N'~ 'N: no I',
                                  ref == 'G' ~ 'G',
                                  ref == 'A' ~ 'A')) %>%
  group_by(color, cond) %>% 
   # summarise(se = sd(diff*100)/n()) %>%arrange(color, -se) %>%
   # group_by(color) %>% summarise(se=median(se))
  summarise(med_diff = median(diff*100),
            min_diff = min(diff*100),
            max_diff = max(diff*100)) %>% arrange(color, -med_diff) %>% print(n=500)

 df %>% 
   dplyr::mutate(treat = ifelse(treat=='CE2', 'Treated','Untreated')) %>%
   dplyr::mutate(name = paste(ref,treat)) %>% 
   dplyr::filter(ref %in% c('G','A','N'), nuc == 'T',!grepl('TGIRT', cond)) %>%
   dplyr::mutate(se = sqrt((prop*(1-prop))/depth),
                 sd = sqrt(prop * (1-prop))) %>%
   dplyr::select(end, ref, cond, nuc, prop, treat, se) %>%
   dplyr::mutate(color = case_when(cond == 'I at A site' & end == 21  ~ 'N: I',
                                   cond %in% c('10% I at A site', '50% I at A site', '100% I at A site') & end == 21 ~ 'N: I',
                                   cond == 'Multiple Inosines: I at A site' & end == 21 ~ 'N: I',
                                   cond == 'Multiple Inosines: I at B site' & end == 23 ~ 'N: I',
                                   cond == 'Multiple Inosines: I at A,B sites' & end %in% c(21,23) ~ 'N: I',
                                   cond == 'Multiple Inosines: I at A,B,C sites' & end %in% c(21,23,28) ~ 'N: I',
                                   cond == 'Multiple Inosines: I at all sites' & end %in% c(21,23,27,28,33) ~ 'N: I',
                                   ref == 'N'~ 'N: no I',
                                   ref == 'G' ~ 'G',
                                   ref == 'A' ~ 'A')) %>%
   group_by(color, cond) %>% 
   summarise(med_se = median(se*100),
             min_se = min(se*100),
             max_se = max(se*100)) %>% arrange(color, -med_se) %>% print(n=500) 
 
 

