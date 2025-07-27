setwd('~/ICEMaP/processing/piledriver_output/')
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(vroom)
library(patchwork)
library(BuenColors)
library(sys)

#functions
getcov <- function(df){
  print(max(df$depth))
  df <- df %>% 
    dplyr::mutate("cov" = 0.5 + 10000*df$depth / max(df$depth))
  print(max(df$cov))
  df <-df %>% 
    dplyr::mutate("norm_cov" = df$cov / max(df$cov))
  print(max(df$norm_cov))
  return(df)
}

plot_cov <- function(f1, f2, cond1, cond2, outname, outdir_p, outdir_t){
  # file readin
  df1 <- vroom::vroom(f1)
  df2 <- vroom::vroom(f2)
  
  # add normalized coverage
  df1 <- getcov(df1)
  df2 <- getcov(df2)
  
  # plot and save
  to_plot <- rbind(df1 %>% dplyr::select("start","end","ref","depth","norm_cov") %>% dplyr::mutate(source = cond1),
                   df2 %>% dplyr::select("start","end","ref","depth","norm_cov") %>% dplyr::mutate(source = cond2))
  
  to_plot$source <- factor(to_plot$source, levels = c(cond1, cond2), labels = c(cond1, cond2))
  
  
  mpre1 <- mean(df1 %>% dplyr::filter(start <= 44) %>% .$norm_cov)
  mpost1 <- mean(df1 %>% dplyr::filter(start > 44) %>% .$norm_cov)
  mpre2 <- mean(df2 %>% dplyr::filter(start <= 44) %>% .$norm_cov)
  mpost2 <- mean(df2 %>% dplyr::filter(start > 44) %>% .$norm_cov)
  
  print(c(mpre1, mpost1, mpre2, mpost2))
  
  p <- ggplot(data = to_plot, aes(x=start, y = norm_cov, fill=source))+
    geom_area(position = "identity")+
    scale_fill_manual(values=c(jdb_palette("brewer_spectra")[7], jdb_palette("brewer_spectra")[1]))+
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(aspect.ratio=1,
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())+
    xlab('Position')+
    ylab('Normalized coverage') +
    geom_segment(aes(x = 0, xend = 45, y = mpre1, yend = mpre1, color = cond1))+
    geom_segment(aes(x = 46, xend = 217, y = mpost1, yend = mpost1, color = cond1))+
    geom_segment(aes(x = 0, xend = 45, y = mpre2, yend = mpre2, color = cond2))+
    geom_segment(aes(x = 46, xend = 217, y = mpost2, yend = mpost2, color = cond2))+
    scale_color_manual(values=c('#4e437f', '#c05d25'))
  
  p <- p + geom_text(data = data.frame(xcoord =c(22.5,131.5,22.5,131.5),
                                  ycoord = c(mpre1, mpost1, mpre2, mpost2),
                                  labs = c(paste0(as.character(round(mpre1*100,1)),"%"),
                                           paste0(as.character(round(mpost1*100,1)),"%"),
                                           paste0(as.character(round(mpre2*100,1)),"%"),
                                           paste0(as.character(round(mpost2*100,1)),"%")),
                                  cols = c(cond1, cond1, cond2, cond2)), aes(x=xcoord, y=ycoord, label = labs, fill =  cols, vjust = -0.5))

plt_combined <- p + plot_layout(ncol = 1, nrow = 1, heights = c(6))
cowplot::save_plot(paste0(outdir_p, Sys.Date(), '.coverage.', outname, '.pdf'),
                   plt_combined,
                   base_height = 6,
                   base_width = 12,
                   device = cairo_pdf
)

# save off table

rbind(df1 %>% dplyr::select("start","end","ref","depth", "cov","norm_cov") %>% dplyr::mutate(source = cond1),
      df2 %>% dplyr::select("start","end","ref","depth","cov","norm_cov") %>% dplyr::mutate(source = cond2)) %>%
  vroom_write(paste0(outdir_t, Sys.Date(), '.coverage.', outname, '.txt'))
print('saved!')

  return(p)
}

admix_cov <- function(perc, df1, cond){
  mpre <- mean(df1 %>% dplyr::filter(start <= 45) %>% .$norm_cov)
  mpost <- mean(df1 %>% dplyr::filter(start > 45) %>% .$norm_cov)
  return(c(perc, mpost, mpre, cond))
}

# get equation from lm
lm_eqn <- function(df){
  m <- lm(cov.loss ~ perc.I, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


#file and variable names

df = data.frame(f1 = c(#'170609-Layla.A_unt_RP1v2.piledriver.txt', 
                # '170609-Layla.A_unt_RP1v2.piledriver.txt', 
                # '170609-Layla.A_unt_RP1v3.piledriver.txt',
                # '170609-Layla.A_unt_RP1v3.piledriver.txt',
                # '170609-Layla.I_unt_RPIv2.piledriver.txt',
                # '170609-Layla.I_unt_RPIv2.piledriver.txt',
                # '170609-Layla.I_unt_RPIv3.piledriver.txt',
                # '170609-Layla.I_unt_RPIv3.piledriver.txt',
                '170622-Layla.0percent-untreated.piledriver.txt',
                '170622-Layla.100percent-untreated.piledriver.txt'
                ),
                f2 = c(#'170609-Layla.A_CE1_RPIv2.piledriver.txt',
                       # '170609-Layla.A_CE2_RPIv2.piledriver.txt',
                       # '170609-Layla.A_CE1_RPIv3.piledriver.txt',
                       # '170609-Layla.A_CE2_RPIv3.piledriver.txt',
                       # '170609-Layla.I_CE1_RPIv2.piledriver.txt',
                       # '170609-Layla.I_CE2_RPIv2.piledriver.txt',
                       # '170609-Layla.I_CE1_RPIv3.piledriver.txt',
                       # '170609-Layla.I_CE2_RPIv3.piledriver.txt',
                       '170622-Layla.0percent-CE2.piledriver.txt',
                       '170622-Layla.100percent-CE2.piledriver.txt'),
                cond1 = c(#'A untreated RPIv2',
                          # 'A untreated RPIv2',
                          # 'A untreated RPIv3',
                          # 'A untreated RPIv3',
                          # 'I untreated RPIv2',
                          # 'I untreated RPIv2',
                          # 'I untreated RPIv3',
                          # 'I untreated RPIv3',
                          'A untreated 70622',
                          'I untreated 170622'),
                cond2 = c(#'A CE1 RPIv2',
                          # 'A CE2 RPIv2',
                          # 'A CE1 RPIv3',
                          # 'A CE2 RPIv3',
                          # 'I CE1 RPIv2',
                          # 'I CE2 RPIv2',
                          # 'I CE1 RPIv3',
                          # 'I CE2 RPIv3',
                          'A CE2 170622',
                          'I CE2 170622'),
                outname = c(#'adenosine.ce1.rpiv2',
                            # 'adenosine.ce2.rpiv2',
                            # 'adenosine.ce1.rpiv3',
                            # 'adenosine.ce2.rpiv3',
                            # 'inosine.ce1.rpiv2',
                            # 'inosine.ce2.rpiv2',
                            # 'inosine.ce1.rpiv3',
                            # 'inosine.ce2.rpiv3',
                            '170622.0percent.ce2',
                            '170622.100percent.ce2'))

for (row in 1:nrow(df)){
  print(row)
  f1 <- df$f1[row]
  f2 <- df$f2[row]
  cond1 <- df$cond1[row]
  cond2 <- df$cond2[row]
  outname <- df$outname[row]
  outdir_p <- "~/ICEMaP/figures/"
  outdir_t <- "~/ICEMaP/tables/"

  p <- plot_cov(f1,f2, cond1, cond2, outname, outdir_p, outdir_t)
  p
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
outdir_p <- '~/ICEMaP/figures/'
outdir_t <- '~/ICEMaP/tables/'
outname <- '170622.admixture.covloss'

#file read-in
df.z1 <- vroom::vroom(z1) %>% getcov(.)
df.z2 <- vroom::vroom(z2) %>% getcov(.)
df.h1 <- vroom::vroom(h1) %>% getcov(.)
df.h2 <- vroom::vroom(h2) %>% getcov(.)
df.t1 <- vroom::vroom(t1) %>% getcov(.)
df.t2 <- vroom::vroom(t2) %>% getcov(.)
df.f1 <- vroom::vroom(f1) %>% getcov(.)
df.f2 <- vroom::vroom(f2) %>% getcov(.)

admix_out <- data.frame(admix_cov(0, df.z1, '-ACN'),
                        admix_cov(0, df.z2, '+ACN'),
                        admix_cov(50, df.f1, '-ACN'),
                        admix_cov(50, df.f2, '+ACN'),
                        admix_cov(10,df.t1, '-ACN'),
                        admix_cov(10,df.t2, '+ACN'),
                        admix_cov(100,df.h1, '-ACN'),
                        admix_cov(100,df.h2, '+ACN'),
           row.names = c('perc.I', 'mpost','mpre', 'cond'))
admix_out <- data.frame(t(admix_out)) %>% 
  tibble::rownames_to_column('source') %>% 
  dplyr::select(-source) 
admix_out$perc.I <- as.numeric(admix_out$perc.I)
admix_out$mpost <- as.numeric(admix_out$mpost)*100
admix_out$mpre <- as.numeric(admix_out$mpre)*100
admix_out <- admix_out %>% dplyr::mutate(cov.loss = mpost - mpre)

admix_out %>% vroom_write(paste0(outdir_t, Sys.Date(),'.', outname, '.txt'))

p2<- ggplot(admix_out, aes(x=perc.I, y=cov.loss, group=cond, color = cond)) +
  geom_point()+
  geom_smooth(method = "lm", alpha = .15, aes(fill = cond))+
  scale_color_manual(values=c('#c05d25','#4e437f'))+
  theme(aspect.ratio=1.5,
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  ylim(0,100)

 p2 <- p2 + geom_text(x = 50, y = 35, label = lm_eqn(admix_out %>% dplyr::filter(cond == '-ACN')), parse = TRUE, color = '#c05d25' )+
   geom_text(x = 50, y = 67.5, label = lm_eqn(admix_out %>% dplyr::filter(cond == '+ACN')), parse = TRUE, color ='#4e437f')+
   theme(legend.position = "none")+
   xlab('% I at Input HTR2C A site') +
   ylab('% Coverage Loss')
 
 p2
plt_combined <- p2 + plot_layout(ncol = 1, nrow = 1, heights = c(6))
cowplot::save_plot(paste0(outdir_p, Sys.Date(), '.',outname, '.pdf'),
                   plt_combined,
                   base_height = 6,
                   base_width = 12,
                   device = cairo_pdf)
#anova
lm1 <- lm(cov.loss ~ perc.I + cond, data = admix_out)
lm2 <- lm(cov.loss ~ perc.I + cond + perc.I*cond, data = admix_out)

lm1
lm2

anova(lm1, lm2)

one.way <- aov(cov.loss ~ perc.I + cond + perc.I*cond, data = admix_out)
summary(one.way)
 