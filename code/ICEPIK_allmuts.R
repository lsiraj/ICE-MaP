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
library(circlize)

CE0 <- vroom::vroom('~/ICEMaP/tables/2024-10-01.phasing.CE0.ABCDE.txt')
CE1 <- vroom::vroom('~/ICEMaP/tables/2024-10-01.phasing.CE1.ABCDE.txt')


# duos <- CE1 %>% 
#   dplyr::mutate(AB = ifelse(muts == 'All mutated' | (grepl('A', muts) & grepl('B',muts)), T, F),
#                 AE = ifelse(muts == 'All mutated' | (grepl('A', muts) & grepl('E',muts)), T, F),
#                 AC = ifelse(muts == 'All mutated' | (grepl('A', muts) & grepl('C',muts)), T, F),
#                 AD = ifelse(muts == 'All mutated' | (grepl('A', muts) & grepl('D',muts)), T, F),
#                 BE = ifelse(muts == 'All mutated' | (grepl('B', muts) & grepl('E',muts)), T, F),
#                 BC = ifelse(muts == 'All mutated' | (grepl('B', muts) & grepl('C',muts)), T, F),
#                 BD = ifelse(muts == 'All mutated' | (grepl('B', muts) & grepl('D',muts)), T, F),
#                 EC = ifelse(muts == 'All mutated' | (grepl('E', muts) & grepl('C',muts)), T, F),
#                 ED = ifelse(muts == 'All mutated' | (grepl('E', muts) & grepl('D',muts)), T, F),
#                 CD = ifelse(muts == 'All mutated' | (grepl('C', muts) & grepl('D',muts)), T, F)) %>%
#   dplyr::select(-sites, -prop, -muts,-n) %>%
#   colSums()/sum(CE1$n)
# 
# duos <- data.frame(duos) %>% 
#   dplyr::rename(prop = duos) %>% 
#   tibble::rownames_to_column(., "duos") %>%
#   dplyr::mutate(name1 = gsub('.{1}$', '', duos)) %>% 
#   dplyr::mutate(name2 = substr(duos, 2, 2)) %>%
#   dplyr::select(-duos)%>%
#   dplyr::select(name1, name2,prop)
# colors <- jdb_palette("corona",5) %>% as.character 
# names <- c('A','B','C','D','E')
# colors.df <- data.frame(colors,names)
# 
# chordDiagram(duos, 
#              preAllocateTracks = list(track.height = 0.25),
#              grid.col = setNames(colors.df$color, colors.df$names))
# 
# 
# get_trios <- function(nuc){
#   trios <- 
#     CE1 %>% 
#     dplyr::filter(grepl(nuc,muts) | muts == 'All mutated' ) %>%
#     dplyr::mutate(AB = ifelse(muts == 'All mutated' | (grepl('A', muts) & grepl('B',muts)), T, F),
#                   AE = ifelse(muts == 'All mutated' | (grepl('A', muts) & grepl('E',muts)), T, F),
#                   AC = ifelse(muts == 'All mutated' | (grepl('A', muts) & grepl('C',muts)), T, F),
#                   AD = ifelse(muts == 'All mutated' | (grepl('A', muts) & grepl('D',muts)), T, F),
#                   BE = ifelse(muts == 'All mutated' | (grepl('B', muts) & grepl('E',muts)), T, F),
#                   BC = ifelse(muts == 'All mutated' | (grepl('B', muts) & grepl('C',muts)), T, F),
#                   BD = ifelse(muts == 'All mutated' | (grepl('B', muts) & grepl('D',muts)), T, F),
#                   EC = ifelse(muts == 'All mutated' | (grepl('E', muts) & grepl('C',muts)), T, F),
#                   ED = ifelse(muts == 'All mutated' | (grepl('E', muts) & grepl('D',muts)), T, F),
#                   CD = ifelse(muts == 'All mutated' | (grepl('C', muts) & grepl('D',muts)), T, F)) %>%
#     dplyr::select(-sites, -prop, -muts,-n) %>%
#     colSums()/sum(CE1$n)
#   
#   trios <- data.frame(trios) %>% 
#     dplyr::rename(prop = trios) %>% 
#     tibble::rownames_to_column(., "trios") %>%
#     dplyr::mutate(name1 = gsub('.{1}$', '', trios)) %>% 
#     dplyr::mutate(name2 = substr(trios, 2, 2)) %>%
#     dplyr::select(-trios)%>%
#     dplyr::select(name1, name2,prop) %>%
#     dplyr::filter(name1 != nuc, name2 != nuc)
#   
#   
#   return(trios)
# }
# 
# plot_trios <- function(trios){
#   colors <- jdb_palette("corona",5) %>% as.character 
#   names <- c('A','B','C','D','E')
#   colors.df <- data.frame(colors,names)
#   
#   chordDiagram(trios, 
#                preAllocateTracks = list(track.height = 0.25),
#                grid.col = setNames(colors.df$color, colors.df$names))
#   title(nuc)
# }


props <- CE1 %>% 
  group_by(muts) %>% 
  dplyr::summarise(totprop = sum(prop)) %>% 
  arrange(totprop) %>% 
  dplyr::mutate(A = ifelse(muts == 'All mutated' | grepl('A', muts), 1, 0),
                B = ifelse(muts == 'All mutated' | grepl('B', muts), 1, 0),
                E = ifelse(muts == 'All mutated' | grepl('E', muts), 1, 0),
                C = ifelse(muts == 'All mutated' | grepl('C', muts), 1, 0),
                D = ifelse(muts == 'All mutated' | grepl('D', muts), 1, 0))
col= colorRampPalette(c("white", "white", "maroon"))(256)

jpeg(file="~/ICEMaP/figures/allsites.heatmap.jpg",quality=100, height=1800, pointsize = 14, width=3000, res=600)
heatmap(props %>% dplyr::select(c('A','B','C','D','E')) %>% as.matrix(),
            scale = 'none',
            Colv=NA, Rowv=NA,
            col = col)
dev.off()



jpeg(file="~/ICEMaP/figures/allsites.barplot.jpg",quality=100, height=1800, pointsize = 14, width=3000, res=600)
g2 <- ggplot(props %>% mutate(muts = fct_reorder(muts, totprop)), aes(x = muts, y = totprop, label = paste0(round(totprop*100,1),'%')))+
  geom_col(aes(fill='a'))+
  coord_flip()+
  theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line())+
  scale_fill_manual(values=c('maroon'))+
  geom_text(hjust = -0.1)+
  ylim(0,0.3)

g1 <- props %>% pivot_longer(cols = c('A','B','C','D','E'), names_to = 'site') %>%
  mutate(muts = fct_reorder(muts, totprop)) %>%
  mutate(site = factor(site, levels = c('A','B','E','C','D'))) %>%
  ggplot(.,aes(x=site,y=muts,fill=value)) +
  geom_tile(color = 'black')+
  scale_fill_gradient(low='white',high='maroon')+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank())

plt_combined <- g1 + g2 + plot_layout(ncol = 2, nrow = 1, heights = c(6))
plt_combined
cowplot::save_plot(paste0('~/ICEMaP/figures/', Sys.Date(), '.allsites.pdf'),
                   plt_combined,
                   base_height = 6,
                   base_width = 12,
                   device = cairo_pdf
)

#barplot(mtcars[hm$rowInd,"mpg"],horiz=T,names.arg=row.names(mtcars)[hm$rowInd],las=2,cex.names=0.7,col="purple",2)

#dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none'
