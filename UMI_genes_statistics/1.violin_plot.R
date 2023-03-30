setwd('D:/postdoc/peilab/project1/liver/check')

library('ggplot2')


pdf("nGene per Pixel.pdf")
aa <- read.table('gene_count_wanted',header=T,sep = '\t')
ggplot(aa,aes(x=sample,y=genes,fill=sample))+geom_violin()+scale_y_log10()+ylab("nGene per Pixel")+
  theme_bw()+ theme(panel.background = element_rect(fill = 'white', color = 'black'))

ggplot(aa,aes(x=density,y=genes,fill=density))+geom_violin()+scale_y_log10()+ ylab("nGene per Pixel")+
  theme_bw()+ theme(panel.background = element_rect(fill = 'white', color = 'black'))
dev.off()

pdf("nUMI per Pixel.pdf")
aa <- read.table('UMI_count_wanted',header=T,sep = '\t')
ggplot(aa,aes(x=sample,y=UMIs,fill=sample))+geom_violin()+scale_y_log10()+ylab("nUMI per Pixel")+
  theme_bw()+ theme(panel.background = element_rect(fill = 'white', color = 'black'))

ggplot(aa,aes(x=density,y=UMIs,fill=density))+geom_violin()+scale_y_log10()+ ylab("nUMI per Pixel")+
  theme_bw()+ theme(panel.background = element_rect(fill = 'white', color = 'black'))
dev.off()

pdf("nGene per ¦Ìm2.pdf")
aa <- read.table('30pixel2_gene_count_wanted',header=T,sep = '\t')
ggplot(aa,aes(x=sample,y=genes,fill=sample))+geom_violin()+scale_y_log10()+ylab("nGene per ¦Ìm2")+
  theme_bw()+ theme(panel.background = element_rect(fill = 'white', color = 'black'))
dev.off()

pdf("nUMI per ¦Ìm2.pdf") 
aa <- read.table('30pixel2_UMI_count_wanted',header=T,sep = '\t')
ggplot(aa,aes(x=sample,y=UMIs,fill=sample))+geom_violin()+scale_y_log10()+ylab("nUMI per ¦Ìm2")+
  theme_bw()+ theme(panel.background = element_rect(fill = 'white', color = 'black'))
dev.off()

pdf("pixel density per ¦Ìm2.pdf")
aa <- read.table('30pixel2_pixel_density',header=T,sep = '\t')
ggplot(aa,aes(x=sample,y=counts,fill=sample))+geom_violin()+scale_y_log10()+ylab(expression(paste("nUMI per ",¦Ì, m^2)))+
  theme_bw()+ theme(panel.background = element_rect(fill = 'white', color = 'black'))
dev.off()

