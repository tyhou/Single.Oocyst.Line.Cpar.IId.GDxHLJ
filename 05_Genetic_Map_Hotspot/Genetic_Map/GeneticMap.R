
all.map <- read.csv("map-15-Lines.csv",header = 1)
chr_info <- read.table("GeneticMap_Info_Chr.txt",header = 1)
library("tidyr")
library("ggplot2")
library("ggchicklet")

all.map.2 <- all.map
all.map.2$genetic_distance <- all.map$genetic_distance *4

# wider to longer
all.df <- gather(all.map.2, key = "group",
                  value = "pos",-`id`,-`Chr`,
                  na.rm = FALSE, convert = FALSE, factor_key = FALSE)

chr.df <- gather(chr_info, key = "group",
                 value = "pos",-`Chr`,
                 na.rm = FALSE, convert = FALSE, factor_key = FALSE)
chr.df$pos <- ifelse(chr.df$group %in% "genetic_distance", chr.df$pos *4, chr.df$pos)

g <- ggplot()+
  geom_line(data = all.df, aes(x=group,y=pos,group=id), color="gray" ,position = position_dodge(0))+
  geom_col(data = chr.df, aes(x=group,y=pos,fill=group),position = 'dodge',width = 0.2)+
  facet_wrap( ~ Chr,strip.position = "bottom",ncol = 8)+
  labs(x = '', y = 'Genetic distance (cM)')  + 
  scale_y_continuous(limits = c(0,1300),
                     breaks = c(0,500,1000),
                     labels = c(0,125,200),
                     expand = c(0.1,0.1),
                     sec.axis = sec_axis(~.,
                                         name = 'Physical distance (kb)',
                                         breaks = c(0,500,1000)))+
  scale_fill_manual(values=c("#B71C1C", "#3949AB"))+
  theme(panel.border = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank (),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) + 
  theme(legend.position = 'none')
g
ggsave("physical_and_genetic_maps.pdf",g,width = 6,height = 4,dpi = 300 )
ggsave("physical_and_genetic_maps.png",g,width = 6,height = 4,dpi = 600 )




