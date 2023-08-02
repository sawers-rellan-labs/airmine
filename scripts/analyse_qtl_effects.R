library(tidyr)
library(dplyr)
library(qtl)
library(ggplot2)


mineral_names <-  c("Al", "B",  "Ca", "Fe",
                    "K",  "Mg", "Mn", "Mo",
                    "Ni", "P",  "S",  "Zn",
                    "Na", "Cu")

traits <- c("weight",mineral_names)

reliable <- c( "Ca", "Fe",  "K",
               "Mg", "Mn", "P", 
               "S",  "Zn")

load(file="../run/AIR_cross.Rdata")
cross <- BC2S3
qtls  <- read.csv("../run/mim_qtl_explained_var.csv")

reliable_peaks <- qtls %>%
  # filter(lod>thresh_adj) %>%
  filter(trait %in% reliable) %>%
  select(trait,name:lod,var,everything()) %>%
  arrange(-lod)


p<- reliable_peaks  %>%
  ggplot(aes(x =var, fill =trait, group=trait,label =trait)) +
  geom_dotplot(stackgroups = TRUE,method = "histodot",binwidth = 0.2)+
  ggtitle("Distribution of QTL effects")+
  xlab("Variance Explained %")+
  ylab("Count")+
  theme(legend.position = "none") 

to_plot <- ggplot_build(p)[[1]][[1]]


label_at <- function(n) function(x) ifelse(x %% n == 0, x, "")


x_max <- ceiling(max(to_plot$x))
x_min <- floor(max(to_plot$x))

y_max <- ceiling(max(to_plot$stackpos))
y_min <- floor(max(to_plot$stackpos))

quartz(width=3.2, height = 10.15)

to_plot %>%
  ggplot(aes(x = x, y = stackpos, label = label,group=label, fill =label)) +
  geom_vline(xintercept = 2.1, linetype =2)+
  geom_point(pch =21, size =7) +
  #  scale_fill_brewer(palette="Set1") +
  ggtitle("QTL effects")+
  xlab("Variance Explained %")+
  ylab("QTL Count")+
  scale_x_continuous(breaks = seq(1, x_max, 0.2),  
                     labels = label_at(1),
                     limits = c(1,x_max)) +
  scale_y_continuous(breaks = 0:y_max,  
                     labels = label_at(1),
                     limits = c(0,y_max)) +
  geom_text(size=3)+
  ggpubr::theme_classic2(base_size = 15)+
  coord_flip() +
  theme(legend.position = "none") 


####### Effect plots #####

new_peaks <- qtls %>%
  filter(lod>thresh_adj) %>%
  group_by(trait) %>%
  arrange(-lod) %>%
  dplyr::slice(1) %>%
  arrange(-lod)


npeaks <- nrow(new_peaks)

effect <- lapply(1:npeaks, FUN = function(peak_idx) {
  this_peak <-  new_peaks[peak_idx,]
  data.frame( qtl_name = this_peak$name,
              marker = this_peak$marker, 
              genotype = factor(pull.geno(cross,chr=this_peak$chr)[,this_peak$marker]),
              trait = this_peak$trait,
              mass = cross$pheno[,this_peak$trait])
}
) %>% dplyr::bind_rows()




levels(effect$genotype) <- list(CML530=1,Het=2,LR=3)

plots <- list()


for(t in  unique(effect$trait)){
  t_data <- effect[effect$trait==t,]
  plots[[paste0("effect_",t)]] <- t_data %>%
    filter(!is.na(genotype)) %>%
    mutate(qtl_name= paste(trait,qtl_name)) %>%
    ggplot2::ggplot(aes(x = genotype, y = mass, col =genotype)) +
    ggplot2::ylab(t) +
    ylab(paste(t,"ppm [mg/kg]")) +
    xlab("Genotype") +
    ggbeeswarm::geom_quasirandom() + 
    #ggplot2::geom_text( x= 1, y= 1, label = paste("R2 =", lms$r.squared)) +
    stat_summary( fun.data = mean_cl_normal,col ="black") +  
    #ggtitle(paste("R2", lms$r.squared) +
    ggplot2::facet_wrap(~ qtl_name) +
    ggpubr::theme_classic2(base_size = 20) +
    ggplot2::theme(legend.position = "none",
                   strip.background = element_rect(colour = "white", fill = "white"))
}


# summary(lm(data= effect[effect$marker == pqm[2],], P ~ genotype))


# Phosphorus boxplots ####
plots$p_donor <- cross$pheno %>%
  ggplot(aes(x= donor, y =P, col =range)) +
  ggbeeswarm::geom_quasirandom() +
  stat_summary( fun.data = mean_cl_normal,col ="black") +  
  ggpubr::theme_classic2(base_size = 15) +
  ggplot2::theme(legend.position = "top")




##### Print Plots to pdf ######

pdf(file="../run/ICPMS_effects.pdf")
lapply(plots, cowplot::ggdraw)
dev.off()



quartz(height =12, width =6.5)

ggpubr::ggarrange(cowplot::ggdraw(plots$effect_Mg),
                  cowplot::ggdraw(plots$effect_P),
                  ncol=2)

#####
# Range effect


ionome <- read.csv("../data/AIR_ICPMS_processed.csv") %>%
  dplyr::rename(id="line")

ionome[duplicated(ionome$id),] %>%
  dplyr::pull(id) %>% table


#  Should I add this info?
donor_range <-read.csv("../data/donor_range.csv") %>%
  arrange(range, Code)

# discard outliers :/

no_outliers <- apply(ionome[,mineral_names],2,FUN=function(x){
  Q1 <- quantile(x, .25)
  Q3 <- quantile(x, .75)
  IQR <- IQR(x)
  x[x < Q1 -3*IQR | x>Q3+3*IQR] <- NA
  x
}) %>% as.matrix()

m <- no_outliers[,reliable]

ionome$donor <- factor(ionome$donor, levels = donor_range$Code )
levels(ionome$donor)


# nope this one does not work
# I needd the CML530 for comparison load the raw data
# add the CML530s

range_ionome <- rbind(
ionome %>%
  filter(donor == "CML530"),
ionome %>% 
  filter(donor %in% levels(cross$pheno$donor))
)


plots$K_range <- 
  range_ionome  %>%
  ggplot(aes(x=range, y =K, col =range, alpha=1)) +
  ylab("K ppm [mg/kg]") +
  ggbeeswarm::geom_quasirandom() +
  stat_summary( fun.data = mean_cl_normal,col ="black") +  
  ggpubr::theme_classic2(base_size = 15) +
  ggplot2::theme(legend.position = "none",
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

plots$P_range <- 
  range_ionome  %>%
  ggplot(aes(x= range, y =P, col=range, alpha=1)) +
  ylab("P ppm [mg/kg]") +
  ggbeeswarm::geom_quasirandom() +
  stat_summary( fun.data = mean_cl_normal,col ="black") +  
  ggpubr::theme_classic2(base_size = 15) +
  ggplot2::theme(legend.position = "none",
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

plots$Mg_range <- 
  range_ionome  %>%
  ggplot(aes(x=range, y =Mg, col =range, alpha=1)) +
  ylab("Mg ppm [mg/kg]") +
  ggbeeswarm::geom_quasirandom() +
  stat_summary( fun.data = mean_cl_normal,col ="black") +  
  ggpubr::theme_classic2(base_size = 15) +
  ggplot2::theme(legend.position = "none",
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )


quartz(width=5,height=10)
ggpubr::ggarrange(cowplot::ggdraw(plots$K_range),
                  cowplot::ggdraw(plots$P_range),
                  cowplot::ggdraw(plots$Mg_range),
                  ncol=3)

range_ionome
names(plots)



lm(data = range_ionome, K ~ range) %>% summary()
lm(data = range_ionome, P ~ range) %>% summary()
lm(data = range_ionome, Mg ~ range) %>% summary()
class(range_ionome[,reliable])
mm <- manova(cbind(Ca,Fe,K,Mg,Mn,P,S,Zn) ~ range, data=range_ionome)
summary.aov(mm)
mm
effectsize::eta_squared(mm)
summary.aov(mm)

mvs_test <- goft::mvshapiro_test(range_ionome[,reliable] %>% as.matrix())


library(MASS)

lda <- lda(range ~ Ca+Fe+K+Mg+Mn+P+S+Zn,data=range_ionome, CV = F)
abs(lda$scaling[,"LD1"]) %>% sort(decreasing = TRUE) 
abs(lda$scaling[,"LD2"]) %>% sort(decreasing = TRUE) 
quartz()

range_ionome <- cbind(range_ionome,predict(lda)$x)

range_ionome$range <- factor(range_ionome$range, levels=c("TMVB","SMCG","CML530"))
quartz()
range_ionome %>% arrange(range) %>%
ggplot() +
  geom_point(aes(x = LD1, y = LD2, color = range)) +
  theme_classic()


