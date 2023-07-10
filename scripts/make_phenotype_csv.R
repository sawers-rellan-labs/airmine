library(tidyr)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)

### Mineral PCA ####

mineral_names <-  c("Al", "B",  "Ca", "Fe",
                    "K",  "Mg", "Mn", "Mo",
                    "Ni", "P",  "S",  "Zn",
                    "Na", "Cu")


reliable <- c( "Ca", "Fe",  "K",
               "Mg", "Mn", "P", 
               "S",  "Zn")

unreilable <- c("Al", "B", "Mo","Ni","Na","Cu")


ionome <- read.csv("../data/AIR_ICPMS_processed.csv") %>%
  dplyr::rename(id="line")

ionome[duplicated(ionome$id),] %>%
  dplyr::pull(id) %>% table


#  Should I add this info?
# donor_range <-read.csv("../data/donor_range.csv") %>%
#   arrange(range, Code)

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

# initialize plot output list
plots <- list()

# Make PCA
pca <- PCA(m, ncp=10, graph = FALSE)

# quartz()
# fviz_pca_ind(
#   pca,axes = c(1,2),
#   col.ind = factor(ionome$donor), 
#   geom ="point",   addEllipses = TRUE) +
#   ggtitle("ionome PCA")


plots$pc12 <- fviz_pca_ind(
  pca,axes = c(1,2),
  col.ind = factor(ionome$range), 
  geom ="point",   addEllipses = TRUE) +
  ggtitle("ionome PCA")


plots$pc34 <- fviz_pca_ind(
  pca,axes = c(3,4), 
  col.ind = factor(ionome$range), 
  geom ="point", addEllipses = TRUE) +
  ggtitle("ionome PCA")


plots$cor_pc12 <- fviz_pca_var(pca,axes = c(1,2))


plots$cor_pc34 <- fviz_pca_var(pca,axes = c(3,4))


plot_con <- function(pca, axes = c(1,2)){
  l <- sqrt(length(axes))
  to_grid <- lapply( axes, 
                     function(pc){
                       fviz_contrib(pca, choice="var", axes = pc, top = ncol(pca$var$coord)) +
                         ggplot2::ggtitle(paste0("PC",pc))
                     }
  )
  ggpubr::ggarrange(plotlist = to_grid, 
                    ncol = ceiling(l),
                    nrow = floor(l))
}

plots$con <- plot_con(pca,1:6)



ionome <- cbind(
  ionome %>%
  select(id,range:ICPMS_sample,weight),
  no_outliers)


distro <- ionome  %>%
  gather(key="var", value="value", all_of(reliable)) %>% # Convert to key-value pairs
  dplyr::mutate(donor = factor(donor), var =factor(var)) 

plots$density <- distro  %>%
  ggplot(aes(value, group = donor)) +  # Plot the values
  xlab("Mass") +
  facet_wrap(~ var, scales = "free") +     # In separate panels
  geom_density( aes(col = donor)) + 
  ggpubr::theme_classic2() 


plots$box <- distro %>%# Convert to key-value pairs
  ggplot(aes(x=value, y = donor, group = donor, col= range)) + # Plot the values
  facet_wrap(~ var, scales = "free_x") +    # In separate panels
  geom_vline(data = distro %>% 
               group_by(var) %>%
               summarise(mean = mean(value, na.rm = TRUE)),
             aes(xintercept = mean)) +
  geom_boxplot(aes(alpha =1)) +
  scale_alpha(guide = 'none') +
  ggpubr::theme_classic2(base_size = 10) +
  ggplot2::theme(legend.position = "top")


# Hierarchical clustering of features ######

# Get centered and scaled data
# imputing NAs with column means


X <- scale(pca$call$X) %>% as.matrix()/sqrt(2) 

# correlation plot ####

cols.cor <- cor(X)

library(corrplot)
require(gridGraphics)
var_corplot <- function(x){
  quartz()
  corrplot(x)
  pl <- recordPlot()
  dev.off()
  pl
}

plots$cor <- cowplot::ggdraw(var_corplot(cols.cor))
class(plots$cor)


plots$cor

D <- as.dist(sqrt(2*(1-cols.cor)))


hc <- hclust(D, method = "ward.D")
plot(hc)

library(pvclust)

abscor.pv <- pvclust(X, method.hclust =  "ward.D", method.dist =  "abscor", nboot =1000)


vclust <- function(){
  plot(abscor.pv, print.num=FALSE, hang = -15)
  pvrect(abscor.pv)
}  



plots$vclust <- cowplot::ggdraw(vclust) 


##### Print Plots to pdf ######

pdf(file="../run/ICPMS_pca_qtl.pdf")
lapply(plots, cowplot::ggdraw)
dev.off()

# Make phenotype table ######
# Field traits (maybe I should change this name to whole plant_traits)

field_traits <-  read.csv("../data/RR-19-Field-Cycle_book_Clayton_A5C.csv")

field_traits[duplicated(field_traits$id),]

# there are 15 CML530 obserbations 
# and one for AN40, AN28, AN10, AN11AN12 each

# Get introgression lines 

field_traits <- field_traits %>%
  dplyr::select(id,plant_height:buggy_whip) %>%
  dplyr::filter(grepl("LANAP",id))

field_traits$donor <- field_traits$id
field_traits$donor <- paste0("AN",substr(field_traits$donor,6,7))


range_donor <- ionome %>% 
  dplyr::select(range:donor) %>%
  unique()


pheno <- field_traits %>% 
  dplyr::inner_join(range_donor)  %>% 
  dplyr::full_join(ionome %>% filter(grepl("LANAP",id))) %>%
  dplyr::select(id,range,collection,AN,donor, everything())

write.csv(pheno, "../run/AIR_PHENO.csv", quote = FALSE, row.names = FALSE)

