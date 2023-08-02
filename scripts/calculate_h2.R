library(qtl)
library(asreml)
library(ASRgenomics)
library(dplyr)
library(lme4)
library(ggplot2)



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
  dplyr::rename( id= "line") %>%
  dplyr::filter(donor != "CML530")  # select just introgression lines

ionome[duplicated(ionome$id),] %>%
  dplyr::pull(id) %>% table


## Using ANOVA


library(lme4)

# Family H2 -------

# ANOVA

var_comp <- data.frame()
for (t in reliable) { 
  t
  formula <- paste0(t," ~ range + donor")
  mod_fam <- lm(formula = formula,
                  data= ionome)
  aov_tab <- anova(mod_fam) %>% as.data.frame() %>%
    tibble::rownames_to_column("grp") %>%
    dplyr::mutate(trait =t) %>%
    dplyr::select(trait,grp, everything())
  var_comp <- rbind(var_comp, aov_tab)
}



fam_hered_aov <- janitor::clean_names(var_comp)  %>%
  dplyr::select(trait:grp,sum_sq) %>%
  tidyr::pivot_wider( names_from = grp, values_from = "sum_sq") %>%
  dplyr::mutate(H2 = donor/(donor+range+Residuals)) %>%
  dplyr::arrange(-H2)


# mixed model variance components


# LMM
var_comp <- data.frame()
for (t in reliable) {
formula <- as.formula(paste0(t," ~ 1 + (1| donor)"))
mod_fam <- lmer(formula,
                data= ionome , 
                control = lmerControl(sparseX = TRUE),
                REML=TRUE)
trait_vc <- VarCorr(mod_fam) %>% as.data.frame()
trait_vc$trait <- t
var_comp <- rbind(var_comp, trait_vc)
}

fam_hered_lme4 <- var_comp %>%
  dplyr::select(-var1,-var2,-sdcor) %>%
  tidyr::pivot_wider( names_from = grp, values_from = vcov) %>%
  dplyr::mutate(H2 = donor/(donor+Residual)) %>%
  arrange(-H2)

fam_hered_lme4 

# ASREML

load("AIR_cross.Rdata")
cross<-BC2S3
BC2S3 <- NULL
mineral_pheno <- cross$pheno


var_comp <- data.frame()
for (t in reliable) {
  fix_fx <- as.formula(paste0(t," ~ 1 + range"))
  mod_fam <- asreml(
    fixed = fix_fx, random = ~ donor,
    residual = ~ idv(units),
    data = mineral_pheno
  )
  
  
  trait_vc <- mod_fam$vparameters %>% t() %>% as.data.frame()
  trait_vc$trait <- t
  var_comp <- rbind(var_comp, trait_vc)
}

colnames(var_comp)[1:2] <- c("donor","Residual")

options(scipen = 0)

fam_hered_asreml <- var_comp %>%
  dplyr::mutate(H2 = donor/(donor+Residual)) %>%
  dplyr::select(trait,donor,Residual,H2) %>%
  arrange(-H2) %>% tibble::tibble()


fam_hered_lme4 
fam_hered_asreml 


# Inbred Line H2 ------

# I can't calculate inbred line H2 with anova because  I do not have replicates in this experiment
# perfect fit with line id, no residuals
lm( formula= P ~ range + id, data = mineral_pheno) %>%  anova() 

# But I can estimate h2 with the genotype data

snp <- qtl::pull.geno(cross,1:10)

snp[1:10,1:10]
dim(snp)
row.names(snp) <- cross$pheno$id

snp <- qtl::pull.geno(cross = cross, chr= 1:10) 

G <- AGHmatrix::Gmatrix( snp-1)
hist(G)
rownames(G) <- cross$pheno$id
colnames(G) <- cross$pheno$id

# G <- G.matrix(
#   M = snp-1,
#   method = "VanRaden",
#   na.string = "NA",
#   sparseform = FALSE,
#   digits = 8,
#   message = TRUE
# )



G.inv <- G.inverse(G = G, bend = FALSE, blend = FALSE, align = FALSE, rcn.thr=1e-12, sparseform = TRUE)



var_comp <- data.frame()

for (t in reliable) {
  fix_fx <- as.formula(paste0(t," ~ 1 + range"))
  mm <- asreml(
    fixed = fix_fx, random = ~ vm(id, G.inv$Ginv), ai.sing = TRUE,
    residual = ~ idv(units),
    data = cross$pheno
  )
  vc <- summary(mm)$varcomp
  trait_vc <- vc$component[1:2] %>% t() %>% as.data.frame()
  colnames(trait_vc)<- c("donor","Residual")
  trait_vc$trait <- t
  var_comp <- rbind(var_comp, trait_vc)
}


line_h2_asreml <- var_comp %>%
  dplyr::mutate(h2 = donor/(donor+Residual)) %>%
  dplyr::select(trait,donor,Residual,h2) %>%
  arrange(-h2) %>% tibble::tibble()


heritability <- fam_hered_aov %>%
  select(trait, H2_anova = H2) %>%
  dplyr::left_join(
  fam_hered_lme4 %>%
  select(trait, H2_lme4 = H2)) %>%
  
  dplyr::left_join(
    fam_hered_asreml %>%
      select(trait, H2_asreml = H2)
    ) %>% 
  
  dplyr::left_join(
  line_h2_asreml %>%
    select(trait, h2_asreml = h2)
) %>% arrange(-h2_asreml)

heritability$trait <- factor(heritability$trait)
heritability$trait <- forcats::fct_reorder(heritability$trait, -heritability$h2_asreml)

toplot <- heritability %>% 
  tidyr::pivot_longer(cols= c("H2_anova","H2_lme4","H2_asreml","h2_asreml"),
                      names_to = "type", values_to = "Heritability")

toplot$type <- factor(toplot$type, 
                      levels= c("H2_anova","H2_lme4","H2_asreml","h2_asreml" ))
quartz()
toplot %>%
  ggplot(aes(x=trait, y=Heritability, group=type, fill=type)) +
  geom_bar(stat="identity", position=position_dodge())+
  ggpubr::theme_classic2() +
  theme(legend.position = "top")
    

quartz()
toplot %>%
  filter(type %in% c("H2_anova","h2_asreml"))%>%
  ggplot(aes(x=trait, y=Heritability, group=type, fill=type)) +
  ylab("Heritability") +
  scale_fill_discrete(name = "Dose", labels = c("Wide sense", "Narrow sense")) +
  geom_bar(stat="identity", position=position_dodge())+
  ggpubr::theme_classic2() +
  theme(legend.title = element_blank())


# 
# #### old attempt
# 
# # Family level effect  ------- 
# # Maybe include CML530?
# 
# pheno <- read.csv("../run/AIR_PHENO.csv")
# 
# quartz()
# pheno %>%
#   ggplot(aes(x=donor,y=plant_height, group= donor, alpha=0)) +
#   geom_hline(yintercept = mean(pheno$plant_height, na.rm = TRUE), 
#              linetype = "dashed") +
#   ggbeeswarm::geom_quasirandom() +
#   geom_boxplot()
# 
# anova(lm(data= pheno, plant_height  ~ donor))
# # quartz()
# # hist(pheno$dta)
# 
# 
# pheno %>%
#   ggplot(aes(x=donor,y=dts, group= donor, alpha=0)) +
#   geom_hline(yintercept = mean(pheno$dta, na.rm = TRUE), 
#              linetype = "dashed") +
#   ggbeeswarm::geom_quasirandom() +
#   geom_boxplot()
# 
# anova(lm(data= pheno, ear_height  ~ donor))
# 
# quartz()
# pheno %>%
#   ggplot(aes(x=donor,y=dta, group= donor, alpha=0)) +
#   geom_hline(yintercept = mean(pheno$dta, na.rm = TRUE), 
#              linetype = "dashed") +
#   ggbeeswarm::geom_quasirandom() +
#   geom_boxplot()
# 
# lm(data= pheno, dta  ~ donor)
# anova(lm(data= pheno, dta  ~ donor))
# 
# #

