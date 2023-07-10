library(dplyr)
library(ggplot2)

mineral_names <-  c("Al", "B",  "Ca", "Fe",
                    "K",  "Mg", "Mn", "Mo",
                    "Ni", "P",  "S",  "Zn",
                    "Na", "Cu")


reliable <- c( "Ca", "Fe",  "K",
               "Mg", "Mn", "P", 
               "S",  "Zn")

unreilable <- c("Al", "B", "Mo","Ni","Na","Cu")


pheno <- read.csv("../run/AIR_PHENO.csv")


# Family level effect  ------- 
# Maybe include CML530?


quartz()
pheno %>%
  ggplot(aes(x=donor,y=plant_height, group= donor, alpha=0)) +
  geom_hline(yintercept = mean(pheno$plant_height, na.rm = TRUE), 
             linetype = "dashed") +
  ggbeeswarm::geom_quasirandom() +
  geom_boxplot()

anova(lm(data= pheno, plant_height  ~ donor))
# quartz()
# hist(pheno$dta)


pheno %>%
  ggplot(aes(x=donor,y=dts, group= donor, alpha=0)) +
  geom_hline(yintercept = mean(pheno$dta, na.rm = TRUE), 
             linetype = "dashed") +
  ggbeeswarm::geom_quasirandom() +
  geom_boxplot()

anova(lm(data= pheno, ear_height  ~ donor))

quartz()
pheno %>%
  ggplot(aes(x=donor,y=dta, group= donor, alpha=0)) +
  geom_hline(yintercept = mean(pheno$dta, na.rm = TRUE), 
             linetype = "dashed") +
  ggbeeswarm::geom_quasirandom() +
  geom_boxplot()

lm(data= pheno, dta  ~ donor)
anova(lm(data= pheno, dta  ~ donor))
