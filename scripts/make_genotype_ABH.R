library(dplyr)
library(vcfR)
library(ggplot2)

# Add dummy parents so tassel can convert to ABH ----

vcf_input <- read.vcfR("../run/AIR_B73_imputed_1.vcf")
fixed <- getFIX(vcf_input)

# get biallelic sites                         & exclude indels
biallelic <- vcf_input [is.biallelic(vcf_input) & !is.na(fixed[,"ALT"]) & !is.na(fixed[,"REF"])]

# Get taxa by parental group either CML530 or ANDOSOL
recur_name <- "CML530"
is_recur <- grepl("CML530", colnames(biallelic@gt))
sum(is_recur)
donor_name <- "ANDOSOL"
is_donor <- grepl("-AN", colnames(biallelic@gt))
sum(is_donor)

# Get FORMAT column
is_recur[1] <- TRUE
is_donor[1] <- TRUE

# biallelic[,is_recur]@gt
# just 33 ref  sites for CML530
# grep("1/1",biallelic[,is_recur]@gt)

# biallelic[,is_donor]@gt
# Get CML530 (5) samples genotypes

# Build consensus genotype for dummy parent
get_consensus <- function(gt, name = "consensus"){
  cs <- apply(gt,1,function(x) names(which.max(table(x)))) 
  # add missing data 
  cs[cs=="NULL"] <- "./."
  cs <- cs %>% t() %>% t
  colnames(cs)[1] <- name
  cs
  #cs %>% t() %>% t() %>% paste0(":.:.:.:.:.:.:.")
}

# Switch genotypes
# I use this function to obtain the opposite genotype to CML530

switch_gt <- function(gt){
  cnames <-colnames(gt)
  rnames <- rownames(gt)
  #Switch homozygous ref to alt
  switched<- apply(gt,2,function(x) gsub("1/1", "0|0", x))
  #Switch homozygous alt to ref
  switched  <- gsub("0/0", "1|1", switched)
  # remove phasing separator
  switched <- apply(switched,2,function(x) gsub("|", "/",x, fixed = TRUE) )
  rownames(switched) <- rnames
  colnames(switched) <- cnames
  switched
}

# make dummy parents based on consensus

recur_gt_in <- extract.gt(biallelic[,is_recur])

donor_gt_in <- extract.gt(biallelic[,is_donor])
donor_cs <- get_consensus(donor_gt_in,  name= donor_name)

recur_gt_out <- get_consensus(recur_gt_in, name= recur_name)
donor_gt_out <- switch_gt(recur_gt_out)

colnames(donor_gt_out)[1] <- donor_name

# which(recur_cs == donor_cs)
# which(switch_gt(donor_cs) != recur_cs)
# which(switch_gt(recur_cs) != donor_cs)
nrow(biallelic)

out_gt <- cbind(
  biallelic@gt[,"FORMAT"], 
  recur_gt_out, 
  donor_gt_out,
  biallelic@gt[,-1])

colnames(out_gt)[1]<-"FORMAT"

biallelic@gt <- out_gt

write.vcf(biallelic, file ="../run/AIR_dummy_parents.vcf.gz")

# gunzip to produce AIR_dummy_parents.vcf
R.utils::gunzip("../run/AIR_dummy_parents.vcf.gz", remove=FALSE)


# Make MDS from genetic data -------

rTASSEL::startLogger(fullPath = "../run/", fileName = "MDS.log")


# Load VCF data
genoPathVCF <- file.path(
  "..","run",
  "AIR_dummy_parents.vcf"
)
genoPathVCF 

# Load in VCF file
tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(
  path = genoPathVCF
)


tasDist <- rTASSEL::distanceMatrix(tasObj = tasGenoVCF)
nrow(tasGenoVCF@jGenotypeTable)

# Exclude dummy parents
MDS <- mds(tasDist)  %>% 
  dplyr::filter(!grepl("^ANDOSOL$|^CML530$",Taxa)) %>% # remove dummy parents
  dplyr::arrange(-PC1)



# Extract line, and MDS plot labels info from the sample name (Taxa)
MDS$label <- gsub("-Parent.*", "",MDS$Taxa, perl = TRUE)
MDS$label <- gsub("^.*-", "",MDS$label, perl = TRUE)
MDS$label <- gsub("RIL", "",MDS$label, perl = TRUE)
MDS$line <- MDS$label
MDS$label <- gsub("^.*LANAP", "",MDS$label, perl = TRUE)
MDS$label <- gsub("\\d\\d\\d$", "",MDS$label, perl = TRUE)
MDS$is_parent <- !grepl("^\\d",MDS$label)
MDS$label <- gsub("CML", "CML530",MDS$label, perl = TRUE)
MDS$label[MDS$is_parent]
MDS$line
MDS$label[!MDS$is_parent] <- paste0("AN",MDS$label[!MDS$is_parent])


# These duplicated lines are the ones that have a RIL attached to their name
# and LANAP36022 :/

dup_lines <- MDS$line[duplicated(MDS$line)] %>% sort %>% unique()

MDS %>%
  dplyr::arrange(line) %>%
  filter(line %in% dup_lines)

# remove RILs and duplicates
MDS <- MDS %>%
  dplyr::filter(!grepl("RIL|P05-H07-LANAP36022", Taxa))  %>%
  dplyr::group_by(line)%>%
  dplyr::slice(1) %>%
  dplyr::arrange(line)



llab <- MDS %>% 
  group_by(line) %>%
  dplyr::filter(PC1>0.5) %>%
  dplyr::filter(is_parent) %>%
  dplyr::slice_max(PC2, n = 1, with_ties = FALSE)

rlab <- MDS %>% group_by(line) %>%
  dplyr::filter(PC1<0.5) %>%
  dplyr::filter(is_parent) %>%
  dplyr::slice_min(PC2, n = 1, with_ties = FALSE)



pdf(file ="../run/genetic_MDS.pdf")

MDS %>%
  dplyr::filter(!is_parent) %>%
  #dplyr::filter(!line %in% gene_outliers) %>%
  ggplot(aes(x= PC1, y=PC2, col= label, label= label)) +
  xlim(-0.5,1)+
  geom_text(size =2)  +
  ggrepel::geom_label_repel( size =2,
                             data = llab, max.overlaps = 100,# xlim = c(0.7,NA),
                             mapping = aes(x= PC1, y=PC2, label =label)) +
  ggrepel::geom_label_repel(size =2,
                            data = rlab, max.overlaps = 100, xlim = c(NA,-0.3),
                            mapping = aes(x= PC1, y=PC2, label =label)) +
  ggpubr::theme_classic2() +
  ggplot2::theme(legend.position = "none") 


# genetic outliers:

gene_outliers <-MDS %>%
  dplyr::filter(!is_parent) %>%
  dplyr::filter(PC1>0.2) %>%
  dplyr::select(line, everything(),-Taxa) %>% pull(line)


MDS %>%
  dplyr::filter(!is_parent) %>%
  dplyr::filter(!line %in% gene_outliers) %>%
  ggplot(aes(x= PC1, y=PC2, col= label, label= label)) +
  xlim(-0.5,1)+
  geom_text(size =2)  +
  ggrepel::geom_label_repel( size =2,
                             data = llab, max.overlaps = 100,# xlim = c(0.7,NA),
                             mapping = aes(x= PC1, y=PC2, label =label)) +
  ggrepel::geom_label_repel(size =2,
                            data = rlab, max.overlaps = 100, xlim = c(NA,-0.3),
                            mapping = aes(x= PC1, y=PC2, label =label)) +
  ggpubr::theme_classic2() +
  ggplot2::theme(legend.position = "none") 

dev.off()

cat(gene_outliers, sep = "\n", file = "../run/gene_outliers.list")

cat(MDS$Taxa[! MDS$Taxa %in% gene_outliers], sep = "\n", file = "../run/AIR_ICPMS_taxa.list")


with_gt <- MDS %>%
  dplyr::filter(is_parent==FALSE) %>%
  dplyr::filter(!line %in% gene_outliers) %>%
  dplyr::select(id=Taxa,line) %>%
  dplyr::arrange(line)
nrow(with_gt)



cat(with_gt$id, sep ="\n", file="../run/AIR_taxa_id_for_qtl.tab")

# Make ABH ---------

toABH_cmd <- paste(
 "'/Applications/TASSEL 5/run_pipeline.pl'",
 "-vcf","../run/AIR_dummy_parents.vcf",
 "-GenosToABHPlugin",
 "-o","../run/AIR_B73_ABH.csv",
 "-parentA","../data/A.txt",
 "-parentB","../data/B.txt",
 "-outputFormat c",
 "-endPlugin")


system(toABH_cmd)


# filter ABH ----
abh_in <- read.csv("../run/AIR_B73_ABH.csv") 
abh_in[1:5,1:5]
# sorted according to line (LANAP NUMBER)
abh_out <- rbind( abh_in[1,],
  abh_in %>% 
  inner_join(with_gt) %>% select(line,id,everything()) %>%
  arrange(line) %>%
  select(-line)
)

abh_out[1:5,1:5]

nrow(abh_in)
nrow(abh_out)

write.csv(abh_out, "../run/AIR_GENO_ABH.csv", quote = FALSE, row.names = FALSE) 




