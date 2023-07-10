library(dplyr)
library(qtl)

# Mendelian segregation calculations
# AA Aa aa
#AA
#Aa
#aa
AA <- c(1, 1/2, 0,
        0, 1/2, 1,
        0,   0, 0) %>% 
  matrix(nrow = 3, byrow = TRUE)

aa <- c(0,   0, 0,
        1, 1/2, 0,
        0, 1/2, 1) %>% 
  matrix(nrow = 3, byrow = TRUE)

Aa <- c(1/2, 1/4,   0,
        1/2, 1/2, 1/2,
        0,   1/4, 1/2) %>% 
  matrix(nrow = 3, byrow = TRUE)

S  <- c(1, 1/4, 0,
        0, 1/2, 0,
        0, 1/4, 1) %>% 
  matrix(nrow = 3, byrow = TRUE)

f1 <- c(0,1,0)

bc2  <- AA %*% AA %*% f1

bc2s3 <- (S %*% S %*% S %*% bc2)[,1]
#bc2s4 <- (S %*% bc2s3)[,1]

AIR_pheno <- read.csv("../run/AIR_PHENO.csv")

AIR_geno  <- read.csv("../run/AIR_GENO_ABH.csv", header = TRUE)
AIR_geno$id <- gsub(".*-","",AIR_geno$id)

chr <- AIR_geno[1,]

pheno_id <- AIR_pheno %>%
  select(id)

pheno_id$pheno <- TRUE


geno_id <- AIR_geno %>%
  select(id)

geno_id$geno <- TRUE

geno_pheno_id <-  geno_id %>%
  inner_join(pheno_id) %>%
  arrange(id) %>% select(id)



write.csv(
  geno_pheno_id %>% inner_join(AIR_pheno), 
  file = "../run/AIR_PHENO_sorted.csv",
  quote = FALSE, row.names = FALSE)

geno_out <- rbind(
  chr,
  geno_pheno_id %>% inner_join(AIR_geno)
)

write.csv(
  geno_out, 
  file = "../run/AIR_GENO_sorted.csv",
  quote = FALSE, row.names = FALSE)



BC2S3<-read.cross(
  format="csvs",
  genfile = "../run/AIR_GENO_sorted.csv", 
  phefile= "../run/AIR_PHENO_sorted.csv",
  BC.gen=2, F.gen=3)

summary(BC2S3)

gt <- pull.geno(BC2S3)
rownames(gt) <- BC2S3$pheno$id

# quartz()
# heatmap(gt, Colv = NA, col=c("red","green","purple"))



# Calcultae consensus map --------

# Remove distorted markers > 95% quantile
# Add this to maptools.R

for (chr in names(BC2S3$geno)){ 
  # due to excess heterozyocity
  # I had to renormalize the chi square distribution
  df <- BC2S3$geno[[chr]]$data %>% as.data.frame()
  gt_count <- apply(df,2,function(x) tabulate(x)[1:3])
  gt_count[is.na(gt_count)] <- 0
  chi_test <- apply(gt_count,2, 
                    function(x) {                          #A B H
                      chi <-  chisq.test(x = x, p = bc2s3[c(1,3,2)]) # ABH not AHB
                      chi$statistic
                    }
  )
  # make a function for this
  dfn <- ecdf(chi_test)
  fit <- smooth.spline(chi_test, dfn(chi_test))
  u <- predict(fit, sort(chi_test))$y
  u[u>0.999999] <- 1-1/sum(nmar(BC2S3))  # can't be > 1
  z <- qnorm(u)
  
  distorted <- z[z > 1.96] # 95% standard normal quantile
  length(distorted)/nrow(df)
  distorted <- distorted[!is.na(distorted)] %>% names
  chi_test[distorted]
  BC2S3 <- drop.markers(BC2S3, distorted )
}



map_ini <- est.map(BC2S3)



BC2S3 <- replace.map(BC2S3, map_ini)

# d <- lapply(map_ini, FUN = function(x){x - dplyr::lag(x)}) %>% unlist()
# quartz()
# hist(d, breaks= 100)


library(igraph)
# add this to maptools.R
get_quirky <- function (map, thresh = 10, c.size =1 ){
  lapply(map, FUN = function(x){
    m <- dist(x) %>% as.matrix()
    m[m>= thresh] <- 0
    diag(m) <- 0
    g <- igraph::graph_from_adjacency_matrix(m, mode ="undirected")
    c <- igraph::components(g)
    excluded <-which(c$membership %in% which(c$csize <= c.size))
    names(c$membership)[excluded]
  }
  ) %>% unlist()
}


quirky <- get_quirky(map_ini, thresh = 10, c.size =1)
length(quirky)
BC2S3 <- drop.markers(BC2S3, quirky)
length(map_ini %>% unlist())

consensus_map <- est.map(BC2S3)

quartz()
plot(consensus_map)


d <- lapply(consensus_map, FUN = function(x){x - dplyr::lag(x)}) %>% unlist()
mean(d, na.rm = TRUE)
median(d, na.rm = TRUE)
d[d>18]
# quartz()
# hist(d, breaks= 100)

# A little bit inflated compare to mapsline predictions
# ~10  gaps 10 to 25 cM in size


# checking if I need to shift the map
BC2S3
bc2s3
pull.map(BC2S3, as.table = TRUE) %>%
  group_by(chr) %>%
  dplyr::slice(1) 

# origins are 0 so I don't need to shift it manually

# Move this section to map_tools.R
# with detect_peaks
# and map_spline

# shifting map (adjusting origin to 0) for plotting

shifted_map <- pull.map(BC2S3, as.table = TRUE)
shifted <- shifted_map[shifted_map$chr==1,"pos"]
shift <- shifted[1]
shifted_map[shifted_map$chr==1,"pos"] <- shifted_map[shifted_map$chr==1,"pos"] - shift 

#fix maps for plotting
shifted <- shifted_map[shifted_map$chr==10,"pos"]
shift <- shifted[1]
shifted_map[shifted_map$chr==10,"pos"] <- shifted_map[shifted_map$chr==10,"pos"] - shift 

# I am happy!!
shifted_map %>%
  group_by(chr) %>%
  dplyr::slice(1) 



BC2S3 <- replace.map(BC2S3, consensus_map)
BC2S3
bc2s3


# If I removeI got an strange gentotype pattern. maybe it is because there no heterocygotes in this map
# they are coded as missing data

# Warning messages:
#   1: In summary.cross(x) : Strange genotype pattern.
#   2: In summary.cross(x) : Strange genotype pattern.


save(BC2S3, file="../run/AIR_cross.Rdata")

######################### THE END ##########################################


