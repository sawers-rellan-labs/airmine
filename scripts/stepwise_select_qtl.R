library(tidyr)
library(dplyr)
#library(purrr)
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

load(file="../run/QTL_intervals_joint.Rdata")

print("simulating genotypes for Multiple Interval Mapping")
cross <- BC2S3
cross.prob <- calc.genoprob(BC2S3)
cross.sim <- sim.geno(BC2S3)
BC2S3 <-NULL

n_perm <- 10
n_pheno <-  nphe(cross)
n_threads <- 4

qtl_count <- 
  as.data.frame(qtl_peaks) %>%
  dplyr::left_join(
    as.data.frame(qtl_peaks) %>%
      dplyr::group_by(trait) %>%
      dplyr::summarise(n = length(trait)) %>%
      dplyr::arrange(n) %>% 
      dplyr::mutate(trait_n = factor(paste0(trait,"_",n), levels = factor(paste0(trait,"_",n))))
  ) %>% 
  dplyr::arrange(trait_n,lod) %>% as.data.frame()


mim_qtl <- data.frame()
mim_fit <- list()
mim_rqtl <- list()
tidx <- 1

with_qtl <-  gsub("_.*","", levels(qtl_count$trait_n) ,perl =TRUE)
with_qtl


source("../run/detect_peaks.R")

for(t in with_qtl){
nt <- which(colnames(cross$pheno) == t)

# Calculate the two marker scan for estimating likelihood penalties.


# print(paste("calculating scantwo for:", t))


# perm2 <- scantwopermhk(cross, pheno.col= nt,
#                n.perm =1, 
#                verbose=TRUE)

# # scantwo(cross, pheno.col= nt,
#                  method="hk",
#                  n.perm =n_perm, 
#                  n.cluster=n_threads,
#                  verbose=TRUE)
#lpen <- calc.penalties(perm2)

# Calculate likelihood  penalties

# print(paste("calculatingLOD penalties for:", t))

# Use LOD penalties
print(paste("using LOD penalties for:", t))

# I have chosen this arbitrary penalty because computing time is
# way too much
lpen <- c(4,4,4)
names(lpen) <-c("main","int.light","int.heavy")

# Getting initial qtl_peaks from single marker analysis

chr <- qtl_count %>% 
  dplyr::filter(trait ==t) %>% 
  dplyr::arrange(trait_n,-lod) %>%
  dplyr::pull("chr")
pos <- qtl_count %>% 
  dplyr::filter(trait ==t) %>% 
  dplyr::arrange(trait_n,-peak_marker_lod) %>%
  dplyr::pull("peak_pos")

# for some reason cross object was not enough
qtl <- makeqtl(cross.prob, chr, pos, what = "prob") 

# Building an initial additive model

add_str <- paste(qtl$altname, collapse = " + ")
add_formula <- paste("y ~ ",add_str) %>% as.formula()

# # Fit the additive model
# fitqtl(cross, pheno.col=nt, qtl=qtl, formula= add_formula )


# Finaly stepwise selection including forward and  backward selection

stepout1 <- stepwiseqtl(cross.prob, pheno.col=nt,
                        method="hk",
                        additive.only=TRUE, 
                        qtl = qtl,
                        penalties = lpen,
                        verbose=TRUE)

if(attr(stepout1,"pLOD") == 0){ 
  print("All qtl_peaks dropped in multiqtl model!")
  next
}

# Initial additive multilocus model fit
# I need it in order to add interations for epistasis analysis

add_str <- paste(stepout1$altname, collapse = " + ")
add_formula <- paste("y ~ ",add_str) %>% as.formula()

additive <- fitqtl(cross.prob,
                   pheno.col=nt,
                   method="hk",
                   qtl=stepout1,
                   formula= add_formula )

# Check that I retained peaks after the multilocus analysis

if(!is.null(additive$result.drop)){
  
  by_var <- additive$result.drop %>% as.data.frame() %>% arrange(-`%var`) %>% row.names()
  
  qtl_step <- makeqtl(cross.prob,
                      what = "prob",
                      stepout1$chr[match(by_var,stepout1$name)], 
                      stepout1$pos[match(by_var,stepout1$name)]
  )
  
} else {
  qtl_step <- stepout1
}

# Fit again the full additive model with the selected markers

full.fit <- fitqtl(
   calc.genoprob(cross.prob), 
   pheno.col=nt, 
   qtl=qtl_step, 
   method = "hk",
   formula = paste("y ~ ", paste(qtl_step$altname, collapse = " + ")),
   get.ests = TRUE )


mim_fit[[t]] <- full.fit

# Get the MIM scan for plotting !
rqtl <- refineqtl(cross.prob, method = "hk", pheno.col=nt, qtl=qtl_step)

mim_rqtl[[t]] <- rqtl

# Get the peaks for visualization  with lodint
# My custom function refine_peaks frm the file detect_peaks.R fails with the MIM scans


mim_peaks <- 
  lapply(1:rqtl$n.qtl,function(x) { 
    lodpeak <- lodint(rqtl, qtl.index = x, drop=1.5, expandtomarkers = TRUE)
    lastidx <- nrow(lodpeak)
    data.frame(
      #  Dummy coordinates
      IRanges(start=lodpeak$pos[1]*1e6 , end = lodpeak$pos[lastidx]*1e6),
      chr = lodpeak$chr[1],
      trait = t,
      lod = lodpeak$lod[2],
      peak_pos = lodpeak$pos[2]*1e6,
      peak_marker = row.names(lodpeak)[2],
      peak_marker_lod = lodpeak$lod[2]) 
  }) %>% dplyr::bind_rows()

mim_qtl <- rbind(mim_peaks, mim_qtl)
tidx <- tidx + 1
}


mim_qtl %>% arrange(-peak_marker_lod)


single_scan <- scanone(
  cross.prob, 
  pheno.col = traits,
  method = "hk",
  verbose = TRUE)

cross_perms <- scanone(cross,  pheno.col = traits,
                       n.perm = 1000, method = "hk", 
                       n.cluster = 8, verbose = TRUE)

cross_thresh <- summary(cross_perms, controlAcrossCol=TRUE, alpha = 0.05)




# summary(single_scan, perms=cross_perms, alpha=0.05,
#         format = "tabByChr", pvalues=TRUE, ci.function="lodint",
#         expandtomarkers = TRUE) 

# get_peak_table(
#   cross = cross.prob,
#   perms = cross_perms, 
#   single_scan = single_scan,
#   controlAcrossCol=FALSE)

summary(single_scan)
  

mim_peaks <- data.frame()
mim_refined_peaks <- data.frame()

for(t in names(mim_rqtl)){

  print(t)
# r<-attr(mim_rqtl[[1]],"lodprofile")
# qtl_names<- names(r)
# chr_names<- gsub("@.*","",qtl_names)
# names(r) <- chr_names
 
result <- mim_rqtl[[t]]
result_sum <- summary(result) %>%as.data.frame()


# when result class is `qtl`
t
peak_count <- table(result$chr)

peak_table <- lapply(1:nrow(result_sum), FUN= function(x){
   peak <- lodint(result,
          result$chr[x],
          x, # just the qtl index
          drop=1.5,
          expandtomarkers=TRUE) %>%
     as.data.frame()
   peak$pos
   peak$marker <- rownames(peak)
   peak$n_qtl <- x
   peak$point <-c("ci_left","peak","ci_right")
   peak$name <-result$name[x]
   peak
}) %>% dplyr::bind_rows()

peak_pos <- peak_table %>%
  tidyr::pivot_wider(
    id_cols = c("chr","n_qtl","name"), 
    values_from = "pos", names_from = "point") %>%
  dplyr::rename(ci.low="ci_left",ci.high = "ci_right", pos="peak")


peak_markers <- peak_table %>%
  tidyr::pivot_wider(
    id_cols = c("chr","n_qtl"), 
    values_from = "marker", names_from = "point") %>%
  dplyr::rename(peak_marker="peak")


peak_lods <- peak_table %>%
  tidyr::pivot_wider(
    id_cols = c("chr","n_qtl"),
    values_from = "lod", names_from = "point") %>%
  dplyr::rename(lod = "peak",) %>%
  dplyr::select(-starts_with("ci"))

peaks <- peak_pos %>%
  inner_join(peak_lods) %>%
  inner_join(peak_markers)

peaks$marker_lod <- peaks$lod
peaks$ci_left  <-  gsub(".*_","", peaks$ci_left, perl =TRUE) %>% as.integer()
peaks$ci_right <-  gsub(".*_","", peaks$ci_right, perl =TRUE) %>% as.integer()
peaks$width <- (peaks$ci_right - peaks$ci_left)/1000000
peaks$trait <- t

lod_sum <- summary(cross_perms, alpha = 0.05)
lod_adj_sum <- summary(cross_perms , controlAcrossCol=TRUE, alpha = 0.05)

peaks$thresh <-  lod_sum[1,t]
peaks$thresh_adj <- lod_adj_sum[1,t]

single_scan <- attr(mim_rqtl[[t]],"lodprofile") 
names(single_scan)

# pass<- peaks$lod > peaks$thresh_adj

# if(any(pass)){
# solution for overlapping peaks
qnames <- names(single_scan)


narrow_peaks <- lapply(qnames, FUN = function(x){
    lp <-single_scan[[x]]
    names(lp)[3]<-t
    refine_peaks(
    peaks = peaks[peaks$name==x,], #assumes single row
    single_scan =  lp, 
    perms =cross_perms, cross = cross.prob) 
}) %>% dplyr::bind_rows()
  
# MIM peaks might have lower single marker lod than scanone trait lod thresh
# refine peaks should remove them

narrow_peaks  <- narrow_peaks[narrow_peaks$lod > narrow_peaks$thresh,]

narrow_peaks  <- narrow_peaks[narrow_peaks$name==narrow_peaks$new_name,]


mim_peaks <- rbind(mim_peaks,peaks)
mim_refined_peaks <- rbind(mim_refined_peaks,narrow_peaks)
} 



mim_refined_peaks[mim_refined_peaks$lod >mim_refined_peaks$thresh_adj,]

mim_names <- NULL
  
mim_names <-lapply(names(mim_rqtl), FUN = function(x){
  r <- mim_rqtl[[x]]
  data.frame(
    trait = x,
    name = r$name,
    chr = r$chr,
    pos = r$pos
  )
}) 

mim_names$marker <- find.marker(cross, chr = mim_names$chr, pos = mim_names$pos)


explained <- lapply(names(mim_fit), FUN = function(x){
  # This model combines additive AND dominant effects
  # as seen in  r$ests$ests
  r <- mim_fit[[x]]
  v <- r$result.drop[,4]
  if(is.null(v)){
    v <- r$result.full[1,5]
    names(v) <- gsub("a$","",names(r$ests$ests[2]))
  }
  data.frame(trait=x,name= names(v),var =v) 
}) %>% dplyr::bind_rows() 



write.csv(mim_refined_peaks %>%
              left_join(explained) , row.names=FALSE,
           file = paste0("../run/mim_qtl_explained_var.csv"))

# Selected QTLs  

new_peaks <-
  mim_refined_peaks %>%
  left_join(explained) %>%
  # filter(lod>thresh_adj) %>%
  filter(trait %in% reliable) %>%
  select(trait,name:lod,var,everything()) %>%
  arrange(-lod)

  
qtl_peaks <- GRanges(seqnames= new_peaks$chr, 
                     IRanges(start=new_peaks$ci_left, end = new_peaks$ci_right),
                     chr = new_peaks$chr,
                     trait = new_peaks$trait,
                     lod = new_peaks$lod,
                     peak_pos = new_peaks$pos,
                     peak_marker = new_peaks$marker,
                     peak_marker_lod = new_peaks$marker_lod) 



onmap <- as.data.frame(qtl_peaks)
onmap$start_bp <- onmap$start
onmap$end_bp <- onmap$start
onmap$peak_bp <- gsub("S.*_","",onmap$peak_marker, perl =TRUE) %>% as.integer()

mapdf <- pull.map(cross,as.table = TRUE) %>% as.data.frame()
mapdf$bp <- gsub("S.*_","",rownames(mapdf), perl =TRUE) %>% as.integer()




onmap <- onmap %>%
  dplyr::inner_join(mapdf, by= c(start="bp",chr="chr")) %>%
  dplyr::mutate(start =pos) %>%
  dplyr::select(-pos) %>%
  dplyr::inner_join(mapdf, by= c(end="bp",chr="chr")) %>%
  dplyr::mutate(end =pos) %>%
  dplyr::select(-pos) 


#calcCis(BC2S3, s1.output = joint_result, perm.output = BC2S3_perms)



# Move this section to map_tools.R
# with detect_peaks
# and map_spline

# shifting map (adjusting origin to 0) for plotting

shifted_map <- pull.map(cross, as.table = TRUE)


#fix maps for plotting
shifted <- shifted_map[shifted_map$chr==10,"pos"]
shift <- shifted[1]
shifted_map[shifted_map$chr==10,"pos"] <- shifted_map[shifted_map$chr==10,"pos"] - shift 

# I am happy!!
shifted_map %>%
  group_by(chr) %>%
  dplyr::slice(1) 


library(qtlTools)


trait_names <- sort(onmap$trait) %>% unique() 
nshow <- length(trait_names)



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#pal <- c(RColorBrewer::brewer.pal(9,"Set1")[-6],"black")
#pal <- RColorBrewer::brewer.pal( nshow,"Set1")
#pal <- RColorBrewer::brewer.pal( nshow,"Set2")
pal <- gg_color_hue (nshow)
names(pal) <-trait_names
# names(pal) <- c("P","Mo","Mg","Mn","S","Cu","K","Zn","Fe")
scales::show_col(pal)



cross <-replace.map(cross,table2map(shifted_map))

qtl_segments <-  function(cross,x){
  nshow <- sort(x$trait) %>% unique() %>% length()
  legend_names <- sort(onmap$trait) %>% unique()
  pal[legend_names]
  quartz()
  segmentsOnMap(
    cross,
    phe = x$trait,
    chr = x$chr,
    peaklod = x$peak_marker_lod,
    l= x$start,
    h = x$end,
    lwd = "byLod",
    leg.lwd =6,
    max.lwd =6,
    min.lwd =3,
    orderBy = "lod",
    col= pal[legend_names],
    chrBuffer = c(0.3, -1),
    shift = FALSE,
    usr= c(0,350,1,12),
    xaxs = "e", 
    legendPosition = list(x=9,y=260),
    legendCex = 0.7,
    cex.lab=3,
    cex.axis=3, 
    ex.main=3,
    cex.sub=3
    #legendPosition = list(x=12,y=100)
  )
  pl <- recordPlot()
  dev.off()
  pl
}
plots <- list()
plots$qtl_joint <-NULL

plots$qtl_joint <-  cowplot::ggdraw(qtl_segments(cross,onmap))  +
  cowplot::draw_label("Joint Family Analysis", y= 0.9, x = 0.2, 
                      hjust = 0,size = 20)








