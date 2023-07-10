library(qtl)
library(dplyr)
library(ggplot2)



mineral_names <-  c("Al", "B",  "Ca", "Fe",
                    "K",  "Mg", "Mn", "Mo",
                    "Ni", "P",  "S",  "Zn",
                    "Na", "Cu")

ionome_traits <- c("weight",mineral_names)

reliable <- c( "Ca", "Fe",  "K",
               "Mg", "Mn", "P", 
               "S",  "Zn")

reliable_traits <- c("weight",reliable)

unreilable <- c("Al", "B", "Mo","Ni","Na","Cu")


# load BC2S3
load(file="../run/AIR_cross.Rdata")
class(BC2S3)

# run scanone ------
n_threads <- 6
n_perm <- 1000
seed <- 1234567890
set.seed(seed)


BC2S3 <- calc.genoprob(BC2S3, step=1)

results <- list()


joint_result <- scanone(BC2S3,  pheno.col= ionome_traits, method="hk")


BC2S3_perms <- scanone(BC2S3,  pheno.col = ionome_traits,
                       n.perm = n_perm, method = "hk", 
                       n.cluster = n_threads)

BC2S3_thresh <- summary(BC2S3_perms, alpha = 0.05)


pdf(file = paste0("../run/ICPMS_QTL_joint_scans.pdf"), 
    height = 5, width = 15)
for( trait in mineral_names){
  in_col <- colnames(joint_result) %in% trait
  lod_col <- which(in_col) - 2
  thresh <- BC2S3_thresh[lod_col]
  ulim <- max(thresh,joint_result[[trait]])
  title <- paste("Joint Family Analysis", trait)
  plot(joint_result, lodcolumn=lod_col,
       main = title, 
       ylim = c(0, ulim))
  abline(h=thresh, col = "red")
}
dev.off()


# Process Peak QTLs
source("/Volumes/GoogleDrive/My\ Drive/repos/PTxB73_lipid_QTL_GxE/detect_peaks.R")
# had to modify these functminerals from the source file commit and update on orginal file
# move this fuction to maptools.R

get_peak_table <- function(cross = NULL, perms = NULL, 
                           single_scan = NULL,
                           controlAcrossCol=FALSE){
  
  perm_sum = summary(single_scan, perms=perms, alpha=0.05,
                     format = "tabByChr", pvalues=TRUE, ci.function="lodint",
                     expandtomarkers = TRUE) 
  
  lod_sum <- summary(perms, alpha = 0.05)
  
  lod_thresh <- data.frame(
    trait = dimnames(lod_sum)[[2]],
    thresh = as.vector(lod_sum))
  
  lod_adj_sum <- summary(perms , controlAcrossCol=controlAcrossCol, alpha = 0.05)
  
  #lod_sum <- summary(perms , controlAcrossCol=TRUE, alpha = 0.05)
  
  lod_thresh$thresh_adj <- as.vector(lod_adj_sum)
  
  lod_thresh
  empty <- sapply(perm_sum, is.null) %>% unlist
  if(is.null(empty)){
    perm_sum <- perm_sum[-which(empty)]  #remove empty elements
  }
  class(perm_sum) <- "list" # couldn't bind rows with a  list subclass.
  # has to stick to a list class
  peaks <- perm_sum %>% 
    dplyr::bind_rows() %>%
    mutate(id = unlist(lapply(perm_sum, rownames))) %>%
    tidyr::separate(id, c("trait","pseudom"), sep = " : ") %>%
    left_join(lod_thresh) %>%
    mutate(trait = factor(trait))
  
  peaks$marker <- find.marker(cross,chr = peaks$chr, pos = peaks$pos)
  peaks$marker_lod <-  get_peak_marker_lod(peaks,single_scan)
  peaks$ci_left  <- find.marker(cross,chr = peaks$chr, pos = peaks$ci.low) 
  peaks$ci_left  <-  gsub(".*_","", peaks$ci_left, perl =TRUE) %>% as.integer()
  peaks$ci_right <- find.marker(cross,chr = peaks$chr, pos = peaks$ci.high)
  peaks$ci_right <-  gsub(".*_","", peaks$ci_right, perl =TRUE) %>% as.integer()
  peaks$width <- (peaks$ci_right - peaks$ci_left)/1000000
  peaks
}

peaks <- get_peak_table(cross = BC2S3, perms = BC2S3_perms, single_scan = joint_result)

# Peak refinement for later #####

refine_peaks<- function(peaks = NULL, single_scan =NULL, perms =NULL, cross = NULL){
  # cross <-BC2S3 %>%
  #   #convert2riself() %>%
  #   jittermap() %>%
  #  calc.genoprob( step = 1)
  
  out_peaks <- peaks[0,]
  for(idx in 1:nrow(peaks)) {
    row <- peaks[idx,]
    trait_name <- as.character(row$trait)
    data <- as.data.frame(single_scan)[single_scan$chr == row$chr,c("pos",trait_name)]
    
    intensity <- 10
    akima_lod <- akimaInterp(data$pos,data[,trait_name], 
                             seq(from = min(data$pos), to = max(data$pos), by = 1/intensity))

    subpeaks <- pracma::findpeaks(akima_lod, 
                                  nups =3,
                                  minpeakdistance = 5*intensity, 
                                  minpeakheight = row$thresh,
                                  npeaks = 5 )
   
    if(length(subpeaks)==0){
      trio <- floor(intensity*(c(row$pos, row$ci.low, row$ci.high)))
      subpeaks <- matrix(c(row$lod, trio), nrow = 1)
    }
    
    # get the markers inside each subpeak with a drop of at least 1.5 lod
    
    in_drop <- t(apply(subpeaks,1,function(x) {
      d <- round(which((x[1] - akima_lod) > 1.5) - x[2],0)
      y <-  c( x[2] + max(d[d < 0]), x[2] + min(d[d > 0]))
      if(xor(y[1] == -Inf,y[2] == Inf)) {
        inf_idx <- which(y == Inf | y == -Inf)
        y[[inf_idx]] <- c(1,length(akima_lod))[inf_idx]} 
      else if (y[1] == -Inf & y[2] == Inf){
        inf_idx <- which( y == -Inf)
        y <- c(1,length(akima_lod))} 
      y
    }))
    row
    subpeaks
    in_drop 
    # Consider multiple subpeaks with a lod drop of 1.5
    if(nrow(subpeaks) > 1 ) {
      edge_list <- as.matrix(IRanges::findOverlaps(IRanges::IRanges(in_drop[,1], in_drop[,2])))
      
      g <- graph_from_edgelist(edge_list ,directed  = FALSE)
      # Iam asuumig small number of subpeaks
      # if subpeaks >10 this algorithm may never end
      
      ivs_list <- ivs(g, min=2)
      if(length(ivs_list) == 0){
        ivs_list <- ivs(g, min=1)
      }
      # Consider multiple solutions with the same vertex size
      vsize <- sapply(ivs_list, function(x) length(x)) 
      
      vs <- ivs_list[vsize == max(vsize)]
      
      is_max_lod <- sapply(ivs_list, function(x) max(subpeaks[x,1]) == max(subpeaks[,1]))
      
      # Select the solution containning the peak with max lod
      
      max_peak <- vs[which(is_max_lod)]    
    } else {
      max_peak <-list()
      max_peak[1] <- 1 
    }
    length(max_peak[1]) 
    length(max_peak[[1]]) 
    
    # Consider  multiple solutions with the same vertex size and including max lod
    # get the solution with the greatest lod sum.
    # Hideous!!!!
    
    if (length(max_peak) > 1 ){
      sum_lod <- sapply(max_peak, function(x) sum(subpeaks[x,1]))
      new_peaks <- cbind(subpeaks, in_drop)[ max_peak[[which(sum_lod == max(sum_lod))]],]
    }else if(length(max_peak)== 1 & length(max_peak[[1]]) > 1 ){
      new_peaks <- cbind(subpeaks, in_drop)[ max_peak[[1]],]
    } else if((length(max_peak) == 1) & (length(max_peak[[1]]) ==  1 )){
      # Remember this is an array of multiple vertices
      new_peaks <- matrix(cbind(subpeaks, in_drop)[max_peak[[1]],], nrow=1)
    }
    
    #finding the markers just outside
    chr_map <-  cross$geno[[row$chr]]$map
    
    bound_m <- t(apply(new_peaks,1, function(x) {
      left_pos <- max(chr_map[ chr_map  <= x[5]/intensity]) 
      right_pos <- min(chr_map[ chr_map >= x[6]/intensity]) 
      find.marker(cross, chr =row$chr, pos= c(left_pos,right_pos))
    }))
    
    out_row <- data.frame(
      chr = row$chr,
      pos = new_peaks[,2]/intensity,
      ci.low  = find.markerpos(cross, bound_m[,1])[,2],
      ci.high = find.markerpos(cross, bound_m[,2])[,2],
      lod = new_peaks[,1],
      pval =  sum(perms[,trait_name] > new_peaks[,1]) / length(perms[,trait_name]),
      trait = trait_name,
      pseudom = find.pseudomarker(cross, chr = row$chr, 
                                  pos = new_peaks[,2]/intensity,
                                  where = "prob"),
      thresh = row$thresh,
      thresh_adj = row$thresh_adj,
      marker = find.marker(cross, chr = row$chr, pos = new_peaks[,2]/intensity)
    )
    
    np_marker <- as.character(out_row$marker)
    np_trait <- as.character(row$trait)
    
    out_row$marker_lod <- single_scan[np_marker,np_trait]
    out_row$ci_left <- gsub(".*_","", bound_m[,1], perl =TRUE) %>% as.integer()
    out_row$ci_right <- gsub(".*_","", bound_m[,2], perl =TRUE) %>% as.integer()
    out_row$width <- (out_row$ci_right - out_row$ci_left) /1000000
    out_peaks <- rbind(out_peaks,out_row)
    
    # quartz()
    # plot(akima_lod, main = row$trait, type="l",
    #      xlab = paste("Chromosome", row$chr),
    #      ylab = "Akima Interpolated LOD")
    # points(subpeaks[, 2], akima_lod[subpeaks[, 2]], pch =2)
    # points(subpeaks[, 3], akima_lod[subpeaks[, 3]], pch = 0)
    # points(subpeaks[, 4], akima_lod[subpeaks[, 4]], pch = 0)
    # points(floor(intensity*out_row$ci.low), akima_lod[floor(intensity*out_row$ci.low)], pch = 15)
    # points(floor(intensity*out_row$ci.high), akima_lod[floor(intensity*out_row$ci.high)], pch = 15)
  }
  
  
  
  out_peaks <- out_peaks[out_peaks$lod > out_peaks$thresh,]
  out_peaks[out_peaks$lod > out_peaks$thresh_adj,]
  nrow(out_peaks)
  out_peaks
  
}




#consensus_map <- consensus.map
BC2S3 <- calc.genoprob(BC2S3)


new_peaks <- refine_peaks(peaks, 
                          perms = BC2S3_perms, single_scan = joint_result, cross=BC2S3)



library(GenomicRanges)

qtls <- GenomicRanges::GRanges(seqnames= new_peaks$chr, 
                               IRanges::IRanges(start=new_peaks$ci_left, end = new_peaks$ci_right),
                               chr = new_peaks$chr,
                               trait = new_peaks$trait,
                               lod = new_peaks$lod,
                               peak_marker = new_peaks$marker,
                               peak_marker_lod = new_peaks$marker_lod)

class(qtls)

# save(qtls,  file = paste0("QTL_intervals_joint.Rdata"))


gr <- GenomeInfoDb::sortSeqlevels(qtls)
onmap <- as.data.frame(gr)
onmap$start_bp <- onmap$start
onmap$end_bp <- onmap$start
onmap$peak_bp <- gsub("S.*_","",onmap$peak_marker, perl =TRUE) %>% as.integer()

mapdf <- pull.map(BC2S3,as.table = TRUE) %>% as.data.frame()
mapdf$bp <- gsub("S.*_","",rownames(mapdf), perl =TRUE) %>% as.integer()




onmap <- onmap %>%
  dplyr::filter(trait %in% reliable_traits) %>%
  dplyr::inner_join(mapdf, by= c(start="bp",chr="chr")) %>%
  dplyr::mutate(start =pos) %>%
  dplyr::select(-pos) %>%
  dplyr::inner_join(mapdf, by= c(end="bp",chr="chr")) %>%
  dplyr::mutate(end =pos) %>%
  dplyr::select(-pos) 



library(qtlTools)
trait_names <- sort(onmap$trait) %>% unique() 
nshow <- length(trait_names)


#pal <- c(RColorBrewer::brewer.pal(9,"Set1")[-6],"black")
#pal <- RColorBrewer::brewer.pal( nshow,"Set1")
pal <- RColorBrewer::brewer.pal( nshow,"Set3")
names(pal) <-trait_names
# names(pal) <- c("P","Mo","Mg","Mn","S","Cu","K","Zn","Fe")
scales::show_col(pal)
onmap

#calcCis(BC2S3, s1.output = joint_result, perm.output = BC2S3_perms)



# Move this section to map_tools.R
# with detect_peaks
# and map_spline

# shifting map (adjusting origin to 0) for plotting

shifted_map <- pull.map(BC2S3, as.table = TRUE)
pull.map(BC2S3)

#fix maps for plotting
shifted <- shifted_map[shifted_map$chr==10,"pos"]
shift <- shifted[1]
shifted_map[shifted_map$chr==10,"pos"] <- shifted_map[shifted_map$chr==10,"pos"] - shift 

# I am happy!!
shifted_map %>%
  group_by(chr) %>%
  dplyr::slice(1) 


BC2S3 <-replace.map(BC2S3,table2map(shifted_map))

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
    leg.lwd =5,
    max.lwd =5,
    min.lwd =2,
    orderBy = "lod",
    col= pal[legend_names],
    chrBuffer = c(0.3, -1),
    shift = FALSE,
    usr= c(0,350,1,12),
    xaxs = "e", 
    legendPosition = list(x=9.5,y=220),
    legendCex = 0.7,
    #legendPosition = list(x=12,y=100)
  )
  pl <- recordPlot()
  dev.off()
  pl
}
plots <- list()

reliable

plots$qtl_joint <-  cowplot::ggdraw(qtl_segments(BC2S3,onmap))  +
  cowplot::draw_label("Joint Family Analysis", y= 0.9, x = 0.35, hjust = 0)
pdf(file = "../run/ICPMS_joint_qtl_map.pdf")
plots$qtl_joint 
dev.off()

quartz()
plots$qtl_joint 
