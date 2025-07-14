#########################################################################################################
#########################################################################################################
#########################################################################################################
#### Functions to use on markers prior to uploading to R/qtl ####
#Convert SNP from c(0|0, 1|1) to AA or BB 
convert_calls <- function(data_in){
  data_in_parents <- data_in[, c(162:163)]
  data_in <- data_in[,-c(162:163)]
  for(i in 1:nrow(data_in)){
    data_in[i, which(data_in[i,] == data_in_parents[i,'HA467'])] <- 'AA'
    data_in[i, which(data_in[i,] == data_in_parents[i,'PI_170415'])] <- 'BB'
  }
  return(data_in)
}

#compare sequential SNP triplets, if a SNP at the ends has greater similarity greater than threshold (0.98), keep the middle SNP
marker_filter2 <- function(data_in){
  mtd <- c()
  for(i in unique(data_in$chr)){
    temp_df <- data_in[data_in$chr == i,]
    temp_ncol <- ncol(temp_df)-2
    for(j in seq(1, (nrow(temp_df)-2),3)){
      pair1 <- length(which(temp_df[j,3:ncol(temp_df)] == temp_df[j+1,3:ncol(temp_df)]))/temp_ncol
      pair2 <- length(which(temp_df[j+1, 3:ncol(temp_df)] == temp_df[j+2, 3:ncol(temp_df)]))/temp_ncol
      if(pair1 > 0.98 | pair2 > 0.98){
        mtd <- c(mtd, row.names(temp_df[j+1,])) #mtd = markers to drop
      }
    }
  }
  if(length(mtd) > 0){
    mtd <- as.numeric(mtd)
    data_in <- data_in[-which(row.names(data_in) %in% mtd),]
    row.names(data_in) <- seq(1, nrow(data_in),1)
  }
  return(data_in)
}

#compare sequential SNP triplets, if a SNP at the ends is less similar to the middle SNP than threshold (0.25), drop end SNP
marker_filter3 <- function(data_in){
  flag <- TRUE
  while(flag == TRUE){
    mtd <- c()
    print('start loop')
    for(i in unique(data_in$chr)){
      temp_df <- data_in[data_in$chr == i,]
      temp_ncol <- ncol(temp_df)-2
      for(j in seq(1, (nrow(temp_df)-2),3)){
        pair1 <- length(which(temp_df[j,3:ncol(temp_df)] == temp_df[j+1,3:ncol(temp_df)]))/temp_ncol
        pair2 <- length(which(temp_df[j+1, 3:ncol(temp_df)] == temp_df[j+2, 3:ncol(temp_df)]))/temp_ncol
        if(pair1 < 0.25 | pair2 < 0.25){
          if(pair1 < 0.25 & pair2<0.25){
            mtd <- c(mtd, row.names(temp_df[j+1,])) #if both pairs suck drop middle marker
          } else { # else drop j if pair1 is bad or j+2 for pair2
            if(pair1 < 0.25){
              mtd <- c(mtd, row.names(temp_df[j,])) #mtd = markers to drop
            } else {
              mtd <- c(mtd, row.names(temp_df[j+2,]))
            }
          }
        }
      }
    }
    if(length(mtd) < 50){
      flag <- FALSE
    }
    if(length(mtd) != 0){
      print(length(mtd))
      mtd <- as.numeric(mtd)
      data_in <- data_in[-which(row.names(data_in) %in% mtd),]
      row.names(data_in) <- seq(1, nrow(data_in),1)
    } else {
      flag <- FALSE
    }
  }
  return(data_in)
}

#compare physical position of groups of SNP, if the physical distance of the first and last SNP are less than user defined vaue (dist_in), keep 1 SNP selected at random.
marker_filt_dist <- function(data_in, dist_in){
  mtk <- c()
  for(i in unique(data_in$chr)){
    temp_df <- data_in[data_in$chr == i,]
    # temp_ncol <- ncol(temp_df)-2
    for(j in seq(1, (nrow(temp_df)-4),5)){
      td2 <- temp_df[j:(j+4),]
      temp_dist <- td2[nrow(td2), 'pos'] - td2[1, 'pos']
      if(temp_dist < dist_in){
        mtk <- c(mtk,sample(row.names(td2),1))
      } else {
        mtk <- c(mtk, row.names(td2))
      }
    }
  }
  mtk <- as.numeric(mtk)
  mtk <- unique(mtk)
  data_in <- data_in[which(row.names(data_in) %in% mtk),]
  return(data_in)
}

library(qtl)
setwd("enter/your/directory")

gstart <- read.csv('genotypes_by_markers.csv')

#format data
temp_df <- data.frame('chr' = sapply(gstart$ID, function(x) strsplit(x, ':')[[1]][1]), 'pos' = sapply(gstart$ID, function(x) strsplit(x, ':')[[1]][2]))
temp_df$chr <- gsub('Ha412HO', '',temp_df$chr)
temp_df$chr <- gsub('Chr','S',temp_df$chr)
gstart <- cbind.data.frame(temp_df, gstart)
row.names(gstart) <- seq(1, nrow(gstart),1)
gstart <- gstart[,-3]
rm(temp_df)
gstart$pos <- as.numeric(gstart$pos)

#start filtering data. The order of filtering functions, number of times a filter, distances and thresholds is determined by the data and user
gstart2 <- marker_filt_dist(gstart,100000)
gstart3 <- marker_filt_dist(gstart2,100000)
gstart4 <- marker_filt_dist(gstart3,30000)
gstart5 <- marker_filt_dist(gstart4,30000)
gstart6 <- marker_filt_dist(gstart5,30000)
gstart7 <- marker_filt_dist(gstart6,10000)
gstart8 <- marker_filt_dist(gstart7,10000)


gstart9 <- marker_filter2(gstart8)
gstart10 <- marker_filt_dist(gstart9,50000)

gstart11 <- marker_filter3(gstart10)

#write.csv(gstart11, 'xxx.csv', row.names = F)


temp_pos <- gstart11[,c(1:2)]
# for input in to R/qtl, for each chr change physical position to sequence of 1 to number of rows by 1. 
for(i in unique(temp_pos$chr)){
  tp2 <- temp_pos[temp_pos$chr == i,]
  tp2$gp <- seq(1, nrow(tp2), 1)
  if(i == unique(temp_pos$chr)[1]){
    out <- tp2
  } else {
    out <- rbind.data.frame(out, tp2)
  }
}
# gstart9[1:10,1:10]
gstart14_1 <- convert_calls(gstart11)
gstart14_1[1:10, 1:10]
gstart14_1[1,]
gstart14_1_backup <- gstart14_1
#convert heterozygous calls to AB
gstart14_1[gstart14_1 == '0|1'] <- 'AB'
gstart14_1[gstart14_1 == '1|0'] <- 'AB'
#write.csv(gstart14_1, 'fileName.csv', row.names=F)
minput <- gstart14_1
out <- minput[,1:2]

#input into R/QTL requires a row containing chromosome number for corresponding SNP
for(i in unique(out$chr)){
  td <- out[out$chr == i,]
  td$pos2 <- seq(1, nrow(td),1)
  if(i == unique(out$chr)[1]){
    out2 <- td
  } else {
    out2 <- rbind.data.frame(out2, td)
  }
}
#formatting
minput$pos <- out2$pos2
minput <- as.data.frame(t(minput))
#minput[1:10, 1:10]
marker_names <- paste(out2$chr, out2$pos, sep = '_')
#marker_names[1:10]
colnames(minput) <-marker_names

minput <- cbind.data.frame(row.names(minput), minput)
row.names(minput) <- seq(1, nrow(minput),1)

minput[1:2,1] <- ''
colnames(minput)[1] <- 'x'
#minput[1:5, 1:10]

# write.csv(minput, 'input_to_rqtl.csv', row.names = F)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#### Functions to use on linkage map ####
#used genetic map and user defined distance. Data within function is passed to the function 'calc_distance'
marker_info <- function(map_in, dist_in){
  temp_out <- numeric()
  gt <- geno.table(map_in, scanone.output = T)
  for(i in unique(gt$chr)){
    temp_chr <- gt[gt$chr == i, ]
    chr_out <- calc_distance(temp_chr, dist_in)
    if(is.data.frame(chr_out) == TRUE){
      temp_out <- rbind.data.frame(temp_out, chr_out)
    }
  }
  return(temp_out)
}
# calculates genetic distance of SNP triplets. If distance of markers at the ends are a greater distance than user defined distance (dist_in), 
# keep name of middle SNP to be dropped
calc_distance <- function(chr_in, dist_in){
  temp_out <- numeric()
  temp_length <- nrow(chr_in)
  for(i in seq(2, (temp_length-1), 1)){
    temp_dist_pre <- chr_in[i, 'pos'] - chr_in[i-1, 'pos']
    temp_dist_post <- chr_in[i+1, 'pos'] - chr_in[i, 'pos']
    if(temp_dist_post > dist_in & temp_dist_pre > dist_in){
      temp_out2 <- chr_in[i,]
      temp_out <- rbind.data.frame(temp_out, temp_out2)
    }
  }
  return(temp_out)
}

# compares genetic position of neighboring SNP. If neighboring SNP are at the same genetic position, compare their chi-square test for Mendelian segregation.
# Drop SNP with greater deviation from expected segregation ratio
thinning_loop <- function(map_in){
  gt<-geno.table(map_in, scanone.output=TRUE)
  snp<-as.character(rownames(gt))
  out<-data.frame(id=as.character(), dist=as.numeric(), pos=as.numeric(), neglog10P=as.numeric(), missing=as.numeric())
  #out1<-data.frame(pos=as.numeric(), neglog10P=as.numeric(), missing=as.numeric())
  len<-as.numeric(length(snp))
  for(i in snp[2:len]){
    t<-which(rownames(gt)==i)
    pre<-t-1
    if(round(gt$pos[t],digits=0)==round(gt$pos[pre], digits=0)){
      dist<-round(gt$pos[t],digits=0) - round(gt$pos[pre], digits=0)
      if(gt$neglog10P[t]>gt$neglog10P[pre]){
        ugh<-gt[t, c(2:4)]
      } else{
        ugh<-gt[pre,c(2:4)]
      }
      ugh[,"id"]<-paste(rownames(ugh))
      out<-rbind(out, cbind(id=ugh[[4]], dist=dist, pos=ugh[[1]], neglog10P=ugh[[2]],missing=ugh[[3]]))}
    #rm(ugh)
  }
  out
  # length(unique(out$id))
  todrop<-unique(out$id)
  return(todrop)
}

# Similar to function 'calc_distance', but looks for blocks of SNP that are greater than a distance from the neighboring SNP
# dist_in: distance threshold for dropping block of SNP
# num_mar_in: block size of SNP to look for
# inner_dist: distance threshold for first and last SNP in block
calc_distance2 <- function(chr_in, dist_in, num_mar_in, inner_dist){
  temp_out <- numeric()
  temp_length <- nrow(chr_in)
  window_width <- num_mar_in-1
  for(i in seq(2, (temp_length-num_mar_in), 1)){
    if(i+num_mar_in < nrow(chr_in)){
      if(chr_in[(i+window_width), 'pos'] - chr_in[i, 'pos'] < inner_dist){ #if the marker distance of 1st and last SNP in window is less than...
        temp_dist1 <- chr_in[i, 'pos'] - chr_in[i-1, 'pos'] #difference from start of window and marker just before (i - 1) the window
        temp_dist2 <- chr_in[i + num_mar_in, 'pos'] - chr_in[i+window_width, 'pos'] #difference from end of window (i + window_width) and marker just after window (i + num_mar_in)
        if(temp_dist1 > dist_in & temp_dist2 > dist_in){
          temp_out2 <- chr_in[i:(i+(num_mar_in-1)),]
          temp_out <- rbind.data.frame(temp_out, temp_out2)
        }
      }
    }
  }
  return(temp_out)
}

#formats data for input into 'calc_distance2' and returns output
marker_info2 <- function(map_in, dist_in, num_mar_in, inner_dist){
  temp_out <- numeric()
  gt <- geno.table(map_in, scanone.output = T)
  for(i in unique(gt$chr)){
    temp_chr <- gt[gt$chr == i, ]
    if(nrow(temp_chr) > num_mar_in){
      chr_out <- calc_distance2(temp_chr, dist_in, num_mar_in, inner_dist)
      if(is.data.frame(chr_out) == TRUE){
        # print(i)
        temp_out <- rbind.data.frame(temp_out, chr_out, inner_dist)
      }
    }
  }
  if(sum(is.na(temp_out)) != 0){
    temp_out <- na.omit(temp_out)
  }
  return(temp_out)
}

#calculates genetic distance between neighboring markers
distbetweenMarkers <- function(data_in){
  data_in$dist <- NA
  colnames(data_in)[1] <- 'pos'
  for(i in 2:nrow(data_in)){
    data_in[i, 'dist'] <- round(data_in[i,'pos'] - data_in[i-1,'pos'],2)
  }
  return(data_in)
}
#########################################################################################################
#########################################################################################################
#########################################################################################################

# This code with estimate linkage map. Depending on number of markers, can take a very long time to run.
library(qtl)
setwd("your/directory")
b <- read.cross('csv', "/your/directory/",
                "input_to_rqtl.csv", estimate.map = F, genotypes = c('AA', 'AB', 'BB'))
b<-convert2riself(b)
nm<-est.map(b, error.prob=0.0025, map.function="kosambi")

b<-replace.map(b,nm)
# write.cross(b, 'csv', 'new2_map_v1')

#########################################################################################################
#########################################################################################################
#########################################################################################################
# Code below provides examples of how filters are run. 
# User defined variables, which filters are used, when they are used, and how many times they are used is unique for each dataset.  
 library(qtl)
 setwd("your/directory")
 b <- read.cross('csv', "your/directory",
                 "new2_map_v1.csv", estimate.map = F, genotypes = c('AA', 'AB', 'BB'))
b<-convert2riself(b)
nm<-est.map(b, error.prob=0.0025, map.function="kosambi")
b<-replace.map(b,nm)
plot.map(b)
b <- est.rf(b)
plotRF(b, col.scheme = 'redblue', mark.diagonal = T)
badMar1 <- marker_info(b, 5)
badMar1 <- badMar1[badMar1$chr != 17,]
mtd <- row.names(badMar1)
b <- drop.markers(b, mtd)
nm<-est.map(b, error.prob=0.0025, map.function="kosambi")
b<-replace.map(b,nm)
# plot.map(b, nm)

badMar1 <- marker_info(b, 5)
badMar1 <- badMar1[badMar1$chr != 17,]
mtd <- row.names(badMar1)
badMar2 <- marker_info2(b,50,2,10)
badMar3 <- marker_info2(b,50,3,10)
badMar4 <- marker_info2(b,50,4,10)
mtd <- c(row.names(badMar1), row.names(badMar2), row.names(badMar3), row.names(badMar4))
b <- drop.markers(b, mtd)
nm<-est.map(b, error.prob=0.0025, map.function="kosambi")
# plot.map(b, nm)
b<-replace.map(b,nm)

thin1 <- thinning_loop(b)
b2 <- drop.markers(b, thin1)
nm<-est.map(b2, error.prob=0.0025, map.function="kosambi")
b2 <- replace.map(b2, nm)
plotRF(b2, col.scheme = 'redblue', mark.diagonal = T)

badMar1 <- marker_info(b2, 5)
badMar1 <- badMar1[badMar1$chr != 17,]
b2 <- drop.markers(b2, mtd)
nm<-est.map(b2, error.prob=0.0025, map.function="kosambi")
# plot.map(b, nm)
b2 <-replace.map(b2,nm)

thin2 <- thinning_loop(b2)
b2 <- drop.markers(b2, thin2)
nm<-est.map(b2, error.prob=0.0025, map.function="kosambi")
b2 <- replace.map(b2, nm)

badMar2 <- marker_info2(b2,6,2,3)
badMar3 <- marker_info2(b2,6,3,4)
badMar4 <- marker_info2(b2,6,4,5)
mtd <- c(row.names(badMar2), row.names(badMar3), row.names(badMar4))

b2 <- drop.markers(b2, mtd)
nm<-est.map(b2, error.prob=0.0025, map.function="kosambi")
b3 <-replace.map(b2,nm)

write.cross(b3, 'csv', 'final_map')
