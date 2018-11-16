
#! /usr/bin/env Rscript
# author: qinmao
# date: 2017-10-30

options(warn=-1);
if(!require("optparse"))
{
  install.packages("optparse", repos = "https://cloud.r-project.org");
  if(!library("optparse", logical.return = T))
  {
    stop("could not load library \"optparse\"\n");
  }
}
cat("load optparse done!\n");
library("base");
option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "file for 1)nucmer program output or 2)cleaned by perl script", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "out.png", help = "output file, default out.png", metavar = "character"),
  make_option(c("-f", "--filter"), type = "integer", default = 500, help = "minimum contig to be plot, default = 500 bp", metavar = "number"),
  make_option(c("-x", "--xsize"), type = "double", default = 800, help = "xlab size in px", metavar = "number"),
  make_option(c("-y", "--ysize"), type = "double", default = 800, help = "ylab size in px", metavar = "number"),
  make_option(c("-l", "--layout"), action = "store_true", default = FALSE,  help = "layout the plot"),
  make_option(c("-c", "--clean"), action = "store_true", default = FALSE, help = "data after by perl script"),
  make_option(c("-s", "--linesize"), action = "double", default = 3, help = "adjust the line size of plot, default = 3", metavar = "number")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$input))
{
  print_help(opt_parser);
  stop("The input file is requried!");
}

out <- opt$out;
xsize <- opt$xsize;
ysize <- opt$ysize;
filter <- opt$filter;
layout <- opt$layout;
clean <- opt$clean;
linesize <- opt$linesize;
cat("Parameters:\n");
cat("Out file :\t", out, "\n");
cat("x size in px :\t", xsize, "\n");
cat("y size in px :\t", ysize, "\n");
cat("filter size of contig :\t", filter, "\n");
cat("line size of plot :\t", linesize, "\n");
cat("is layout: \t", layout, "\n");

if(clean)
{
  plot.dat <- read.table(opt$input, sep = "\t", stringsAsFactors = FALSE);
  colnames(plot.dat) <-
    c("RefId",
      "QryId",
      "RefLen",
      "QryLen",
      "RefStart",
      "RefEnd",
      "QryStart",
      "QryEnd", 
      "RefOLLen",
      "BeRefStart",
      "BeRef");
  
  refs.level <- unique(plot.dat$RefId);
  refs.level <- refs.level[order(nchar(refs.level), refs.level)];
  qrys.level <- unique(plot.dat$QryId);
  refs.len <- 0
  refs.pos <- c()
  x_start <- 0
  
  # compute refs levels and their start and end position
  dat.first <- plot.dat[match(refs.level, plot.dat$RefId),];
  #head(dat.first);
  #cat(nrow(dat.first));
  for (i in 1:nrow(dat.first))
  {
    ref.len <- dat.first[i,3];
    refs.len <- refs.len + ref.len
    record <-
      c(as.character(refs.level[i]),
        x_start,
        x_start + ref.len,
        ref.len)
    refs.pos <- rbind(refs.pos, record)
    x_start <- x_start + ref.len
  }
  refs.pos <- as.data.frame(refs.pos, row.names = 1:nrow(refs.pos));
  colnames(refs.pos) <- c("RefId", "XStart", "XEnd", "RefLen");
  # refs.pos;
  # for(i in 2:4)
  # {
  #   refs.pos[,i] <- as.numeric(as.character(refs.pos[,i]));
  # }
  refs.pos[,-1] <- apply(refs.pos[,-1], 2, as.character);
  refs.pos[,-1] <- apply(refs.pos[,-1], 2, as.numeric);
  # compue qrys.len
  qrys.len <- 0;
  for(i in 1:length(qrys.level))
  {
    qrys <- subset(plot.dat, plot.dat$QryId == qrys.level[i]);
    qrys.len <- qrys.len +  qrys[1, 4];
  }
  rm(dat.first);
} else
{
  start.time <- Sys.time();
  dat <- readLines(opt$input);
  dat <- dat[sapply(strsplit(dat,"\\s", perl = TRUE), length) > 1][-1];
  cat("read the delta file done!", " Elapsed time:", Sys.time() - start.time, "s\n");
  
  start.time <- Sys.time();
  alns.start <- which(startsWith(dat, ">"));
  # adding the length of dat to alns.start
  alns.start <- c(alns.start, length(dat) + 1);
  groups <- rep(1:(length(alns.start) - 1), times = diff(alns.start));
  inner.func <- function(x){
    ids = strsplit(x[1], "\\s", perl = TRUE)[[1]];
    refid = strsplit(ids[1], ">")[[1]][2];
    qryid = ids[2];
    reflen = ids[3];
    qrylen = ids[4];
    alns.records.num <- length(x);
    records <- c();
    for(i in 2:alns.records.num)
    {
      alns = strsplit(x[i], "\\s", perl = TRUE)[[1]];
      record <- c(refid,
                  qryid,
                  reflen,
                  qrylen,
                  alns[1],
                  alns[2],
                  alns[3],
                  alns[4],
                  as.numeric(alns[2]) - as.numeric(alns[1]),
                  0,
                  NA);
      records <- rbind(records, record);
    }
    return(records);
  }
  plot.dat <- lapply(split(dat, groups), inner.func);
  #plot.dat <- do.call(rbind.data.frame, plot.dat);
  plot.dat <- data.frame(Reduce(rbind, plot.dat), row.names =  NULL);
  colnames(plot.dat) <-
    c("RefId",
      "QryId",
      "RefLen",
      "QryLen",
      "RefStart",
      "RefEnd",
      "QryStart",
      "QryEnd", 
      "RefOLLen",
      "BeRefStart",
      "BeRef");
  #library("ggplot2");
  #plot.dat <- as.data.frame(plot.dat, row.names = 1:nrow(plot.dat));
  # for (i in 3:8)
  # {
  #   plot.dat[, i] <- as.numeric(as.character(plot.dat[, i]))
  # }
  plot.dat[,c(1,2,11)] <- apply(plot.dat[,c(1,2,11)], 2, as.character);
  plot.dat[,3:10] <- apply(plot.dat[,3:10], 2, as.character);
  plot.dat[,3:10] <- apply(plot.dat[,3:10], 2, as.numeric);
  rm(dat);
  cat("#4 method, format delta file done!", "Elapse time:", Sys.time() - start.time, "s\n");
  #head(plot.dat);
  start.time <- Sys.time();
  plot.dat <- subset(plot.dat, plot.dat$QryLen >= filter);
  if(nrow(plot.dat) == 0)
  {
    stop("Error: Empty set after filtering ", filter, " bp.");
  }
  refs.level <- unique(plot.dat$RefId);
  refs.level <- refs.level[order(nchar(refs.level), refs.level)];
  qrys.level <- unique(plot.dat$QryId);
  refs.len <- 0
  refs.pos <- c()
  x_start <- 0
  
  # compute refs levels and their start and end position
  dat.first <- plot.dat[match(refs.level, plot.dat$RefId),];
  #head(dat.first);
  #cat(nrow(dat.first));
  for (i in 1:nrow(dat.first))
  {
    ref.len <- dat.first[i,3];
    refs.len <- refs.len + ref.len
    record <-
      c(as.character(refs.level[i]),
        x_start,
        x_start + ref.len,
        ref.len)
    refs.pos <- rbind(refs.pos, record)
    x_start <- x_start + ref.len
  }
  refs.pos <- as.data.frame(refs.pos, row.names = 1:nrow(refs.pos));
  colnames(refs.pos) <- c("RefId", "XStart", "XEnd", "RefLen");
  # refs.pos;
  # for(i in 2:4)
  # {
  #   refs.pos[,i] <- as.numeric(as.character(refs.pos[,i]));
  # }
  refs.pos[,-1] <- apply(refs.pos[,-1], 2, as.character);
  refs.pos[,-1] <- apply(refs.pos[,-1], 2, as.numeric);
  qrys.len <- 0;
  rm(dat.first);
  temp.dat <- c()
  
  ## compute each qrys belong to which reference and it's start position
  for(i in 1:length(qrys.level))
  {
    qrys <- subset(plot.dat, plot.dat$QryId == qrys.level[i]);
    qrys.sum = aggregate(qrys$RefOLLen, by = list(qrys$QryId, qrys$RefId), FUN = sum);
    belongref = qrys.sum[which(qrys.sum$x == max(qrys.sum$x))[1],2];
    temp <- qrys[as.character(qrys$RefId) == belongref, ];
    refstart <- temp[which(temp$RefOLLen == max(temp$RefOLLen))[1],5];
    qrys.len <- qrys.len + qrys[1, 4];
    plot.dat[which(plot.dat$QryId == qrys.level[i]),c("BeRef")] <-as.character(belongref);
    plot.dat[which(plot.dat$QryId == qrys.level[i]),c("BeRefStart")] <-refstart;
  }
  cat("compute belong and start position of reference done!", "Elapsed time:", Sys.time() - start.time, "s\n");
}

# plot a null graph
start.time <- Sys.time();
png(out, width = xsize, height = ysize);
xy_coords = matrix(
  data = c(-refs.len * 0.3, -qrys.len * 0.3, refs.len, qrys.len),
  nrow = 2,
  ncol = 2,
  byrow = T
)
colnames(xy_coords) <- c("x", "y")
#old.pars <- par();
#par(mar=c(5,4,1,1));
plot(
  xy_coords,
  type = "n",
  xlab = "Reference",
  ylab = "Assemblied",
  axes = F
)
#define box
{
  # bottom
xy_coords[1, 1] = 0 # x1
xy_coords[1, 2] = 0 # y1
xy_coords[2, 1] = refs.len # x2
xy_coords[2, 2] = 0 # y2
lines(xy_coords, col = "black")

# left
xy_coords[1, 1] = 0
xy_coords[1, 2] = 0
xy_coords[2, 1] = 0
xy_coords[2, 2] = qrys.len
lines(xy_coords, col = "black")

# right
xy_coords[1, 1] = refs.len
xy_coords[1, 2] = 0
xy_coords[2, 1] = refs.len
xy_coords[2, 2] = qrys.len
lines(xy_coords, col = "black")

# top
xy_coords[1, 1] = 0
xy_coords[1, 2] = qrys.len
xy_coords[2, 1] = refs.len
xy_coords[2, 2] = qrys.len
lines(xy_coords, col = "black")
}

# define a data.frame to store each chr
chr_start = 0
chr_length = 0
chrs_num = length(refs.level)

# define a dataframe to store each chr start and end position
chrs_position <-
  as.data.frame(matrix(NA, nrow = chrs_num, ncol = 4))

colnames(chrs_position) <- c("target", "start", "end", "length")

# plot basic chrs
y_start = 0
cat("Forward alignment is represented in Red colour!\n");
cat("Reverse alignment is represented in blue colour!\n");

for (i in 1:chrs_num)
{
  # subdata is the data which should be plot in the same reference.
  subdata = subset(plot.dat, as.character(plot.dat$BeRef) == refs.level[i]);
  chr_length = as.numeric(refs.pos[refs.pos$RefId == refs.level[i],][4]);
  # plot reference id;
  text((chr_start + chr_length / 2),
       0,
       labels = refs.level[i],
       adj = c(1, 1), 
       srt = 90);
  #new start position and plot reference border
  chr_start = chr_start + chr_length;
  if (i < chrs_num)
  {
    # vertical
    xy_coords[1, 1] = chr_start;
    # x1
    xy_coords[1, 2] = 0;
    # y1
    xy_coords[2, 1] = chr_start;
    # x2
    xy_coords[2, 2] = qrys.len;
    # y2
    lines(xy_coords, col = "gray60", lty = linesize);
  }
  # sort subdata
  subdata <- subdata[with(subdata, order(BeRefStart)), ];
  while (nrow(subdata) != 0)
  {
    qryid <- subdata[1, 2]
    qrylen <- subdata[1, 4]
    qrys <- subset(subdata, subdata$QryId == qryid)
    for (k in 1:nrow(qrys))
    {
      # according to original ref id to identify x start position;
      qrysrefid <- qrys[k, 1];
      x_start <- as.numeric(refs.pos[refs.pos$RefId == qrysrefid,][2]);
      col = "red"
      xy_coords[1, 1] = qrys[k, 5] + x_start;
      #x1
      xy_coords[2, 1] = qrys[k, 6] + x_start;
      #x2
      if (qrys[k, 7] > qrys[k, 8])
      {
        col = "blue"
        if(layout)
        {
          xy_coords[1, 2] = y_start + qrylen - qrys[k, 7] #y1
          xy_coords[2, 2] = y_start + qrylen - qrys[k, 8] #y2
        } else
        {
          xy_coords[1, 2] = y_start + qrys[k, 8] #y1
          xy_coords[2, 2] = y_start + qrys[k, 7] #y2
        }
      }
      else
      {
        xy_coords[1, 2] = y_start + qrys[k, 7] #y1
        xy_coords[2, 2] = y_start + qrys[k, 8] #y2
      }
      lines(xy_coords, col = col, lty = linesize);
    }
    # plot query id;
    text(
      0,
      y_start + qrylen / 2,
      labels = qryid,
      cex = 0.6,
      pos = 2
    );
    # plot query border
    if ((y_start + qrylen) != qrys.len)
    {
      xy_coords[1, 1] = 0;
      #x1
      xy_coords[2, 1] = refs.len #x2
      xy_coords[1, 2] = y_start + qrylen;
      #y1
      xy_coords[2, 2] = y_start + qrylen;
      #y2
      lines(xy_coords, col = "gray60", lty = linesize)
    }
    y_start = y_start + qrylen;
    subdata <- subdata[subdata$QryId != qryid, ];
  }
}
cat("plot done!", "Elapsed time:", Sys.time() - start.time, "s\n");
dev.off();
options(warn=0);
