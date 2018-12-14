#! /usr/bin/env Rscript
# author: qinmao
# date: 2018-12-13

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
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "a list for minimap alignment files", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "out.png", help = "output file, default out.png", metavar = "character"),
  make_option(c("-t", "--title"), type = "character", default = "AssemblyCompare", help = "The title for this plot, default = \"AssemblyCompare\".", metavar = "character"),
  make_option(c("-x", "--xsize"), type = "double", default = 800, help = "xlab size in px, default = 800.", metavar = "number"),
  make_option(c("-y", "--ysize"), type = "double", default = 800, help = "ylab size in px, default = 800.", metavar = "number"),
  make_option(c("-r", "--olratio"), type = "double", default = 0.1, help = "The ratio of overlap, default = 0.1", metavar = "number"),
  make_option(c("-g", "--gap"), type = "integer", default = 100000, help = "merge the alignme within the gap, default = 100000 bp", metavar = "number"),
  make_option(c("-p", "--legendposition"), type = "character", default = "topright", help = "The legend position, default = \"topright\".", metavar = "character"),
  make_option(c("-l", "--labels"), type = "character", default = NULL, help = "The Labels of data, according to alignment.", metavar = "character"),
  make_option(c("-c", "--cex"), type = "double", default = 0.5, help = "The cex parameter, default = 0.5.", metavar = "number"),
  make_option(c("-b", "--breakmarks"), type = "integer", default = 10, help = "Break reference to setting marks, default = 10.", metavar = "number")
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
olratio <- opt$olratio;
title <- opt$title;
legendposition <- opt$legendposition;
cex <- opt$cex;
gap <- opt$gap;
breakmarks <- opt$breakmarks;
labels <- opt$labels;
cat("Parameters:\n");
cat("Out file :\t", out, "\n");
cat("x size in px :\t", xsize, "\n");
cat("y size in px :\t", ysize, "\n");
cat("Title :\t", title, "\n");
cat("legend position:\t", legendposition, "\n");
cat("cex: \t", cex, "\n");
cat("gap: \t", gap, "\n");
cat("olratio: \t", olratio, "\n");
cat("break marks: \t", breakmarks, "\n");
cat("labels: \t", labels, "\n");

cat("Reading alignment file starting...\n");
col.names <- c("CId","CLen", "CStart", "CEnd", "COrient", "RId", "RLen", "RStart", "REnd");
files <- readLines("./files.txt", n = -1);
if(is.null(labels))
{
	labels <- files;
} else {
	tmp <- strsplit(labels, ",");
	labels <- tmp[[1]];
	if(length(labels) != length(files))
	{
		stop("The labels is not equal to the alignment files size!");
	}
}
dat <- list();
for(i in 1:length(files))
{
   tmp <- readLines(files[i], n = -1);
   tmp <- strsplit(tmp, "\t");
   x <- t(sapply(tmp, FUN = function(x){x[1:9]}, simplify = "array"));
   x <- as.data.frame(x);
   colnames(x) <- col.names;
   x[, c(1,5,6)] <- apply(x[,c(1,5,6)], 2, as.character);
   x[, -c(1,5,6)] <- apply(x[,-c(1,5,6)], 2, as.character);
   x[, -c(1,5,6)] <- apply(x[,-c(1,5,6)], 2, as.numeric);
   dat[[i]] <- x;
}
cat("reading alignment file done...\n");

# source findBinCov function
findBinCov <- function(dat)
{
  fEnd = 1;
  output <- c();
  for(i in 1:nrow(dat))
  {
    cStart = dat[i, "RStart"];
    cEnd = dat[i, "REnd"]
    if(cStart >= fEnd)
    {
      output <- rbind(output, c(fEnd, cStart, 0));
      output <- rbind(output, c(cStart, cEnd, 1));
      fEnd <- cEnd;
    } else
    {
      search.depth <- c();
      locus <- c();
      for(j in nrow(output):1)
      {
        if(output[j,2] > cStart)
        {
          search.depth <- c(search.depth, j);
          if(!(output[j,1] %in% locus))
          {
            locus <- c(locus, output[j,1]);
          }
          if(!(output[j,2] %in% locus))
          {
            locus <- c(locus, output[j,2]);
          }
        } else
        {
          break;
        }
      }
      if(!(cStart %in% locus))
      {
        locus <- c(locus, cStart);
      }
      if(!(cEnd %in% locus))
      {
        locus <- c(locus, cEnd);
      }
      locus <- sort(locus);
      bins <- c();
      for(j in 1:(length(locus) - 1))
      {
        bins <- rbind(bins, c(locus[j], locus[j + 1], 0));
      }
      for(j in 1:nrow(bins))
      {
        for(k in search.depth)
        {
          if(bins[j,2] <= output[k,2])
          {
            bins[j,3] <- bins[j,3] + 1;
          }
        }
        if((bins[j,2] > cStart) && (bins[j,2] <= cEnd))
        {
          bins[j,3] <- bins[j,3] + 1;
        }
      }
      output <- output[-search.depth,];
      output <- rbind(output, bins);
      fEnd <- locus[length(locus)];
    }
  }
  output <- as.data.frame(output);
  for(j in 1:ncol(output))
  {
    output[,j] <- as.numeric(as.character(output[,j]));
  }
  output;
}

mergeAlignments <- function(dat, gap = 1000, ol = 0.1)
{
  c.levels <- unique(dat$CId);
  if(length(c.levels) == 0)
  {
    return("Dat is empty!\n");
  }
  dat$CALen <- dat$CEnd - dat$CStart;
  dat$RALen <- dat$REnd - dat$RStart;
  output <- c();
  ref.levels <- unique(dat$RId);
  ref.lens <- c();
  for(i in 1:length(ref.levels))
  {
    ref.lens <- rbind(ref.lens, c(ref.levels[i], dat[which(dat$RId == ref.levels[i]),"RLen"][1]));
  }
  ref.lens <- as.data.frame(ref.lens);
  colnames(ref.lens) <- c("RId", "Len");
  ref.lens[,2] <- as.numeric(as.character(ref.lens[,2]));
  
  for(i in 1:length(c.levels))
  {
    c.alns <- dat[which(dat$CId == c.levels[i]), ];
    c.len <- c.alns$CLen[1];
    r.levels <- unique(c.alns$RId);
    c.output <- list();
    for(j in 1:length(r.levels))
    {
      temp.dat <- subset(c.alns, c.alns$RId == r.levels[j]);
      while(TRUE)
      {
        outcome <- inner.iteration.function(temp.dat, ol = ol, gap = gap);
        c.output[[r.levels[j]]] <- rbind(c.output[[r.levels[j]]], outcome[["ALNS"]]);
        if(outcome[["DONE"]])
        {
          break;
        }
        temp.dat <- outcome[["REMAINS"]];
      }
    }
    maxlen <- 0;
    max.loc <- NA;
    rid <- NA;
    format.output <- c();
    for(j in 1:length(c.output))
    {
      if(j == 1)
      {
        maxlen <- sum(c.output[[j]][,2] - c.output[[j]][,1]);
        max.loc <- c.output[[j]];
        rid <- names(c.output)[j];
      } else
      {
        if(sum(c.output[[j]][,2] - c.output[[j]][,1]) > maxlen)
        {
          maxlen <- sum(c.output[[j]][,2] - c.output[[j]][,1]);
          max.loc <- c.output[[j]];
          rid <- names(c.output)[j];
        }
      }
    }
    for(j in 1:nrow(max.loc))
    {
      if((max.loc[j,2] - max.loc[j,1]) / c.len > 0.1)
      {
        format.output <- rbind(format.output, c(c.levels[i],  c.len, max.loc[j,1], max.loc[j,2], 
                                                max.loc[j,2] - max.loc[j,1], rid, ref.lens[ref.lens$RId == rid, "Len"],
                                                max.loc[j,3], max.loc[j,4],
                                                max.loc[j,4] - max.loc[j,3]));
      }
    }
    output <- rbind(output, format.output); 
  }
  output <- as.data.frame(output);
  colnames(output) <- c("CId", "CLen", "CStart", "CEnd", "CALen", "RId", "RLen", "RStart", "REnd", "RALen");
  for(i in c(2:5, 7:10))
  {
    output[,i] <- as.numeric(as.character(output[,i]));
  }
  output;
}

inner.iteration.function <- function(temp.dat, ol = 0.1, gap = 1000)
{
  temp.dat <- temp.dat[with(temp.dat, order(RStart)),];
  nrecords <- nrow(temp.dat);
  maxRARecordIndex <- which(temp.dat$RALen == max(temp.dat$RALen))[1];
  maxRARecord <- temp.dat[maxRARecordIndex,];
  orient <- maxRARecord$COrient;
  aln.loc <- c(maxRARecord$CStart, maxRARecord$CEnd, maxRARecord$RStart, maxRARecord$REnd);
  output <- list();
  reverse.set <- c();
  is.reverseable <- FALSE;
  if(maxRARecordIndex >= 2)
  {
    reverse.set <- (maxRARecordIndex - 1) : 1;
    is.reverseable <- TRUE;
  }
  forward.set <- c();
  is.forwardable <- FALSE;
  if(maxRARecordIndex < nrecords)
  {
    forward.set <- (maxRARecordIndex + 1):nrecords;
    is.forwardable <- TRUE;
  }
 
  # reverse search
  visted.reverse.set <- c();
  if(is.reverseable)
  {
    for(j in reverse.set)
    {
      if(temp.dat[j,"COrient"] == orient)
      {
        if(orient == "-")
        {
          rdif <- aln.loc[3] - temp.dat[j,"REnd"];
          rolr <- abs(rdif) / temp.dat[j,"RALen"];
          if(rdif >= 0 && rdif > gap)
          {
            # do nothing now
            # break;
          } else
          {
            if((rdif >= 0 && rdif <= gap) || (rdif < 0 && rolr <= ol))
            {
              cdif <- temp.dat[j,"CStart"] - aln.loc[2];
              colr <- abs(cdif) / temp.dat[j,"CALen"];
              if((cdif >= 0 && cdif <= gap) || (cdif < 0 && colr <= ol))
              {
                aln.loc[2] <- temp.dat[j, "CEnd"];
                aln.loc[3] <- temp.dat[j, "RStart"];
              }
            }
            visted.reverse.set <- c(visted.reverse.set, j);
          }
        } else
        {
          rdif <- aln.loc[3] - temp.dat[j,"REnd"];
          rolr <- abs(rdif) / temp.dat[j,"RALen"];
          if(rdif >= 0 && rdif > gap)
          {
            # do nothing now;
            # break;
          } else
          {
            if((rdif >= 0 && rdif <= gap) || (rdif < 0 && rolr <= ol))
            {
              cdif <- aln.loc[1] - temp.dat[j,"CEnd"];
              colr <- abs(cdif) / temp.dat[j,"CALen"];
              if((cdif >= 0 && cdif <= gap) || (cdif < 0 && colr <= ol))
              {
                aln.loc[1] <- temp.dat[j, "CStart"];
                aln.loc[3] <- temp.dat[j, "RStart"];
              }
            }
            visted.reverse.set <- c(visted.reverse.set, j);
          }
        }
      } else
      {
        visted.reverse.set <- c(visted.reverse.set, j);
      }
    }
  }
  # forward search
  visted.forward.set <- c();
  if(is.forwardable)
  {
    for(j in forward.set)
    {
      if(temp.dat[j,"COrient"] == orient)
      {
        if(orient == "-")
        {
          rdif <- temp.dat[j,"RStart"] - aln.loc[4];
          rolr <- abs(rdif) / temp.dat[j,"RALen"];
          if(rdif > 0 && rdif > gap)
          {
            # do nothing now;
            # break;
          } else
          {
            if((rdif >= 0 && rdif <= gap) || (rdif < 0 && rolr <= ol))
            {
              cdif <- aln.loc[1] -  temp.dat[j,"CEnd"];
              colr <- abs(cdif) / temp.dat[j,"CALen"];
              if((cdif >= 0 && cdif <= gap) || (cdif < 0 && colr <= ol))
              {
                aln.loc[1] <- temp.dat[j, "CStart"];
                aln.loc[4] <- temp.dat[j, "REnd"];
              }
            }
            visted.forward.set <- c(visted.forward.set, j);
          }
        } else
        {
          rdif <- temp.dat[j,"RStart"] - aln.loc[4];
          rolr <- abs(rdif) / temp.dat[j,"RALen"];
          if(rdif > 0 && rdif > gap)
          {
            # do nothing now;
            # break;
          } else
          {
            if((rdif >= 0 && rdif <= gap) ||(rdif < 0 && rolr <= ol))
            {
              cdif <- temp.dat[j,"CStart"] - aln.loc[2];
              colr <- abs(cdif) / temp.dat[j,"CALen"];
              if((cdif >= 0 && cdif <= gap )|| (cdif < 0 && colr <= ol))
              {
                aln.loc[2] <- temp.dat[j, "CEnd"];
                aln.loc[4] <- temp.dat[j, "REnd"];
              }
            }
            visted.forward.set <- c(visted.forward.set, j);
          }
        }
      } else
      {
        visted.forward.set <- c(visted.forward.set, j);
      }
    }
  }
  output[["ALNS"]] <- aln.loc;
  output[["REMAINS"]] <- temp.dat[-c(visted.reverse.set, maxRARecordIndex, visted.forward.set),];
  if(nrow(output[["REMAINS"]]) != 0)
  {
    output[["DONE"]] <- FALSE;
  } else
  {
    output[["DONE"]] <- TRUE;
  }
  output[["VISTED"]] <- c(visted.reverse.set, maxRARecordIndex, visted.forward.set);
  output;
}

# plot function
rlen.maxs <- NULL;
ref.levels <- NULL;
for(i in 1:length(dat))
{
  rlen.maxs <- c(rlen.maxs, max(dat[[i]]$RLen));
  ref.levels <- c(ref.levels, unique(dat[[i]]$RId));
}
x.max <- max(rlen.maxs);
ref.levels <- sort(unique(ref.levels));
png(filename = out, height = ysize, width = xsize);
plot(0, type="n", xlab="", ylab="", xlim=c(0, x.max), 
     ylim=c(0, round(x.max) * (1 + length(ref.levels))), axes = "F");
y = 0;
count_index = 0;
step = round(x.max * 0.2);
mark.range = round(x.max * 0.05);
for(i in 1:length(ref.levels))
{
  count_index <- count_index + 1;
  # ref line;
  ch1.len <- dat[[1]][dat[[1]]$RId == ref.levels[i],"RLen"][1];
  y <- count_index * step;
  segments(x0 = 0, x1 = ch1.len, y0 = y, y1 = y, lwd = 3);
  text(x= ch1.len/2, y = y - 4 * mark.range , labels = ref.levels[i], cex = cex);
  marks.loc <- c();
 # if(ref.levels[i] == "ecoli")
  #{
  #  marks.loc <- seq(from = 0, to = ch1.len, length.out = 5);
  #} else{
  #  marks.loc <- seq(from = 0, to = ch1.len, length.out = 10);
  #}
  marks.loc <- seq(from = 0, to = ch1.len, length.out = breakmarks);
  for(j in 1:length(marks.loc))
  {
    segments(x0 = marks.loc[j], x1 = marks.loc[j], y0 = y - mark.range, y1 = y + mark.range, lwd = 2);
    text(x = marks.loc[j], y = y - 2 * mark.range, labels=paste(round(marks.loc[j]/1000000, digits = 1), "Mb"), cex = cex);
  }
  for(j in 1:length(dat))
  {
    w2r <- dat[[j]];
    w2r.format <- mergeAlignments(w2r, gap = 100000);
    w2r.format.ch1 <- w2r.format[which(w2r.format$RId == ref.levels[i]),];
    w2r.format.ch1 <- w2r.format.ch1[with(w2r.format.ch1, order(RStart)),];
    ch1.len <- w2r.format.ch1$RLen[1];
    
    count_index <- count_index + 1;
    # wtdbg line
    y <- count_index * step;
    segments(x0 = 0, x1 = ch1.len, y0 = y, y1 = y, lwd = 3, col = "gray");
    bins <- findBinCov(w2r.format.ch1);
    for(k in 1:nrow(bins))
    {
      if(bins[k,3] == 0)
      {
        rect(bins[k,1], y - mark.range, bins[k,2], y , col = "blue", border = "blue");
      } else if(bins[k,3] > 1)
      {
        rect(bins[k,1],  y , bins[k,2], y + mark.range, col = "red", border = "red");
      }
    }
    text(x = 0, y = y + 1.5*mark.range, labels = labels[j], adj = c(0, 0), cex = cex);
    
  }
  count_index <- count_index + 1;
}
legend(legendposition, legend =c(">1", "0"), fill = c("red", "blue"), title = "Coverages", cex = 2);
dev.off();

