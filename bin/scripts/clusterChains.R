## required library
library(genomeIntervals)

myClusterChains <- function(in.filename, out.filename=""){

  ## read data
  chn <- read.table(in.filename, sep="\t", stringsAsFactors=FALSE)

  ## create intervals
  strand <- rep("+", nrow(chn))
  idx <- chn[,9] > chn[,10]
  st <- chn[,9]
  end <- chn[,10]

  if (any(idx == TRUE)){
    strand[idx] <- "-"                                                         
    st[idx] <- chn[idx,10]
    end[idx] <- chn[idx,9]
  }
  if (any(idx) && any(!idx)){
    class <- "Genome_intervals_stranded"
  }
  else {
    class <- "Genome_intervals"
  } 

  chn.ivals <- new(class, as.matrix(data.frame(st,end)),
                   closed = matrix(rep(c(TRUE, TRUE), each=nrow(chn)), ncol=2),
                   annotation = data.frame(seq_name = chn[,2], strand = strand, inter_base = FALSE,
                     score = chn[,12]))
  
  ## cluster chains
  chn.union <- close_intervals(interval_union(chn.ivals))
  chn.ovl <- interval_overlap(chn.union, chn.ivals)

  ## extract best-scoring chain within each cluster
  chn.idx <- unlist(lapply(chn.ovl, function(x){intersect(x, which(chn[,12] == max(chn[x,12])))[1]}))
  chn.idx <- unique(chn.idx)

  ## write best to file
  if (out.filename != ""){
    write.table(chn[chn.idx,], file=out.filename,
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  }
  return(chn[chn.idx,])
}



## COMMAND LINE OPTIONS
args <- commandArgs(TRUE)

if (length(args) > 0){
  if (length(args) != 2){
    stop("usage: Rscript $0 in.filename out.filename")
  }
  else {
    tmp <- myClusterChains(as.character(args[1]), as.character(args[2]))
  }  
}

