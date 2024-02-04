library(Rsamtools)
library(DuffyTools)


if (!is.list(targets)) {
  if (file.exists(targets))
    targets <- readTargets(targets)
  else
    stopwrap("targets must be either the result of readTargets ",
             "function or a valid targets file!")
}

# The BAM files
bams <- unlist(targets$files,use.names=FALSE)
nams <- unlist(targets$samples,use.names=FALSE)


getbasequal <- function(bam){
  
  #scan BAM file
  bamfile <- BamFile(bam)
  aln <- scanBam(bamfile)
  
  #extract quality scores for th first 50 bases/read
  bqual <- aln[[1]]$qual
  bqual50 <- substr(bqual, 1 , 50)
  
  #randomly sampled 100.000 reads from file
  set.seed(10)
  bqualt <- sample (bqual50, size = 100000)
  
  #convert score
  bqualt_num <- phredScoreStringToInt(bqualt, scoreType = "Phred33")
  
  print("file complete")
  return(bqualt_num)
  
}


qual_list <- lapply(bam,getbasequal)

#write.table(bqualt_num, "/media/galadriel/fleming/semeli_spanou_work/sra_data/testsamplequal.txt")