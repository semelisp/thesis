library(Rsamtools)
library(DuffyTools)

# The BAM files
targets <- "yourbamlistfile.txt"
bam <- unlist(targets,use.names=FALSE)

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

write.table(bqualt_num, "sample_qual.txt")
