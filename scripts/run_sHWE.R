.libPaths(c("/usr/local/lib/R/site-library", .libPaths())) # This is necessary so that singularity doesn't look in bound home directory for libraries.  The pre-pended directory ought to be the one WITHIN the singularity image.

suppressPackageStartupMessages(library(lfa))
suppressPackageStartupMessages(library(optparse))

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-p", "--plink_prefix"),
              help="Full path of plink prefix"),
  make_option(c("-o", "--output"), help="output_file"),
  make_option(c("-d", "--num_latent_factors"), help="Number of latent factors to infer, including the intercept"))

opt <- parse_args(OptionParser(option_list=option_list))

######################### END: Read in arguments #########################

#generates signed integer to genotypes conversion matrix
#columns of matrix coorespond to unsigned integers PLUS ONE
binary.genotype.map <- function(){
    combinations <- as.matrix(expand.grid(0:3,0:3,0:3,0:3))
    snp.map <- matrix(0, 4, 256)
    colnames(combinations) <- NULL

    #again offset by 1
    bitstring <- list()
    bitstring[[1]] <- "00"
    bitstring[[2]] <- "01"
    bitstring[[3]] <- "10"
    bitstring[[4]] <- "11"

    #generate the indices
    indices <- apply(combinations, 1, function(x){
            strtoi(paste(bitstring[[x[1]+1]],
                         bitstring[[x[2]+1]],
                         bitstring[[x[3]+1]],
                         bitstring[[x[4]+1]], sep=""), base=2)})

    indices <- indices+1

    #this order matters...
    combinations[combinations==1] <- NA #PLINK IS BACKWARDS
    combinations[combinations==2] <- 1  #PLINK IS BACKWARDS
    combinations[combinations==0] <- 2  #PLINK IS BACKWARDS
    combinations[combinations==3] <- 0  #PLINK IS BACKWARDS

    snp.map <- apply(combinations, 1, rev)
    snp.map[,indices] <- snp.map

    snp.map
}

read.bed <- function(bed.prefix){
    bed.filename <- paste(bed.prefix, ".bed", sep="")
    bim.filename <- paste(bed.prefix, ".bim", sep="")
    fam.filename <- paste(bed.prefix, ".fam", sep="")

    if(!file.exists(bed.filename))
        stop("need .bed file")
    if(!file.exists(bim.filename))
        stop("need .bim file")
    if(!file.exists(fam.filename))
        stop("need .fam file")

    #figure out number of individuals by counting newlines in fam
    buffer <- read.table(fam.filename, stringsAsFactors=FALSE, colClasses="character", comment.char = "")
    n <- nrow(buffer)

    print(paste("reading in", n, "individuals"))

    #figure out number of SNPs by counting newlines in bim
    buffer <- read.table(bim.filename, stringsAsFactors=FALSE, colClasses="character")
    m <- nrow(buffer)
    rm(buffer)

    print(paste("reading in", m, "snps"))

    #initialize genotype matrix
    X <- matrix(0, m, n)
    snp.map <- binary.genotype.map()

    #open stream
    bed <- file(bed.filename, "rb")

    #check that beginning of bim satisfies magic numbers
    if(readBin(bed, what="integer", n=1, size=1) != 108)
        stop("not valid bed file (magic number fail)")
    if(readBin(bed, what="integer", n=1, size=1) != 27)
        stop("not valid bed file (magic number fail)")
    buffer <- readBin(bed, what="integer", n=1, size=1)
    if(buffer == 0) {
        stop("individual major mode not yet supported")
    } else if(buffer==1) {
        print("snp major mode")
    } else{
        stop("bed mode problem")
    }

    #calculate the blocksize
    numbytes <- ceiling(n/4)


    #read in SNPs!
    for(i in 1:m){
        indices <- readBin(bed, what="int", n=numbytes, size=1,
	  signed=FALSE)+1
        snp.in <- snp.map[,indices]
        X[i,] <- as.vector(snp.in[1:n])

        if(i %% 20000 == 0)
            print(paste("reading snp", i))
    }

    close(bed)
    X
}

######################### Main Body #########################

dat = read.bed(opt$p)
print(dim(dat))
print(opt$num_latent_factors)
print(opt$o)
dat = dat[complete.cases(dat),]
LF = lfa(dat, as.numeric(opt$num_latent_factors))
p = sHWE(dat, LF, 1)
write.table(p, opt$o, quote=F,row.names=F,col.names = F)
