args <- commandArgs(trailingOnly = TRUE) #Gets the arguments put in at the command line. Need to input output_file=, numeratorR=, and permutations= 

#Formats the arguments from the command line properly.
if(length(args>0)){
  arg.matrix <- do.call("rbind", strsplit(args, "="))
  options(warn=-1); arg.char <- which(is.na(as.numeric(arg.matrix[,2]))); options(warn=0)
  if(length(arg.char>0)) arg.matrix[arg.char, 2] <- paste("'", arg.matrix[arg.char, 2], "'", sep="")
  eval(parse(text=apply(arg.matrix, 1, paste, collapse="=")))
  }

#Set default values for the arguments if they're not passed to the function in the command line.
if(!exists("output_file")) output_file <- "my_out"
if(!exists("numeratorR")) numeratorR <- 5
if(!exists("permutations")) permutations <- 1000

library(MatrixEQTL)

#The first parameter to pass is basically the SNPs' weird coord names from the GTEx data. Easiest to just do rownames(snps) to get this.
SNPpositions <- function(SNPnames){
  snpspos <- data.frame(snpid = character(), chr=integer(), pos=integer(), stringsAsFactors = FALSE) #initializes the blank data frame with the appropriate column names.
  for(SNP in 1:length(SNPnames)){ #Iterates through every SNP passed to the function.
    snpspos[SNP,] <- c(SNPnames[SNP], unlist(strsplit(SNPnames[SNP], "_"))[1:2])
    #The first part gets the rsID, the second part extracts the chromosome and position of the SNP based off the formatting of the GTEx data in the snps sliced data, chromo_position_ref_alt_build.
  }
  snpspos[,2] <- as.numeric(snpspos[,2]) #Converts the second column to numeric type.
  snpspos[,3] <- as.numeric(snpspos[,3]) #Converts the second column to numeric type.
  return(snpspos)
}

genespos <- data.frame(geneid = character(), chr=integer(), left=integer(), right=integer(), stringsAsFactors = FALSE)
genespos[1,] <- c("ENSG00000268903.1", 1, 135141, 135895)
genespos[2,] <- c("ENSG00000231709.1", 1, 521369, 523833)
genespos[3,] <- c("ENSG00000229344.1", 1, 568137, 568818)
genespos[4,] <- c("ENSG00000224956.5", 1, 661611, 663527)
genespos[5,] <- c("ENSG00000228327.2", 1, 700237, 714006)
genespos[6,] <- c("ENSG00000237491.4", 1, 721180, 723052)
genespos[7,] <- c("ENSG00000177757.1", 1, 752751, 755214)
genespos[8,] <- c("ENSG00000240453.1", 1, 745489, 753092)
genespos[9,] <- c("ENSG00000225880.4", 1, 761586, 762902)
genespos[10,] <- c("ENSG00000228794.4", 1, 762988, 794826)
genespos[11,] <- c("ENSG00000223764.2", 1, 852250, 855072)
genespos[,2] <- as.numeric(genespos[,2])
genespos[,3] <- as.numeric(genespos[,3])
genespos[,4] <- as.numeric(genespos[,4])

setwd("/mnt/gluster/home/ittai/inputs/data") #Sets the wd to the appropriate place

gene.wise.p.vals <- function(output_file, numeratorR, permutations, genepos = genespos, SNP_file = "SNP.txt", expression_file = "GE.txt", covariates_file = "Covariates.txt", dataWD = "/mnt/gluster/home/ittai/inputs/data", outputWD = "/mnt/gluster/home/ittai/outputs", useModel = modelLINEAR, output_file_name = tempfile(), pvOutputThreshold = 0, pvOutputThreshold.cis = 1, cisDist=1e5, errorCovariance = numeric()){
  
  ## First, define the working directory where SNP_file, expression_file, and covariate_file can   be found.
  setwd(dataWD)
  
  ## Now, define the SNP, gene expression, and covariate dataframes, to later run both normal     (real data) and null (permuted data) analyses:
  
  #Loads SNP genotype data.
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file);
  
  #Loads gene expression data.
  expr = SlicedData$new();
  expr$fileDelimiter = "\t";      # the TAB character
  expr$fileOmitCharacters = "NA"; # denote missing values;
  expr$fileSkipRows = 1;          # one row of column labels
  expr$fileSkipColumns = 1;       # one column of row labels
  expr$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  expr$LoadFile(expression_file);
  
  #Loads covariates data.
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file)>0) {
    cvrt$LoadFile(covariates_file);
  }
  
  ## Runs a normal matrix_eQTL_main analysis for cis SNP eQTLs on the main dataset, finding the SNPs with the minimum p-value for each gene, in order to store the "real" data analysis.
  real.data <- Matrix_eQTL_main(snps, expr, cvrt, output_file_name = output_file_name,
                                output_file_name.cis = output_file_name, pvOutputThreshold = pvOutputThreshold, pvOutputThreshold.cis = pvOutputThreshold.cis, snpspos = SNPpositions(rownames(snps)), genepos = genepos, cisDist = cisDist, useModel = useModel, errorCovariance = errorCovariance, verbose = FALSE, pvalue.hist = FALSE, min.pv.by.genesnp = TRUE, noFDRsaveMemory = FALSE)
  
  ## Uses the lowest.pvals function on this real data, storing the most significant SNPs per each gene in the new matrix of real.pvals.
  real.pvals <- as.matrix(real.data$cis$min.pv.gene)
  
  ## Initializes a list to store all the null distribution p-values for each gene, in addition to each of their real p-values.
  pvals.list <- list()
  pvals.list$RealValues <- real.pvals #Creates the first object in the pvals list, a matrix of each gene's most relevant (i.e. lowest p-value) SNP, and the associated p-value.
  
  ## Begins iterating through every gene, creating null distributions of the data by using permuted genotypes. Matrix_eQTL_main is then run on each permuted dataset, and the minimum       p-value is stored for the gene in question. This continues until there are 5 instances of the    null data creating more significant p-values for each gene than the real data. The number of     permutations done, as well as all the p-values, are saved for each gene.
  for(gene in rownames(real.pvals)){ #Iterates through all the genes.
    real.pval <- as.numeric(real.pvals[gene,]) #Pulls out the gene's real p-value.
    null.pvals <- vector() #Initializes a blank vector to store null p-values.
    
    ## Keeps permuting new data until "numeratorR" or more p-values from null distributions are more significant than the real p-value. Could take a while.
    while(((sum(real.pval >= null.pvals)) < numeratorR) & length(null.pvals) <permutations){
      snps$ColumnSubsample(sample(seq(1:298))) #Permutes the genotype data.
      
      ## Matrix_eQTL_main call to create the null data.
      null.data <- Matrix_eQTL_main(snps, expr, cvrt, output_file_name = output_file_name,
                                    output_file_name.cis = output_file_name, pvOutputThreshold = pvOutputThreshold, pvOutputThreshold.cis = pvOutputThreshold.cis, snpspos = SNPpositions(rownames(snps)), genepos = genepos, cisDist = cisDist, useModel = useModel, errorCovariance = errorCovariance, verbose = FALSE, pvalue.hist = FALSE, min.pv.by.genesnp = TRUE, noFDRsaveMemory = FALSE)
      
      ## Pulls out the minimal null p-value for the gene, and adds it to the null.pvals vector.
      null.pvals <- c(null.pvals, as.numeric(as.matrix(null.data$cis$min.pv.gene)[gene,]))
      
      ## So I can know about the progress of the function.
      cat(paste("Still working on ", gene, ", currently on row ", length(null.pvals), "\n", sep=""))
    }
    
    ## Once there are at least numeratorR number of null distribution p-values more significant than the real one, assigns the vector of these null p-values to the gene's value in the pvals.list.
    pvals.list[[gene]] <- null.pvals
    print(paste(gene, "is complete."))
  }
  
  ## Once all of the genes have been through the process of getting at least numeratorR p-values from null distributions that are more significant than the p-value from the real data, or have exceeded the maximum number of permutations, calculates the new p-value for each gene and puts it into a new matrix, final_result.
  final_result <- matrix(,nrow = length(real.pvals), ncol = 1)
  rownames(final_result) <- rownames(real.pvals)
  for(gene in rownames(real.pvals)){
    final_result[gene,] <- (((sum(pvals.list$RealValues[gene,] > pvals.list[[gene]]))+1)/(length(pvals.list[[gene]])+1))
  }
  pvals.list$final.result <- final_result
  
  setwd(outputWD)
  saveRDS(pvals.list, file = paste(output_file, ".rds", sep = ""))
}

gene.wise.p.vals(output_file=output_file, numeratorR=numeratorR, permutations=permutations)
