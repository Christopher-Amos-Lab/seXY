### Enter directory containing downloaded source files
source_dir="/home/source_dir/"

### Enter directory containing SNP array genotype data
array_dir="/home/array_dir/"

### Is the genotype data in PLINK *.ped format or GenomeStudio matrix *.txt format?
### Enter "PLINK" or "GS":
format="PLINK"

### Enter filename of X chromosome genotype data
X_raw="Prostate_X_plink.ped"

### Enter filename of Y chromosome genotype data
Y_raw="Prostate_Y_plink.ped"

### Enter desired filename of sex inference output file
output="sex_results.txt"


##################################
##################################
##################################

### 1. Fit logistic regression model to males and females from prostate cancer and ovarian cancer GWAS TRAINING SET

ref <- read.delim(paste(c(source_dir, "XY_data.txt"), collapse=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
fit <- glm(factor(true_sex) ~ XH + YM, family="binomial", data=ref)


### 2. Import TEST SET genotype data

if (format == "GS") {
  X <- read.delim(paste(c(array_dir, X_raw), collapse=""), sep="\t", skip=9, stringsAsFactors=FALSE, header=TRUE)
  X <- as.matrix(X[,(2:ncol(X))])
  
  Y <- read.delim(paste(c(array_dir, Y_raw), collapse=""), sep="\t", skip=9, stringsAsFactors=FALSE, header=TRUE)
  Y <- as.matrix(Y[,(2:ncol(Y))])
} else {
  X <- read.delim(paste(c(array_dir, X_raw), collapse=""), sep=" ", stringsAsFactors=FALSE, header=FALSE)
  X <- as.matrix(X[,(7:ncol(X))])
  
  Y <- read.delim(paste(c(array_dir, Y_raw), collapse=""), sep=" ", stringsAsFactors=FALSE, header=FALSE)
  Y <- as.matrix(Y[,(7:ncol(Y))])
}


### 3. Compute X chromosome heterzygosity for TEST SET data

if (format == "GS") {
  XH <- rep(0, ncol(X))
  scorecard <- matrix(0, nrow=length(XH), ncol=2)
  
  HZ <- function(x) {
    if (nchar(x) == 2 & substr(x, 1, 1) != "-" & substr(x, 2, 2) != "-") {
      if (substr(x, 1, 1) != substr(x, 2, 2)) {
        return(c(1, 0))
      } else {
        return(c(0, 0))
      }
    } else {
      return(c(0, 1))
    }
  }
  
  pb <- txtProgressBar(min = 0, max = length(XH), style = 3)
  for (i in 1:length(XH)) {
    Sys.sleep(0.001)
    setTxtProgressBar(pb, i)
    
    prelim <- lapply(X[,i], HZ)
    scorecard <- do.call(rbind, prelim)
    
    XH[i] <- sum(scorecard[,1])/(nrow(scorecard) - sum(scorecard[,2]))
  }
  close(pb)
  rm(pb)
} else {
  XH <- rep(0, nrow(X))
  
  pb <- txtProgressBar(min = 0, max = length(XH), style = 3)
  for (i in 1:length(XH)) {
    Sys.sleep(0.001)
    setTxtProgressBar(pb, i)
    
    allele1 <- X[i, seq(1, (ncol(X)-1), 2)]
    allele2 <- X[i, seq(2, ncol(X), 2)]
    missing <- union(which(allele1 == "0"), which(allele2 == "0"))
    
    XH[i] <- length(which(allele1 != allele2))/(ncol(X)/2 - length(missing))
  }
  close(pb)
  rm(pb)
}


### 4. Compute Y chromosome missingness for TEST SET data

if (format == "GS") {
  YM <- rep(0, ncol(Y))
  
  pb <- txtProgressBar(min = 0, max = length(YM), style = 3)
  for (i in 1:length(YM)) {
    Sys.sleep(0.001)
    setTxtProgressBar(pb, i)
    
    YM[i] <- length(which(Y[,i] == "--"))/nrow(Y)
  }
  close(pb)
  rm(pb)
} else {
  YM <- rep(0, nrow(Y))
  
  pb <- txtProgressBar(min = 0, max = length(YM), style = 3)
  for (i in 1:length(YM)) {
    allele1 <- Y[i, seq(1, (ncol(Y)-1), 2)]
    allele2 <- Y[i, seq(2, ncol(Y), 2)]
    YM[i] <- length(union(which(allele1 == "0"), which(allele2 == "0")))/(ncol(Y)/2)
  }
  close(pb)
  rm(pb)
}


### 5. Infer sex from XH and YM
### Male = 1, Female = 2

testset <- data.frame(Subject=1:length(XH), XH, YM)
probs <- as.numeric(predict(fit, testset, type="response"))
testset$Predicted <- (probs > 0.5)*1+1


### 6. Print results

write.table(testset, paste(c(array_dir, output), collapse=""), sep="\t", row.name=FALSE, col.names=TRUE, quote=FALSE)
