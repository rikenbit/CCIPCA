\name{CCIPCA}
\alias{CCIPCA}
\title{A function to perform CCIPCA}
\description{
Candid Covariance-Free Incremental PCA (CCIPCA) without estimation of covariance matrix}
\usage{
CCIPCA(data = NA, runmode = "OnMemory", filelist = NA, dim, param)
}
\arguments{
\item{data}{data matrix.}
\item{runmode}{"OnMemory" or "FileLoading" is available.}
\item{filelist}{Local path for specifying data file, when runmode is "FileLoading".}
\item{dim}{Dimention to be compressed (No. of principal component)}
\item{param}{Amnesia parameter. 2 to 4 is recommended.}
}
\details{Unlike with the standard PCA (principal component analysis) functions such as \code{prcomp}, \code{eigen}, or \code{svd}, this function loads the data from each row of the data matrix and estimates eigen vectors and eigen values with incremental calculation. Standard PCA algorithms become very slow, when the variables we want to analyze at once is huge (e.g., 10^4 to 10^6), caused by expanded covariance matrix. CCIPCA estimates the eigen vectors and eigen values directlly from data, without calculation of covariance matrix.}
\value{
A list with fields
\item{values}{Eigen values}
\item{vectors}{Eigen vectors}
}
\references{
Weng J. et al. (2003). Candid covariance-free incremental principal component analysis.  \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}
}
\author{Koki Tsuyuzaki}

\examples{
# Standart PCA
data <- t(USArrests)
pca1 <- svd(var(data))
pca2 <- eigen(var(data))
pca3 <- prcomp(data)

# Onmemory mode
pca4 <- CCIPCA(data=data, dim=3, param=2)

# FileLoading mode
\dontrun{
opendir <- function(dir = getwd()){
  if(.Platform['OS.type'] == "windows"){
    shell.exec(dir)
  }else{
    system(paste(Sys.getenv("R_BROWSER"), dir))
  }
}}
\dontrun{opendir(tempdir())}
\dontrun{filelist <- paste0(tempdir(), "/", rownames(data))}
\dontrun{
for(i in 1:nrow(data)){
  write.table(data[i,], file=filelist[i], sep=",", row.names=F, col.names=F)
}}
\dontrun{pca5 <- CCIPCA(runmode="FileLoading", filelist=filelist, dim=3, param=2)}
\dontrun{library("rgl")}
\dontrun{plot3d(pca1$u[,1], pca1$u[,2], pca1$u[,3])}
\dontrun{plot3d(pca2$vectors[,1], pca2$vectors[,2], pca2$vectors[,3])}
\dontrun{plot3d(pca3$rotation[,1], pca3$rotation[,2], pca3$rotation[,3])}
\dontrun{plot3d(pca4$vectors[,1], pca4$vectors[,2], pca4$vectors[,3])}
\dontrun{plot3d(pca5$vectors[,1], pca5$vectors[,2], pca5$vectors[,3])}
}







