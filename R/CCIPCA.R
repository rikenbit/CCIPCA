CCIPCA <- function(data=NA, runmode="OnMemory", filelist=NA, dim, param){
    # OnMemory Mode
    if(runmode == "OnMemory"){
        if(!is.matrix(data)){
            stop("data is assumed to be matrix!")
        }
        # v : Initial Eigen vectors
        v <- t(data[1:dim, ])
        colnames(v) <- paste0("PA", 1:dim)
        # u : Initial observed data fro calculation of eigen vectors
        u <- v
        u[, ] <- NA
        u <- cbind(u, NA)
        ################### CCIPCA ###################
        pb <- txtProgressBar(min = 1, max = nrow(data), style = 3)
        for(n in 1:nrow(data)){
            # Evaluate Progress Bar
            setTxtProgressBar(pb, n)
            u[, 1] <- t(data[n, ])
            for(i in 1:min(dim, n)){
                if(i == n){
                    v[, i] <- u[, i] 
                }else{
                    w1 <- (n-1-param)/n
                    w2 <- (1+param)/n
                    # Update of the i-th vector, which has same direction of i-th eigen vector
                    v[, i] <- w1*v[, i, drop=F] + w2*u[, i, drop=F]%*%t(u[, i, drop=F])%*%(v[, i, drop=F]/as.numeric(sqrt(crossprod(v[, i]))))
                    # Update of the (i+1)-th data for next (i+1)-th eigen vector 
                    u[, (i+1)] <- u[, i, drop=F] - as.numeric(t(u[, i, drop=F])%*%(v[, i, drop=F]/as.numeric(sqrt(crossprod(v[, i])))))*(v[, i, drop=F]/as.numeric(sqrt(crossprod(v[, i]))))
                }
            }
        }
    }
    # Fileloading mode, when using huge matrix
    else if(runmode == "FileLoading"){
        if(length(filelist) < dim){
            stop("'filelist' must be longer than at least number of principle components")
        }
        # v : Initial Eigen vector
        # u : Initial Observed data for calculation of Eigen vectors
        for(i in 1:dim){
            if(i == 1){
                pre_data <- read.big.matrix(filelist[i], type="double")
                v <- big.matrix(nrow=length(pre_data), ncol=dim)
                u <- big.matrix(nrow=length(pre_data), ncol=(dim+1))
                v[,1] <- pre_data[1,]
                rm(pre_data)
            }else{
                v[, i] <- read.big.matrix(filelist[i], type="double")[,]
            }
        }
        ################### CCIPCA ###################
        # Prepare Progress Bar
        cat('    "Patience is bitter, but its fruit is sweet."\n'); flush.console()
        pb <- txtProgressBar(min = 1, max = length(filelist), style = 3)
        for(n in 1:length(filelist)){
            # Evaluate Progress Bar
            setTxtProgressBar(pb, n)
            u[, 1] <- read.big.matrix(filelist[n], type="double")[, ]
            for(i in 1:min(dim, n)){
                if(i == n){
                    v[, i] <- u[, i]
                }else{
                    w1 <- (n-1-param)/n
                    w2 <- (1+param)/n
                    # Update of the i-th vector, which has same direction of i-th eigen vector
                    v[, i] <- w1*v[, i, drop=F] + w2*u[, i, drop=F]%*%(t(u[, i, drop=F])%*%(v[, i, drop=F]/as.numeric(sqrt(crossprod(v[, i])))))
                    # Update of the (i+1)-th data for next (i+1)-th eigen vector 
                    u[, (i+1)] <- u[, i, drop=F] - as.numeric(t(u[, i, drop=F])%*%(v[, i, drop=F]/as.numeric(sqrt(crossprod(v[, i])))))*(v[, i, drop=F]/as.numeric(sqrt(crossprod(v[, i]))))
                }
            }
        }
    }else{
        stop("'runmode' must be specified, 'OnMemory' or 'FileLoading'")
    }
    # Output eigen values and eigen vectors
    values <- rep(NA, length=dim)
    vectors <- matrix(NA, nrow=nrow(v), ncol=dim)
    for(i in 1:ncol(v)){
        values[i] <- as.numeric(sqrt(crossprod(v[, i])))
        vectors[, i] <- v[, i, drop=F]/as.numeric(sqrt(crossprod(v[, i])))
    }
    return(list(values=values, vectors=vectors))
}