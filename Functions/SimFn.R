
####*************************************************************************************
# "Simulations
####*************************************************************************************

###################################
# generate arbitrary \Sigma
###################################
# Fix the covariance matrix across studies and observations

multiStudySim <- function(sim.num = 130,
                          iter = c(),
                          sampSize = 400,
                          covariate.var = 10, # scales the variance of the MVNormal that generates the true means of the covaraites
                          beta.var = 1, # variance of betas
                          clust.num.mu = 80, # number of clusters of means (if this is set to num.studies, than each true mean of covariate is different)
                          clust.num.beta = NULL, # number of clusters of betas (if this is set to num.studies, than each true mean of covariate is different)
                          num.covariates = 20,
                          zero_covs = 11:20, # indices of covariates which are 0
                          fixed.effects = c(), # indices of fixed effects
                          fixed.effectsX = c(), # indices of "fixed covariates" across studies
                          studyNoise = c(1, 1), # range of noises for the different studies 
                          num.studies = 80,
                          beta.mean.range = 10, # true means of hyperdistribution of beta are drawn from a unif(-beta.mean.range, beta.mean.range)
                          perturb = beta.var * 0.1 / 2, # perturb = 0 means all clusters are identical. otherwise perturnance of betas are elementwise drawn from a unif(-perturb, perturb)
                          SB = 1,
                          params = TRUE,
                          sigmaDiag = FALSE, # if true then the covariance matrix of covariates (X and Z) is diagonal
                          sigmaIdentity = FALSE, # if true then covariance matrix is the identity matrix
                          Xmeans_constant = FALSE, # if TRUE then the means of all covariates in a study are shifted by the same amount
                          Xmean0 = FALSE # IF TRUE then the marginal distribution of the Xs is mean 0
                          ){
    
    library(MASS) # for multivariate normal
    if( length(iter) > 0){
        filename <- paste0("Sim ", sim.num, "_",iter)
    }else{
        filename <- paste0("Sim ", sim.num)
    }
    
    if( is.null(clust.num.beta) ){
        clust.num.beta <- clust.num.mu
    }
    
    n <- num.covariates
    p <- qr.Q(qr(matrix(rnorm(n^2), n)))
    sig.vec <- abs(rnorm(n)) # diagonal elements so variances of the covariates
    Sigma <- crossprod(p, p*(sig.vec))
    
    if(sigmaDiag)    Sigma <- diag( diag( Sigma )  ) # remove off diagonal terms
    
    if(sigmaIdentity)    Sigma <- diag( ncol(Sigma) ) # make identity matrix
    
    #### Varies between Studies
    if( length(beta.var) != num.covariates + 1){
        # if the length of the vector does not match, just use the first element
        message("Length of beta.var not equal to number of model coefficients: only the first element used for all coefficients")
        beta.var <- rep(beta.var[1], num.covariates + 1)
    }
    
    # had mean 5 until 7/28/20
    if(Xmean0){
        mean.vec <- rnorm(num.covariates, 0, sqrt(10)) # the true means of the MVRN used to generate the Betas
    }else{
        mean.vec <- rep(0, num.covariates)
    }
    
    # keep degree of variability across studies constant across covariates
    mu.sigma.vec <- rep(1, num.covariates)#runif(num.covariates, 0, 10) # a vector where each element is the constant by which we multiple 
    
    # number of observations per study
    obs.vec <- vector(length = num.studies)
    
    if( length(sampSize) == 1){
        obs.vec <- rep(sampSize, num.studies) # number of observations per study is constant
    }else if(length(sampSize) !=  num.studies){
        
        message("Sample Sizes provided not equal to total studies: Using first sample size provided for all studies")
        obs.vec <- rep(sampSize[1], num.studies) # number of observations per study is constant
        
    }else if (length(sampSize) ==  num.studies ){
        # if equal then just use that
        obs.vec <- sampSize
        
    }
    
    # otherwise use the vector of sample sizes provided if it is the correct length (i.e., as long as num.studies)
    
    # row corresponds to study and column corresponds to True Beta for corresponding covaraite
    beta.matrix <- matrix(NA, nrow = num.studies, ncol = num.covariates + 1)
    
    # fixed across covariates and used as a dial for between study variability
    
    # fix the covariance matrix and mean of the hyperdistribution for the Betas constant
    p <- qr.Q(qr(matrix(rnorm((clust.num.beta)^2), (clust.num.beta))))
    sig.vec <- abs(rnorm(clust.num.beta)) # diagonal elements so variances of the covariates
    beta.Sigma <- crossprod(p, p*(sig.vec)) #* SB
    diag(beta.Sigma) <- diag(beta.Sigma) * SB # ONLY SCALE DIAGONAL TERMS
    
    # generate different mean vector to generate betas from for MVN
    beta.mean.vec <- runif(num.covariates + 1, -beta.mean.range, beta.mean.range) # a vector where each element is the constant by which we multiple 
    
    
    for (i in 1:(num.covariates + 1)){
        
        # generate matrix of true betas
        beta.matrix[, i] <-  rnorm(clust.num.beta, beta.mean.vec[i], sqrt(beta.var[i]))
        # perturb matrix: add matrix of noise from unif() 
        
    }
    
    elems <- nrow(beta.matrix) * ncol(beta.matrix) # elements in beta.mat
    perturb.mat <- matrix(runif(elems, -perturb, perturb), ncol = ncol(beta.matrix))
    beta.matrix <- beta.matrix + perturb.mat
    
    if (length(zero_covs) > 0){
        # if there are zero_covs then set them to 0
        beta.matrix[, zero_covs + 1 ] <- 0 # + 1 is because of intercept--intercept is never set to 0
    }
    
    if(length(fixed.effects) > 0){
        # fix betas across studies for fixed effects
        # if fixed.effects = 0, that is intercept
        fixed.effects <- fixed.effects + 1
        
        indxS <- sample.int(nrow(beta.matrix), 1)
        beta.matrix[, c(fixed.effects)] <- t( replicate(nrow(beta.matrix), beta.matrix[indxS, c(fixed.effects)]) ) # just arbitrarily choose one of the betas
    }
    
    # setwd("/n/home12/gloewinger/sims_true_betas")
    # write.csv(beta.matrix, paste0(filename, "_true_betas"))
    # write.csv(beta.Sigma, paste0(filename, "Sigma_beta"))
    # write.csv(beta.mean.vec, paste0(filename, "_beta_mean_vec"))
    # setwd("/n/home12/gloewinger")
    
    
    # generate true covariate means across studies within covariate
    mu.matrix <- perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
    # Rows: Study Number
    # Columns: Covaraite
    
    
    
    perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
    # 
    for (i in 1:num.covariates){
        # means used in the multivariate normal
        # true means use the mu covariance matrix multiplied by some constant that differs across covariates
        mu.matrix[, i] <- rnorm(clust.num.mu, mean.vec[i], mu.sigma.vec[i] * covariate.var)
        cov.perturb <- mu.sigma.vec[i] * covariate.var * 0.1 / 2 
        # add noise to all true means. scale elements in the noise matrix according to degree of within-covariate across study variability
        # perturb.mat[,i] <- runif(clust.num.mu, -cov.perturb, cov.perturb)
        perturb.mat[,i] <- runif(num.studies, -cov.perturb, cov.perturb)
    }
    
    if(Xmeans_constant){
        # each covariate in a study is shifted by the same amount
        
        perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
        # 
        d <- vector(length = num.studies)
        d[1:num.studies] <- rnorm(clust.num.mu, 0, mu.sigma.vec[1] * covariate.var)
        
        for (i in 1:num.covariates){
            # means used in the multivariate normal
            # true means use the mu covariance matrix multiplied by some constant that differs across covariates
            mu.matrix[, i] <- mean.vec + d # mu.sigma.vec + d
            #cov.perturb <- mu.sigma.vec[i] * covariate.var * 0.1 / 2 
            # add noise to all true means. scale elements in the noise matrix according to degree of within-covariate across study variability
            # perturb.mat[,i] <- runif(clust.num.mu, -cov.perturb, cov.perturb)
            perturb.mat[,i] <- runif(num.studies, -cov.perturb, cov.perturb)
        }
        
        # for(i in 1:num.studies)       mu.matrix[, i] <- mu.matrix[,i ] + d
        
    }else{
        d <- NA
    }
    
    mu.matrix <- mu.matrix + perturb.mat
    
    if(length(fixed.effectsX) > 0){
        # fix betas across studies for fixed effects
        # if fixed.effects = 1, that is intercept
        
        indxS <- sample.int(nrow(mu.matrix), 1)
        mu.matrix[, c(fixed.effectsX)] <- t( replicate(nrow(mu.matrix), mu.matrix[indxS, c(fixed.effectsX)]) ) # just arbitrarily choose one of the betas
    }
    
    # different variance levels of \epsilon for different studies
    noiseVec <- vector(length = num.studies)
    for(j in 1:num.studies){
        
        noiseVec[j] <- runif(1, studyNoise[1], studyNoise[2])
        
    }
    
    noiseVec <- sqrt( noiseVec ) # make it so it is on standard deviation scale to sample from rnorm() below
    

    
    # Simulate data with above parameters
    
    for(y in 1:num.studies){
        data.sim <- matrix(NA, ncol = num.covariates, nrow = obs.vec[y])
        
        
        ##############################
        # Generate the design matrix
        ##############################
        data.sim <- mvrnorm(obs.vec[y], mu.matrix[y , ], Sigma)
        
        # generate a vector Y and add noise
        Y <- cbind(1, data.sim) %*% beta.matrix[y, ] + rnorm(obs.vec[y], 0, noiseVec[y]) # noise is mean 0 with study specific noise levels
        
        # bind it to data
        data.sim <- cbind(y, Y, data.sim) # first column is study number, then Y outcome, then design matroix
        colnames(data.sim) <- c("Study", "Y", paste0("V_", 1:num.covariates))
        assign( paste0(filename, y), data.sim )
        
    }
    
    for(i in 1:num.studies){
        if (i == 1){
            final.results <- get(paste0(filename, i))
            #file.remove(paste0(filename, i))
            #rm(get(paste0(filename, i)))
        }
        else{
            final.results <- rbind(final.results, get(paste0(filename, i)))
            #file.remove(paste0(filename, i))
        }
    }
    
    colnames(final.results) <- c("Study", "Y", paste0("V_", 1:num.covariates))
    
    # write.csv(final.results, paste0(filename, "_Combined"))
    message(paste0("Data Simulation Complete: ", filename))
    
    if(params){
        return( list(data = final.results, 
                     betaMeans = beta.mean.vec, 
                     betas = beta.matrix, 
                     xMean = mean.vec,
                     Xs = mu.matrix,
                     Sigma = Sigma,
                     d = d))
    }else{
        return( final.results )
    }
    
    
    
    
}



###################
# Binary Outcomes
###################

multiStudySimBin <- function(sim.num = 130,
                          iter = c(),
                          sampSize = 400,
                          covariate.var = 10, # scales the variance of the MVNormal that generates the true means of the covaraites
                          beta.var = 1, # variance of betas
                          clust.num.mu = 80, # number of clusters of means (if this is set to num.studies, than each true mean of covariate is different)
                          clust.num.beta = NULL, # number of clusters of betas (if this is set to num.studies, than each true mean of covariate is different)
                          num.covariates = 20,
                          zero_covs = 11:20, # indices of covariates which are 0
                          num.studies = 80,
                          beta.mean.range = 10, # true means of hyperdistribution of beta are drawn from a unif(-beta.mean.range, beta.mean.range)
                          perturb = beta.var * 0.1 / 2, # perturb = 0 means all clusters are identical. otherwise perturnance of betas are elementwise drawn from a unif(-perturb, perturb)
                          SB = 1){
    
    library(MASS) # for multivariate normal
    if( length(iter) > 0){
        filename <- paste0("Sim ", sim.num, "_",iter)
    }else{
        filename <- paste0("Sim ", sim.num)
    }
    
    if( is.null(clust.num.beta) ){
        clust.num.beta <- clust.num.mu
    }
    
    n <- num.covariates
    p <- qr.Q(qr(matrix(rnorm(n^2), n)))
    sig.vec <- abs(rnorm(n)) # diagonal elements so variances of the covariates
    Sigma <- crossprod(p, p*(sig.vec))
    
    #### Varies between Studies
    
    mean.vec <- rnorm(num.covariates, 0, sqrt(10)) # the true means of the MVRN used to generate the Betas
    
    # keep degree of variability across studies constant across covariates
    mu.sigma.vec <- rep(1, num.covariates)#runif(num.covariates, 0, 10) # a vector where each element is the constant by which we multiple 
    
    # number of observations per study
    obs.vec <- vector(length = num.studies)
    
    for(i in 1:num.studies){
        obs.vec[i] <- sampSize # number of observations per study is constant
    }
    
    # row corresponds to study and column corresponds to True Beta for corresponding covaraite
    beta.matrix <- matrix(NA, nrow = num.studies, ncol = num.covariates + 1)
    
    # fixed across covariates and used as a dial for between study variability
    
    # fix the covariance matrix and mean of the hyperdistribution for the Betas constant
    p <- qr.Q(qr(matrix(rnorm((clust.num.beta)^2), (clust.num.beta))))
    sig.vec <- abs(rnorm(clust.num.beta)) # diagonal elements so variances of the covariates
    beta.Sigma <- crossprod(p, p*(sig.vec)) #* SB
    diag(beta.Sigma) <- diag(beta.Sigma) * SB # ONLY SCALE DIAGONAL TERMS
    
    # generate different mean vector to generate betas from for MVN
    beta.mean.vec <- runif(num.covariates + 1, -beta.mean.range, beta.mean.range) # a vector where each element is the constant by which we multiple 
    
    
    for (i in 1:(num.covariates + 1)){
        
        # generate matrix of true betas
        beta.matrix[, i] <-  rnorm(clust.num.beta, beta.mean.vec[i], sqrt(beta.var))
        # perturb matrix: add matrix of noise from unif() 
        
    }
    
    elems <- nrow(beta.matrix) * ncol(beta.matrix) # elements in beta.mat
    perturb.mat <- matrix(runif(elems, -perturb, perturb), ncol = ncol(beta.matrix))
    beta.matrix <- beta.matrix + perturb.mat
    
    if (length(zero_covs) > 0){
        # if there are zero_covs then set them to 0
        beta.matrix[, zero_covs + 1 ] <- 0 # + 1 is because of intercept--intercept is never set to 0
    }
    
    
    # setwd("/n/home12/gloewinger/sims_true_betas")
    # write.csv(beta.matrix, paste0(filename, "_true_betas"))
    # write.csv(beta.Sigma, paste0(filename, "Sigma_beta"))
    # write.csv(beta.mean.vec, paste0(filename, "_beta_mean_vec"))
    # setwd("/n/home12/gloewinger")
    
    
    # generate true covariate means across studies within covariate
    mu.matrix <- perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
    # Rows: Study Number
    # Columns: Covaraite
    
    
    
    perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
    # 
    for (i in 1:num.covariates){
        # means used in the multivariate normal
        # true means use the mu covariance matrix multiplied by some constant that differs across covariates
        mu.matrix[, i] <- rnorm(clust.num.mu, mean.vec[i], mu.sigma.vec[i] * covariate.var)
        cov.perturb <- mu.sigma.vec[i] * covariate.var * 0.1 / 2 
        # add noise to all true means. scale elements in the noise matrix according to degree of within-covariate across study variability
        # perturb.mat[,i] <- runif(clust.num.mu, -cov.perturb, cov.perturb)
        perturb.mat[,i] <- runif(num.studies, -cov.perturb, cov.perturb)
    }
    
    mu.matrix <- mu.matrix + perturb.mat
    
    
    
    # Simulate data with above parameters
    
    for(y in 1:num.studies){
        data.sim <- matrix(NA, ncol = num.covariates, nrow = obs.vec[y])
        
        
        ##############################
        # Generate the design matrix
        ##############################
        data.sim <- mvrnorm(obs.vec[y], mu.matrix[y , ], Sigma)
        
        # generate a vector Y and add noise
        systematicVec <- cbind(1, data.sim) %*% beta.matrix[y, ] #+ rnorm(obs.vec[y])
        pVec <- as.vector( exp(systematicVec) / (1 + exp(systematicVec)) )
        # pVec <- max(pVec, 0)
        # pVec <- min(pVec, 1)
        nSize <- as.integer(length(pVec))
        tots <- as.integer(1)
        Y <- rbinom(nSize, tots, pVec)
        
        # bind it to data
        data.sim <- cbind(y, Y, data.sim) # first column is study number, then Y outcome, then design matroix
        colnames(data.sim) <- c("Study", "Y", paste0("V_", 1:num.covariates))
        assign( paste0(filename, y), data.sim )
        
    }
    
    for(i in 1:num.studies){
        if (i == 1){
            final.results <- get(paste0(filename, i))
            #file.remove(paste0(filename, i))
            #rm(get(paste0(filename, i)))
        }
        else{
            final.results <- rbind(final.results, get(paste0(filename, i)))
            #file.remove(paste0(filename, i))
        }
    }
    
    colnames(final.results) <- c("Study", "Y", paste0("V_", 1:num.covariates))
    
    # write.csv(final.results, paste0(filename, "_Combined"))
    message(paste0("Data Simulation Complete: ", filename))
    return( final.results )
    
}


###################
# Binary Outcomes
###################

multiStudySimNL <- function(sim.num = 130,
                             iter = c(),
                             sampSize = 400,
                             covariate.var = 10, # scales the variance of the MVNormal that generates the true means of the covaraites
                             beta.var = 1, # variance of betas
                             clust.num.mu = 80, # number of clusters of means (if this is set to num.studies, than each true mean of covariate is different)
                             clust.num.beta = NULL, # number of clusters of betas (if this is set to num.studies, than each true mean of covariate is different)
                             num.covariates = 20,
                             zero_covs = 11:20, # indices of covariates which are 0
                             num.studies = 80,
                             beta.mean.range = 10, # true means of hyperdistribution of beta are drawn from a unif(-beta.mean.range, beta.mean.range)
                             perturb = beta.var * 0.1 / 2, # perturb = 0 means all clusters are identical. otherwise perturnance of betas are elementwise drawn from a unif(-perturb, perturb)
                             SB = 1){
    
    library(MASS) # for multivariate normal
    if( length(iter) > 0){
        filename <- paste0("Sim ", sim.num, "_",iter)
    }else{
        filename <- paste0("Sim ", sim.num)
    }
    
    if( is.null(clust.num.beta) ){
        clust.num.beta <- clust.num.mu
    }
    
    n <- num.covariates
    p <- qr.Q(qr(matrix(rnorm(n^2), n)))
    sig.vec <- abs(rnorm(n)) # diagonal elements so variances of the covariates
    Sigma <- crossprod(p, p*(sig.vec))
    
    #### Varies between Studies
    
    mean.vec <- rnorm(num.covariates, 0, sqrt(10)) # the true means of the MVRN used to generate the Betas
    
    # keep degree of variability across studies constant across covariates
    mu.sigma.vec <- rep(1, num.covariates)#runif(num.covariates, 0, 10) # a vector where each element is the constant by which we multiple 
    
    # number of observations per study
    obs.vec <- vector(length = num.studies)
    
    for(i in 1:num.studies){
        obs.vec[i] <- sampSize # number of observations per study is constant
    }
    
    # row corresponds to study and column corresponds to True Beta for corresponding covaraite
    beta.matrix <- matrix(NA, nrow = num.studies, ncol = num.covariates + 1)
    
    # fixed across covariates and used as a dial for between study variability
    
    # fix the covariance matrix and mean of the hyperdistribution for the Betas constant
    p <- qr.Q(qr(matrix(rnorm((clust.num.beta)^2), (clust.num.beta))))
    sig.vec <- abs(rnorm(clust.num.beta)) # diagonal elements so variances of the covariates
    beta.Sigma <- crossprod(p, p*(sig.vec)) #* SB
    diag(beta.Sigma) <- diag(beta.Sigma) * SB # ONLY SCALE DIAGONAL TERMS
    
    # generate different mean vector to generate betas from for MVN
    beta.mean.vec <- runif(num.covariates + 1, -beta.mean.range, beta.mean.range) # a vector where each element is the constant by which we multiple 
    
    
    for (i in 1:(num.covariates + 1)){
        
        # generate matrix of true betas
        beta.matrix[, i] <-  rnorm(clust.num.beta, beta.mean.vec[i], sqrt(beta.var))
        # perturb matrix: add matrix of noise from unif() 
        
    }
    
    elems <- nrow(beta.matrix) * ncol(beta.matrix) # elements in beta.mat
    perturb.mat <- matrix(runif(elems, -perturb, perturb), ncol = ncol(beta.matrix))
    beta.matrix <- beta.matrix + perturb.mat
    
    if (length(zero_covs) > 0){
        # if there are zero_covs then set them to 0
        beta.matrix[, zero_covs + 1 ] <- 0 # + 1 is because of intercept--intercept is never set to 0
    }
    
    
    # setwd("/n/home12/gloewinger/sims_true_betas")
    # write.csv(beta.matrix, paste0(filename, "_true_betas"))
    # write.csv(beta.Sigma, paste0(filename, "Sigma_beta"))
    # write.csv(beta.mean.vec, paste0(filename, "_beta_mean_vec"))
    # setwd("/n/home12/gloewinger")
    
    
    # generate true covariate means across studies within covariate
    mu.matrix <- perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
    # Rows: Study Number
    # Columns: Covaraite
    
    
    
    perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
    # 
    for (i in 1:num.covariates){
        # means used in the multivariate normal
        # true means use the mu covariance matrix multiplied by some constant that differs across covariates
        mu.matrix[, i] <- rnorm(clust.num.mu, mean.vec[i], mu.sigma.vec[i] * covariate.var)
        cov.perturb <- mu.sigma.vec[i] * covariate.var * 0.1 / 2 
        # add noise to all true means. scale elements in the noise matrix according to degree of within-covariate across study variability
        # perturb.mat[,i] <- runif(clust.num.mu, -cov.perturb, cov.perturb)
        perturb.mat[,i] <- runif(num.studies, -cov.perturb, cov.perturb)
    }
    
    mu.matrix <- mu.matrix + perturb.mat
    
    
    
    # Simulate data with above parameters
    
    for(y in 1:num.studies){
        data.sim <- matrix(NA, ncol = num.covariates, nrow = obs.vec[y])
        
        
        ##############################
        # Generate the design matrix
        ##############################
        data.sim <- mvrnorm(obs.vec[y], mu.matrix[y , ], Sigma)
        
        # generate a vector Y and add noise
        quadTerms <- 1:3 # which covaraites are also added as quadratic terms
        systematicVec <- cbind(1, data.sim) %*% beta.matrix[y, ] + data.sim[,quadTerms]^2 %*% beta.matrix[y, quadTerms + 1] # add one because of intercept
        
        Y <- systematicVec + rnorm(obs.vec[y])
        
        # bind it to data
        data.sim <- cbind(y, Y, data.sim) # first column is study number, then Y outcome, then design matroix
        colnames(data.sim) <- c("Study", "Y", paste0("V_", 1:num.covariates))
        assign( paste0(filename, y), data.sim )
        
    }
    
    for(i in 1:num.studies){
        if (i == 1){
            final.results <- get(paste0(filename, i))
            #file.remove(paste0(filename, i))
            #rm(get(paste0(filename, i)))
        }
        else{
            final.results <- rbind(final.results, get(paste0(filename, i)))
            #file.remove(paste0(filename, i))
        }
    }
    
    colnames(final.results) <- c("Study", "Y", paste0("V_", 1:num.covariates))
    
    # write.csv(final.results, paste0(filename, "_Combined"))
    message(paste0("Data Simulation Complete: ", filename))
    return( final.results )
    
    
    
}






###################
# Binary Outcomes - TEST
###################

multiStudySimBinTest <- function(sim.num = 130,
                             iter = c(),
                             sampSize = 400,
                             covariate.var = 10, # scales the variance of the MVNormal that generates the true means of the covaraites
                             beta.var = 1, # variance of betas
                             clust.num.mu = 80, # number of clusters of means (if this is set to num.studies, than each true mean of covariate is different)
                             clust.num.beta = NULL, # number of clusters of betas (if this is set to num.studies, than each true mean of covariate is different)
                             num.covariates = 20,
                             zero_covs = 11:20, # indices of covariates which are 0
                             num.studies = 80,
                             beta.mean.range = 10, # true means of hyperdistribution of beta are drawn from a unif(-beta.mean.range, beta.mean.range)
                             perturb = beta.var * 0.1 / 2, # perturb = 0 means all clusters are identical. otherwise perturnance of betas are elementwise drawn from a unif(-perturb, perturb)
                             SB = 1){
    
    flag <- TRUE
    
    while(flag){
        # while loop is ensure there are proper amount of outcomes (i.e., Y = 1)
        
        library(MASS) # for multivariate normal
        if( length(iter) > 0){
            filename <- paste0("Sim ", sim.num, "_",iter)
        }else{
            filename <- paste0("Sim ", sim.num)
        }
        
        if( is.null(clust.num.beta) ){
            clust.num.beta <- clust.num.mu
        }
        
        n <- num.covariates
        p <- qr.Q(qr(matrix(rnorm(n^2), n)))
        sig.vec <- abs(rnorm(n)) # diagonal elements so variances of the covariates
        Sigma <- crossprod(p, p*(sig.vec)) #/ 10 #/ 20
        
        #### Varies between Studies
        
        mean.vec <- rnorm(num.covariates, 0, 1 / 2) # the true means of the MVRN used to generate the Betas
        
        # keep degree of variability across studies constant across covariates
        mu.sigma.vec <- rep(1, num.covariates)#runif(num.covariates, 0, 10) # a vector where each element is the constant by which we multiple 
        
        # number of observations per study
        obs.vec <- vector(length = num.studies)
        
        for(i in 1:num.studies){
            obs.vec[i] <- sampSize # number of observations per study is constant
        }
        
        # row corresponds to study and column corresponds to True Beta for corresponding covaraite
        beta.matrix <- matrix(NA, nrow = num.studies, ncol = num.covariates + 1)
        
        # fixed across covariates and used as a dial for between study variability
        
        # fix the covariance matrix and mean of the hyperdistribution for the Betas constant
        p <- qr.Q(qr(matrix(rnorm((clust.num.beta)^2), (clust.num.beta))))
        sig.vec <- abs(rnorm(clust.num.beta)) # diagonal elements so variances of the covariates
        beta.Sigma <- crossprod(p, p*(sig.vec)) #* SB
        diag(beta.Sigma) <- diag(beta.Sigma) * SB # ONLY SCALE DIAGONAL TERMS
        
        # generate different mean vector to generate betas from for MVN
        beta.mean.vec <- runif(num.covariates + 1, -beta.mean.range, beta.mean.range) # a vector where each element is the constant by which we multiple 
        
        
        for (i in 1:(num.covariates + 1)){
            
            # generate matrix of true betas
            beta.matrix[, i] <-  rnorm(clust.num.beta, beta.mean.vec[i], sqrt(beta.var))
            # perturb matrix: add matrix of noise from unif() 
            
        }
        
        elems <- nrow(beta.matrix) * ncol(beta.matrix) # elements in beta.mat
        perturb.mat <- matrix(runif(elems, -perturb, perturb), ncol = ncol(beta.matrix))
        beta.matrix <- beta.matrix + perturb.mat
        
        ######## FIXED
        # beta.matrix[,1] <- 3 # all studies have same intercept
        # beta.matrix[,2] <- sample(c(1,-1), nrow(beta.matrix), replace = TRUE)
        # beta.matrix[,3] <- sample(c(1,-1), nrow(beta.matrix), replace = TRUE)
        ######## FIXED
        
        if (length(zero_covs) > 0){
            #message("Zero Covariates")
            #print("zero_covs + 1")
            # if there are zero_covs then set them to 0
            beta.matrix[, zero_covs + 1 ] <- 0 # + 1 is because of intercept--intercept is never set to 0
        }
        
        #print(beta.matrix)
        # setwd("/n/home12/gloewinger/sims_true_betas")
        # write.csv(beta.matrix, paste0(filename, "_true_betas"))
        # write.csv(beta.Sigma, paste0(filename, "Sigma_beta"))
        # write.csv(beta.mean.vec, paste0(filename, "_beta_mean_vec"))
        # setwd("/n/home12/gloewinger")
        
        
        # generate true covariate means across studies within covariate
        mu.matrix <- perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
        # Rows: Study Number
        # Columns: Covaraite
        
        
        
        perturb.mat <- matrix(NA, nrow = num.studies, ncol = num.covariates)
        # 
        for (i in 1:num.covariates){
            # means used in the multivariate normal
            # true means use the mu covariance matrix multiplied by some constant that differs across covariates
            mu.matrix[, i] <- rnorm(clust.num.mu, mean.vec[i], mu.sigma.vec[i] * covariate.var)
            cov.perturb <- mu.sigma.vec[i] * covariate.var * 0.1 / 2 
            # add noise to all true means. scale elements in the noise matrix according to degree of within-covariate across study variability
            # perturb.mat[,i] <- runif(clust.num.mu, -cov.perturb, cov.perturb)
            perturb.mat[,i] <- runif(num.studies, -cov.perturb, cov.perturb)
        }
        
        mu.matrix <- mu.matrix + perturb.mat
        #print("mu matrix")
        #print( mu.matrix )
        
        out_flag <- TRUE # indicates that no study has too few or too many outcomes (set default)
        
        
        # Simulate data with above parameters
        
        for(y in 1:num.studies){
            data.sim <- matrix(NA, ncol = num.covariates, nrow = obs.vec[y])
            
            
            ##############################
            # Generate the design matrix
            ##############################
            data.sim <- mvrnorm(obs.vec[y], mu.matrix[y , ], Sigma)
            
            # generate a vector Y and add noise
            systematicVec <- cbind(1, scale( data.sim) ) %*% beta.matrix[y, ] #+ rnorm(obs.vec[y])
            pVec <- as.vector(    1 / (1 + exp(-systematicVec) )    )
            
            # pVec <- max(pVec, 0)
            # pVec <- min(pVec, 1)
            nSize <- as.integer(length(pVec))
            tots <- as.integer(1)
            Y <- rbinom(nSize, tots, pVec)
            

            
            # bind it to data
            data.sim <- cbind(y, Y, data.sim) # first column is study number, then Y outcome, then design matroix
            colnames(data.sim) <- c("Study", "Y", paste0("V_", 1:num.covariates))
            assign( paste0(filename, y), data.sim )
            
            print(paste("p: ", sum(Y) / length(Y)))
            
            # test to make sure there are enough outcomes in ever study
            if(
                sum(Y) / length(Y) < 0.1 || 
                    sum(Y) / length(Y) > 0.9 ||
                        anyNA(Y) ){
                
                out_flag <- FALSE
                print(paste0("Repeat: ", sum(data.sim[2, ]) / nrow(data.sim) ) )
                
            }
            
        }
        
        
        
                
        # iterate through studies and make sure there are proper amount of outcomes (Y = 1)
        # for(study_number in 1:num.studies){
        #     rowsStudy <- which(data.sim[, 1] == study_number)
        #     print(length(rowsStudy))
        #     if(
        #         sum(data.sim[2, rowsStudy]) / length(rowsStudy) < 0.1 || sum(data.sim[2, rowsStudy]) / length(rowsStudy) > 0.9){
        #         
        #         out_flag <- FALSE
        #         print("Repeat")
        #         
        #     }
        # }
        # 
        # if no study has too few or too many outcomes set flag to true to break while loop and return data
        if(out_flag){
            # if there are proper amount of outcomes, then break while loop and return data
            flag <- FALSE
        }
        
        
    } # end of whil eloop
        
    for(i in 1:num.studies){
        if (i == 1){
            final.results <- get(paste0(filename, i))
                #file.remove(paste0(filename, i))
                #rm(get(paste0(filename, i)))
        }
        else{
                final.results <- rbind(final.results, get(paste0(filename, i)))
                #file.remove(paste0(filename, i))
        }
    }
        
    colnames(final.results) <- c("Study", "Y", paste0("V_", 1:num.covariates))
    
    # write.csv(final.results, paste0(filename, "_Combined"))
    message(paste0("Data Simulation Complete: ", filename))
    
    return( final.results )
    
}


