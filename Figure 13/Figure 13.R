library(ggplot2)
library(latex2exp)
library(grid)
library(gridExtra)
library(RColorBrewer)

setwd("~/Desktop/Research")
full <- read.csv("combined_subSamp")
full <- full[,-c(1,2,4,5)] # remove unnecessary columns


K <- 15 # studies
p <- 1000 # covariates
corr.mat <- matrix(ncol = 1000, nrow = K)
row.names(corr.mat) <- c(paste0("study_", 1:K))
colnames(corr.mat) <- paste0("covariate_", 1:p)


# iterate through studies
for (study in 1:K){
    # iterate through covariates
    for (cov in 1:p){
        corr.mat[study, cov] <- cor( full[full$Electrode == study, cov + 2], full$DA[full$Electrode == study] )
    }
}

# prepare for ggplot2
corr.vec <- as.vector(t(corr.mat)) #vectorize
df <- cbind(corr.vec, rep(1:1000, K), as.factor(rep(1:K, each = p)))
df <- as.data.frame(df)
df[,2] <- as.factor(df[,2])
df[,3] <- as.factor(df[,3])
colnames(df) <- c("Corr", "Cov", "Study")


##################
# plot
##################
    
p1 <- ggplot(data = df[df$Study == c(1,2,3,4),], aes(x=factor(Cov), y=Corr, fill = (Study))) + 
     xlab("Covariate Number") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, 
                                                              color="black", size=12, face="bold")) +
    scale_x_discrete(breaks = c(1,200,400,600,800,1000)) +
    theme(axis.text.x = element_text(face = "bold", 
                                     # angle=45, 
                                     size=rel(2.5),
                                     vjust = 0.6 ),
          axis.text.y = element_text(face = "bold", 
                                     size=rel(2.5)
          ),
          axis.title.x = element_text(face = "bold", 
                                      size=rel(2.25),
                                      vjust = -2
          ),
          axis.title.y = element_text(face = "bold", 
                                      size=rel(2.25),
                                      vjust = 3
          ) ) +
    ylab("Marginal Correlation") +
    annotate("text", x=850, y=1.2,
             label='bold("Electrode 1")',
             size=7, color="#F8766D", parse=TRUE) +
    annotate("text", x=850, y=1.05,
             label='bold("Electrode 2")',
             size=7, color="#7CAE00", parse=TRUE) +
    annotate("text", x=850, y=0.9,
             label='bold("Electrode 3")',
             size=7, color="#00BFC4", parse=TRUE) +
    annotate("text", x=850, y=0.75,
             label='bold("Electrode 4")',
             size=7, color="#C77CFF", parse=TRUE) + 
    geom_line(aes(group = (Study), color = Study), size = 1.5) +
    theme( plot.margin = margin(r = 30, b = 20, l = 15) )  # add margins

#*************************************************
##########################
# derivative
##########################
K <- 15 # studies
p <- 999 # covariates
corr.mat <- matrix(ncol = 999, nrow = K)
row.names(corr.mat) <- c(paste0("study_", 1:K))
colnames(corr.mat) <- paste0("covariate_", 1:p)
full.deriv <- cbind( full[,c(1,2)], t( diff( t( full[,-c(1,2)] ))) ) # derivative of CVs

# iterate through studies
for (study in 1:K){
    # iterate through covariates
    for (cov in 1:p){
        corr.mat[study, cov] <- cor( as.numeric(full.deriv[full.deriv$Electrode == study, cov + 2]), 
                                                as.numeric(full.deriv$DA[full.deriv$Electrode == study] ) )
    }
}


# prepare for ggplot2
corr.vec <- as.vector(t(corr.mat)) #vectorize
df.diff <- cbind(corr.vec, rep(1:p, K), rep(1:K, each = p))
df.diff <- as.data.frame(df.diff)
df.diff[,2] <- as.factor(df.diff[,2])
df.diff[,3] <- as.factor(df.diff[,3])
colnames(df.diff) <- c("Corr", "Cov", "Study")


##################
# Derivative plot
##################

subrows <- c()
for (i in 1:4){
    subrows <- c(subrows, which(df.diff$Study == i))
}

df.diffSub <- df.diff[subrows,]

p2 <- ggplot(data = df.diffSub, aes(x=factor(Cov), y=Corr, fill = (Study))) + 
    xlab("Covariate Number") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, 
                                                              color="black", size=12, face="bold")) +
    scale_x_discrete(breaks = c(1,200,400,600,800, 999)) +
    theme(axis.text.x = element_text(face = "bold", 
                                     # angle=45, 
                                     size=rel(2.5),
                                     vjust = 0.6 ),
          axis.text.y = element_text(face = "bold", 
                                     size=rel(2.5)
          ),
          axis.title.x = element_text(face = "bold", 
                                      size=rel(2.25),
                                      vjust = -2
          ),
          axis.title.y = element_text(face = "bold", 
                                      size=rel(2.25),
                                      vjust = 3
          ) ) +
    ylab("Marginal Correlation") +
    geom_line(aes(group = (Study), color = Study), size = 1.5) +
    theme( plot.margin = margin(r = 30, b = 20, l = 15) )  # add margins
    # add margins

############
# save
############
setwd("~/Desktop/Research Final/Methods Paper/Draft/Figures/fscv_corr") # setwd for saving file
ggsave("fscv_corr.png", 
       plot = p1,
       width = 7, height = 6.5)
ggsave("fscv_corr_deriv.png",
       plot = p2,
       width = 7, height = 6.5)

