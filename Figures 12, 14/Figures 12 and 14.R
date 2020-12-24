##################
# Figure 18
#################
library(tidyverse)
setwd("~/Desktop/apmth221-project-copy")
load("rdas/cv-average.rda")
load("rdas/operating-characteristics.rda")
dslabs::ds_theme_set(new = "theme_minimal")

# gabe added for derivative
A <- matrix(NA, ncol = 999, nrow = 15)
for ( i in 1:15){
    A[i,] <- diff(as.matrix(colMeans(full[full$Electrode == i,(ncol(full) - 999):ncol(full)])))
}

cv.average.gather2 <- matrix(NA, ncol = 4, nrow = 15 * 999)
diff.vec <- as.vector(t(A))
cv.num <- c()
for (i in 1:15){
    cv.num <- c(cv.num, paste0("CV.",rep(i, 999)))
}
cv.average.gather2[,1] <- cv.num
cv.average.gather2[,2] <- diff.vec
cv.average.gather2[,3] <- 1:999
cv.average.gather2 <- tibble(cv.average.gather2)
#################

# -- Viz of similarity metric
f1 <- cv.average.gather %>%
    filter(CV %in% c("CV.1", "CV.11")) %>%
    ggplot(aes(Time.Step, Voltage, color = CV)) +
    geom_line(size = 1) +
    geom_curve(aes(x = 100, y = 1150, xend = 105, yend = 740), lty=2, curvature = 0, color = "black") +
    geom_curve(aes(x = 314, y = 700, xend = 316, yend = 500), lty=2, curvature = 0, color = "black") +
    geom_curve(aes(x = 424, y = 1000, xend = 426, yend = 640), lty=2, curvature = 0, color = "black") +
    geom_curve(aes(x = 844, y = -800, xend = 846, yend = -1230), lty=2, curvature = 0, color = "black") +
    ylab("Current (nA)") +
    xlab("Voltage Potential Index") +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#525252", "#E69F00")) + 
    theme(axis.text.x=element_text(face="bold", color="black", size = rel(2)),
          axis.text.y=element_text(face="bold",color="black", size = rel(2)),
          strip.text = element_text(size = 10,face="bold",color="black"), 
          axis.title.x = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = -2),
          axis.title.y = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = 1.1),
          legend.position = "none") +
    theme( plot.margin = margin(r = 20, b = 20, l = 15, t = 5) ) +  # add margins

# -- Viz of similarity metric 2
f2 <- cv.average.gather %>%
    filter(CV %in% c("CV.1", "CV.11")) %>%
    ggplot(aes(Time.Step, Voltage, color = CV)) +
    geom_line(size = 1) +
    
    geom_point(aes(x = 100, y = 740), color = "black", size = 4, pch = 1) +
    geom_point(aes(x = 314, y = 480), color = "black", size = 4, pch = 1) +
    geom_point(aes(x = 426, y = 660), color = "black", size = 4, pch = 1) +
    geom_point(aes(x = 850, y = -790),color = "black", size = 4, pch = 1) +
    
    geom_point(aes(x = 102, y = 1150), color = "black", size = 4, pch = 1) +
    geom_point(aes(x = 320, y = 720), color = "black", size = 4, pch = 1) +
    geom_point(aes(x = 426, y = 1000), color = "black", size = 4, pch = 1) +
    geom_point(aes(x = 850, y = -1220),color = "black", size = 4, pch = 1) +
    
    ylab("Current (nA)") +
    xlab("Voltage Potential Index") +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#525252", "#E69F00")) + 
    theme(axis.text.x=element_text(face="bold", color="black", size = rel(2)),
          axis.text.y=element_text(face="bold",color="black", size = rel(2)),
          strip.text = element_text(size = 10,face="bold",color="black"), 
          axis.title.x = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = -2),
          axis.title.y = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = 1.1),
          legend.position = "none") +
    theme( plot.margin = margin(r = 20, b = 20, l = 15, t = 5) ) +  # add margins

# save figures
ggsave("similarity-metric-2.png", 
       plot = f1,
       width = 6, height = 4.25)
ggsave("similarity-metric.png", 
       plot = f2,
       width = 6, height = 4.25)

###################
# Figure 16
###################
setwd("~/Desktop/Research")
full <- read.csv("combined_subSamp")

library(tidyverse)
# add working directory here
load("rdas/cv-average.rda")
load("rdas/operating-characteristics.rda")
dslabs::ds_theme_set(new = "theme_minimal")

# gabe added for derivative
A2 <- matrix(NA, ncol = 999, nrow = 15)
A <- matrix(NA, ncol = 1000, nrow = 15)
for ( i in 1:15){
    A[i,] <- (as.matrix(colMeans(full[full$Electrode == i,(ncol(full) - 999):ncol(full)])))
    A2[i,] <- diff(as.matrix(colMeans(full[full$Electrode == i,(ncol(full) - 999):ncol(full)])))
    
}

# derivative
cv.average.gather2 <- data.frame(matrix(NA, ncol = 3, nrow = 15 * 999))
diff.vec <- as.vector(t(A2))
cv.num <- c()
for (i in 1:15){
    cv.num <- c(cv.num, paste0("CV.",rep(i, 999)))
}
cv.average.gather2[,1] <- cv.num
cv.average.gather2[,2] <- as.numeric(diff.vec)
cv.average.gather2[,3] <- as.numeric(1:999)
colnames(cv.average.gather2) <- c("CV", "Voltage", "Time.Step")
cv.average.gather2 <- tibble(cv.average.gather2)


# derivative
cv.average.gather2 <- data.frame(matrix(NA, ncol = 3, nrow = 15 * 999))
diff.vec <- as.vector(t(A2))
cv.num <- c()
for (i in 1:15){
    cv.num <- c(cv.num, paste0("CV.",rep(i, 999)))
}
cv.average.gather2[,1] <- cv.num
cv.average.gather2[,2] <- as.numeric(diff.vec)
cv.average.gather2[,3] <- as.numeric(1:999)
colnames(cv.average.gather2) <- c("CV", "Voltage", "Time.Step")
cv.average.gather2 <- tibble(cv.average.gather2)
#################

# plot figure derivative
f2 <- cv.average.gather2 %>%
    filter(CV %in% c("CV.1", "CV.5", "CV.10", "CV.15")) %>%
    ggplot(aes(Time.Step, Voltage, color = CV)) +
    geom_line(size = 1) +
    ylab("Current (nA)") +
    xlab("Voltage Potential Index") +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#525252", "#E69F00")) + 
    theme(axis.text.x=element_text(face="bold", color="black", size = rel(2)),
          axis.text.y=element_text(face="bold",color="black", size = rel(2)),
          strip.text = element_text(size = 10,face="bold",color="black"), 
          axis.title.x = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = -2),
          axis.title.y = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = 1.1),
          legend.position = "none") +
    theme( plot.margin = margin(r = 20, b = 20, l = 15, t = 5) )  # add margins

# raw

library(tidyverse)
setwd("~/Desktop/apmth221-project-copy")
load("rdas/cv-average.rda")
load("rdas/operating-characteristics.rda")
dslabs::ds_theme_set(new = "theme_minimal")

# -- Viz of sample of raw electrodes 
f1 <- cv.average.gather %>%
    filter(CV %in% c("CV.1", "CV.5", "CV.10", "CV.15")) %>%
    ggplot(aes(Time.Step, Voltage, color = CV)) +
    geom_line(size = 1) +
    ylab("Current") +
    xlab("Voltage Potential Index") +
    scale_color_manual(values = c("#ca0020", "#0868ac", "#525252", "#E69F00")) + 
    theme(axis.text.x=element_text(face="bold", color="black", size = rel(2)),
          axis.text.y=element_text(face="bold",color="black", size = rel(2)),
          strip.text = element_text(size = 10,face="bold",color="black"), 
          axis.title.x = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = -2),
          axis.title.y = element_text(face="bold", color="black", 
                                      size = rel(2),
                                      vjust = 1.1),
          legend.position = "none") +
    theme( plot.margin = margin(r = 20, b = 20, l = 15, t = 5) ) +  # add margins
annotate("text", x=700, y=1000,
         label='bold("Electrode 1")',
         size=6, color="#ca0020", parse=TRUE) +
annotate("text", x=700, y=800,
         label='bold("Electrode 5")',
         size=6, color="#0868ac", parse=TRUE) +
annotate("text", x=716, y=600,
         label='bold("Electrode 10")',
         size=6, color="#525252", parse=TRUE) +
annotate("text", x=716, y=400,
         label='bold("Electrode 15")',
         size=6, color="#E69F00", parse=TRUE) +
  ylim(-1300, 1300)



# save figures
ggsave("sample-electrodes.png", 
       plot = f1,
       width = 6, height = 4.25)
ggsave("sample-electrodes-diff.png", 
       plot = f2,
       width = 6, height = 4.25)

####################################################
