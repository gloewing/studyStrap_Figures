library(ggplot2)
library(gridExtra)
library(viridis)
library(knitr)
library(kableExtra)
library(grid)
library(latex2exp)
library(tidyverse)
setwd("~/Desktop/Research/sims5")

sim.name <- c(160:191) # order to correspond with figure # add 138 later
num.sims <- length(sim.name)


sims <- cbind(160:183, expand.grid(c(0.05, 0.25, 1, 3), c(0.05, 5, 20)^2, c(0,4)))
sims2 <- cbind( 184:191, expand.grid(c(0.05,0.25,1,3),c(10,15), c(0)))
colnames(sims) <- c("simNum", "betaVar", "XVar", "Clusts")
colnames(sims2) <- c("simNum", "betaVar", "XVar", "Clusts")
sims <- rbind(sims,sims2)
sim.matrix <- sims

rm <- c(17,19) # remove these indices for visual purposes
bagSize.vec <- c(1:3, seq(4,18, by = 2), seq(20,100, by = 10), 200, 500, 1000)[-rm] # 
labs <- c("OSE", "OSE-Stack", "OSE-CPS", "SSE", "SSE-Stack", "SSE-CPS", "AR", "AR-CPS", "AR-Stack")
bags <- bagSize.vec

#Our transformation function to ensure rounding is same decimal places for all
scaleFUN <- function(x) sprintf("%.2f", x)

for(y in 1:num.sims){
    
    num <- sim.num <- sim.name[y]
    
    clusts <- sim.matrix[y,4]
    sigma <- sim.matrix[y,2]
    sigX <- sim.matrix[y,3]
   
    dat <- read.csv(paste0("LASSO_Sims5_Results_",num)) 
    merged <- dat[,1]
    naive <- mean( dat[,2] /merged) # sse
    merged.mean <- mean(merged)
    dat <- dat[,-1] # removed /merged from here since we do it below for sse.stack and sse.cps
    sse.stack <- mean( dat[,2] / merged )
    sse.cps <- mean( dat[,3] / merged )
    
    res <- colMeans(read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Avg")) / merged ) 
    res.cpsw <- colMeans(read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_CPS")) / merged ) 
    res.stack <- colMeans(read.csv(paste0("LASSO SS BagSize_Test_Sims5_",num,"_Stack")) / merged ) 
    
    # remove 2 columns for visual purposes as these are superfluous
    res <- res[-rm] #[-c(22,23)]
    res.cpsw <- res.cpsw[-rm] #[-c(22,23)]
    res.stack <- res.stack[-rm] #[-c(22,23)]

    bags <- bagSzs <- bagSize.vec
    bags <- c(bags)
    df.cpsw <- data.frame(RMSE = res.cpsw, Bag.Size = bags, Method = rep(2, length(bags)))
    df.stack <- data.frame(RMSE = res.stack, Bag.Size = bags, Method = rep(3, length(res.stack)))
    
    
    bags <- factor(bags, levels = c(bagSzs))
    df <- data.frame(RMSE = res, Bag.Size = bags, Method = 1)
    df <- rbind(df,df.cpsw,df.stack)

    p <- ggplot(data = df, aes(x=factor(Bag.Size), y=log(RMSE), color = factor(Method), group = factor(Method))) + 
        geom_line( size = 2) + 
        theme_classic(base_size = 12) + 
        geom_point( aes(shape = factor(Method) ), size = 3 ) + 
        xlab("Bag Size") + 
        ggtitle(paste0(sim.num)) + 
        # geom_point(colour = "black")+
        scale_y_continuous(labels=scaleFUN) +
        scale_color_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + # "#525252",
        scale_fill_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + # "#525252",
        theme(legend.position = "none", 
              plot.title = element_text(hjust = 0.5, size=33, face="bold")) +
        scale_x_discrete(breaks = bagSize.vec[seq(1,21,4)] ) + #c(1, 4, 10, 15, 20, 70, 200)) +
        theme(axis.text.x = element_text(face="bold", 
                                         angle=45, 
                                         size=rel(3),
                                         vjust = 0.6 ),  
              axis.text.y = element_text(face="bold", 
                                         size=rel(3) ),  
              axis.title.y =  element_blank(),
              axis.title.x =  element_blank() ) +   
        
        ylab(TeX('$RMSE/RMSE_{TOM}$')) +
        ggtitle( TeX(sprintf(paste0('$\\mathbf{\\sigma_{\\beta}^2 =$', sigma, 
                                    ', $\\sigma_{\\mathbf{X}}^2 =$', sigX ,
                                    "$}$")))) + 
        geom_hline(yintercept=0, 
                   linetype="dashed", 
                   color = "black", 
                   size = rel(0.5),
                   alpha = 0.7)  + 
        theme( plot.margin = margin(r = 20, b = 15) )  # add margins
    
    if(y == 1 ){
        p <- p + 
            annotate("text", x=18, y=0.035,  
                     label='bold("Stacking")',
                     size=10, color="#0868ac", parse=TRUE) + #  "#619CFF",
            annotate("text", x=18, y=0.02,  
                     label='bold("CPS")', 
                     size=10, color="#E69F00", parse=TRUE) +   # "#00BA38"  "#ca0020", , 
            annotate("text", x=18, y=0.05,  
                     label='bold("Average")', 
                     size=10, color="#ca0020", parse=TRUE) + #+     "#F8766D"   
            coord_cartesian(ylim = c(-0.03, 0.06))
        
    }else if( sim.num == 172){
        p <- p + 
            annotate("text", x=18, y=-0.1,  
                     label='bold("Stacking")',
                     size=10, color="#0868ac", parse=TRUE) + #  "#619CFF",
            annotate("text", x=18, y=-0.15,  
                     label='bold("CPS")', 
                     size=10, color="#E69F00", parse=TRUE) +   # "#00BA38"  "#ca0020", , 
            annotate("text", x=18, y=-0.05,  
                     label='bold("Average")', 
                     size=10, color="#ca0020", parse=TRUE) + #+     "#F8766D"   
        coord_cartesian(ylim = c(-0.2, 0.05))
    }
    
    assign(paste0("p",sim.num), p) 
    
    
    #### only avg
    
    pa <- ggplot(data = df[df$Method == 1, ], aes(x=factor(Bag.Size), y=RMSE, color = factor(Method), group = factor(Method))) + 
        geom_line( size = 2) + 
        geom_point( aes(shape = factor(Method) ), size = 3 ) + 
        theme_classic(base_size = 12) + 
        xlab("Bag Size") + 
        ggtitle(paste0(sim.num)) + 
        # geom_point(colour = "black")+
        scale_y_continuous(labels=scaleFUN) +
        scale_color_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + # "#525252",
        scale_fill_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + # "#525252",
        theme(legend.position = "none", 
              plot.title = element_text(hjust = 0.5, size=30, face="bold")) +
        scale_x_discrete(breaks = bagSize.vec[seq(1,21,4)] ) + #c(1, 4, 10, 15, 20, 70, 200)) +
        theme(axis.text.x = element_text(face="bold",
                                         angle=45, size=rel(3),
                                         vjust = 0.6 ),  
              axis.text.y = element_text(face="bold", 
                                         size=rel(3) ),  
              axis.title.y =  element_blank(),
              axis.title.x =  element_blank() ) +   
        
        ylab(TeX('$RMSE/RMSE_{TOM}$')) +
        ggtitle( TeX(sprintf(paste0('$\\mathbf{\\sigma_{\\beta}^2 =$', sigma, 
                                    ', $\\sigma_{\\mathbf{X}}^2 =$', sigX ,
                                    "$}$")))) + 
        geom_hline(yintercept=0, 
                   linetype="dashed", 
                   color = "black", 
                   size = rel(0.5),
                   alpha = 0.7)  + 
        theme( plot.margin = margin(r = 20, b = 15) )  # add margins
    
    
    ######################
    
    # stacking and average
    
    ps <- ggplot(data = df[df$Method !=2, ], aes(x=factor(Bag.Size), y = log(RMSE), color = factor(Method), group = factor(Method))) + 
        geom_line( size = 2) + 
        geom_point( aes(shape = factor(Method) ), size = 3 ) + 
        theme_classic(base_size = 12) + 
        scale_color_manual(values = c("#ca0020", "#0868ac")) + # "#525252",
        scale_fill_manual(values = c("#ca0020", "#0868ac")) + # "#525252",
        xlab("Bag Size") + 
        ggtitle(paste0(sim.num)) + 
        # geom_point(colour = "black")+
        scale_y_continuous(labels=scaleFUN) +
        theme(legend.position = "none", 
              plot.title = element_text(hjust = 0.5, 
                                        size=35, face="bold")) +
        scale_x_discrete(breaks = bagSize.vec[seq(1,21,4)] ) + #c(1, 4, 10, 15, 20, 70, 200)) +
        theme(axis.text.x = element_text(face="bold", 
                                         angle=45, 
                                         size=rel(3.25),
                                         vjust = 0.6 ),  
              axis.text.y = element_text(face="bold", 
                                         size=rel(3.25) ),  
              axis.title.y =  element_blank(),
              axis.title.x =  element_blank() ) + 
        scale_y_continuous(labels=scaleFUN) +
        ylab(TeX('$log(RMSE/RMSE_{TOM})$')) +
        ggtitle( TeX(sprintf(paste0('$\\mathbf{\\sigma_{\\beta}^2 =$', sigma, 
                                    ', $\\sigma_{\\mathbf{X}}^2 =$', sigX ,
                                    "$}$")))) + 
        geom_hline(yintercept=0, 
                   linetype="dashed", 
                   color = "black", 
                   size = rel(1),
                   alpha = 0.7) + 
        theme( plot.margin = margin(r = 20, b = 15) )  # add margins
    
    if(num == 172 | num == 175 | num == 180 | num == 183){
        # remove bag size labels and alter tick marks / ylimits
        ps <- ps + theme( plot.margin = margin(l = 30, b = 20, r = 15 ) ) # )  # add margins
        
        #######################
        #  tick marks on y axis
        #######################
        if(num == 172){
            ps <- ps +
                scale_y_continuous(breaks = c(-0.01, 0, 0.01), limits = c(-.018, 0.013) ) +
                theme(axis.text.x = element_blank() ) + 
                theme( plot.margin = margin(l = 10, b = 100, r = 25 ) ) # compensate for extra space on bottom in other plots
        }

        if(num == 175){
            ps <- ps +
                scale_y_continuous(breaks = c(-0.01, 0, 0.01), limits = c(-.018, 0.013) ) +
                theme(axis.text.x = element_blank() ) +
                theme(axis.text.y = element_blank() ) + 
                theme( plot.margin = margin(l = 50, b = 100, r = 30) ) # compensate for extra space on bottom in other plots
        }

        if(num == 180){
            ps <- ps +
                scale_y_continuous(breaks = round( c(-1.25, 0, 1.25), 2), 
                                   limits = c(-1.25, 2.25) ) + 
                theme( plot.margin = margin(l = 10, b = 30, r = 25 ) ) # )  # add margins
        }

        if(num == 183){
            ps <- ps +
                scale_y_continuous(breaks = round(c(-1.00, 0, 1.00), 2), 
                                   limits = c(-1.25, 2.25) ) +
                theme(axis.text.y = element_blank() ) + 
                theme( plot.margin = margin(l = 50, b = 30, r = 30 ) )  # )  # add margins
        }
    }
    
    if(y == 1 ){
        ps <- ps + 
            annotate("text", x=16, y=  0.0075,  
                     label='bold("Stacking")', 
                     size=13, color= "#0868ac", parse=TRUE) + 
            annotate("text", x=16, y= -0.0075, 
                     label='bold("Average")', 
                     size=13, color="#ca0020" , parse=TRUE) 
    }else if( sim.num == 172){
        ps <- ps + 
            annotate("text", x=16, y= 0.0075, 
                     label='bold("Stacking")', 
                     size=13, color= "#0868ac", parse=TRUE) + 
            annotate("text", x=16, y= -0.0075,  
                     label='bold("Average")', 
                     size=13, color="#ca0020", parse=TRUE) 
    }
    
    assign(paste0("ps",sim.num), ps) 
    ###################
    
    assign(paste0("pa",sim.num), pa) 
    assign(paste0("paz", num), pa + coord_cartesian(ylim = c(0.8, 8.5)) ) # zoom
    assign(paste0("p0", num), pa + coord_cartesian(ylim = c(0.9, 1.1)) )
    assign(paste0("pi", num), pa + coord_cartesian(ylim = c(1, 2.75)) )
    }
#}

ggsave("Sims5_bagSize_SSTest_Hybrid2_LineLog.png", 
       plot = grid.arrange(p160, p161,  p163, #p162,
                           p164, p165,  p167, #p166,
                           p168, p169,  p171, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3),
                                           hjust = 0.45,
                                           vjust = 0.4),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3.5))),
       width = 20, height = 15)


ggsave("Sims5_bagSize_SSTest_Clusters2_LineLog.png", 
       plot = grid.arrange(p172, p173,  p175, #p162,
                           p176, p177,  p179, #p166,
                           p180, p182,  p183, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3),
                                           hjust = 0.45,
                                           vjust = 0.4),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3.5))),
       width = 20, height = 15)


## only average
ggsave("Sims5Avg_bagSize_SSTest_Hybrid2_LineLog.png", 
       plot = grid.arrange(p0160, p0161,  p0163, #p162,
                           pi164, pi165,  pi167, #p166,
                           paz168, paz169,  paz171, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3)),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3))),
       width = 20, height = 15)


ggsave("Sims5Avg_bagSize_SSTest_Clusters2_LineLog.png", 
       plot = grid.arrange(p0172, p0173,  p0175, #p162,
                           paz176, paz177,  paz179, #p166,
                           paz180, paz182,  paz183, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{RMSE/RMSE_{TOM}}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3)),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3))),
       width = 20, height = 15)

###########
 # avg + stack

ggsave("Sims5_bagSize_SSTest_Hybrid2AvgStack_LineLog.png", 
       plot = grid.arrange(ps160, ps163, #p162,
                           ps168, ps171, nrow = 2,
                           left = textGrob(TeX('$\\mathbf{RMSE/RMSE_{TOM}}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3)),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3))),
       width = 15, height = 10)


ggsave("Sims5_bagSize_SSTest_Clusters2AvgStack_LineLog.pdf", 
       plot = grid.arrange(ps172, ps175, #p162,
                           ps180, ps183, nrow = 2,
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3.25))),
       width = 15, height = 14)
