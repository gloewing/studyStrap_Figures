library(ggplot2)
library(gridExtra)
library(viridis)
library(knitr)
library(kableExtra)
library(grid)
library(latex2exp)
library(tidyverse)
dslabs::ds_theme_set(new = "theme_minimal")
setwd("~/Desktop/Research/sims5")

sim.name <- c(160:191) # order to correspond with figure # add 138 later
num.sims <- length(sim.name)

sims <- cbind(160:183, expand.grid(c(0.05, 0.25, 1, 3), c(0.05, 5, 20)^2, c(0,4)))
sims2 <- cbind( 184:191, expand.grid(c(0.05,0.25,1,3),c(10,15)^2, c(0)))
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
    
    res <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Avg")) / merged ) 
    res.cpsw <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_CPS")) / merged ) 
    res.stack <- colMeans(read.csv(paste0("LASSO AR BagSize_Test_Sims5_",num,"_Stack")) / merged ) 
    
    # remove 2 columns for visual purposes as these are superfluous
    res <- res[-rm]
    res.cpsw <- res.cpsw[-rm]
    res.stack <- res.stack[-rm]


    minBSz <- which.min(res)

    bags <- bagSzs <- bagSize.vec
    bags <- c(bags)
    df.cpsw <- data.frame(RMSE = res.cpsw, Bag.Size = bags, Method = rep(2, length(bags)))
    df.stack <- data.frame(RMSE = res.stack, Bag.Size = bags, Method = rep(3, length(res.stack)))
    
    
    bags <- factor(bags, levels = c(bagSzs))
    df <- data.frame(RMSE = res, Bag.Size = bags, Method = 1)
    df <- rbind(df,df.cpsw,df.stack)

    p <- ggplot(data = df, aes(x=factor(Bag.Size), y=log(RMSE), color = factor(Method), group = factor(Method))) + 
        geom_line( size = 2) +  # aes(linetype = factor(Method)), 
        geom_point( aes(shape = factor(Method) ), size = 3 ) + 
        xlab("Bag Size") + 
        ggtitle(paste0(sim.num)) + 
        # geom_point(colour = "black")+
        scale_y_continuous(labels=scaleFUN) +
        theme_classic(base_size = 12) + 
        theme(legend.position = "none", 
              plot.title = element_text(hjust = 0.5, size=33, face="bold")) +
        scale_x_discrete(breaks = bagSize.vec[seq(1,21,4)] ) + #c(1, 4, 10, 15, 20, 70, 200)) +
        theme(axis.text.x = element_text(face="bold", #color="#993333", 
                                         angle=45, size=rel(3.1),
                                         vjust = 0.6),  
              axis.text.y = element_text(face="bold", 
                                         size=rel(3.1) ),  
              axis.title.y =  element_blank(),
              axis.title.x =  element_blank(),
              plot.margin = margin(10, 15, 10, 10) # make space for 200 on RHS
        ) +    
        scale_color_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + # "#525252",
        scale_fill_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + # "#525252",
        ylab(TeX('$RMSE/RMSE_{TOM}$')) +
        ggtitle( TeX(sprintf(paste0('$\\mathbf{\\sigma_{\\beta}^2 =$', sigma, 
                                    ', $\\sigma_{\\mathbf{X}}^2 =$', sigX ,
                                    "$}$")))) + 
        geom_hline(yintercept=0, 
                   linetype="dashed", 
                   color = "black", 
                   size = rel(1),
                   alpha = 0.7) #   
    
    p <- p + theme( plot.margin = margin(r = 30, b = 15, l = 20) )  # add margins
    
    if(y == 1 ){
        p <- p + 
            annotate("text", x=18, y=0.075,  
                     label='bold("CPS")', 
                     size=10, color="#E69F00", parse=TRUE) + #  "#619CFF",
            annotate("text", x=18, y=0.1,  
                     label='bold("Stack")', 
                     size=10, color="#0868ac", parse=TRUE) +   # "#00BA38"  "#ca0020", , 
            annotate("text", x=18, y=0.125,  
                     label='bold("Average")', 
                     size=10, color="#ca0020"  , parse=TRUE)  #+     "#F8766D"
        # coord_cartesian(ylim = c(0.6, 1.4))
        
    }else if( sim.num == 172){
        p <- p + 
            annotate("text", x=18, y=-0.2,  
                     label='bold("CPS")', 
                     size=10, color="#E69F00", parse=TRUE) + #  "#619CFF",
            annotate("text", x=18, y=-0.15,  
                     label='bold("Stack")', 
                     size=10, color="#0868ac", parse=TRUE) +   # "#00BA38"  "#ca0020", , 
            annotate("text", x=18, y=-0.1,  
                     label='bold("Average")', 
                     size=10, color="#ca0020"  , parse=TRUE)  #+     "#F8766D"
           # coord_cartesian(ylim = c(0.6, 1.4))
    }
    
    assign(paste0("p",sim.num), p) 
    
    
    #### only avg
    
    pa <- ggplot(data = df[df$Method == 1, ], aes(x=factor(Bag.Size), 
                    y=log(RMSE), color = factor(Method), group = factor(Method))) + 
        geom_line( aes(linetype = factor(Method)), size = 2) + 
        geom_point( aes(shape = factor(Method) ), size = 3 ) + 
        xlab("Bag Size") + 
        ggtitle(paste0(sim.num)) + 
        # geom_point(colour = "black")+
        scale_y_continuous(labels=scaleFUN) +
        # theme_void(base_size = 12) + 
        theme_classic(base_size = 12) + 
        theme(legend.position = "none", 
              plot.title = element_text(hjust = 0.5, size=33, face="bold")) +
        scale_x_discrete(breaks = bagSize.vec[seq(1,21,4)], #c(1, 4, 10, 15, 20, 70, 200)) +
                         limits = factor(bags)) +
        theme(axis.text.x = element_text(face="bold", #color="#993333", 
                                         angle=45, size=rel(3.1),
                                         vjust = 0.6),  
              axis.text.y = element_text(face="bold", 
                                         size=rel(3.1) ),  
              axis.title.y =  element_blank(),
              axis.title.x =  element_blank(),
              plot.margin = margin(10, 15, 10, 10) # make space for 200 on RHS
              ) +   
        scale_color_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + # "#525252",
        scale_fill_manual(values = c("#ca0020", "#E69F00", "#0868ac")) + # "#525252",
        ylab(TeX('$log(RMSE/RMSE_{TOM})$')) +
        ggtitle( TeX(sprintf(paste0('$\\mathbf{\\sigma_{\\beta}^2 =$', sigma, 
                                    ', $\\sigma_{\\mathbf{X}}^2 =$', sigX ,
                                    "$}$")))) + 
        geom_hline(yintercept=0, 
                   linetype="dashed", 
                   color = "black", 
                   size = rel(0.5),
                   alpha = 0.7) #  
    
    

    if(num == 172 | num == 176 | num == 180){
        # remove bag size labels and alter tick marks / ylimits
        pa <- pa + theme( plot.margin = margin(l = 30, r = 20, b = 15) )  # add margins
        
        #######################
        #  tick marks on y axis
        #######################
        if(num == 172){
            pa <- pa +
                scale_y_continuous(breaks = c(-0.2, -0.1, 0), limits = c(-.25, 0.025) ) + 
                theme(axis.text.x = element_blank() ) + 
                theme( plot.margin = margin(l = 30, r = 20, b = 70) )  # add margins
        }
        
        if(num == 180){
            pa <- pa +
                scale_y_continuous(breaks = c(-0.2, 0.4, 1.0), limits = c(-.2, 1.1) )
        }
        
        if(num == 176){
            pa <- pa +
                scale_y_continuous(breaks = c(-0.3, 0, 0.3), limits = c(-.30, 0.5) ) + 
                theme(axis.text.x = element_blank() ) + 
                theme( plot.margin = margin(l = 30, r = 20, b = 70) )  # add margins
        }
    }
    
  
    
    assign(paste0("pa",sim.num), pa) 
    assign(paste0("paz", num), pa + #coord_cartesian(ylim = c(0.7, 2))  + 
        geom_vline(xintercept=minBSz, color = "black") ) # zoom and add vertical line
    
    assign(paste0("p0", num), pa )#+ coord_cartesian(ylim = c(0.9, 1.2)) )
    assign(paste0("pi", num), pa ) #+ coord_cartesian(ylim = c(1, 3)) )
    assign(paste0("pl", num), pa ) #+ coord_cartesian(ylim = c(1, 6)) )
    
}

ggsave("Sims5_bagSize_ARTest_Hybrid5Log.png", 
       plot = grid.arrange(p160, p161,  p163, #p162,
                           p164, p165,  p167, #p166,
                           p168, p169,  p171, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3),
                                           hjust = 0.45,
                                           vjust = 0.4),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3))),
       width = 20, height = 15)


ggsave("Sims5_bagSize_ARTest_Clusters5Log.png", 
       plot = grid.arrange(p172, p173,  p175, #p162,
                           p176, p177,  p179, #p166,
                           p180, p182,  p183, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3),
                                           hjust = 0.45,
                                           vjust = 0.4),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3))),
       width = 20, height = 15)


## only average
ggsave("Sims5Avg_bagSize_ARTest_Hybrid5Log.png", 
       plot = grid.arrange(p0160, p0161,  p0163, #p162,
                           pi164, pi165,  pi167, #p166,
                           pl168, pl169,  pl171, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3)),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3))),
       width = 20, height = 15)


ggsave("Sims5Avg_bagSize_ARTest_Clusters5Log.png", 
       plot = grid.arrange(paz172, paz173,  paz175, #p162,
                           paz176, paz177,  paz179, #p166,
                           paz180, paz182,  paz183, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3)),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3))),
       width = 20, height = 15)
#

ggsave("Sims5Avg_bagSize_ARTest_Clusters5LogReduced.pdf", 
       plot = grid.arrange(paz172, #paz173,  paz175, #p162,
                           paz176, #paz177,  paz179, #p166,
                           paz180, #paz182,  paz183, nrow = 3,
                           left = textGrob(TeX('$\\mathbf{log(RMSE/RMSE_{TOM})}$'), rot = 90,
                                           gp = gpar(fontface = "bold", cex = 3)),
                           bottom = textGrob('Bag Size', vjust = 0,
                                             gp = gpar(fontface = "bold", cex = 3.03))),
       width = 7, height = 13)
