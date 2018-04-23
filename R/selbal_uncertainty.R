#------------------------------------------------------------------------------#
#             WORKING with SELBAL acknowledging the VARIABILITY
#------------------------------------------------------------------------------#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# IMPORTANT: the MCMCpack package is required to run this code

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


#------------------------------------------------------------------------------#
# FUNCTION: selbal.dir
# 
# INPUT:
#         
#
#         x: a matrix object with the information of variables (\emph{columns}) 
#            for each sample (\emph{rows}).
#         y: the response variable, either continuous or dichotomous.
#    n.fold: number of folds in which to divide the whole data set.
#    n.iter: number of iterations for the cross - validation process.
#      seed: a seed to make the results reproducible.
#    th.imp: the minimum increment needed when adding a new variable into
#            the balance, in order to consider an improvement.
#     covar: a data.frame with the variables to adjust for (columns)
#       col: vector of two colours for differentiate the variables
#            appearing in the numerator and in the denominator of the balances.
#      col2: vector of three colours for the lines of the barplot with
#            the aim of identifying if each variable appears in the Global 
#            balance, in the CV - balance or in both of them.
# logit.acc: the measure to compute for the correlation between \code{Y}
#            and the proposed \emph{balance} adjusting for covariates. 
#            One of the following values: "Rsq" (default), "AUC" or "Tjur".
#      maxV: numeric value defining the maximum number of variables
#            composing the balance. Default 1e10 to give prevalence to th.imp
#            parameter.
#  zero.rep: "bayes" for BM-replacement or "one" to add one read tho each cell
#            of the matrix.
#    mc.sam: Monte-Carlo samples to generate
#------------------------------------------------------------------------------#


  selbal.dir <- function(x, y, n.fold = 5, n.iter = 10, seed = 31415,
                         covar = NULL, col = c("steelblue1", "tomato1"),
                         col2 = c("darkgreen", "steelblue4","tan1"),
                         logit.acc = "AUC", maxV = 20,
                         zero.rep = "bayes", mc.sam = 50){
  
  # Load library
    library(selbal)
      
    
  #----------------------------------------------------------------------------#
  # STEP 1: run selbal.cv 
  #----------------------------------------------------------------------------#
    
    global <- selbal.cv(x,y, n.fold = n.fold, n.iter = n.iter, seed = seed,
                        covar = covar, col = col, col2 = col2, 
                        logit.acc = logit.acc, maxV = maxV, 
                        zero.rep = zero.rep)
    
  # Define variables in NUM and DEN
    NUM <- global$balance$Taxa[global$balance$Group=="NUM"]
    DEN <- global$balance$Taxa[global$balance$Group=="DEN"]
    
    
  # Extract the optimal number of variables in the balance
    optV <- global[[length(global)]]

    
  #----------------------------------------------------------------------------#
  # STEP2 : Dirichlet Monte-Carlo sampling and selbal for each one
  #----------------------------------------------------------------------------#
    
  # The matrix without zeros
    x.non0 <- cmultRepl2(x, zero.rep = zero.rep)
    # Compute sampling depth for each sample
      s.depth <- apply(x.non0,1,sum)
  # Load library
    library(MCMCpack)
    
  # Generate a list of  "mc.sam" elements with different values
    
  #----------------------------------------------------------------------------#
  # AUXILIAR FUNCTION
  #----------------------------------------------------------------------------#
      
    sim.mc <- function(x, s.depth){
      # Generate the matrix using the Dirichlet distribution
        M <- t(apply(x,1,function(x) rdirichlet(1,x)))*s.depth
        colnames(M) <- colnames(x)
      # Return M  
        return(M)
    }
    
    
    # Build a list to save the new matrices
      L.mat <- list()
      Res <- list()
    
    # Set seed
      set.seed(seed)
    # Build the "mc.sam" matrices  
      for (i in 1:mc.sam){ 
        print(i)
        L.mat[[i]] <- sim.mc(x,s.depth)
        Res[[i]] <- selbal.aux(L.mat[[i]],y,maxV = optV)
      }

  #----------------------------------------------------------------------------#
  # STEP 3: extract the information of all the balances
  #----------------------------------------------------------------------------#
      
  # Unlist balances
    Balances <- Res
  # Build a matrix to save the information of each balance  
    BR<-matrix(0,nrow =mc.sam , ncol = ncol(x))
    colnames(BR) <- colnames(x)
    # Codification of each balance
      for (i in 1:length(Balances)){
        BR[ i, colnames(BR) %in% Balances[[i]][Balances[[i]][,"Group"]=="NUM","Taxa"]] <- 1
        BR[ i, colnames(BR) %in% Balances[[i]][Balances[[i]][,"Group"]=="DEN","Taxa"]] <- 2
    }

  # Individual information of each variable
    # ROB.TAB
      ROB.TAB <- matrix(0, nrow = ncol(x), ncol=3)
      colnames(ROB.TAB) <- c("Prop. Included", "Prop_Numerator",
                             "Prop_Denominator")
      row.names(ROB.TAB) <- colnames(x)
    
  # Complete the table
    ROB.TAB[,1] <- apply(BR!=0,2,function(x) 100*mean(x))
    ROB.TAB[,2] <- apply(BR==1,2,function(x) 100*mean(x))
    ROB.TAB[,3] <- apply(BR==2,2,function(x) 100*mean(x))
      
  # Variables with at least included in a balance
    fil1 <- which(ROB.TAB[,1]!=0)
    ord.ROB.TAB <- ROB.TAB[fil1,]
  # Order by the first row
    sel.ord <- order(ord.ROB.TAB[,1], decreasing = F)
    ord.ROB.TAB <- ord.ROB.TAB[sel.ord,]
      # Define a data.frame to plot the results
        BAL.SEL.TAB <- data.frame(name = row.names(ord.ROB.TAB),
                                  sel = ord.ROB.TAB[,1])
      # Define the levels of the $name
        BAL.SEL.TAB$name <- factor(BAL.SEL.TAB$name, levels = BAL.SEL.TAB$name)
      # Define the color for variables in the numerator (as overall)
        COLOR.BAL <- rep(col[2],nrow(BAL.SEL.TAB))
      # If the variable is in the denominator modify the color:
        # Variables in the denominator
          vDEN <- row.names(ord.ROB.TAB)[which(ord.ROB.TAB[,"Prop_Denominator"]!=0)]
          COLOR.BAL[row.names(BAL.SEL.TAB) %in% vDEN] <- col[1]
      # Add COLOR.BAL to BAL.SEL.TAB
        BAL.SEL.TAB$COLOR.BAL <- factor(COLOR.BAL, levels = col, labels = col)
      
      
      
      
      #------------------------------------------#
      # Barplot representation
      #------------------------------------------#
      
      # Load library ggplot2
        suppressMessages(library(ggplot2))
      
      # IMP.plot
        IMP.plot <- ggplot(BAL.SEL.TAB, aes(x=factor(name), y=sel)) +
          geom_bar(stat="identity", aes(fill = COLOR.BAL),
                 size=1) +
          guides(size = FALSE) + # Not to show the legend of the size
          scale_fill_manual(name = "Group of . . .",
                            values = c(col[1], col[2]),
                            breaks = c(col[1], col[2]),
                            labels=c("DEN","NUM")) +
          scale_color_manual(name ="Variables \n appearing in . . .",
                             values = c(col2,"white"),
                             breaks = c(col2,"white"),
                             labels = c("Both", "Global", "CV" ,"NONE"),
                             drop=F,
                             guide=guide_legend(
                              override.aes = list(fill="gray90"))) +
          ylab("% of times included in a balance") +
          xlab("") + theme_bw() +
          coord_flip() +
          ggtitle("Cross validation in balance selection") +
          theme(strip.text.x = element_text(size=12, angle=0,
                                            face="bold",colour="white"),
                strip.text.y = element_text(size=12, face="bold"),
                strip.background = element_rect(colour="black",
                                                fill="black"),
                plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                          face = "bold"),
                legend.title = element_text(face="bold"),
                legend.text = element_text(face="bold"))      
      
      
      
    #-----------------------------------#
    # Most repeated balances
    #-----------------------------------#
        
    # Balances' strings
      BAL.str <- apply(BR,1, function(x) paste(x, collapse=""))
    # Resume the information
      BAL.tab <- prop.table(table(BAL.str))
    # Names of appearing balances
      nam.str <- names(BAL.tab)
    # Values
      nam.A <- t(sapply(nam.str, FUN = function(x) unlist(strsplit(x,""))))
      x.nam <- colnames(x)
      
    # Variables included in the most abundant balances
      INC <- apply(nam.A, 1, function(x) x.nam[x!=0])
    # Variables included in the numerator of each selected balance
      INC.NUM <- alply(nam.A, 1, function(x) x.nam[x==1])
    # Variables included in the denominator of each selected balance
      INC.DEN <- alply(nam.A, 1, function(x) x.nam[x==2])
        
      # Variables selected
        UNIQUE.VAR <- unique(c(as.vector(unlist(INC)),NUM, DEN))
      # Build a data.frame to represent
        RESUME.BAL <- as.data.frame(matrix(0, nrow =length(UNIQUE.VAR),
                                           ncol = length(BAL.tab)))
        row.names(RESUME.BAL) <- UNIQUE.VAR
        
      # Put "NUM" if the variable is in the numerator of the balance
        RESUME.BAL[sapply(INC.NUM, function(x) UNIQUE.VAR %in% x)] <- "NUM"
      # Put "DEN" if the variable is in the denominator of the balance
        RESUME.BAL[sapply(INC.DEN, function(x) UNIQUE.VAR %in% x)] <- "DEN"
        
      # Add the relative frequency of the balances
        RESUME.BAL <- rbind(RESUME.BAL, FREQ=as.numeric(BAL.tab))
        
        
      # Add two new columns (one for Global Balance and another for percentages)
        RESUME.BAL <- cbind(RESUME.BAL, 0, 0)
        
      # Add the information of the global balance
        RESUME.BAL[row.names(RESUME.BAL)%in%NUM,ncol(RESUME.BAL)] <- "NUM"
        RESUME.BAL[row.names(RESUME.BAL)%in%DEN,ncol(RESUME.BAL)] <- "DEN"
        
      # NEW
        RESUME.BAL[-nrow(RESUME.BAL) ,ncol(RESUME.BAL) - 1] <-
          ROB.TAB[row.names(RESUME.BAL)[-nrow(RESUME.BAL)],1]
        
      # Order RESUME.BAL by FREQ
        RESUME.BAL <- RESUME.BAL[,c(ncol(RESUME.BAL), ncol(RESUME.BAL)-1,
                                    order(RESUME.BAL[nrow(RESUME.BAL),
                                                     -c(ncol(RESUME.BAL),
                                                        ncol(RESUME.BAL)-1)],
                                          decreasing = T))]
        
      # No frequency for the Global balance and the CV.Balance
        RESUME.BAL[nrow(RESUME.BAL),1:2]<- "-"
        
        
      # Data to plot (maximum 5 different balances)
        dat <- RESUME.BAL[,c(1,2:(min(5,ncol(RESUME.BAL))))]
        W <- which(apply(dat[,-2]==0,1,mean)==1)
        if(length(W)!=0){ dat <- dat[-as.numeric(W),]}
        # Change the order of first and second colum
          dat <- dat[,c(2,1,3:ncol(dat))]
          colnames(dat)[1:2]<-c("%","Global")
        
        # Order dat (rows ordered by their presence percentage)
          dat<-dat[c(order(as.numeric(dat[-nrow(dat),1]),decreasing=T),nrow(dat)),]
        
  # Return dat for the representation
    return(dat)
      
    
  }
  
#------------------------------------------------------------------------------#
  

#------------------------------------------------------------------------------#
# Load data
#------------------------------------------------------------------------------#

# Load selbal library
  library(selbal)
    
# Define x and y (Crohn data set)
  x <- as.matrix(Crohn[,-49])
  y <- factor(Crohn$y, labels = c(0,1))
  
# Run selbal dir function  
  pc <- proc.time()
  A <- selbal.dir(x,y, mc.sam = 100)
  proc.time() - pc
  
#------------------------------------------------------------------------------#  

# Graphical representation
  plot.tab(A)

