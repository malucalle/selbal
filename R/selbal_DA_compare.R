#..............................................................................#
#                   COMPARE selbal WITH OTHER METHODOLOGIES
#..............................................................................#

#------------------------------------------------------------------------------#
# Functions to consider for the comparison:

  # DESeq2
  # edgeR
  # ANCOM
#------------------------------------------------------------------------------#


################################################################################
# STEP 1: Build auxiliar functions in order to define features differentially
#         significant.
################################################################################

#-----------------------#
#       DESeq.aux
#-----------------------#

  DESeq.aux <- function(x, y, covar = NULL){
    
    #----------------------------------------------------------------------------#
    # STEP 0: load libraries and extract information
    #----------------------------------------------------------------------------#
    
    # Load library
      library(DESeq2)
      library(pROC)
    
    #----------------------------------------------------------------------------#
    # STEP 1: auxiliar elements
    #----------------------------------------------------------------------------#
    
    # Levels of factor y
      ylev <- levels(y)
    # y as numeric  
      numy <- as.numeric(factor(y, labels = c(0,1))) - 1
    
    
    #----------------------------------------------------------------------------#
    # STEP 2: significantly differentiated bacteria
    #----------------------------------------------------------------------------#
    
    # Define input data  
      ds.counts <- t(x)
      if(is.null(covar)){covariates <- data.frame(y)
      }else{
        covariates <- data.frame(y, covar)
      }
    
    # DESeqDATAset object  
      dds<-DESeqDataSetFromMatrix(countData = ds.counts,colData = covariates,
                                  design = ~ y)
    
    #----------------------------------------------------------------------------#
    # NOTE: sometimes estimateSizeFactors() return an error, and to solve it, it 
    #       is usefull to run these lines of code:
    #
    # Modify function not to have problems with the execution (because of zeroes)
         gm_mean = function(x, na.rm=TRUE){
            exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
         }
    #----------------------------------------------------------------------------#    
    
    # Geometric means      
      geoMeans = apply(counts(dds), 1, gm_mean)
    # Estimate size factors
      diagdds = estimateSizeFactors(dds, geoMeans = geoMeans)
  
    
    # DESeq function for the analysis  
      diagdds = DESeq(diagdds, fitType="mean")
    
    # Results
      res.diagdds<-results(diagdds)
    
    # With a FDR of 0.1, DESeq analysis detects differences for these taxa
      Dif.DESeq<-colnames(x)[which(res.diagdds$padj<0.05)]
    
    
    # Define dat.model 
      dat.model <- data.frame(cbind(numy, covar, x[,Dif.DESeq]))
    
    # The regression model
      FIT <- glm(numy ~ .,data=dat.model, family = "binomial")
    
    #----------------------------------------------------------------------------#
    # STEP 3: Values to return
    #----------------------------------------------------------------------------#
    
  
  
    
    # List with elements to return
      L <- list(Dif.DESeq, FIT)
    
    # Return
      return(L)
    
  }

#-----------------------#
#       edgeR.aux
#-----------------------#

  edgeR.aux <- function(x, y, covar = NULL){
    
    #--------------------------------------------------------------------------#
    # STEP 0: load libraries and extract information
    #--------------------------------------------------------------------------#
    
    # Load library
      library(edgeR)
      library(pROC)
    
    
    #--------------------------------------------------------------------------#
    # STEP 1: auxiliar elements
    #--------------------------------------------------------------------------#
    
    # Levels of factor y
      ylev <- levels(y)
    # y as numeric  
      numy <- as.numeric(factor(y, labels = c(0,1))) - 1
      
    # Define input data  
      ds.counts <- t(x)
      if(is.null(covar)){covariates <- data.frame(y)
        }else{
          covariates <- data.frame(y, covar)
      }
    
    #----------------------------------------------------------------------------#
    # STEP 2: comparison
    #----------------------------------------------------------------------------#
    
    # Build dgList object
      dgList <- DGEList(counts = ds.counts, group = y)
    # Normalization
      dgList <- calcNormFactors(dgList, method="TMM")
    # Desing matrix
      designMat <- model.matrix(~ y)
      
    # Estimate dispersions
      dgList <- estimateGLMCommonDisp(dgList, design = designMat)
    # Differential Expression
      fit <- glmFit(dgList, designMat)
      lrt <- glmLRT(fit, coef=2)
      
    # Differentiated taxa (adjusting p-values)  
      Dif.edgeR <- colnames(x)[which(p.adjust(lrt$table$PValue, method = "BH")<0.05)]
      
    # Define dat.model 
      dat.model <- data.frame(cbind(numy, covar, x[,Dif.edgeR]))
      
    # The regression model
      FIT <- glm(numy ~ .,data=dat.model, family = "binomial")
      
    #--------------------------------------------------------------------------#
    # STEP 3: Values to return
    #--------------------------------------------------------------------------#

      
    # List with elements to return
      L <- list(Dif.edgeR, FIT)
      
    # Return
      return(L)
      
      
    
  }


#-----------------------#
#       ancom.aux
#-----------------------#

  ancom.aux <- function(x, y, covar = NULL){
    
    #----------------------------------------------------------------------------#
    # STEP 0: load libraries and extract information
    #----------------------------------------------------------------------------#
    
    # Load library
      library(ancom.R)
      library(pROC)
    
    
    #--------------------------------------------------------------------------#
    # STEP 1: auxiliar elements
    #--------------------------------------------------------------------------#
    
    # Levels of factor y
      ylev <- levels(y)
    # y as numeric  
      numy <- as.numeric(factor(y, labels = c(0,1))) - 1
    
    #--------------------------------------------------------------------------#
    # STEP 2: comparison
    #--------------------------------------------------------------------------#
      
    # Table for ANCOM analysis (including a column with the factor variable)
      ANCOM.tab <- as.data.frame(cbind(x, Group = factor(y)))
    # Run ANCOM function
      run.ANCOM <- ANCOM(OTUdat = ANCOM.tab, repeated=FALSE, sig = 0.05)
    # Results
      # Names of different detected taxa
        Dif.ANCOM <- run.ANCOM$detected
        
    # Define dat.model 
      dat.model <- data.frame(cbind(numy, covar, x[,Dif.ANCOM]))
        
      # The regression model
        FIT <- glm(numy ~ .,data=dat.model, family = "binomial")
        
    #--------------------------------------------------------------------------#
    # STEP 3: Values to return
    #--------------------------------------------------------------------------#
        
        
    # List with elements to return
      L <- list(Dif.ANCOM, FIT)
        
    # Return
      return(L)
        
    
  }
  
  
#-----------------------#
#      ALDEx2.aux
#-----------------------#
  
  ALDEx2.aux <- function(x, y, covar = NULL){
    
    #----------------------------------------------------------------------------#
    # STEP 0: load libraries and extract information
    #----------------------------------------------------------------------------#
    
    # Load library
      library(compositions)
      library(ALDEx2)
      library(pROC)
    # Clr transformed data
      x.clr <- clr(x + .5)
    
    
    #--------------------------------------------------------------------------#
    # STEP 1: auxiliar elements
    #--------------------------------------------------------------------------#
    
    # Levels of factor y
      ylev <- levels(y)
    # y as numeric  
      numy <- as.numeric(factor(y, labels = c(0,1))) - 1
    
    #--------------------------------------------------------------------------#
    # STEP 2: comparison
    #--------------------------------------------------------------------------#
    
    # Aldex.clr 
      ald.clr <- aldex.clr(reads = as.data.frame(t(x)),
                           conds = as.character(y), mc.samples = 128)
    # GLM  
      ald.glm <- aldex.glm(clr = ald.clr, conditions = y)
    
    # Differentiated taxa
      Dif.ALDEx2 <- colnames(x)[which(ald.glm$glm.eBH<0.05)]
    
    #--------------------------------------------------------------------------#
    # STEP 3: Values to return
    #--------------------------------------------------------------------------#
    
    # Build data.frame with the variables to include  
      dat.model <- data.frame(cbind(numy, covar, x.clr[,Dif.ALDEx2]))
    
    # The regression model
      FIT <- glm(numy ~ .,data=dat.model, family = "binomial")  
      
    
    # List with elements to return
      L <- list(Dif.ALDEx2, FIT)
    
    # Return
      return(L)
    
    
  }  



################################################################################
# STEP 2: Build a general function in order to evaluate the AUC
################################################################################

  DA.method.eval <- function(x, y, n.fold = 5, n.iter = 10, seed = 31415,
                          covar = NULL){
  
  # Load package plyr
    suppressMessages(library(plyr))
  
  
  #----------------------------------------------------------------------------#
  
  #----------------------------------------------------------------------------#
  # STEP 0: load libraries and extract information
  #----------------------------------------------------------------------------#
  
  # Load library
    library(compositions)
    library(DESeq2)
    library(edgeR)
    library(ancom.R)
    library(ALDEx2)
    library(pROC)
    library(ggplot2)
    library(CMA)
  
  #----------------------------------------------------------------------------#
  # STEP 1: auxiliar elements
  #----------------------------------------------------------------------------#
  
  # Levels of factor y
    ylev <- levels(y)
  # y as numeric  
    numy <- as.numeric(y) - 1
  
  # Names of x
    x.nam <- colnames(x)
  
  # clr-transformed x
    x.clr <- clr(x + 0.05)
  
  #---------------------------------------#
  # 1.1: Define cross - validation groups
  #---------------------------------------#
  
  # CROSS - VALIDATION groups
    # Load library
      suppressMessages(library(CMA))
  
    # Fix a seed
      set.seed(seed)
    # Matrix where rows represent the individuals taking part for the
    # learningsets
      CV.groups<-GenerateLearningsets(y = y, fold = n.fold, niter = n.iter,
                                      method = "CV",
                                      strat = ifelse(class(y)!="factor",
                                                     F,T))@learnmatrix
    # Unload CMA package not to have issues with pROC package
      suppressMessages(detach("package:CMA", unload=TRUE))
      suppressMessages(library(pROC))
  
  
  #--------------------------------------------#
  # 1.2: Define cross - validation function
  #--------------------------------------------#
  
      cv.DA <- function(k){
        
        #-------------------------------------#
        # Necessary matrices for the function
        #-------------------------------------#
        
        # Folds associated to a "complete FOLD"
          CV.group<-CV.groups[((k-1)*n.fold + 1):(k*n.fold),]
        
        #---------------#  
        # Build objects
        #---------------#  
        # List all the variables included
          Var.List.DESeq <- Var.List.edgeR <-list()
          Var.List.ancom <- Var.List.ALDEx2 <-list()
        # Vector for AUC  values
          AUC.test.DESeq <- AUC.test.edgeR <- vector() 
          AUC.test.ancom <- AUC.test.ALDEx2 <- vector()
        # Number of differentially abundant features
          n.diff.DESeq <- n.diff.edgeR <- vector()
          n.diff.ancom <- n.diff.ALDEx2 <- vector()
        
        # For each fold
          for (i in 1:nrow(CV.group)){
          
          # Define training data.set
            train.idx<-CV.group[i,]
          # Training dataset (x, y and covar)
            x.train<-x[train.idx,]
            x.test<-x[-train.idx,]
            y.train<-numy[train.idx]
            y.test<-numy[-train.idx]
            covar.train<-covar[train.idx,]
            covar.test<-covar[-train.idx,]
          
          # Information for training data set
            DES <- DESeq.aux(x = x.train, y = y.train, covar = covar.train)
            EDG <- edgeR.aux(x = x.train, y = y.train, covar = covar.train)
            ANC <- ancom.aux(x = x.train, y = y.train, covar = covar.train)
            ALD <- ALDEx2.aux(x = x.train, y = y.train, covar = covar.train)
          
          # Variables included in the balance (as NUMERATOR | DENOMINATOR)
            Var.List.DESeq[[i]] <- DES[[1]]
            Var.List.edgeR[[i]] <- EDG[[1]]
            Var.List.ancom[[i]] <- ANC[[1]]
            Var.List.ALDEx2[[i]] <- ALD[[1]]
          
          # Number of significant features
            n.diff.DESeq[i] <- length(DES[[1]])
            n.diff.edgeR[i] <- length(EDG[[1]])
            n.diff.ancom[i] <- length(ANC[[1]])
            n.diff.ALDEx2[i] <- length(ALD[[1]])
          
          #--------------------------------------------------------------------------#
          # Compute the AUC for the predicted values
          #--------------------------------------------------------------------------#
  
          # Build test data.frame 
            
            # Define x.test for ALDEx2 (from x.clr)
              x.test.clr <- x.clr[-train.idx,]
            
            if(!is.null(covar)){
              df.test.DESeq <- data.frame(cbind(covar.test, x.test[,DES[[1]]]))
              df.test.edgeR <- data.frame(cbind(covar.test, x.test[,EDG[[1]]]))
              df.test.ancom <- data.frame(cbind(covar.test, x.test[,ANC[[1]]]))
              df.test.ALDEx2 <- data.frame(cbind(covar.test, x.test.clr[,ALD[[1]]]))
            } else{
              df.test.DESeq <- data.frame(x.test[,DES[[1]]])
              df.test.edgeR <- data.frame(x.test[,EDG[[1]]])
              df.test.ancom <- data.frame(x.test[,ANC[[1]]])
              df.test.ALDEx2 <- data.frame(x.test.clr[,ALD[[1]]])
            }

          
          # Predictions
            PRED.DESeq <- predict(object = DES[[2]], newdata = df.test.DESeq)
            PRED.edgeR <- predict(object = EDG[[2]], newdata = df.test.edgeR)
            PRED.ancom <- predict(object = ANC[[2]], newdata = df.test.ancom)
            PRED.ALDEx2 <- predict(object = ALD[[2]], newdata = df.test.ALDEx2)
          
          # AUC for predicted values
            AUC.test.DESeq[i] <- auc(y.test, PRED.DESeq)
            AUC.test.edgeR[i] <- auc(y.test, PRED.edgeR)
            AUC.test.ancom[i] <- auc(y.test, PRED.ancom)
            AUC.test.ALDEx2[i] <- auc(y.test, PRED.ALDEx2)
          
          
        } 
        
        A <- list(Var.List.DESeq, AUC.test.DESeq, n.diff.DESeq,
                  Var.List.edgeR, AUC.test.edgeR, n.diff.edgeR,
                  Var.List.ancom, AUC.test.ancom, n.diff.ancom,
                  Var.List.ALDEx2, AUC.test.ALDEx2, n.diff.ALDEx2)
        
        return(A)
      }
  
  #------------------------------------------------------------------------------#
  
  
  #------------------------------------------------------------------------------#
  # STEP 2: Cross - validation process (parallelized)
  #------------------------------------------------------------------------------#    
  
  ################################################################################
  
  cat(paste("\n\n#-------------------------------------------------------------#",
            "\n# Starting the cross - validation procedure . . ."))
  
  # Build a parallelization scenario
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
  # Number of cores of the computer but one
  no_cores <- detectCores() - 2
  # Register the number of cores
  registerDoParallel(no_cores)
  
  # Define the function comb
  # Define the function comb
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  
################################################################################
  
# CV - procedure computed in parallel
  INTEREST <- foreach(h=1:n.iter,
                      .export=c("DESeq.aux","edgeR.aux","ancom.aux",
                                "ALDEx2.aux", "auc","CV.groups","x","numy",
                                "covar", "n.fold"),
                      .combine="comb",
                      .multicombine=TRUE,
                      #                        .init=list(list(), list()))%dopar%{
                      .init=list(list(), list(), list(), list(),
                                 list(), list(), list(), list(),
                                 list(), list(), list(), list()))%dopar%{
                        cv.DA(h)
                      }
# Stop the parallelization
  stopImplicitCluster()
  
  cat(paste("\n . . . finished.",
            "\n#-------------------------------------------------------------#",
            "\n###############################################################"))
  
#------------------------------------------------------------------------------#
  
#------------------------------------------------------------------------------#     
# STEP 3: Extract the information
#------------------------------------------------------------------------------#
  
# Variables selected by DESeq
  Var.DES <- unlist(INTEREST[[1]])
  Var.EDG <- unlist(INTEREST[[4]])
  Var.ANC <- unlist(INTEREST[[7]])
  Var.ALD <- unlist(INTEREST[[10]])
# AUC values
  AUC.DES <- unlist(INTEREST[[2]])
  AUC.EDG <- unlist(INTEREST[[5]])
  AUC.ANC <- unlist(INTEREST[[8]])
  AUC.ALD <- unlist(INTEREST[[11]])
# Number of optimal variables
  n.opt.DES <- unlist(INTEREST[[3]])
  n.opt.EDG <- unlist(INTEREST[[6]])
  n.opt.ANC <- unlist(INTEREST[[9]])
  n.opt.ALD <- unlist(INTEREST[[12]])
  
#------------------------------------------------------------------------------#
# STEP 4: Graphical representation
#------------------------------------------------------------------------------#
  
# Build a frequency table with the selected variables
  Dif.Ab.Var.DES <- (100/(n.iter*n.fold))*table(Var.DES)
  Dif.Ab.Var.EDG <- (100/(n.iter*n.fold))*table(Var.EDG)
  Dif.Ab.Var.ANC <- (100/(n.iter*n.fold))*table(Var.ANC)
  Dif.Ab.Var.ALD <- (100/(n.iter*n.fold))*table(Var.ALD)
# Order the table
  Dif.Ab.Var.DES <- Dif.Ab.Var.DES[order(Dif.Ab.Var.DES)]
  Dif.Ab.Var.EDG <- Dif.Ab.Var.EDG[order(Dif.Ab.Var.EDG)]
  Dif.Ab.Var.ANC <- Dif.Ab.Var.ANC[order(Dif.Ab.Var.ANC)]
  Dif.Ab.Var.ALD <- Dif.Ab.Var.ALD[order(Dif.Ab.Var.ALD)]
  
  # Build a data.frame
    DF.DES <- data.frame(name = names(Dif.Ab.Var.DES), 
                         sel = as.numeric(Dif.Ab.Var.DES))
    DF.DES$name <- factor(DF.DES$name, levels = names(Dif.Ab.Var.DES))
    
    DF.EDG <- data.frame(name = names(Dif.Ab.Var.EDG), 
                         sel = as.numeric(Dif.Ab.Var.EDG))
    DF.EDG$name <- factor(DF.EDG$name, levels = names(Dif.Ab.Var.EDG))
    
    DF.ANC <- data.frame(name = names(Dif.Ab.Var.ANC), 
                         sel = as.numeric(Dif.Ab.Var.ANC))
    DF.ANC$name <- factor(DF.ANC$name, levels = names(Dif.Ab.Var.ANC))
  
    DF.ALD <- data.frame(name = names(Dif.Ab.Var.ALD),
                         sel = as.numeric(Dif.Ab.Var.ALD))
    DF.ALD$name <- factor(DF.ALD$name, levels = names(Dif.Ab.Var.ALD))
    
    
  #---------------------#  
  # IMP.plot for DESeq
  #---------------------#
    
  IMP.plot.DES <- ggplot(DF.DES, aes(x=factor(name), y=sel)) +
    geom_bar(stat="identity", fill = "darkgreen",
             size=1) +
    guides(size = FALSE) + # Not to show the legend of the size
    ggtitle("DESeq2") +
    ylab("% of times resulted significant") +
    xlab("") + theme_bw() +
    coord_flip() +
    theme(strip.text.x = element_text(size=12, angle=0,
                                      face="bold",colour="white"),
          strip.text.y = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="black",
                                          fill="black"),
          plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                    face = "bold"))
  
  
  #---------------------#  
  # IMP.plot for edgeR
  #---------------------#
  
  IMP.plot.EDG <- ggplot(DF.EDG, aes(x=factor(name), y=sel)) +
  geom_bar(stat="identity", fill = "darkgreen",
           size=1) +
  guides(size = FALSE) + # Not to show the legend of the size
  ggtitle("edgeR") +
  ylab("% of times resulted significant") +
  xlab("") + theme_bw() +
  coord_flip() + 
  theme(strip.text.x = element_text(size=12, angle=0,
                                    face="bold",colour="white"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black",
                                        fill="black"),
        plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                  face = "bold"))
  
  
  #---------------------#  
  # IMP.plot for ANCOM
  #---------------------#  
    
      
  IMP.plot.ANC <- ggplot(DF.ANC, aes(x=factor(name), y=sel)) +
  geom_bar(stat="identity", fill = "darkgreen",
           size=1) +
  guides(size = FALSE) + # Not to show the legend of the size
  ggtitle("ANCOM") +
  ylab("% of times resulted significant") +
  xlab("") + theme_bw() +
  coord_flip() + 
  theme(strip.text.x = element_text(size=12, angle=0,
                                    face="bold",colour="white"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black",
                                        fill="black"),
        plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                  face = "bold"))
  
  
  #---------------------#  
  # IMP.plot for ALDEX2
  #---------------------#  
  
  IMP.plot.ALD <- ggplot(DF.ALD, aes(x=factor(name), y=sel)) +
    geom_bar(stat="identity", fill = "darkgreen",
             size=1) +
    guides(size = FALSE) + # Not to show the legend of the size
    ggtitle("ALDEx2") +
    ylab("% of times resulted significant") +
    xlab("") + theme_bw() +
    coord_flip() + 
    theme(strip.text.x = element_text(size=12, angle=0,
                                      face="bold",colour="white"),
          strip.text.y = element_text(size=12, face="bold"),
          strip.background = element_rect(colour="black",
                                          fill="black"),
          plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                    face = "bold"))
  
  
  
  
# Draw IMP.plot
  IMP.plot.DES; IMP.plot.EDG; IMP.plot.ANC; IMP.plot.ALD
  
#------------------------------------------------------------------------------#  

  
# Save the results into a list of lists:
  
   # List of plots
     L.Plot <- list(IMP.plot.DES, IMP.plot.EDG, IMP.plot.ANC, IMP.plot.ALD)
   # List of AUCs
     L.AUC <- list(DESeq2 = AUC.DES, edgeR = AUC.EDG, 
                   ANCOM = AUC.ANC, ALDEx2 = AUC.ALD)
   # List of number of selected features
     L.num <- list(n.opt.DES, n.opt.EDG, n.opt.ANC, n.opt.ALD)
  
    
  
  # Thins to return
  L <- list(L.Plot, L.AUC, L.num)
  
}

#------------------------------------------------------------------------------#
  
  
################################################################################
# STEP 3: Run the example for Crohn's disease
################################################################################
  
# Define x and y
  x <- as.matrix(Crohn[,-49])
  y <- factor(Crohn$y, labels = c(0,1))

#----------------------------------------------------#    
# Analysis for all the considered methods but selbal
#----------------------------------------------------#  
  # Initial time
    pc <- proc.time()
    # Run DA.method.eval function
      R <-  DA.method.eval(x,y)
  # Elapsed time
    proc.time() - pc
  
#----------------------------------------------------#
# Analysis with selbal
#----------------------------------------------------#
  # Initial time  
    pc <- proc.time()
    # Run selbal
      sel <- selbal.cv(x,y)
  # Elapsed time
    proc.time() - pc
    
#------------------------------------------------------------------------------#
        
    
################################################################################
# STEP 4: Extract and represent the results
################################################################################
    
# AUC results
  AUC.res <- R[[2]]; names(AUC.res) <- c("DESeq", "edgeR", "ANCOM","ALDEx2")
  # For all the methods but selbal  
    lapply(AUC.res,summary)
  # For selbal  
    summary(sel[[5]])
    
# Results for the number of variables
  AUC.var <- R[[3]]; names(AUC.var) <- c("DESeq", "edgeR", "ANCOM","ALDEx2")
  # For all the methods but selbal
    lapply(AUC.var,summary)
  # For selbal
    sel[[8]]
    
#------------------------------------------------------------------------------#
# RESULTS
#------------------------------------------------------------------------------#
    
#  $DESeq
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.7159  0.7545  0.7779  0.7752  0.7983  0.8254 
    
#  $edgeR
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.7139  0.7527  0.7665  0.7721  0.7899  0.8273 
    
#  $ANCOM
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.6053  0.6974  0.7169  0.7125  0.7365  0.7873 
    
#  $ALDEx2
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.7255  0.8033  0.8140  0.8156  0.8327  0.8729 
    
# selbal
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.7389  0.8024  0.8219  0.8196  0.8423  0.8792
    
#------------------------------------------------------------------------------#
    
    
################################################################################  
# GRAPHICAL REPRESENTATION of the RESULTS    
################################################################################
    
#-----------------------------#        
# Boxplot with the AUC values
#-----------------------------#    
    
# Values in a data.frame instead of a list
  AUC.tab <- data.frame(DESeq = R[[2]][[1]],
                        edgeR = R[[2]][[2]],
                        ANCOM = R[[2]][[3]],
                        ALDEx2 = R[[2]][[4]],
                        selbal = sel[[5]])
# Melt object
  library(reshape2)
  AUC.tab.melt <- melt(AUC.tab)
    
# Graphical representation
  library(ggplot2)
  ggplot(AUC.tab.melt, aes(x=factor(variable), y=value)) +
      geom_boxplot() +
      xlab("") + theme_bw() +
      ylab ("AUC") +
      ggtitle("Method comparison") +
      theme(strip.text.x = element_text(size=12, angle=0,
                                        face="bold",colour="white"),
            strip.text.y = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",
                                            fill="black"),
            plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                      face = "bold"))
    
#-----------------------------#        
# Method comparison
#-----------------------------#  
    
# Add the sample ID      
  AUC.tab2 <- cbind(AUC.tab, SampleID = factor(1:nrow(AUC.tab)))  
  # Melt 
    library(reshape2)
    AUC.tab2.melt <- melt(AUC.tab2)
    
  # Graphical representation
    ggplot(data=AUC.tab2.melt ,aes(x=variable,y=value, group=SampleID)) +
      geom_line(size=1, color="black")+theme_bw()+
      ggtitle("AUC values (according to folds)") +
      theme(legend.position="none",
            strip.background=element_blank(),
            strip.text=element_blank())+
      geom_point(color="firebrick1")+
      geom_vline(xintercept=6, linetype="dotted", color="black", size=2)+
      theme(axis.text=element_text(size=12), axis.title.x=element_blank(),
            axis.title.y=element_blank())+
      theme(panel.border=element_rect(linetype="solid", color="black", size=0.5),
            plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                      face = "bold"))
    
    
#--------------------------------------#
#  Represent the number of variables
#--------------------------------------#  
    
# Values in a data.frame instead of a list
  nvar.tab <- data.frame(DESeq = R[[3]][[1]],
                         edgeR = R[[3]][[2]],
                         ANCOM = R[[3]][[3]],
                         ALDEx2 = R[[3]][[4]])
# Melt object
  library(reshape2)
  nvar.tab.melt <- melt(nvar.tab)
    
# Graphical representation
  library(ggplot2)
  ggplot(nvar.tab.melt, aes(x=factor(variable), y=value)) +
      geom_boxplot() +
      xlab("") + theme_bw() +
      ylab ("Significant variables") +
      ggtitle("Number of significant variables") +
      theme(strip.text.x = element_text(size=12, angle=0,
                                        face="bold",colour="white"),
            strip.text.y = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",
                                            fill="black"),
            plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                      face = "bold"))
  
#------------------------------------------------------------------------------#  