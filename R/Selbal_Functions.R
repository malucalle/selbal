#library(selbal)


# Define the function selbal
  selbal <- function(x, y, th.imp = 0, covar = NULL, logit.acc="AUC",
                     logt=T, col = c("steelblue1", "tomato1"), tab=T,
                     draw=T, maxV = 1e10, zero.rep = "bayes"){

      # y1=as.numeric(y)
      # nulldev0<-deviance(glm(y1~1))
      #
  #----------------------------------------------------------------------------#
  # STEP 0: load libraries and extract information
  #----------------------------------------------------------------------------#

    #----------------------------------------------#
    # 0.1: information about the response variable
    #----------------------------------------------#
    # Class of the response variable
      classy <- class(y)
    # Family for the glm (default: gaussian)
      f.class <- "gaussian"
    # numy to be y and it will be modified if y is a factor
      numy <- y
    # If y is a factor, compute the number of levels of y
      if (classy == "factor"){
        ylev <- levels(y)
        numy <- as.numeric(y) - 1
        f.class <- "binomial"
  #      nulldev0<-deviance(glm(numy~1, family=f.class))
      }

    #------------------------------------------------------------------#
    # 0.2: information and transformation of the independent variables
    #------------------------------------------------------------------#

    # Load library
      suppressMessages(library(zCompositions))

    # Variables name
      var.nam <- rem.nam <- colnames(x)

    # Build a table with the response variable and covariates for correction
      if (!is.null(covar)){ dat <- data.frame(cbind(numy, covar))
      } else { dat <-data.frame(numy)}

    # The logCounts (with zero replacement)
      if (logt == F){ logCounts <- x
      } else{
        logCounts <- log(cmultRepl2(x, zero.rep = zero.rep))
      }


    #--------------------------------------------------------------------------#
    # 0.3: auxiliar functions
    #--------------------------------------------------------------------------#

    #--------------------------------------------------------------------------#
    # Auxiliar function to compute the first balance
    #--------------------------------------------------------------------------#

    first.bal<- function(logCounts, Y, covar=NULL){

      #------------------------------------------------------------------------#
      # STEP 0: extract information
      #------------------------------------------------------------------------#

      # Number and name of variables
        n <- ncol(logCounts)
        nam <- colnames(logCounts)


      #------------------------------------------------------------------------#
      # STEP 1: build the output matrix
      #------------------------------------------------------------------------#

      # The matrix
        if (classy=="factor"){ M<-matrix(0, nrow=n, ncol=n)
        }else{ M<-matrix(1e99, nrow=n, ncol=n)}
      # Row.names and colnames
        row.names(M)<-colnames(M)<-nam

      #------------------------------------------------------------------------#
      # STEP 2: complete the matrix
      #------------------------------------------------------------------------#

      if(classy=="factor"){

        # Solve the problem with libraries
          suppressWarnings(suppressMessages(library("CMA")))
          suppressMessages(detach("package:CMA", unload=TRUE))
          suppressMessages(library(pROC))

          for (i in 2:n){
            for (j in 1:(i-1)){
            # Build a table with the information
              TAB <- data.frame(cbind(Y,logCounts[,i]-logCounts[,j],covar))
            # Fit the regression model
              FIT <- glm(Y ~ .,data=TAB, family = f.class)

            # Add the value into the matrix
              ifelse(FIT$coefficients[2]>0,
                     M[i,j] <- logit.cor(FIT, y = y, covar = covar, logit.acc),  #CANVI Y
                     M[j,i] <- logit.cor(FIT, y = y, covar = covar, logit.acc)) #CANVI Y

            } # End j
          } # End i



        # Indices for the highest logit.cor value
          r <- which(M == max(M), arr.ind = T)

        } else {
          for (i in 2:n){
            for (j in 1:(i-1)){
            # Build a table with the information
              TAB <- data.frame(cbind(Y,logCounts[,i]-logCounts[,j],covar))
            # Fit the regression model
              FIT <- glm(Y ~ .,data=TAB, family = f.class)
            # Complete the matrix
              ifelse(FIT$coefficients[2]>0,
                     M[i,j] <- mean(FIT$residuals^2),
                     M[j,i] <- mean(FIT$residuals^2))
            } # End j
          } # End i



        # Indices for the lowest MSE value
          r <- which(M == min(M), arr.ind = T)
        }

#browser()
      # Return the row and column of the maximum value
        return(r)
      }

    #--------------------------------------------------------------------------#

    #--------------------------------------------------------------------------#
    # Auxiliar function to compute the "association value" when adding a new
    # variable into the balance*
    #--------------------------------------------------------------------------#

    balance.adding.cor <- function(x, LogCounts, POS, NEG, numy, covar=NULL){

      #----------------------------------------#
      # If x added into the numerator, . .
      #----------------------------------------#

      # The "numerator"
        S1.pos <- rowM(LogCounts[,c(POS,x)]); s1 <- length(POS) + 1
      # The "denominator"
        S2.pos <- rowM(LogCounts[,NEG])     ; s2 <- length(NEG)
      # The balance
        BAL <- sqrt((s1*s2)/(s1+s2))*(S1.pos - S2.pos)

      # Data.frame with the variables
        D.pos <- data.frame(cbind(numy, BAL, covar))

      # Regression model
        FIT.pos <- glm(numy~., data=D.pos, family=f.class)
      # The MSE or the corresponding value for dichotomous responses
        if(classy=="numeric"){ C.pos <- mean(FIT.pos$residuals^2)
     #   }else{ C.pos <- logit.cor(FIT.pos,numy,covar = covar, logit.acc)} #diferencies
        }else{ C.pos <- logit.cor(FIT.pos,y,covar = covar, logit.acc)} #diferencies
      #----------------------------------------#
      # If x added into the numerator, . .
      #----------------------------------------#

      # The numerator
        S1.neg <- rowM(LogCounts[,POS])       ; s1 <- length(POS)
      # The denominator
        S2.neg <- rowM(LogCounts[,c(NEG,x)])  ; s2 <- length(NEG) + 1
      # The balance
        BAL <- sqrt((s1*s2)/(s1+s2))*(S1.neg - S2.neg)

      # Data.frame with the variables
        D.neg <- data.frame(cbind(numy, BAL, covar))

      # Regression model
        FIT.neg <- glm(numy~., data=D.neg, family=f.class)
      # The MSE or the corresponding value for dichotomous responses
        if(classy=="numeric"){ C.neg <- mean(FIT.neg$residuals^2)
        #}else{ C.neg <- logit.cor(FIT.neg,numy,covar = covar, logit.acc)}  #diferencies
        }else{ C.neg <- logit.cor(FIT.neg,y,covar = covar, logit.acc)}  #diferencies

      # Correlation values
        COR <- c(C.pos, C.neg)
      # Return the values
        return(COR)
      }

#------------------------------------------------------------------------------#


    #--------------------------------------------------------------------------#
    # STEP 1: depending on the response variable class, . . .
    #--------------------------------------------------------------------------#

    # Define the first balance
      A1 <- first.bal(logCounts, Y = numy, covar=covar)
    # Variables taking parti into the first balance
      POS <- colnames(x)[A1[1,1]]
      NEG <- colnames(x)[A1[1,2]]

    # Included variables in the model
      INC.VAR <- c(POS, NEG)
    # Delete these variables from rem.nam
      rem.nam <- setdiff(rem.nam, INC.VAR)

    # Define the initial balance (B_1)
      S1 <- logCounts[,POS]
      S2 <- logCounts[,NEG]
    # Assign the values to B
      B <- sqrt(1/2)*(S1 - S2)

    #--------------------------------------------------------------------------#
    # Information about the ACC for the Balance values
    #--------------------------------------------------------------------------#
    # Build a new data.frame
      dat.ini <- cbind(dat, B)
    # Fit the regression model
      FIT.initial <- glm(numy ~ .,data=dat.ini, family = f.class)

    # Solve the problem with libraries
      suppressWarnings(suppressMessages(library("CMA")))
      suppressMessages(detach("package:CMA", unload=TRUE))
      suppressMessages(library(pROC))

    # Define the initial "accuracy" or "association" value
      if(classy=="numeric"){ ACC.Bal <- mean(FIT.initial$residuals^2)
      }else{ ACC.Bal <- logit.cor(FIT.initial, y, covar = covar, logit.acc)}   #diferencies  numy
    #  }else{ ACC.Bal <- logit.cor(FIT.initial, numy, covar = covar, logit.acc)}   #new diferencies  numy

#browser()

  #----------------------------------------------------------------------------#

    # ACC reference
    ACC.ref <- ACC.Bal

    #------------------------------------------------------------------------#
    # Improve the balances
    #------------------------------------------------------------------------#
    # Define some parameters
    # The p.value to compare 2 balances (one of them with an additional
    # variable)
      ACC.set <- ACC.ref

    # Generate a data.frame with the relevant information
      if(tab){
      # Build a data.frame with the information
        EVOL <- data.frame(POS,NEG,ACC.ref, ACC.ref)
        colnames(EVOL) <- c("NUMERATOR","DENOMINATOR", "ACC", "Increase")

      # Index of factor (character) columns
        f.indx <- sapply(EVOL, is.factor)
      }
    # Index of the number of variables for the balance
      nV <- 2


    #------------------------------#
    # For numeric responses, . . .
    #------------------------------#

      if (classy=="numeric"){

      # While there is an improvement and the maximun number of variables
      # has not been reached, . . .
        while (ACC.set <= ACC.ref && length(rem.nam)!=0 && nV<maxV){


        # The new p.bal.ref is the p.set of the previous step
          ACC.ref <- ACC.set

        # Function to extract the p-value
          add2bal.ACC <- matrix(0, nrow = 0, ncol = 2)

        # Solve the problem with libraries
          suppressWarnings(suppressMessages(library("CMA")))
          suppressMessages(detach("package:CMA", unload=TRUE))
          suppressMessages(library(pROC))


        # Extract the p-values
          add2bal.ACC <- t(sapply(rem.nam, function(x)
            balance.adding.cor(x, LogCounts = logCounts, POS, NEG, numy = numy,
                               covar = covar)))
        # Add names to the rows
          row.names(add2bal.ACC) <- rem.nam


        # Extract which is the variable (only the first row)
          ACC.opt <- which(add2bal.ACC==min(add2bal.ACC),arr.ind = T)[1,]
        # Modify p.set
          ACC.set <- min(add2bal.ACC)

        # If there is an improvement, . . .
          #if (abs(ACC.set - ACC.ref) > th.imp){
          if ((ACC.set - ACC.ref) < th.imp){
            INC.VAR <- c(INC.VAR, rem.nam[ACC.opt[1]])
            ACC.Bal <- c(ACC.Bal, ACC.set)
            nV <- nV + 1
            if (ACC.opt[2]==1){
              POS <- c(POS, rem.nam[ACC.opt[1]])
              if(tab){
                EVOL[f.indx] <- lapply(EVOL[f.indx], as.character)
                EVOL <- rbind(EVOL, c(rem.nam[ACC.opt[1]], "-", ACC.set,
                                    ACC.set - ACC.ref))}
            } else if (ACC.opt[2]==2){
              NEG <- c(NEG, rem.nam[ACC.opt[1]])
              if(tab){
                EVOL[f.indx] <- lapply(EVOL[f.indx], as.character)
                EVOL <- rbind(EVOL, c("-", rem.nam[ACC.opt[1]], ACC.set,
                                    ACC.set - ACC.ref))}
          }
        } else {ACC.set <- 0 }

        # Remainng variables (possible to add to the balance)
          rem.nam <- rem.nam[-ACC.opt[1]]

      } # End while
    }else{
      #-----------------------------------#
      # For non-numeric responses, . . .
      #-----------------------------------#

      # While there is an improvement and the maximun number of variables has not
      # been reached, . . .
        while (ACC.set >= ACC.ref && length(rem.nam)!=0 && nV<maxV){

        # The new p.bal.ref is the p.set of the previous step
          ACC.ref <- ACC.set

        # Function to extract the p-value
          add2bal.ACC <- matrix(0, nrow = 0, ncol = 2)

        # Solve the problem with libraries
          suppressWarnings(suppressMessages(library("CMA")))
          suppressMessages(detach("package:CMA", unload=TRUE))
          suppressMessages(library(pROC))

        # Extract the p-values
          add2bal.ACC <- t(sapply(rem.nam, function(x)
            balance.adding.cor(x, LogCounts = logCounts, POS, NEG, numy = numy,
                               covar = covar)))
        # Add names to the rows
          row.names(add2bal.ACC) <- rem.nam

        # Extract which is the variable (only the first row)
          ACC.opt <- which(add2bal.ACC==max(add2bal.ACC),arr.ind = T)[1,]
        # Modify p.set
          ACC.set <- max(add2bal.ACC)

        # If there is an improvement, . . .
          if ((ACC.set - ACC.ref) > th.imp){
            # Add the included variable
              INC.VAR <- c(INC.VAR, rem.nam[ACC.opt[1]])
            # Add the Accuracy value
              ACC.Bal <- c(ACC.Bal, ACC.set)
            # A new variable into the balance
              nV <- nV + 1
            # Variable included into NUMERATOR or DENOMINATOR?
              if (ACC.opt[2]==1){
                POS <- c(POS, rem.nam[ACC.opt[1]])
                if(tab){
                  EVOL[f.indx] <- lapply(EVOL[f.indx], as.character)
                  EVOL <- rbind(EVOL, c(rem.nam[ACC.opt[1]], "-", ACC.set,
                                      ACC.set - ACC.ref))}
              } else if (ACC.opt[2]==2){
                NEG <- c(NEG, rem.nam[ACC.opt[1]])
                if(tab){
                  EVOL[f.indx] <- lapply(EVOL[f.indx], as.character)
                  EVOL <- rbind(EVOL, c("-", rem.nam[ACC.opt[1]], ACC.set,
                                        ACC.set - ACC.ref))}
              }
          } else {ACC.set <- 0 }

        # Remainng variables (possible to add to the balance)
          rem.nam <- rem.nam[-ACC.opt[1]]

        }
      }

    # K1 and k2
      k1 <- length(POS); k2 <- length(NEG)

    # The final Balance
      FINAL.BAL <- sqrt((k1*k2)/(k1+k2))*
        (rowM(logCounts[,POS])- rowM(logCounts[,NEG]))

  # Draw the plot if draw == T
  if (draw){
  #----------------------------------------------------------------------------#
  # GRAPHICAL REPRESENTATION
  #----------------------------------------------------------------------------#

    # Load library
      library(ggplot2)

    #-----------------------------------------#
    # FIRST: The names of included variables
    #-----------------------------------------#

    # Variables included
      T1 <- c("NUMERATOR", POS)
      T2 <- c("DENOMINATOR",NEG)

    # Parameter to specify the limits for writting
      yl <- max(length(T1), length(T2)) + .5

    # Empty plot with text
      df <- data.frame()
      Imp.table <- ggplot(df) + xlim(0, 100) + ylim(-0.5, yl) + theme_void() +
        annotate("text",
                 x = 75,
                 y = floor(yl):ceiling(yl-length(T1)),
                 label = T1,
                 colour = c("royalblue1",rep("black",length(T1)-1)),
                 fontface = 2) +
        annotate("text",
                 x = 25,
                 y = floor(yl):ceiling(yl-length(T2)),
                 label = T2,
                 colour = c("royalblue1",rep("black",length(T2)-1)),
                 fontface = 2)

      # Parameter 2 to specify the limits for writting
      T1 <- c(POS,"A")
      T2 <- c(NEG,"B")
      yl2 <- max(length(T1), length(T2));
      yl2 <- yl2*1.05;
      escal <- 3;
      bot <- 5*escal*yl2/100;

      # Empty plot 2 with text
      ndiv = max(3,floor(yl2))+1;
      lineh <- 0 # 0.5*ceiling(yl-length(T1));
      df2 <- data.frame()
      colbalance<-colbalance<-"brown3"
      Imp.table2 <- ggplot(df2) + xlim(0, 100) + ylim(-bot, ndiv) + theme_void() +
        geom_segment(aes(x = 10, y = lineh, xend = 90, yend = lineh),color=colbalance, size=1.3) +
        geom_segment(aes(x = 50-escal, y = lineh-bot, xend = 50+escal, yend = lineh-bot),color=colbalance, size=1) +
        geom_segment(aes(x = 50+escal, y = lineh-bot, xend = 50+escal, yend = lineh),color=colbalance, size=1) +
        geom_segment(aes(x = 50-escal, y = lineh-bot, xend = 50-escal, yend = lineh),color=colbalance, size=1) +
        annotate("text",
                 x = 75,
                 y = seq(bot, floor(yl2), length.out = ndiv),
                 label = c(T1,rep("",ndiv-length(T1))),
                 colour = c(rep("black",length(T1)-1),rep(colbalance,ndiv-length(T1)+1)),
                 fontface = 2) +
        annotate("text",
                 x = 25,
                 y = seq(bot, floor(yl2), length.out = ndiv),
                 label = c(T2,rep("",ndiv-length(T2))),
                 colour = c(rep("black",length(T2)-1),rep(colbalance,ndiv-length(T2)+1)),
                 fontface = 2)


    #-----------------------------------------#
    # SECOND: The representation of the plots
    #-----------------------------------------#

    # Auxiliar data.frame for graphical representation
      U <- data.frame(dat, FINAL.BAL)
      colnames(U)[ncol(U)] <- "V1"
    # Regression model
      FIT.final <- glm(numy~., data=U, family = f.class)

    # The plot depending of the class of the response variable
      if (classy=="factor"){

      #---------------------------------------------------#
      # The composition of the final plot
      #---------------------------------------------------#

        # BOXPLOT 1
        BoxP <-  ggplot(U, aes(x=y, y=V1, fill=y)) +
          geom_boxplot(color="black", size=1) +
          scale_fill_manual(values=col) +
          theme_bw() +
          ylab("Balance") +
          xlab("Factor") +
          theme(legend.position = "none")
        # BOXPLOT 2
        BoxP2 <-  ggplot(U, aes(x=y, y=V1, fill=y)) +
          geom_boxplot(color="black", size=1) +
          scale_fill_manual(values=col) +
          theme_bw() +
          ylab("Balance") +
          xlab("") +
          theme(legend.position = "none")+coord_flip()

        # Density plot 1 for the balance
        ydensity <- ggplot(U, aes(V1, fill=y)) +
                    geom_density(alpha=.5, size=1.25) +
                    scale_fill_manual(values = col) +
                    theme_bw() + xlab("") + ylab("") +
                    theme(legend.position = "none",
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()) +
                    coord_flip()
        # Density plot 2 for the balance
        ydensity2 <- ggplot(U, aes(V1, fill=y)) +
          geom_density(alpha=.5, size=1) +
          scale_fill_manual(values = col) +
          theme_bw() + xlab("") + ylab("") +
          theme(legend.position = "none")

      # ROC - curve
        library(pROC)
      # Build ROC curve
        A<-roc(response = U$numy,predictor = FIT.final$fitted.values)
      # Extract the sensitivity and specificiti values
        ROC.TAB <- data.frame(x=1-A$specificities, y=A$sensitivities)
      # Order them for a correct representation
        ROC.TAB <- ROC.TAB[order(ROC.TAB$y),]
      # AUC value
        #auc.val <- round(logit.cor(FIT.final,y = U$numy, covar = covar, logit.acc = logit.acc),3)
        #auc.val<-round(as.numeric(auc(y, FIT.final$fitted.values)),3)   # diferencies
        auc.val<-round(as.numeric(auc(U$numy, FIT.final$fitted.values)),3)   # new diferencies
      # The plot
        ROC.plot <- ggplot(data=ROC.TAB, aes(x=x, y=y)) +
                    geom_line() +
                    ggtitle("ROC curve") +
                    xlab("FPR") + ylab("TPR") +
                    geom_step() +
                    annotate("text", x = .75, y = .2,
                              label = paste("AUC-ROC \n",auc.val), size = 2.5) +
                    theme_bw() +
                    theme(plot.title = element_text(hjust = 0.5))

      # Load libraries
        library("gridExtra")
        library("grid")
        FINAL.P <- arrangeGrob(Imp.table, ROC.plot, BoxP, ydensity,
                               ncol=2, nrow=2, widths=c(5,1.25), heights=c(2, 5),
                               vp=viewport(width=0.8, height=0.8))
        FINAL.P2 <- arrangeGrob(Imp.table2, BoxP2, ydensity2,
                               ncol=1, nrow=3, widths=4, heights=c(4, 4, 4),
                               vp=viewport(width=0.75, height=1))

      } else {

        # Fit the regression model
          FIT.p <- glm(y ~ FINAL.BAL, family = f.class)

          PLOT.G <- ggplot(U, aes(V1, y)) +
                    geom_point(colour = "black", size = 3) +
                    geom_abline(intercept=FIT.p$coefficients[1],
                                slope=FIT.p$coefficients[2], lwd=3, col="blue") +
                    theme_bw() +
                    xlab("Balance value") + ylab("Response variable") +
                    ggtitle("") +
                    theme(strip.text.x = element_text(size=12, angle=0,
                                                      face="bold",colour="white"),
                          strip.text.y = element_text(size=12, face="bold"),
                          strip.background = element_rect(colour="black",
                                                          fill="black"),
                    plot.title = element_text(size=20, vjust=2.25, hjust= .5,
                                              face = "bold"),
                    legend.title = element_text(face="bold"),
                    legend.text = element_text(face="bold"))

      # Load libraries
        library(grid)
        library(gridExtra)

        FINAL.P <- arrangeGrob(Imp.table, PLOT.G, nrow=2,
                               heights=c(0.2,0.5),vp=viewport(width=0.8,
                                                              height=0.8))
        FINAL.P2 <- arrangeGrob(Imp.table2, PLOT.G, nrow=2,
                                  heights=c(1,2),vp=viewport(width=0.8,
                                                                 height=0.8))
      }

      # Draw the plot if draw == T
      grid.draw(FINAL.P)
	}  #end plot

      if (classy=="numeric") ROC.plot=NULL

    # Round the values
    if(tab){
        EVOL[,3]<-round(as.numeric(EVOL[,3]),5)
        EVOL[,4]<-round(as.numeric(EVOL[,4]),5)
    		if (draw== TRUE){
    		  L <- list(FINAL.BAL, POS, NEG, INC.VAR, ACC.Bal, EVOL, global.plot=FINAL.P,
                      FIT.final,global.plot2 = FINAL.P2, ROC.plot = ROC.plot)
    		} else {
    			L <- list(FINAL.BAL, POS, NEG, INC.VAR, ACC.Bal, EVOL)
       	}
    } else {
    		if (draw== TRUE){
    			L <- list(FINAL.BAL, POS, NEG, INC.VAR, ACC.Bal, global.plot=FINAL.P,
                      FIT.final,global.plot2 = FINAL.P, ROC.plot = ROC.plot)
    		} else {
    			L <- list(FINAL.BAL, POS, NEG, INC.VAR, ACC.Bal)
       	}
    }

      return(L)
  }


#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# FUNCTION: selbal.cv
# INPUT:
#              x: matrix with counts
#              y: response variable
#         n.fold: number of folds for each iteration
#         n.iter: number of iterations
#           seed: seed for results replication (used for build the groups)
#         th.imp: threshold of improvement when adding a variable. If the
#                 improvement is not higher than this quantity, it stops.
#          covar: data.frame with covariates
#      logit.acc: the accuracy value to return when a logit regression is made
#    user_numVar: number of variables the user decides to include

# OUTPUT:
#           $NUM: variables included in the numerator (whole data set)
#           $DEN: variables included in the denominator (whole data set)
#       $ROB.TAB: table of "robustness" with the percentage of times that the
#                 variable appears in the cross - validation model
#    $GLOBAL.ACC: ACC for the whole data set after computing the balance
#------------------------------------------------------------------------------#

#' Cross - validation process for the selection of the optimal number of
#' variables and robustness evaluation
#'
#'
#' @param x a \code{matrix} object with the information of variables
#' (\emph{columns}) for each sample (\emph{rows}).
#' @param y the response variable, either continuous or dichotomous.
#' @param n.fold number of folds in which to divide the whole data set.
#' @param n.iter number of iterations for the cross - validation process.
#' @param seed a seed to make the results reproducible.
#' @param th.imp the minimum increment needed when adding a new variable into
#' the balance in order to consider an improvement.
#' @param covar \code{data.frame} with the variables to adjust for
#' (\emph{columns}).
#' @param col \code{vector} of two colours for differentiate the variables
#' appearing in the numerator and in the denominator of the balances.
#' @param col2 \code{vector} of three colours for the lines of the barplot with
#' the aim of identifying if each variable appears in the Global balance, in the
#' CV - balance or in both of them.
#' @param logit.acc when \code{y} is dichotomous, the measure to compute for
#' the correlation between \code{y} and the proposed \emph{balance}
#' adjusting for covariates. One of the following values: \code{"AUC"} (default),
#'  \code{"Dev"}, \code{"Rsq"} or \code{"Tjur"}.
#' @param maxV \code{numeric} value defining the maximum number of variables
#' composing the balance. Default 1e10 to give prevalence to \code{th.imp}
#' parameter.
#' @param zero.rep a value defining the method to use for zero - replacement.
#' \code{"bayes"} for BM-replacement or \code{"one"} to add one read tho each
#' cell of the matrix.
#' @param opt.cri parameter indicating the method to determine the optimal
#' number of variables. \code{"max"} to define this number as the number of
#' variables which maximizes the association value or \code{"1se"} to take also
#' the standard error into account.
#' @param user_numVar parameter to modify the choosen optimal number of variables.
#' If it is used, it is the final number of variables used in the method.
#'
#'
#' @return A \code{list} with the following objects:
#'
#' \itemize{
#' \item a boxplot with the mean squared errors (numeric responses) or AUC
#' values (dichotomous responses) for the test data sets using the balances
#' resulted in the cross - validation. Branches represent the standard error and
#' the optimal number of components according with the \code{opt.cri} criteria
#' is highlighted with a dashed line.
#' \item barplot with the proportion of times a variable appears in the
#' cross - validation balances.
#' \item a graphical representation of the Global Balance (draw it using
#' \code{grid.draw} function).
#' \item a table with the infromation of Global Balance, CV Balance and the
#' three most repeated balances in the cross - validation process (draw it using
#' \code{plot.tab} function).
#' \item a vector with the accuracy values (MSE for continuous variables and
#' AUC for dichotomous variables) obtained in the cross - validation procedure.
#' \item a table with the variables appearing in the Global Balance in a useful
#' format for \code{bal.value} function in order to get the balance score for
#' new datasets.
#' \item the regression model object where the covariates and the final balance
#' are the explanatory variables and \code{y} the response variable.
#' \item the optimal number of variables estimated in the cross -
#' validation.
#'
#' }
#'
#' @examples
#' # Load data set
#'   load("HIV.rda")
#' # Define x and y
#'   x <- HIV[,1:60]
#'   y <- HIV[,62]
#' # Run the algorithm
#'   CV.Bal <- selbal.cv(x,y)
#' @export selbal.cv




  selbal.cv <- function(x, y, n.fold = 5, n.iter = 10, seed = 31415,
                        covar = NULL, col = c("steelblue1", "tomato1"),
                        col2 = c("darkgreen", "steelblue4","tan1"),
                        logit.acc = "AUC", maxV = 20, zero.rep = "bayes",
                        opt.cri = "1se", user_numVar = NULL){

    # Load package plyr
    suppressMessages(library(plyr))


    #------------------------------------------------------------------------------#

    #----------------------------------------------------------------------------#
    # STEP 0: build the necessary objects
    #----------------------------------------------------------------------------#

    #-----------------------------------------------#
    # 0.1: Build necessary objects for the function
    #-----------------------------------------------#

    # Class of the response
    classy <- class(y)
    # Family for the glm (default: gaussian)
    f.class <- "gaussian"
    # numy to be y and it will be modified if y is a factor
    numy <- y
    # If y is a factor, compute the number of levels of y
    if (classy == "factor"){
      ylev <- levels(y)
      numy <- as.numeric(y) - 1
      f.class <- "binomial"
    }

    # Names of x
    x.nam <- colnames(x)

    # ROB.TAB
    ROB.TAB <- matrix(0, nrow = ncol(x), ncol=3)
    colnames(ROB.TAB) <- c("Prop. Included", "Prop_Numerator",
                           "Prop_Denominator")
    row.names(ROB.TAB) <- x.nam

    # BAL.resume
    BAL.resume <- matrix(0,nrow =n.fold*n.iter , ncol = ncol(x))
    colnames(BAL.resume) <- colnames(x)

    # Build a table with the response variable and covariates for correction
    if (!is.null(covar)){ dat <- data.frame(cbind(numy, covar))
    } else { dat <-data.frame(numy)}

    # Message starting the algorithm
    cat(paste("\n\n###############################################################",
              "\n STARTING selbal.cv FUNCTION",
              "\n###############################################################"))
    # Log-transformed counts (with zero replacement
    cat(paste(
      "\n\n#-------------------------------------------------------------#",
      "\n# ZERO REPLACEMENT . . .\n\n"))

    # Define log-transformed data with the zero-replacement made
    logc <- log(cmultRepl2(x, zero.rep = zero.rep))


    cat(paste("\n, . . . FINISHED.",
              "\n#-------------------------------------------------------------#"))

    #---------------------------------------#
    # 0.2: Define cross - validation groups
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


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Function for cross - validation
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

    cv.MSE <- function(k){

      #-------------------------------------#
      # Necessary matrices for the function
      #-------------------------------------#

      # Folds associated to a "complete FOLD"
      CV.group<-CV.groups[((k-1)*n.fold + 1):(k*n.fold),]

      # Build objects
      # Table with all the variables included
      Bal.List <- list()
      # Matrix for MSE values (or ACC "Accuracy")
      ACC.Mat <- matrix(0,nrow = maxV-1, ncol = nrow(CV.group))


      # For each fold
      for (i in 1:nrow(CV.group)){

        # Define training data.set
        train.idx<-CV.group[i,]
        # Training dataset (x, y and covar)
        x.train<-logc[train.idx,]
        x.test<-logc[-train.idx,]
        y.train<-y[train.idx]
        y.test<-y[-train.idx]
        covar.train<-covar[train.idx,]
        covar.test<-covar[-train.idx,]

        # Compute the balances for the training data set
        BAL <- selbal.aux(x.train, y.train, th.imp = 0, covar = covar.train,
                          logit.acc, logt=F, maxV = maxV)

        # Variables included in the balance (as NUMERATOR | DENOMINATOR)
        Bal.List[[i]]<-BAL

        # A matrix for predictions for Y
        PRED.y <- matrix(0, nrow = length(y.test), ncol = nrow(BAL) - 1)

        # For each number of variables (2:nrow(BAL))
        for (l in 2:min(maxV,nrow(BAL))){
          # Data frame for train data
          df <-data.frame(Y = y.train, B = bal.value(BAL[1:l,],x.train))

          # Data frame for test data
          df.test <- data.frame(Y = y.test, B = bal.value(BAL[1:l,],x.test))

          if(!is.null(covar)){
            df <- cbind(df,cov=covar.train)
            df.test <- cbind(df.test,cov=covar.test)
          }

          # Regression model for test data
          FIT <- glm(Y ~ ., data= df, family = f.class)
          # Predictions
          PRED.y[,l-1] <- predict(FIT,df.test, type="response")
        }

        #------------------------------------------------------------------------------#
        # FUNCTION: Measure the error value
        #------------------------------------------------------------------------------#

        ACC.eval <- function(y, pred, classy, logit.acc=NULL){

          if (classy == "numeric"){
            ACC <- apply(pred, 2, function(x) mean((y-x)^2))
          }else{
            if (logit.acc == "AUC"){
              # Load library
                library(pROC)
                ACC <- apply(pred, 2, function(x) auc(y,x))
            } else if(logit.acc == "Rsq"){
                ACC <- apply(pred, 2, function(x) cor (y, x)^2)
            } else if (logit.acc == "Tjur"){
                ACC <- apply(pred, 2, function(x) mean(x[y==1]) - mean(x[y==0]))
            } else if (logit.acc == "Dev"){
				ACC<-apply(pred, 2, function(x) 1-(deviance(glm(y ~ x, data= df, family = binomial()))/glm(y~1, family=binomial())[[10]]) )  # proportion of explained deviance
			}

          }

        return(ACC)
      }

        #------------------------------------------------------------------------------#

        # Run ACC.eval function
          R <- ACC.eval(numy[-train.idx], PRED.y, classy=classy, logit.acc)
        # Add the information
          ACC.Mat[,i] <- c(R, rep(R[length(R)], maxV-length(R)-1))
      } # End of i

      return(list(Bal.List, ACC.Mat))

    }

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
    comb <- function(x, ...) {
      lapply(seq_along(x),
             function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }

    ################################################################################

    # CV - procedure computed in parallel
    INTEREST <- foreach(h=1:n.iter,
                        .export=c("logit.cor", "rowM","selbal.aux", "bal.value",
                                  "logit.acc", "cmultRepl","cmultRepl2"),
                        .combine='comb',
                        .multicombine=TRUE,
                        .init=list(list(), list())) %dopar% {
                          cv.MSE(h)
                        }
    # Stop the parallelization
    stopImplicitCluster()

    cat(paste("\n . . . finished.",
              "\n#-------------------------------------------------------------#",
              "\n###############################################################"))

    #------------------------------------------------------------------------------#

    # Rebuild the objects from INTEREST
      Balances <- unlist(INTEREST[[1]], recursive = F)
      ACC.Matrix <- do.call(cbind,INTEREST[[2]])

    # ACC mean values for each number of variables
      ACC.mean <- apply(ACC.Matrix,1,mean)
      ACC.se <- apply(ACC.Matrix,1,function(x) sd(x)/sqrt(length(x)))

    if(classy == "numeric"){
      # Define the minimum mean value
        m <- which.min(ACC.mean)
      # The minimum value under ACC.mean[m] + SE
        if(length(which((ACC.mean<(ACC.mean[m] + ACC.se[m]))==T))>0){
        mv <- min(which((ACC.mean<(ACC.mean[m] + ACC.se[m]))==T))
        } else {mv<-m}

      # Depending on "opt.cri":
        if (opt.cri == "1se"){opt.M <- mv + 1
        }else {opt.M <- m + 1}
      }else{
        # Define the maximum ACC value
          m <- which.max(ACC.mean)
        # The minimum value whith ACC.mean over ACC.mean[m] - ACC.sd[m]
          if(length(which((ACC.mean>(ACC.mean[m] - ACC.se[m]))==T))>0){
          mv <- min(which((ACC.mean>(ACC.mean[m] - ACC.se[m]))==T))
          } else {mv<-m}

      # Depending on "opt.cri":
      if (opt.cri == "1se"){opt.M <- mv + 1
      }else { opt.M <- m + 1}
    }

  # Print a message indicating the number of optimal variables
    cat(paste("\n\n The optimal number of variables is:", opt.M, "\n\n"))

	if (!is.null(user_numVar)){
		opt.M <- user_numVar;
	}


    # Define NUM and DEN according to opt.M
    suppressMessages(BAL <- selbal.aux(x, y, th.imp = 0, covar = covar,
                                       logit.acc, logt=T, maxV = opt.M)) #diferencies logit.acc=logit.acc
  # Variables in the NUMERATOR and the DENOMINATOR
    NUM <- BAL[BAL[,2]=="NUM","Taxa"]
    DEN <- BAL[BAL[,2]=="DEN","Taxa"]

    # Information about the GLM
    # Final balance (number of components)
    k1 <- length(NUM); k2 <- length(DEN)

    # The final Balance
    FINAL.BAL <- sqrt((k1*k2)/(k1+k2))*(rowM(logc[,NUM])- rowM(logc[,DEN]))

    # Auxiliar data.frame for graphical representation
    U <- data.frame(dat, FINAL.BAL)
    colnames(U)[ncol(U)] <- "V1"
    # Regression model
    FIT.final <- glm(numy~., data=U, family = f.class)


    #------------------------------------------------------------------------------#
    #                           GRAPHICAL REPRESENTATION
    #------------------------------------------------------------------------------#

    # Build a data.frame with the information
    if(classy=="numeric"){
      df.boxplot <- data.frame(mean = ACC.mean, se = ACC.se, n =2:maxV)
      # Load library
      library(ggplot2)
      # The plot
      MSE.Boxplot <- ggplot(df.boxplot, aes(x=n, y=mean)) +
        geom_errorbar(aes(ymin=mean - se, ymax= mean + se),
                      width = 0.2, col = "gray") +
        geom_vline(xintercept = opt.M, linetype = "dotted",
                   col = "blue") +
        geom_point(color = "red", lwd=2) +
        theme_bw() +
        xlab("Number of variables") +
        ylab("Mean-Squared Error") +
        scale_x_continuous(breaks=seq(2,maxV,1)) +
        theme(strip.text.x = element_text(size=12, angle=0,
                                          face="bold",colour="white"),
              strip.text.y = element_text(size=12, face="bold"),
              strip.background = element_rect(colour="black",
                                              fill="black"),
              plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                        face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

    }else{
      df.boxplot <- data.frame(mean = ACC.mean, se = ACC.se, n =2:maxV)
	  ylabelName="Accuracy (AUC)";
	  if (logit.acc=="Dev"){
		ylabelName="Explained Deviance";
	  }
      # Load library
      library(ggplot2)
      # The plot
      MSE.Boxplot <- ggplot(df.boxplot, aes(x=n, y=mean)) +
        geom_errorbar(aes(ymin=mean - se, ymax= mean + se),
                      width = 0.2, col = "gray") +
        geom_vline(xintercept = opt.M, linetype = "dotted",
                   col = "blue") +
        geom_point(color = "red", lwd=2) +
        theme_bw() +
        xlab("Number of variables") +
        ylab(ylabelName) +
        scale_x_continuous(breaks=seq(2,maxV,1)) +
        theme(strip.text.x = element_text(size=12, angle=0,
                                          face="bold",colour="white"),
              strip.text.y = element_text(size=12, face="bold"),
              strip.background = element_rect(colour="black",
                                              fill="black"),
              plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                        face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }



    # Extract information about the variables selected in the balances
    Sub.Balances <- lapply(Balances, function(x) x[1:min(nrow(x),opt.M),])
    # Complete the matrix with the variables appearing
    # Build the matrix
    BR<-matrix(0,nrow =n.fold*n.iter , ncol = length(x.nam))
    colnames(BR) <- x.nam
    # Complete the row of BAL.resume
    for (i in 1:length(Sub.Balances)){
      BR[ i, colnames(BR) %in% Sub.Balances[[i]][Sub.Balances[[i]][,"Group"]=="NUM","Taxa"]] <- 1
      BR[ i, colnames(BR) %in% Sub.Balances[[i]][Sub.Balances[[i]][,"Group"]=="DEN","Taxa"]] <- 2
    }

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


    #------------------------------------------------------------------------------#
    # GRAPHICAL REPRESENTATION OF MAIN BALANCES: global balance and cv - balance
    #------------------------------------------------------------------------------#

    # Define y for the plot
    ifelse(classy %in% c("numeric","integer"),
           y.plot <- y,
           y.plot <- 1:length(y))

    # PLOT GLOBAL
    PLOT.Global <- plot.bal(NUM,DEN,logc,y, covar, col = col, logit.acc)

    # Message starting the algorithm
    cat(paste("\n\n###############################################################",
              "\n . . . FINISHED.",
              "\n###############################################################"))

    # Build a list with the elements of interest
    L <- list(accuracy.nvar = MSE.Boxplot,
              var.barplot = IMP.plot,
              global.plot = PLOT.Global$Global.plot,
              global.plot2 = PLOT.Global$Global.plot2,
              ROC.plot = PLOT.Global$ROC.plot,
              cv.tab = dat,
              cv.accuracy = ACC.Matrix[(opt.M - 1),],
              global.balance = BAL,
              glm = FIT.final,
              opt.nvar = opt.M)

    return(L)

  }


#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# FUNCTION: selbal.aux
# INPUT:
#              x: matrix with counts
#              y: response variable
#         th.imp: threshold of improvement when adding a variable. If the
#                 improvement is not higher than this quantity, it stops.
#          covar: data.frame with covariates
#      logit.acc: the accuracy value to return when a logit regression is made
#           logt: should data be log-transformed?
#           maxV: number of maximum variables for the balance

# OUTPUT:
#        Tab.var: a table including the variables included in the balance in the
# order they were added.
#------------------------------------------------------------------------------#

#' The balance defined in the order its variables were selected
#'
#'
#' @param x a \code{matrix} object with the information of variables
#' (\emph{columns}) for each sample (\emph{rows}).
#' @param y the response variable, either continuous or dichotomous.
#' @param th.imp the minimum increment needed when adding a new variable into
#' the balance, in order to consider an improvement.
#' @param covar \code{data.frame} with the variables to adjust for
#' (\emph{columns}).
#' @param logit.acc when \code{y} is dichotomous, the measure to compute for
#' the correlation between \code{y} and the proposed \emph{balance}
#' adjusting for covariates. One of the following values: \code{"Rsq"} (default),
#'  \code{"AUC"} or \code{"Tjur"}.
#' @param logt a logical value indicating if the data needs to be
#' log-transformed.
#' @param maxV the number maximum of variables defining the balance.
#'
#' @return A \code{data.frame} with the following objects:
#'
#' \itemize{
#' \item \code{Taxa} the name of the taxa added into the balance.
#' \item \code{Group} \emph{NUM} if the variable is added into the numerator
#' and \emph{DEN} if it is added in the denominator.
#'
#' }
#'
#' @examples
#' # Load data set
#'   load("HIV.rda")
#' # Define x and y
#'   x <- HIV[,1:60]
#'   y <- HIV[,62]
#' # Run the algorithm
#'   Bal <- selbal.aux(x,y)
#' @export selbal.aux



  selbal.aux <- function(x, y, th.imp = 0, covar = NULL, logit.acc="AUC",
                         logt=T, maxV = 1e10, zero.rep = "bayes"){

    #--------------------------------------------------------------------------#
    # STEP 0: load libraries and extract information
    #--------------------------------------------------------------------------#

    #----------------------------------------------#
    # 0.1: information about the response variable
    #----------------------------------------------#

    # Class of the response variable
    classy <- class(y)
    # Family for the glm (default: gaussian)
    f.class <- "gaussian"
    # numy to be y and it will be modified if y is a factor
    numy <- y
    # If y is a factor, compute the number of levels of y
    if (classy == "factor"){
      ylev <- levels(y)
      numy <- as.numeric(y) - 1
      f.class <- "binomial"
    }

    #------------------------------------------------------------------#
    # 0.2: information and transformation of the independent variables
    #------------------------------------------------------------------#

    # Load library
    suppressMessages(library(zCompositions))

    # Variables name
    var.nam <- rem.nam <- colnames(x)

    # Build a table with the response variable and covariates for correction
    if (!is.null(covar)){ dat <- data.frame(cbind(numy, covar))
    } else { dat <-data.frame(numy)}



    # The logCounts (with zero replacement)
    if (logt == F){ logCounts <- x
    } else{
      logCounts <- log(cmultRepl2(x, zero.rep = zero.rep))
    }


    #--------------------------------------------------------------------------#
    # 0.3: auxiliar functions
    #--------------------------------------------------------------------------#

    #--------------------------------------------------------------------------#
    # Auxiliar function to compute the first balance
    #--------------------------------------------------------------------------#

    first.bal<- function(logCounts, Y, covar=NULL){

      #------------------------------------------------------------------------#
      # STEP 0: extract information
      #------------------------------------------------------------------------#

      # Number and name of variables
      n <- ncol(logCounts)
      nam <- colnames(logCounts)


      #------------------------------------------------------------------------#
      # STEP 1: build the output matrix
      #------------------------------------------------------------------------#

      # The matrix
        if (classy=="factor"){ M<-matrix(0, nrow=n, ncol=n)
       }else{ M<-matrix(1e99, nrow=n, ncol=n)}
      # Row.names and colnames
      row.names(M)<-colnames(M)<-nam

      #------------------------------------------------------------------------#
      # STEP 2: complete the matrix
      #------------------------------------------------------------------------#

      if(classy=="factor"){

        # Solve the problem with libraries
        suppressWarnings(suppressMessages(library("CMA")))
        suppressMessages(detach("package:CMA", unload=TRUE))
        suppressMessages(library(pROC))

        for (i in 2:n){
          for (j in 1:(i-1)){
            # Build a table with the information
            TAB <- data.frame(cbind(Y,logCounts[,i]-logCounts[,j],covar))
            # Fit the regression model
            FIT <- glm(Y ~ .,data=TAB, family = f.class)

            # Add the value into the matrix
            ifelse(FIT$coefficients[2]>0,
                   M[i,j] <- logit.cor(FIT, y = y, covar = covar, logit.acc),#diferencies Y  covar
                   M[j,i] <- logit.cor(FIT, y = y, covar = covar, logit.acc)) #diferencies Y  covar

          } # End j
        } # End i

        # Indices for the highest logit.cor value
        r <- which(M == max(M), arr.ind = T)


      } else {
        for (i in 2:n){
          for (j in 1:(i-1)){
            # Build a table with the information
            TAB <- data.frame(cbind(Y,logCounts[,i]-logCounts[,j],covar))
            # Fit the regression model
            FIT <- glm(Y ~ .,data=TAB, family = f.class)
            # Complete the matrix
            ifelse(FIT$coefficients[2]>0,
                   M[i,j] <- mean(FIT$residuals^2),
                   M[j,i] <- mean(FIT$residuals^2))
          } # End j
        } # End i

        # Indices for the lowest MSE value
        r <- which(M == min(M), arr.ind = T)
      }


      # Return the row and column of the maximum value
      return(r)
    }

    #--------------------------------------------------------------------------#

    #--------------------------------------------------------------------------#
    # Auxiliar function to compute the "association value" when adding a new
    # variable into the balance
    #--------------------------------------------------------------------------#

    balance.adding.cor <- function(x, LogCounts, POS, NEG, numy, covar=NULL){

      #----------------------------------------#
      # If x added into the numerator, . .
      #----------------------------------------#

      # The "numerator"
      S1.pos <- rowM(LogCounts[,c(POS,x)]); s1 <- length(POS) + 1
      # The "denominator"
      S2.pos <- rowM(LogCounts[,NEG])     ; s2 <- length(NEG)
      # The balance
      BAL <- sqrt((s1*s2)/(s1+s2))*(S1.pos - S2.pos)

      # Data.frame with the variables
      D.pos <- data.frame(cbind(numy, BAL, covar))

      # Regression model
      FIT.pos <- glm(numy~., data=D.pos, family=f.class)
      # The MSE or the corresponding value for dichotomous responses
      if(classy=="numeric"){ C.pos <- mean(FIT.pos$residuals^2)
      #}else{ C.pos <- logit.cor(FIT.pos,numy,covar = covar, logit.acc)}#diferencies covar
      }else{ C.pos <- logit.cor(FIT.pos,y,covar = covar, logit.acc)}#diferencies covar

      #----------------------------------------#
      # If x added into the numerator, . .
      #----------------------------------------#

      # The numerator
      S1.neg <- rowM(LogCounts[,POS])       ; s1 <- length(POS)
      # The denominator
      S2.neg <- rowM(LogCounts[,c(NEG,x)])  ; s2 <- length(NEG) + 1
      # The balance
      BAL <- sqrt((s1*s2)/(s1+s2))*(S1.neg - S2.neg)

      # Data.frame with the variables
      D.neg <- data.frame(cbind(numy, BAL, covar))

      # Regression model
      FIT.neg <- glm(numy~., data=D.neg, family=f.class)
      # The MSE or the corresponding value for dichotomous responses
      if(classy=="numeric"){ C.neg <- mean(FIT.neg$residuals^2)
      #}else{ C.neg <- logit.cor(FIT.neg,numy,covar = covar, logit.acc)}  #diferencies covar
      }else{ C.neg <- logit.cor(FIT.neg,y,covar = covar, logit.acc)}  #diferencies covar
      # Correlation values
      COR <- c(C.pos, C.neg)
      # Return the values
      return(COR)
    }

#------------------------------------------------------------------------------#



    #--------------------------------------------------------------------------#
    # STEP 1: depending on the response variable class, . . .
    #--------------------------------------------------------------------------#

    # Define the first balance
    A1 <- first.bal(logCounts, Y = numy, covar=covar)
    # Variables taking parti into the first balance
    POS <- colnames(x)[A1[1,1]]
    NEG <- colnames(x)[A1[1,2]]

    # Included variables in the model
    INC.VAR <- c(POS, NEG)
    # Delete these variables from rem.nam
    rem.nam <- setdiff(rem.nam, INC.VAR)

    # Define the initial balance (B_1)
    S1 <- logCounts[,POS]
    S2 <- logCounts[,NEG]
    # Assign the values to B
    B <- sqrt(1/2)*(S1 - S2)

    #--------------------------------------------------------------------------#
    # Information about the ACC for the Balance values
    #--------------------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    # NEW: A table with the included variables and the group
    #--------------------------------------------------------------------------#

    Tab.var<- data.frame(Taxa = c(POS,NEG), Group = c("NUM", "DEN"))
    Tab.var[,1]<-as.character(Tab.var[,1])


    # Build a new data.frame
    dat.ini <- cbind(dat, B)
    # Fit the regression model
    FIT.initial <- glm(numy ~ .,data=dat.ini, family = f.class)

    # Solve the problem with libraries
    suppressWarnings(suppressMessages(library("CMA")))
    suppressMessages(detach("package:CMA", unload=TRUE))
    suppressMessages(library(pROC))

    # Define the initial "accuracy" or "association" value
    if(classy=="numeric"){ ACC.Bal <- mean(FIT.initial$residuals^2)
    #}else{ ACC.Bal <- logit.cor(FIT.initial, numy, covar = covar, logit.acc)}#diferencies covar
    }else{ ACC.Bal <- logit.cor(FIT.initial, y, covar = covar, logit.acc)}#diferencies covar


    #----------------------------------------------------------------------------#

    # ACC reference
    ACC.ref <- ACC.Bal

    #------------------------------------------------------------------------#
    # Improve the balances
    #------------------------------------------------------------------------#
    # Define some parameters
    # The p.value to compare 2 balances (one of them with an additional
    # variable)
    ACC.set <- ACC.ref

    # Index of the number of variables for the balance
    nV <- 2


    #------------------------------#
    # For numeric responses, . . .
    #------------------------------#

    if (classy=="numeric"){

      # While there is an improvement and the maximun number of variables has not
      # been reached, . . .
      while (ACC.set <= ACC.ref && length(rem.nam)!=0 && nV<maxV){

        # The new p.bal.ref is the p.set of the previous step
        ACC.ref <- ACC.set

        # Function to extract the p-value
        add2bal.ACC <- matrix(,nrow = 0, ncol = 2)

        # Solve the problem with libraries
        suppressWarnings(suppressMessages(library("CMA")))
        suppressMessages(detach("package:CMA", unload=TRUE))
        suppressMessages(library(pROC))


        # Extract the p-values
        add2bal.ACC <- t(sapply(rem.nam, function(x)
          balance.adding.cor(x, LogCounts = logCounts, POS, NEG, numy = numy,
                             covar = covar)))
        # Add names to the rows
        row.names(add2bal.ACC) <- rem.nam


        # Extract which is the variable (only the first row)
        ACC.opt <- which(add2bal.ACC==min(add2bal.ACC),arr.ind = T)[1,]
        # Modify p.set
        ACC.set <- min(add2bal.ACC)


        # If there is an improvement, . . .
        #if (abs(ACC.set - ACC.ref) > th.imp){
        if ((ACC.set - ACC.ref) < th.imp){
          INC.VAR <- c(INC.VAR, rem.nam[ACC.opt[1]])
          nV <- nV + 1
          if (ACC.opt[2]==1){
            POS <- c(POS, rem.nam[ACC.opt[1]])
            Tab.var <- rbind(Tab.var, c(rem.nam[ACC.opt[1]], "NUM"))
          } else if (ACC.opt[2]==2){
            NEG <- c(NEG, rem.nam[ACC.opt[1]])
            Tab.var <- rbind(Tab.var, c(rem.nam[ACC.opt[1]], "DEN"))
          } else {ACC.set <- 0 }

          # Remainng variables (possible to add to the balance)
          rem.nam <- rem.nam[-ACC.opt[1]]
        }

        } # End while

      }else{

        #-----------------------------------#
        # For non-numeric responses, . . .
        #-----------------------------------#

        # While there is an improvement and the maximun number of variables has not
        # been reached, . . .
        while (ACC.set >= ACC.ref && length(rem.nam)!=0 && nV<maxV){


          # The new p.bal.ref is the p.set of the previous step
          ACC.ref <- ACC.set

          # Function to extract the p-value
          add2bal.ACC <- matrix(,nrow = 0, ncol = 2)

          # Solve the problem with libraries
          suppressWarnings(suppressMessages(library("CMA")))
          suppressMessages(detach("package:CMA", unload=TRUE))
          suppressMessages(library(pROC))


          # Extract the p-values
          add2bal.ACC <- t(sapply(rem.nam, function(x)
            balance.adding.cor(x, LogCounts = logCounts, POS, NEG, numy = numy,
                               covar = covar)))
          # Add names to the rows
          row.names(add2bal.ACC) <- rem.nam


          # Extract which is the variable (only the first row)
          ACC.opt <- which(add2bal.ACC==max(add2bal.ACC),arr.ind = T)[1,]
          # Modify p.set
          ACC.set <- max(add2bal.ACC)

          # If there is an improvement, . . .
          if ((ACC.set - ACC.ref) > th.imp){
            INC.VAR <- c(INC.VAR, rem.nam[ACC.opt[1]])
            nV <- nV + 1
            if (ACC.opt[2]==1){
              POS <- c(POS, rem.nam[ACC.opt[1]])
              Tab.var <- rbind(Tab.var, c(rem.nam[ACC.opt[1]], "NUM"))
            } else if (ACC.opt[2]==2){
              NEG <- c(NEG, rem.nam[ACC.opt[1]])
              Tab.var <- rbind(Tab.var, c(rem.nam[ACC.opt[1]], "DEN"))
            } else {ACC.set <- 0 }
          } # End if

          # Remainng variables (possible to add to the balance)
          rem.nam <- rem.nam[-ACC.opt[1]]

        }
      }

      # K1 and k2
      k1 <- length(POS); k2 <- length(NEG)

      # The final Balance
      FINAL.BAL <- sqrt((k1*k2)/(k1+k2))*
        (rowM(logCounts[,POS])- rowM(logCounts[,NEG]))


      return(Tab.var)

    }

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# FUNCTION: bal.value
# INPUT:
#
#        bal.tab: a table defining a balance as in the output of selbal.aux
#              x: matrix with counts

# OUTPUT:
#        Tab.var: a table including the variables included in the balance in the
# order they were added.
#------------------------------------------------------------------------------#

#' The balance value for a given matrix of counts using the variables appearing
#' in a \code{data.frame} like the output of \code{selbal.aux}
#'
#'
#' @param bal.tab a \code{data.frame} including the variables defining the
#' balance (like the output of \code{selbal.aux} or the sixth element in the
#' output of \code{selbal.cv} function).
#' @param x the matrix with the log-transformed counts for a given subset of
#' individuals.
#'
#' @return A \code{vector} with the balance values for each subject.
#'
#'
#' @examples
#' # Load data set
#'   load("HIV.rda")
#' # Define x and y
#'   x <- HIV[,1:60]
#'   y <- HIV[,62]
#' # Run the algorithm
#'   Bal <- selbal.aux(x,y)
#' # Balance values for the individuals (log-transformed x values with the corresponding zero-replacement)
#'   bal.value(Bal,log(cmultrepl2(x)))
#'
#' @export bal.value


  bal.value <- function(bal.tab, x){
    # Variables in the numerator and the denominator
      vNUM <- bal.tab[bal.tab[,2]=="NUM",1]; k1 <- length(vNUM)
      vDEN <- bal.tab[bal.tab[,2]=="DEN",1]; k2 <- length(vDEN)
    # Value for the balance
      bal <- sqrt(k1*k2/(k1+k2))*(rowM(x[,vNUM]) - rowM(x[,vDEN]))
    # Return the value
      return(bal)
  }

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# FUNCTION: plot.bal
# INPUT:
#
#        POS: variables taking part in the numerator of the balance
#        NEG: variables taking part in the denominator of the balance
#          x: matrix with counts
#          y: response variable
#        col: vector with two colours, each one for a level when dichotomous
#             factor.
# OUTPUT:
#        A graphical representation of the given balance with the corresponding
#     information.
#
#------------------------------------------------------------------------------#

#' The balance value for a given matrix of counts using the variables appearing
#' in a data.frame like the output of selbal.aux
#'
#'
#' @param POS a vector containing the variables taking part in the numerator of
#' the balance.
#' @param NEG a vector containing the variables taking part in the denominator
#' of the balance.
#' @param x a \code{matrix} object with the log-transforemd counts
#' @param y the response variable, either continuous or dichotomous.
#' @param covar a \code{data.frame} with the covariates to adjust for.
#' @param col a vector with two values for representing the individuals for
#' each of the two considered groups.
#' @param logit.acc when \code{y} is dichotomous, the measure to compute for
#' the correlation between \code{y} and the proposed \emph{balance}
#' adjusting for covariates. One of the following values: \code{"AUC"} (default),
#' \code{"Dev"}, \code{"Rsq"} or \code{"Tjur"}.
#'
#' @return A graphical representation of the balance for the selected
#' individuals and the variables taking part on it.
#'
#'
#' @examples
#' # Load data set
#'   load("HIV.rda")
#' # Define x and y
#'   x <- HIV[,1:60]
#'   y <- HIV[,61]
#' # Run the algorithm
#'   Bal <- selbal.aux(x,y)
#' # Balance values for the individuals
#'   POS <- Bal[Bal[,2]=="NUM",1]
#'   NEG <- Bal[Bal[,2]=="DEN",1]
#' # Log-transformed x for the representation
#'   logx <- log(cmultRepl(x))
#'   A <- plot.bal(POS,NEG,logx,y,classy="factor",col=c("red","blue"),
#'   logit.acc="AUC)
#' # The plot
#'   grid.draw(A)
#'
#' @export plot.bal


  plot.bal <- function(POS, NEG, x, y, covar=NULL, col, logit.acc=NULL){


    # Class of the response
      classy <- class(y)
    # Family for the glm (default: gaussian)
      f.class <- "gaussian"
    # numy to be y and it will be modified if y is a factor
      numy <- y
    # If y is a factor, compute the number of levels of y
      if (classy == "factor"){
        ylev <- levels(y)
        numy <- as.numeric(y) - 1
        f.class <- "binomial"
      }


    #--------------------------------------------------------------------------#
    # STEP 1: represent the "tree" with the variables included
    #--------------------------------------------------------------------------#

    # Variables included
      T1 <- c("NUMERATOR", POS)
      T2 <- c("DENOMINATOR",NEG)
      # Parameter to specify the limits for writting
      yl <- max(length(T1), length(T2)) + .5

    # Empty plot 1 with text
      df <- data.frame()
      Imp.table <- ggplot(df) + xlim(0, 100) + ylim(-0.5, yl) + theme_void() +
        annotate("text",
                 x = 75,
                 y = floor(yl):ceiling(yl-length(T1)),
                 label = T1,
                 colour = c("royalblue1",rep("black",length(T1)-1)),
                            fontface = 2) +
        annotate("text",
                 x = 25,
                 y = floor(yl):ceiling(yl-length(T2)),
                 label = T2,
                 colour = c("royalblue1",rep("black",length(T2)-1)),
                            fontface = 2)
      # Parameter 2 to specify the limits for writting
      T1 <- c(POS,"A")
      T2 <- c(NEG,"B")
      yl2 <- max(length(T1), length(T2));
      yl2 <- yl2*1.05;
      escal <- 3;
      bot <- 5*escal*yl2/100;

      # Empty plot 2 with text
      ndiv = max(3,floor(yl2))+1;
      lineh <- 0 # 0.5*ceiling(yl-length(T1));
      df2 <- data.frame()
      colbalance<-"brown3"
      Imp.table2 <- ggplot(df2) + xlim(0, 100) + ylim(-bot, ndiv) + theme_void() +
        geom_segment(aes(x = 10, y = lineh, xend = 90, yend = lineh),color=colbalance, size=1.3) +
        geom_segment(aes(x = 50-escal, y = lineh-bot, xend = 50+escal, yend = lineh-bot),color=colbalance, size=1) +
        geom_segment(aes(x = 50+escal, y = lineh-bot, xend = 50+escal, yend = lineh),color=colbalance, size=1) +
        geom_segment(aes(x = 50-escal, y = lineh-bot, xend = 50-escal, yend = lineh),color=colbalance, size=1) +
        annotate("text",
                 x = 75,
                 y = seq(bot, floor(yl2), length.out = ndiv),
                 label = c(T1,rep("",ndiv-length(T1))),
                 colour = c(rep("black",length(T1)-1),rep(colbalance,ndiv-length(T1)+1)),
                 fontface = 2) +
        annotate("text",
                 x = 25,
                 y = seq(bot, floor(yl2), length.out = ndiv),
                 label = c(T2,rep("",ndiv-length(T2))),
                 colour = c(rep("black",length(T2)-1),rep(colbalance,ndiv-length(T2)+1)),
                 fontface = 2)



    #--------------------------------------------------------------------------#
    # STEP 2: the rest of the plots depending on the type of variable
    #--------------------------------------------------------------------------#

    # Specify FINAL.BAL
      l1 <- length(POS); l2 <- length(NEG)
      Coef <- sqrt((l1*l2)/(l1+l2))
    # The final balance
      FINAL.BAL <- Coef*(rowM(x[,POS]) - rowM(x[,NEG]))

    # Auxiliar data.frame for graphical representation
      U <- as.data.frame(cbind(numy, covar,FINAL.BAL))
      # Factor
        if(classy=="factor"){U$numy <- factor(U$numy, labels = ylev)}
      colnames(U) <- c("numy", colnames(covar),"V1")
    # Final regression model
      FIT.final <- glm(numy~., data=U, family = f.class)

    # Maximum and minimum values for y and V1 (when numeric)
      if (class(y)%in%c("numeric","integer")){
        a1<-min(U$V1); a2 <- max(U$V1)
        b1<-min(U$numy);  b2 <- max(U$numy)
      }

    # The plot depending of the class of the response variable
      if (classy=="factor"){

      #---------------------------------------------------#
      # The composition of the final plot
      #---------------------------------------------------#

        # BOXPLOT 1
        BoxP <-  ggplot(U, aes(x=numy, y=V1, fill=y)) +
          geom_boxplot(color="black", size=1) +
          scale_fill_manual(values=col) +
          theme_bw() +
          ylab("Balance") +
          xlab("Factor") +
          theme(legend.position = "none")
        # BOXPLOT 2
        BoxP2 <-  ggplot(U, aes(x=y, y=V1, fill=y)) +
          geom_boxplot(color="black", size=1) +
          scale_fill_manual(values=col) +
          theme_bw() +
          ylab("Balance") +
          xlab("") +
          theme(legend.position = "none")+coord_flip()
        # Density plot 1 for the balance
        ydensity <- ggplot(U, aes(V1, fill=y)) +
          geom_density(alpha=.5, size=1.25) +
          scale_fill_manual(values = col) +
          theme_bw() + xlab("") + ylab("") +
          theme(legend.position = "none",
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
          coord_flip()
        # Density plot 2 for the balance
        ydensity2 <- ggplot(U, aes(V1, fill=y)) +
          geom_density(alpha=.5, size=1.) +
          scale_fill_manual(values = col) +
          theme_bw() + xlab("") + ylab("") +
          theme(legend.position = "none")

      # ROC - curve
        library(pROC)
      # Build ROC curve
        A<-roc(response = U$numy,predictor = FIT.final$fitted.values)
      # Extract the sensitivity and specificity value
        ROC.TAB <- data.frame(x=1-A$specificities, y=A$sensitivities)
      # Order them for a correct representation
        ROC.TAB <- ROC.TAB[order(ROC.TAB$y),]
      # AUC value
        #auc.val <- round(logit.cor(FIT.final,y = U$numy, logit.acc = logit.acc),3)
        #auc.val<-round(as.numeric(auc(y, FIT.final$fitted.values)),3)  #diferencies
        if (class(y) == "factor") {numy<-as.numeric(y)-1} else {numy<-y} #new diferences
        auc.val<-round(as.numeric(auc(numy, FIT.final$fitted.values)),3)  #new diferencies
      # The plot
        ROC.plot <- ggplot(data=ROC.TAB, aes(x=x, y=y)) +
          geom_line() +
          ggtitle("ROC curve") +
          xlab("FPR") + ylab("TPR") +
          geom_step() +
          annotate("text", x = .5, y = .2,
                   label = paste("AUC-ROC \n",auc.val), size=3) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))

      # Load libraries
        library("gridExtra")
        library("grid")
        FINAL.P <- arrangeGrob(Imp.table, ROC.plot, BoxP, ydensity,
                               ncol=2, nrow=2, widths=c(5,1), heights=c(2, 5),
                               vp=viewport(width=0.8, height=0.8))
        FINAL.P2 <- arrangeGrob(Imp.table2, BoxP2, ydensity2,
                               ncol=1, nrow=3, widths=4, heights=c(4, 4, 4),
                               vp=viewport(width=0.75, height=1))
        # Build a list with the elements of interest
        L <- list(Global.plot = FINAL.P,Global.plot2 = FINAL.P2, ROC.plot = ROC.plot)


      } else {

      # Coefficient of determination
        C.Det <- cor(FIT.final$fitted.values, y)^2

        PLOT.G <- ggplot(U, aes(V1, y)) +
          geom_point(colour = "black", size = 3) +
          geom_abline(intercept=FIT.final$coefficients[1],
                      slope=FIT.final$coefficients[length(FIT.final$coefficients)], lwd=3, col="blue") +
          theme_bw() +
          xlab("Balance value") + ylab("Response variable") +
          ggtitle("") +
          annotate("text", x = a2 - .2*(a2-a1), y = b2 - .1*(b2-b1),
                   label = paste("italic(R) ^ 2 ==", round(C.Det,3)),
                   parse = TRUE) +
          theme(strip.text.x = element_text(size=12, angle=0,
                                            face="bold",colour="white"),
                strip.text.y = element_text(size=12, face="bold"),
                strip.background = element_rect(colour="black",
                                                fill="black"),
                plot.title = element_text(size=20, vjust=2.25, hjust= .5,
                                          face = "bold"),
                legend.title = element_text(face="bold"),
                legend.text = element_text(face="bold"))

      # Load libraries
        library(grid)
        library(gridExtra)

        FINAL.P <- arrangeGrob(Imp.table, PLOT.G, nrow=2,
                               heights=c(0.2,0.5),vp=viewport(width=0.8,
                                                              height=0.8))
        FINAL.P2 <- arrangeGrob(Imp.table2, PLOT.G, nrow=2,
                                heights=c(1,2),vp=viewport(width=0.8,
                                                               height=0.8))

        # Build a list with the elements of interest

        L <- list(Global.plot = FINAL.P,Global.plot2 = FINAL.P2)

      }


      return(L)

#      return(FINAL.P)

    }

#------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  #                           AUXILIARY FUNCTIONS
  #------------------------------------------------------------------------------#

  #------------------------------------------------------------------------------#
  # NAME: rename_OTU
  # FUNCTION: it renames the rows of an phyloseq object not to have repeated
  #           names. The idea is to write at genus level for example,
  #                       f_Lactobacillales_g_unclassified

  # INPUT:    a "phyloseq" object with @tax_table information and the rank
  #           given as a number of the column of @tax_table.
  # OUTPUT:   the "taxonomyTable" object with non-repeated modified names.
  #------------------------------------------------------------------------------#


  #' Rename each taxon
  #'
  #' \code{rename_OTU} assigns a non - repeated name to each bacteria for a given
  #'  taxonomic level.
  #'
  #'
  #' @param phy the phyloseq object where the information is contained.
  #' @param rank a \code{numeric} value indicating the taxonomic level where the
  #' names should be taken. It corresponds with the column number of
  #' \code{phy@tax_table} associated to the desired taxonomic rank.
  #' @param db the bacterial database used to name each OTU of the phyloseq
  #' object: SILVA (\emph{default}) or GreenGenes.
  #'
  #'
  #'
  #' @return A vector with the names for each bacteria. It has not repeated names
  #' and none of them starts with \emph{Incertae_Sedis} or \emph{unclassified}.
  #'
  #' @export rename_OTU


  rename_OTU <- function(phy,rank, db ="SILVA"){

    # Test if the objects are from the expected class
    stopifnot(class(phy) == "phyloseq")
    stopifnot(class(rank) == "numeric")
    # Rank abbreviation (Kingdom, Phylum, Class, Order, Family, Genus, Specie)
    Ranking<-c("k","p","c","o","f","g","s")


    ################################################################################
    # AUXILIAR FUNCTION
    ################################################################################

    replace_rare <- function(j, rk, tax_table, Nam, Ranking, db){

      # The column of tax_table we are working on
      u <- rk[j]
      # The first previous column which is not "unclassified"
      # If the name is unclassified . . .
      if (tax_table[j,u]=="unclassified"){
        while(tax_table[j,u] %in% c("unclassified", "Incertae_Sedis")){
          u <- u - 1}
        # Modify Nam
        if (db=="SILVA"){
          V<-paste(Ranking[u], tax_table[j,u],
                   paste(unlist(strsplit(Nam[j],split="_"))[-c(1:2)],
                         collapse = "_"),
                   sep="_")
        } else {
          v<-paste(tax_table[j,u],
                   paste(unlist(strsplit(Nam[j],split="_"))[-c(1:2)],
                         collapse = "_"),
                   sep="_")}

        # . . . else if the name is Incertae_Sedis
      } else if (tax_table[j,u]=="Incertae_Sedis"){
        while(tax_table[j,u] %in% c("unclassified", "Incertae_Sedis")){
          u <- u - 1}
        if (db=="SILVA"){
          # Modify Nam
          V<-paste(Ranking[u], tax_table[j,u],
                   paste(unlist(strsplit(Nam[j],split="_")), collapse = "_"),
                   sep="_")
        } else {
          V<-paste(tax_table[j,u],
                   paste(unlist(strsplit(Nam[j],split="_")), collapse = "_"),
                   sep="_")
        }
      }

      # Return V
      return(V)
    }

    ################################################################################

    # Vector with initial names
    Nam<-phy@tax_table[,rank]
    if (db=="SILVA"){Nam <- paste(Ranking[rank],Nam,sep="_")}
    # Make a copy to work with (Nam2)
    Nam2<-Nam
    # Initial r value (counting the maximum repeated values for a certain name)
    r<-2
    # Initial value for the number of categories used
    i<-1
    # A vector with the rank of the name (initially rank value)
    rk <- rep(rank, length(Nam))


    # While r>1 (while there are repeated names)
    while (r>1){
      # Repeated names
      Rep.Nam<-names(table(Nam)[(table(Nam)>1)])
      # Indices with repeated names
      Rep.Idx<-which(Nam %in% Rep.Nam)
      # Modify rk
      rk[Rep.Idx]<- rk[Rep.Idx] - 1
      # Modify Nam for Rep.Idx if it is not null
      if (length(Rep.Idx)!=0){
        if(db=="SILVA"){
          Nam[Rep.Idx] <- paste(Ranking[rank-i],phy@tax_table[Rep.Idx,rank-i],
                                Nam2[Rep.Idx],sep="_")
        }else{
          Nam[Rep.Idx] <- paste(phy@tax_table[Rep.Idx,rank-i],
                                Nam2[Rep.Idx],sep="_")
        }
      }
      # Modify r as the number of the maximum repeated name
      r<-max(table(Nam))
      i<-i+1

    }


    # Load library
    library(qdapRegex)
    # Extract the first value for each name
    First.NAM <- unlist(lapply(Nam, function(x)
      unlist(rm_between(x, "_","_", extract = T))[1]))
    # Indices for the names to replace (with the first name as unclassified
    # or Incertae_Sedis)
    IDX.Rep <- which(First.NAM %in% c("unclassified", "Incertae"))
    # If there are unclassified, . . .
    if (length(IDX.Rep) !=0){
      for (i in 1:length(IDX.Rep)){
        Nam[IDX.Rep[i]] <-  replace_rare(IDX.Rep[i],
                                         rk,
                                         phy@tax_table,
                                         Nam,
                                         Ranking,
                                         db)
      }

    }

    return(Nam)
  }

  ################################################################################

  #------------------------------------------------------------------------------#
  # Auxiliar function in order to replace zeros if necesary
  #------------------------------------------------------------------------------#

  #' Zero replacement for compositional data
  #'
  #' \code{cmultRepl2} replaces the zeros for a matrix where each row is
  #' compositional
  #'
  #'
  #' @param x a \code{matrix} object with the information of variables
  #' (\emph{columns}) for each sample (\emph{rows}).
  #' @param zero.rep if \emph{"bayes"} the Bayesian - Multiplicative treatment
  #' implemented in \code{zCompositions} is applied. If \emph{"one"}, a
  #' pseudocount of 1 is added to the whole matrix.
  #'
  #'
  #'
  #' @return The initial matrix after the zero - replacement normalized so that
  #' each sample's composition sums one.
  #'
  #'
  #' @examples
  #'
  #' # Load the count matrix (with zeros)
  #'   x <- HIV[,1:60]
  #' # Zero replacement
  #'   x.non0 <- cmultRepl2(x, zero.rep = "bayes")
  #'
  #'
  #' @export cmultRepl2



  cmultRepl2 <- function(x, zero.rep = "bayes"){

    # Load library
    library(zCompositions)
    # If there are zeros, use cmultRepl
    if(sum(x==0)!=0){
      if (zero.rep =="bayes"){
        new.x <- cmultRepl(x, suppress.print = T)
      } else if (zero.rep =="one"){
        new.x <- x + 1
      }
    }else { new.x <- x}
    # Return new.x
    return(new.x)
  }

  #------------------------------------------------------------------------------#


  #------------------------------------------------------------------------------#
  # NAME: plot.tab
  # INPUT:
  #     dat: the table to represent
  #     col: vector with two colors, each of them referring to the variables
  #          appearing in the numerator and in the denominator respectively.
  # OUTPUT: a plot with the information given in dat (for the first five columns)
  #------------------------------------------------------------------------------#


  #' Plots the cross - validated summary table
  #'
  #'
  #'
  #' @param dat the \code{data.frame} object to draw obtained from
  #'  \code{selbal.cv} function.
  #' @param col vector with two colors, each of them referring to the variables
  #' appearing in the numerator and in the denominator, respectively.
  #'
  #'
  #' @return A colored table with the information of the given \code{data.frame}:
  #'
  #' \itemize{
  #'  \item The first column of the table (%) represents the proportion of times
  #'  a variables has been selected for a balance.
  #'
  #'  \item The second column (\emph{Global}) shows the result for the balance
  #'  obtained using all the available samples.
  #'
  #'  \item The last three columns represent the most repeated balances in the
  #'  cross - validation procedure.
  #'
  #'  }
  #'
  #' @examples
  #' #---------------------------------------------------#
  #' # 1) Compute a cross - validated balance selection
  #' #---------------------------------------------------#
  #'  # Load data set
  #'    load("HIV.rda")
  #'   # Define x and y
  #'     x <- HIV[,1:60]
  #'     y <- HIV[,62]
  #'  # Run the algorithm
  #'     CV.Bal <- selbal.cv(x,y)
  #'
  #' #----------------------------------------#
  #' # 2) Plot the table
  #' #----------------------------------------#
  #'
  #'   plot.tab(CV.Bal[[4]])
  #'
  #' @export plot.tab


  plot.tab <- function(dat, col = c("steelblue1","tomato1")){
    # Load library
    library(gtable)
    # Data.frame dimension
    dim.dat<-dim(dat)
    # Colnames dat
    colnames(dat)<- c("%","Global",paste("BAL", 1:(dim.dat[2]-2), sep=" "))
    # Modify dat
    dat.l<-data.frame(cbind(c(" ",row.names(dat)),rbind(colnames(dat),dat)))
    row.names(dat.l)<-colnames(dat.l)<-NULL

    # Build my.data2 with the row.names(my.data) as a first column
    dat2<-dat.l
    dat2<-apply(dat2,2,function(x) as.character(x))
    # Extract some information
    nc <- ncol(dat2)
    nr <- nrow(dat2)
    n <- nc*nr
    # Change the values for the background colour
    # Legend (NUM = 10, DEN = 9, OTHER <- 0, HEADER/ROW.NAME <- 1)
    dat2[dat2=="NUM"] <- 10
    dat2[dat2=="DEN"] <- 9
    dat2[,1]<-dat2[1,] <- 1
    dat2[-1,2]<- 0
    dat2[nr,]<- 0

    # Information but the last linefor colors
    Letra <- as.character(factor(dat2, labels=c("black", "gray15", col[2],
                                                col[1])))

    # Delete words from dat, only selecting the numbers
    dat1<-dat.l
    dat1[-c(1,nrow(dat1)),-c(1,2)]<-" "
    dat1[is.na(dat1)]<-" "

    # Filling colors
    fill <- as.character(factor(t(dat2),labels=c("white", "gray85", col[2],
                                                 col[1])))
    # Define the background of cells
    fill <- lapply(seq_len(n), function(ii) rectGrob(gp=gpar(fill=fill[ii])))

    # Some calculations for cell sizes
    row_heights <- function(m){
      do.call(unit.c, apply(m, 1, function(l)
        max(do.call(unit.c, lapply(l, grobHeight)))))
    }

    col_widths <- function(m){
      do.call(unit.c, apply(m, 2, function(l)
        max(do.call(unit.c, lapply(l, grobWidth)))))
    }

    # Object as matrix
    label_matrix <- as.matrix(dat1)

    nc <- ncol(label_matrix)
    nr <- nrow(label_matrix)
    n <- nc*nr

    # Text written in the table
    # Auxiliar values for text justification
    pos <- rep(c(0.5, rep(0.96,nr-2),0.5),nc)
    jus <- rep(c(0.5, rep(1,nr-2),0.5),nc)
    # Third column centered
    pos[(nr +1):(2*nr)] <- jus[(nr+1):(2*nr)] <- 0.5

    # Define the text characteristics
    labels <- lapply(seq_len(n), function(ii)
      textGrob(label_matrix[ii],x=pos[ii],
               just=jus[ii],
               gp=gpar(fontface="bold",col=Letra[ii])))
    label_grobs <- matrix(labels, ncol=nc)

    # Place labels in a gtable
    g <- gtable_matrix("table", grobs=label_grobs,
                       widths=col_widths(label_grobs) + unit(8,"mm"),
                       heights=row_heights(label_grobs) + unit(5,"mm"))

    # Add the background
    g <- gtable_add_grob(g, fill, t=rep(seq_len(nr), each=nc),
                         l=rep(seq_len(nc), nr), z=0, name="fill")
    # Graphical representation
    grid.draw(g)

  }

  #------------------------------------------------------------------------------#




  ################################################################################
  # FUNCTION: logit.cor
  ################################################################################

  #' Computes an association value  between a dichotomous variable and a
  #' continuous one.
  #'
  #'
  #'
  #' @param FIT a \code{glm} object referred to the logistic regression of a
  #' dichotomous variable.
  #' @param y the response variable (dichotomous).
  #' @param logit.acc when \code{y} is dichotomous, the measure to compute for
  #' the correlation between \code{y} and the proposed \emph{balance}
  #' adjusting for covariates. One of the following values: \code{"Rsq"} (default),
  #'  \code{"AUC"} or \code{"Tjur"}.
  #'
  #'
  #' @return The association value using the selected method \code{logit.acc}.
  #'
  #'
  #' @export logit.cor

  # Define the function logit.cor
  logit.cor <- function(FIT, y, covar = NULL,logit.acc){
    if (logit.acc == "AUC"){
      #d <- as.numeric(auc(y, FIT$fitted.values)) #
      if (class(y) == "factor") {numy<-as.numeric(y)-1} else {numy<-y} #new diferences
      d <- as.numeric(auc(numy, FIT$fitted.values)) # new diferences
    } else if (logit.acc == "Rsq"){
      d <- cor(as.numeric(y), FIT$fitted.values)^2
    } else if (logit.acc == "Tjur"){
      if (class(y) == "factor") {numy<-as.numeric(y)-1} else {numy<-y}
      d <- mean(FIT$fitted.values[numy==1]) - mean(FIT$fitted.values[numy==0])
    } else if (logit.acc == "Dev"){
      f.class <- ifelse (class(y) == "factor", "binomial", "gaussian")
      if (class(y) == "factor") {numy<-as.numeric(y)-1} else {numy<-y}
      #if (is.null(covar)){
        d<-1-(deviance(FIT)/deviance(glm(numy~1,family=f.class)))  # proportion of explained deviance
      # } else {
      # d<-1-(deviance(FIT)/deviance(glm(numy~., data=cbind(numy,covar),family=f.class)))  # proportion of explained deviance
      # }
      # d<-1-(deviance(FIT)/nulldev0)  # proportion of explained deviance
    }

    # Return the value
    return(d)
  }

  #------------------------------------------------------------------------------#


  ################################################################################
  # FUNCTION: rowM
  ################################################################################

  #' Calculates the mean of each row of a matrix (though having only one column).
  #'
  #'
  #'
  #' @param x a \code{matrix} object.
  #'
  #' @return A \code{vector} with the mean of each row of \code{x}, even if the
  #' matrix only contains one column.
  #'
  #' @examples
  #' # Build a matrix with one column
  #'   M <- matrix(rnorm(10), nrow=1)
  #' # rowM (resulting on the same matrix M)
  #'   rowM(M)
  #'
  #' # Build a matrix
  #'   M <- matrix (runif(100), nrow=10)
  #' # Apply rowM function
  #'   rowM(M)
  #'
  #' @export rowM

  # Define the function rowM
  rowM <- function(x){
    if(is.vector(x)) {u <- x
    } else { u <- rowMeans(x)}
    return(u)
  }

  #------------------------------------------------------------------------------#
