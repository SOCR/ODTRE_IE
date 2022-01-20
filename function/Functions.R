# R functions
#######Utility functions###############

#' Sample treatment A from a matrix of sampling probabilities.
#' Returns a list of sampled treatments.
#'  
#' @param matrix.pi a matrix of sampling probabilities, which could be non-normalized.
A.sim <- function(matrix.pi) {
  # Obtain sample size
  N <- nrow(matrix.pi)
  # Obtain treatment options
  K <- ncol(matrix.pi)
  if (N <= 1 |
      K <= 1)
    stop("Sample size or treatment options are insufficient!")
  if (min(matrix.pi) < 0)
    stop("Treatment probabilities should not be negative!")
  
  # normalize probabilities to add up to 1 and simulate treatment A for each row
  probs <-
    t(apply(matrix.pi, 1, function(x) {
      x / sum(x, na.rm = TRUE)
    }))
  A <- apply(probs, 1, function(x)
    sample(0:(K - 1), 1, prob = x))
  return(A)
}

#' Estimate propensity score by a multinomial model
#' Returns a matrix of with N rows and length(unique(A)) columns, 
#' representing the propensity scores for all records.
#' 
#' @param A treatment vector.
#' @param Xs covariate matrix.
#' @param covars user specified covariate structure for the multinomial model.
M.propen <- function(A, Xs, covars = NULL) {
  if (ncol(as.matrix(A)) != 1)
    stop("Cannot handle multiple stages of treatments together!")
  if (length(A) != nrow(as.matrix(Xs)))
    stop("A and Xs do not match in dimension!")
  if (length(unique(A)) <= 1)
    stop("Treament options are insufficient!")
  class.A <- sort(unique(A))#class.A=Unique treatments
  
  require(nnet)
  s.data <- data.frame(A, Xs)
  # multinomial regression with output suppressed
  if (is.null(covars)) {
    model <- capture.output(mlogit <- multinom(A ~ ., data = s.data))
  } else{
    model <- capture.output(mlogit <-
                              multinom(as.formula(paste(
                                "A ~", paste(covars, collapse = "+")
                              )),
                              data = s.data))
  }
  s.p <- predict(mlogit, s.data, "probs")
  
  if (length(class.A) == 2) {
    s.p <- cbind(1 - s.p, s.p)
  }
  colnames(s.p) <- paste("pi=", class.A, sep = "")
  
  return(s.p)
}

#' Estimate conditional means for multiple stages
#' Returns a list of estimated conditional means and the corresponding model.
#' 
#' @param Y a continuous outcome of interest.
#' @param As (A1, A2, ...) a matrix of treatments at multiple stages.
#' Stage t has treatment K_t options labeled as 0, 1, ..., K_t-1.
#' @param H a matrix of covariates before assigning final treatment, 
#' excluding previous treatment variables.
Reg.mu <- function(Y, As, H) {
  if (nrow(as.matrix(As)) != nrow(as.matrix(H))) {
    stop("Treatment and Covariates do not match in dimension!")
  }
  Ts <- ncol(as.matrix(As))
  N <- nrow(as.matrix(As))
  if (Ts < 0 | Ts > 3)
    stop("Only support 1 to 3 stages!")
  H <- as.matrix(H)
  
  if (Ts == 1L) {
    #one stage
    A1 <- as.matrix(As)[, 1]
    A1 <- as.factor(A1)
    KT <- length(unique(A1)) # treatment options at last stage
    if (KT < 2)
      stop("No multiple treatment options!")
    
    RegModel <- lm(Y ~ H * A1)
    mus.reg <- matrix(NA, N, KT)
    for (k in 1L:KT) {
      mus.reg[, k] <- predict(RegModel,
                              newdata = data.frame(H,
                                                   A1 = factor(rep(
                                                     sort(unique(A1))[k], N
                                                   ))))
    }
  }
  if (Ts == 2L) {
    A1 <- as.matrix(As)[, 1]
    A2 <- as.matrix(As)[, 2]
    A1 <- as.factor(A1)
    A2 <- as.factor(A2)
    KT <- length(unique(A2))
    if (KT < 2)
      stop("No multiple treatment options!")
    
    RegModel <- lm(Y ~ (H + A1) * A2)
    
    mus.reg <- matrix(NA, N, KT)
    for (k in 1L:KT) {
      mus.reg[, k] <- predict(RegModel,
                              newdata = data.frame(H, A1,
                                                   A2 = factor(rep(
                                                     sort(unique(A2))[k], N
                                                   ))))
    }
  }
  if (Ts == 3L) {
    A1 <- as.matrix(As)[, 1]
    A2 <- as.matrix(As)[, 2]
    A3 <- as.matrix(As)[, 3]
    A1 <- as.factor(A1)
    A2 <- as.factor(A2)
    
    A3 <- as.factor(A3)
    KT <- length(unique(A3))
    if (KT < 2)
      stop("No multiple treatment options!")
    
    RegModel <- lm(Y ~ (H + A1 + A2) * A3)
    
    mus.reg <- matrix(NA, N, KT)
    for (k in 1L:KT) {
      mus.reg[, k] <- predict(RegModel,
                              newdata = data.frame(H, A1, A2,
                                                   A3 = factor(rep(
                                                     sort(unique(A3))[k], N
                                                   ))))
    }
  }
  
  output <- list(mus.reg, RegModel)
  names(output) <- c("mus.reg", "RegModel")
  return(output)
}

#' Calculate mus in AIPW 
#' Returns mu that is a weighted average of original Y and the Y* depend on model.
#' 
#' @param Y a continuous outcome of interest.
#' @param A treatment vector. 
#' @param pis.hat estimated propensity matrix 
#' @param mus.reg regression-based conditional means 
mus.AIPW <- function(Y, A, pis.hat, mus.reg) {
  class.A <- sort(unique(A))
  K <- length(class.A)
  N <- length(A)
  if (K < 2 | N < 2)
    stop("No multiple treatments or samples!")
  if (ncol(pis.hat) != K |
      ncol(mus.reg) != K | nrow(pis.hat) != N | nrow(mus.reg) != N) {
    stop("Treatment, propensity or conditional means do not match!")
  }
  #AIPW estimates
  mus.a <- matrix(NA, N, K)
  for (k in 1L:K) {
    mus.a[, k] <- (A == class.A[k]) * Y / pis.hat[, k] +
      (1 - (A == class.A[k]) / pis.hat[, k]) * mus.reg[, k]
  }
  return(mus.a)
}

#' Calculate AIPW adaptive contrasts and working orders
#' Returns a data frame of contrasts and working orders.
#' 
#' @param Y outcome of interest.
#' @param A treatment vector. 
#' @param pis.hat estimated propensity matrix.  
#' @param mus.reg regression-based conditional means.
CL.AIPW <- function(Y, A, pis.hat, mus.reg) {
  class.A <- sort(unique(A))
  K <- length(class.A)
  N <- length(A)
  if (K < 2 | N < 2)
    stop("No multiple treatments or samples!")
  if (ncol(pis.hat) != K |
      ncol(mus.reg) != K | nrow(pis.hat) != N | nrow(mus.reg) != N) {
    stop("Treatment, propensity or conditional means do not match!")
  }
  
  #AIPW estimates
  mus.a <- matrix(NA, N, K)
  for (k in 1:K) {
    mus.a[, k] <- (A == class.A[k]) * Y / pis.hat[, k] +
      (1 - (A == class.A[k]) / pis.hat[, k]) * mus.reg[, k]
  }
  
  # C.a1 and C.a2 are AIPW contrasts; l.a is AIPW working order
  C.a1 <- C.a2 <- l.a <- rep(NA, N)
  for (i in 1:N) {
    # largest vs. second largest
    C.a1[i] <- max(mus.a[i, ]) - sort(mus.a[i, ], decreasing = T)[2]
    # largest vs. smallest
    C.a2[i] <- max(mus.a[i, ]) - min(mus.a[i, ])
    # minus 1 to match A's range of 0,...,K-1
    l.a[i] <- which(mus.a[i, ] == max(mus.a[i, ])) - 1
    #l.a=the index of maximum outcome-1
  }
  output <- data.frame(C.a1, C.a2, l.a)
  return(output)
}

#' Choose the optimal treatment with given estimated mu's matrix 
#' in one specific stage
#' Returns the expected counterfactual outcome under optimal treatment and 
#' the optimal treatment.
#' 
#' @param A the treatment options received by patients in this stage.
#' @param mu.hat the estimated counterfactural outcomes according to the treatments.
Opt.A <- function(A, mus.hat) {
  class.A <- sort(unique(A))
  if (length(class.A) == 1) {
    trt.opt1 <- class.A
    Ey.opt1 <- mean(mus.hat)
  } else{
    if (length(A) != nrow(mus.hat) || length(class.A) != ncol(mus.hat)) {
      stop("Treatment options and mean matrix dimension do not match!")
    }
    # pick a single best treatment for all patients
    c.means <- apply(mus.hat, 2, mean)
    Ey.opt1 <- max(c.means)
    trt.opt1 <- class.A[which(c.means == Ey.opt1)]
  }
  outs <- list(Ey.opt1, trt.opt1)
  names(outs) <- c("Ey.opt1", "trt.opt1")
  return(outs)
}

#' Combine two matrices with different dimensions, 
#' filling in NAs for additional rows/columns
#' Return combined matrix 
#' 
#' @param m1 first matrix.
#' @param m2 second matrix.
#' @param by combining by "column" or by "row".
combine.mat <- function(m1, m2, by = "column") {
  nrow1 <- nrow(m1)
  ncol1 <- ncol(m1)
  nrow2 <- nrow(m2)
  ncol2 <- ncol(m2)
  if (by == "column") {
    combine <- matrix(NA, max(nrow1, nrow2), ncol1 + ncol2)
    combine[1:nrow1, 1:ncol1] <- m1
    combine[1:nrow2, (ncol1 + 1):(ncol1 + ncol2)] <- m2
  }
  if (by == "row") {
    combine <- matrix(NA, nrow1 + nrow2, max(ncol1, ncol2))
    combine[1:nrow1, 1:ncol1] <- m1
    combine[(nrow1 + 1):(nrow1 + nrow2), 1:ncol2] <- m2
  }
  return(combine)
}

#' Split data by X to fit child nodes, calculate new means for each child node
#' Returns a data frame containing the covariate to split, average counterfactual 
#' outcome for the current node, and corresponding treatments. 
#' 
#' @param X one patient covariate.
#' @param A received treatments.
#' @param mus.hat the estimated mean counterfactual outcomes.
#' @param minsplit is the minimum cases in each node.
Split.X <- function(X, A, mus.hat, minsplit = 20) {
  n <- length(X)
  X.val <- unique(X)
  n.X.val <- length(X.val)
  class.A <- sort(unique(A))
  
  if (n < 2 * minsplit ||
      n.X.val < 2L || length(class.A) < 2L)
    return(NULL)
  
  if (is.numeric(X) == TRUE ||
      is.ordered(X) == TRUE || n.X.val == 2L) {
    X.val <- sort(X.val)
    # reduce computation by using quantiles
    if (n.X.val > 100L) {
      X.val <- quantile(X, 1:100 / 100)
      n.X.val <- 100L
    }
    Ey.opt1.X <-
      trt.left.X <-
      trt.right.X <- rep(NA, n.X.val - 1)#initalize E(Y|optX)
    for (i in 1L:(n.X.val - 1)) {
      left <- which(X <= X.val[i])
      if (length(left) >= minsplit && length(left) <= n - minsplit) {
        left.est <- Opt.A(A[left],
                          mus.hat[left, which(class.A %in% unique(A[left]))])
        right.est <- Opt.A(A[-left],
                           mus.hat[-left, which(class.A %in% unique(A[-left]))])
        
        trt.left.X[i] <- left.est$trt.opt1
        trt.right.X[i] <- right.est$trt.opt1
        Ey.opt1.X[i] <- length(left) / n * left.est$Ey.opt1 +
          (1 - length(left) / n) * right.est$Ey.opt1
      }
    }
    # pick the best split of X
    if (sum(!is.na(Ey.opt1.X)) > 0L) {
      mEy.opt1 <- max(Ey.opt1.X, na.rm = T)
      cutoff1 <-
        which(Ey.opt1.X == mEy.opt1)[1]#take the minimum cutoff
      X.cutoff <- X.val[cutoff1]
      trt.L <- trt.left.X[cutoff1]
      trt.R <- trt.right.X[cutoff1]
      
      output <- data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
      names(output) <- c("X.subset", "mEy.opt1", "trt.L", "trt.R")
    } else{
      return(NULL)
    }
    
  }
  
  if (is.numeric(X) == F && is.ordered(X) == F && n.X.val > 2L) {
    n.X.combo <- 2 ^ (n.X.val - 1) - 1
    X.combo <- combn(X.val, 1)
    if (n.X.val > 3L && n.X.val %% 2 == 1L) {
      for (k in 2L:(n.X.val - 1) / 2)
        X.combo <- combine.mat(X.combo, combn(X.val, k))
    }
    if (n.X.val > 3L && n.X.val %% 2 == 0L) {
      for (k in 2L:(n.X.val / 2)) {
        if (k < (n.X.val / 2))
          X.combo <- combine.mat(X.combo, combn(X.val, k))
        if (k == (n.X.val / 2)) {
          temp.mat <- combn(X.val[-1], k - 1)
          first.row <- rep(X.val[1], ncol(temp.mat))
          X.combo <-
            combine.mat(X.combo, rbind(first.row, temp.mat))
        }
      }
    }
    
    Ey.opt1.X <- trt.left.X <- trt.right.X <- rep(NA, n.X.combo)
    for (i in 1L:n.X.combo) {
      left <- which(X %in% X.combo[, i])
      if (length(left) >= minsplit && length(left) <= n - minsplit) {
        left.est <- Opt.A(A[left],
                          mus.hat[left, which(class.A %in% unique(A[left]))])
        right.est <- Opt.A(A[-left],
                           mus.hat[-left, which(class.A %in% unique(A[-left]))])
        
        trt.left.X[i] <- left.est$trt.opt1
        trt.right.X[i] <- right.est$trt.opt1
        Ey.opt1.X[i] <- length(left) / n * left.est$Ey.opt1 +
          (1 - length(left) / n) * right.est$Ey.opt1
      }
    }
    # pick the best split of X
    if (sum(!is.na(Ey.opt1.X)) > 0L) {
      mEy.opt1 <- max(Ey.opt1.X, na.rm = T)
      cutoff1 <- which(Ey.opt1.X == mEy.opt1)[1]
      X.subset <- X.combo[, cutoff1]
      # change a vector into a single string while removing NA's
      X.subset <- paste(X.subset[!is.na(X.subset)], collapse = " ")
      trt.L <- trt.left.X[cutoff1]
      trt.R <- trt.right.X[cutoff1]
      
      output <- data.frame(X.subset, mEy.opt1, trt.L, trt.R)
      names(output) <- c("X.subset", "mEy.opt1", "trt.L", "trt.R")
    } else{
      return(NULL)
    }
  }
  return(output)
}


#' Pick the best X for split
#' Returns the best covariate, cut off, mean expected counterfactual outcome,
#' and treatments for the node split
#' 
#' @param H patient covariate matrix.
#' @param A received treatments.
#' @param mus.hat the estimated mean counterfactual outcomes.
#' @param minsplit the minimum cases in each node.
best.H <- function(H, A, mus.hat, minsplit = 20) {
  p <- ncol(H)
  output <- as.data.frame(matrix(NA, p, 5))
  output[, 1] <- 1:p
  colnames(output) <-
    c("X", "X.subset", "mEy.opt1", "trt.L", "trt.R")
  #output is a p*5 matrix p=number of features.
  
  for (i in 1:p) {
    split.i <-
      Split.X(
        X = H[, i],
        A = A,
        mus.hat = mus.hat,
        minsplit = minsplit
      )
    
    if (!is.null(split.i))
      output[i, -1] <- split.i
    #output column i other than the first cell is replaced by the best split under X
  }
  if (sum(!is.na(output$mEy.opt1)) > 0L) {
    max.p <- which(output$mEy.opt1 == max(output$mEy.opt1, na.rm = T))[1]
    opt.output <- output[max.p, ]
    if (opt.output$trt.L == opt.output$trt.R) {
      return(NULL)
    } else{
      return(opt.output)
    }
  } else{
    return(NULL)
  }
}

#' Look ahead - split X with looking head for 1 more step
#' Returns the best split by looking 1 step ahead with the best covariate, 
#' cut off, mean expected counterfactual outcome, and treatments for the node split.
#' 
#' @param X a patient covariate.
#' @param A received treatments.
#' @param H patient covariate matrix.
#' @param mus.hat the estimated mean counterfactual outcomes.
#' @param minsplit the minimum cases in each node.
Split.X.lh <- function(X, A, H, mus.hat, minsplit = 20) {
  n <- length(X)
  X.val <- unique(X)
  n.X.val <- length(X.val)
  class.A <- sort(unique(A))
  
  if (n < 2 * minsplit ||
      n.X.val < 2L || length(class.A) < 2L)
    return(NULL)
  
  if (is.numeric(X) == T || is.ordered(X) == T || n.X.val == 2L) {
    X.val <- sort(X.val)
    # reduce computation by using quantiles
    if (n.X.val > 100L) {
      X.val <- quantile(X, 1:100 / 100)
      n.X.val <- 100L
    }
    Ey.opt1.X <-
      Ey.opt2.X <- trt.left.X <- trt.right.X <- rep(NA, n.X.val - 1)
    for (i in 1L:(n.X.val - 1)) {
      left <- which(X <= X.val[i])
      if (length(left) >= minsplit && length(left) <= n - minsplit) {
        # further split left and right
        # lookahead one step
        best.H.L <- best.H.R <- NULL
        if (length(left) >= 2 * minsplit) {
          #After the split left still have room for the next split
          best.H.L <- best.H(
            H = H[left, ],
            A = A[left],
            mus.hat = mus.hat[left, which(class.A %in% unique(A[left]))],
            minsplit = minsplit
          )
          #output the best variable to split output cutoff, mean,left, right optimal Y
        }
        if (n - length(left) >= 2 * minsplit) {
          best.H.R <- best.H(
            H = H[-left, ],
            A = A[-left],
            mus.hat = mus.hat[-left, which(class.A %in% unique(A[-left]))],
            minsplit = minsplit
          )
        }
        # if unable to further split, then pick 1 treatment
        left.est <- Opt.A(A[left],
                          mus.hat[left, which(class.A %in% unique(A[left]))])
        #optimal Y in all left area
        right.est <- Opt.A(A[-left],
                           mus.hat[-left, which(class.A %in% unique(A[-left]))])
        
        trt.left.X[i] <- left.est$trt.opt1 #option
        trt.right.X[i] <- right.est$trt.opt1
        
        if (!is.null(best.H.L)) {
          left.Ey <- best.H.L$mEy.opt1
        } else{
          left.Ey <- left.est$Ey.opt1
        }
        if (!is.null(best.H.R)) {
          right.Ey <- best.H.R$mEy.opt1
        } else{
          right.Ey <- right.est$Ey.opt1
        }
        # the max Ey after assigning 2 treatments to two child nodes is used for later splits
        Ey.opt1.X[i] <- length(left) / n * left.est$Ey.opt1 +
          (1 - length(left) / n) * right.est$Ey.opt1
        
        # the max Ey if lookahead is used to choose the optimal X and its threshould
        Ey.opt2.X[i] <-
          length(left) / n * left.Ey + (1 - length(left) / n) * right.Ey
      }
    }
    # pick the best split of X
    if (sum(!is.na(Ey.opt1.X)) > 0L & sum(!is.na(Ey.opt2.X)) > 0L) {
      mEy.opt1 <- max(Ey.opt1.X, na.rm = T)
      mEy.opt2 <- max(Ey.opt2.X, na.rm = T)
      # based on lookahead
      cutoff <- which(Ey.opt2.X == mEy.opt2)[1]
      X.cutoff <- X.val[cutoff]
      trt.L <- trt.left.X[cutoff]
      trt.R <- trt.right.X[cutoff]
      
      output <- data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
      names(output) <- c("X.subset", "mEy.opt1", "trt.L", "trt.R")
    } else{
      return(NULL)
    }
  } else{
    stop("Lookahead currently only supports numerical or ordinal covariates!")
  }
  
  return(output)
}

#' Best covariate X based on looking ahead.
#' Returns the best covariate, cut off, mean expected counterfactual outcome,
#' and treatments for the node split according to looking ahead.
#' 
#' @param H patient covariate matrix.
#' @param A received treatments.
#' @param mus.hat the estimated mean counterfactual outcomes.
#' @param minsplit the minimum cases in each node.
best.H.lh <- function(H, A, mus.hat, minsplit = 20) {
  p <- ncol(H)
  output <- as.data.frame(matrix(NA, p, 5))
  output[, 1] <- 1:p
  colnames(output) <- c("X", "X.subset", "mEy.opt1", "trt.L", "trt.R")
  
  for (i in 1:p) {
    split.i <- Split.X.lh(
      X = H[, i],
      A = A,
      H = H,
      mus.hat = mus.hat,
      minsplit = minsplit
    )
    
    if (!is.null(split.i))
      output[i, -1] <- split.i
  }
  if (sum(!is.na(output$mEy.opt1)) > 0L) {
    max.p <- which(output$mEy.opt1 == max(output$mEy.opt1, na.rm = T))[1]
    opt.output <- output[max.p, ]
    if (opt.output$trt.L == opt.output$trt.R) {
      return(NULL)
    } else{
      return(opt.output)
    }
  } else{
    return(NULL)
  }
}

#######Tree-based DTR##########

#' Estimate optimal DTR for one stage using Tree-based Reinforcement Learning (T-RL)
#' Returns a tree-based DTR with binary cut offs. The returned DTR matrix describes 
#' the tree nodes from top to down, left to right. Each row includes node number, 
#' the covariate for the node split, the corresponding cutoff value for the covariate, 
#' the estimated conditional mean outcome, and optimal treatment. 
#' 
#' @param Y observed outcome.
#' @param A received treatment matrix.
#' @param H patient covariate matrix.
#' @param pis.hat estimated propensity matrix.
#' @param m.method estimating method.  
#' @param mus.reg regression-based conditional means.
#' @param depth maximum tree depth.
#' @param lambda.pct minimum improvement for a split. We consider a split if 
#' improved at least lambda, as a percent of the parent node EY.
#' @param minsplit the minimum cases in each node.
#' @param lookahead a boolean for if or not look ahead when considering the best
#' split.
DTRtree <- function(Y,
           A,
           H,
           pis.hat = NULL,
           m.method = c("AIPW", "randomForest"),
           mus.reg = NULL,
           depth = 5,
           lambda.pct = 0.05,
           minsplit = 20,
           lookahead = F) {
    # initialization
    # indicator for subset data
    n <- length(Y)#number of people
    I.node <- rep(1, n)#indicator of nodes
    class.A <- sort(unique(A))
    output <- matrix(NA, 1, 5)
    colnames(output) <- c("node", "X", "cutoff", "mEy", "trt")
    
    # estimate mus.hat if not given
    if (m.method[1] == "AIPW") {
      # estimate propenstiy matrix if not given, using all data
      # same propensity for all subset data
      if (is.null(pis.hat))
        pis.hat <- M.propen(A = A, Xs = H)
      if (is.null(mus.reg))
        mus.reg <- Reg.mu(Y = Y, As = A, H = H)$mus.reg
      mus.hat <-
        mus.AIPW(
          Y = Y,
          A = A,
          pis.hat = pis.hat,
          mus.reg = mus.reg
        )
    } else if (m.method[1] == "randomForest") {
      require(randomForest)
      RF <- randomForest(Y ~ ., data = data.frame(A, H))
      mus.hat <- matrix(NA, n, length(class.A))
      for (i in 1L:length(class.A)) {
        mus.hat[, i] <-
          predict(RF, newdata = data.frame(A = rep(class.A[i], n), H))
      }
    } else{
      stop("The method for estimating conditional means is not available!")
    }
    
    # expected outcome at root
    root <- Opt.A(A, mus.hat)
    Ey0 <- root$Ey.opt1
    
    # split if improved at least lambda, as a percent of Ey0
    lambda <- abs(Ey0) * lambda.pct
    
    for (k in 1L:depth) {
      output <- rbind(output, matrix(NA, 2 ^ k, 5))
      output[, 1] <- 1L:(2 ^ (k + 1) - 1)
      if (k == 1L) {
        # apply lookahead to the first split, the most important split
        # only to first split so as to save computation time
        # use a larger minsplit for the first split
        if (lookahead) {
          best.H.1 <-
            best.H.lh(
              H = H,
              A = A,
              mus.hat = mus.hat,
              minsplit = 0.15 * n
            )
        } else{
          best.H.1 <-
            best.H(
              H = H,
              A = A,
              mus.hat = mus.hat,
              minsplit = minsplit
            )
        }
        if (is.null(best.H.1) == F &&
            best.H.1$mEy.opt1 > Ey0 + lambda) {
          #meet the split criteria
          output[k, -1] <-
            c(best.H.1$X, best.H.1$X.subset, best.H.1$mEy.opt1, NA)
          I.node[I.node == k &
                   H[, best.H.1$X] <= best.H.1$X.subset] <- 2 * k
          output[2 * k, -1] <- c(NA, NA, NA, best.H.1$trt.L)
          I.node[I.node == k &
                   H[, best.H.1$X] > best.H.1$X.subset] <- 2 * k + 1
          output[2 * k + 1, -1] <- c (NA, NA, NA, best.H.1$trt.R)
        } else{
          output[k, 4:5] <- c(root$Ey.opt1, root$trt.opt1)
          break
        }
      } else{
        for (j in (2 ^ (k - 1)):(2 ^ k - 1)) {
          if (!is.na(output[trunc(j / 2), 2])) {
            best.H.j <- best.H(
              H = H[I.node == j, ],
              A = A[I.node == j],
              mus.hat = mus.hat[I.node == j, ],
              minsplit = minsplit
            )
            if (is.null(best.H.j) == F &&
                best.H.j$mEy.opt1 > output[trunc(j / 2), 4] + lambda) {
              output[j, -1] <-
                c(best.H.j$X,
                  best.H.j$X.subset,
                  best.H.j$mEy.opt1,
                  NA)
              I.node[I.node == j &
                       H[, best.H.j$X] <= best.H.j$X.subset] <- 2 * j
              output[2 * j, -1] <- c(NA, NA, NA, best.H.j$trt.L)
              I.node[I.node == j &
                       H[, best.H.j$X] > best.H.j$X.subset] <- 2 * j + 1
              output[2 * j + 1, -1] <- c(NA, NA, NA, best.H.j$trt.R)
            }
          }
        }
        if (sum(is.na(output[(2 ^ (k - 1)):(2 ^ k - 1), 2])) == 2 ^ (k - 1))
          break
      }
    }
    output <- output[!is.na(output[, 2]) | !is.na(output[, 5]), ]
    return(output)
  }

#' Predict optimal treatment using the output from DTRtree
#' Returns a list of optimal treatments the patients.
#' 
#' @param treeout the output DTR object from DTRtree.
#' @param newdata patient covariate matrix for prediction.
predict.DTR <- function(treeout, newdata) {
  n <- nrow(newdata)
  predicts <- rep(NA, n)
  
  # treeout is supposed to be a matrix
  # if there is no split
  if (length(treeout) == 5) {
    predicts <- rep(treeout[5], n)
  } else{
    # if there are splits
    treeout <- as.data.frame(treeout)
    newdata <- as.data.frame(newdata)
    
    for (i in 1:n) {
      nd <- 1
      while (is.na(treeout$trt[treeout$node == nd])) {
        if (newdata[i, treeout$X[treeout$node == nd]] <= treeout$cutoff[treeout$node == nd]) {
          nd = 2 * nd #yes proceed first
        } else{
          nd = 2 * nd + 1
        }
      }
      predicts[i] <- treeout$trt[treeout$node == nd]
    }
  }
  return(predicts)
}

######Simulation functions#######

#' Simulation for two stages
#' Returns the opt% and expected mean counterfactual outcome for 
#' Case 1 with text, Case 1 without text,Case 2 with text complete,
#' Case 2 with text imputed, Case 2 without text complete, and
#' Case 2 without text imputed, respectively.
#' 
#' @param i seed.
#' @param N number of patients.
#' @param text_data input data containing text information.

sim_text_2stage <- function(i, N, text_data) {
  set.seed(i)
  train_id <- sample(1:nrow(text_data), N)
  N2 <- 1000 # sample size of test data
  rem_id <- c(1:nrow(text_data))[-train_id]
  test_id <- sample(rem_id, N2)
  iter <- 5 # replication
  
  
  select1 <- select2 <- selects <- rep(NA, iter) # percent of optimality
  EYs <- rep(NA, iter) # estimated mean counterfactual outcome
  
  
  set.seed(i * 50)
  x1 <- text_data$smk_extract[train_id]#x1=smk
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  x4 <- text_data$wt_full[train_id]#x4=wt
  x5 <- rnorm(N)
  X0 <- cbind(x1, x2, x3, x4, x5)
  
  #Variables are coded as (names the paper = names in the code):
  #X1 = X2,
  #X2 = X3,
  #X3 = X5,
  #X4 = X4,
  #X5 = X1
  
  ## stage 1 data simulation
  # simulate A1, stage 1 treatment with K1=3
  pi10 <-
    rep(1, N)
  pi11 <- exp(0.005 * x4 + 0.5 * x1)
  pi12 <- exp(0.5 * x5 - 0.5 * x1)
  # weights matrix
  matrix.pi1 <- cbind(pi10, pi11, pi12)
  
  A1 <- A.sim(matrix.pi1)
  class.A1 <- sort(unique(A1))
  
  # propensity stage 1
  pis1.hat <- M.propen(A1, cbind(x1, x4, x5))
  
  # simulate stage 1 optimal g1.opt
  g1.opt <- (x1 < 1) * ((x2 > -0.5) + (x2 > 0.5))
  # stage 1 outcome
  R1 <-
    exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) * (A1 - g1.opt) ^ 2) + rnorm(N, 0, 1)
  
  ## stage 2 data simulation
  # A2, stage 2 treatment with K2=3
  pi20 <- rep(1, N)
  pi21 <- exp(0.2 * R1 - 0.5)
  pi22 <- exp(0.5 * x2)
  matrix.pi2 <- cbind(pi20, pi21, pi22)
  A2 <- A.sim(matrix.pi2)
  class.A2 <- sort(unique(A2))
  
  # propensity stage 2
  pis2.hat <- M.propen(A2, cbind(R1, x2))
  
  # optimal g2.opt
  g2.opt <- (x3 > -1) * ((R1 > 0) + (R1 > 2))
  # stage 2 outcome R2
  R2 <-
    exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) * (A2 - g2.opt) ^ 2) + rnorm(N, 0, 1)
  
  
  ## Different versions of the data
  data_full <- data.frame(X0, A1, A2, R1, R2)
  #data with missing
  text_data_train <- text_data[train_id, ]
  data_miss <- data_full
  data_error <- data_full
  data_miss$x4obs <-
    ifelse(is.na(text_data_train$weightlb), data_miss$x4, NA)
  #- this is the structured only wt(excluding all information in text)
  #Introduce missing: missing is 90% in cells with x2>-0.5 and 10% randomly in other cells
  x4misprob <- sapply(1:nrow(data_full), function(x)
    rbinom(1, 1, prob = 0.1 + 0.5 * (
      data_full$x1[x] < 1 & is.na(data_miss$x4obs[x])
    )))
  
  data_miss$x4obs_sim <- ifelse(x4misprob == 0, data_miss$x4, NA)
  data_miss$x4obstext_sim <-
    ifelse(!is.na(text_data_train$weightlb),
           data_miss$x4,
           data_miss$x4obs_sim)
  
  # Do imputations using missForest
  invisible(capture.output(
    data_miss$x4obs_sim_imp <- missForest(data_miss %>%
                  select(x1, x2, x3, x4obs_sim, A1, A2, R1, R2))$ximp$x4obs_sim
  ))
  invisible(
    capture.output(
      data_miss$x4obstext_sim_imp <- missForest(data_miss %>%
        select(x1, x2, x3, x4obstext_sim, A1, A2, R1, R2))$ximp$x4obstext_sim
  ))
  
  data_miss_str <-
    data_miss %>% select(x1, x2, x3, x4obs_sim_imp, x5, A1, A2, R1, R2)
  names(data_miss_str)[4] = "x4"
  
  data_miss_txt <-
    data_miss %>% select(x1, x2, x3, x4obstext_sim_imp, x5, A1, A2, R1, R2)
  names(data_miss_txt)[4] = "x4"
  
  # Complete cases only.
  data_miss_str_complete <- data_miss %>%
    select(x1, x2, x3, x4obs_sim, x5, A1, A2, R1, R2)
  names(data_miss_str_complete)[4] = "x4"
  data_miss_str_complete = data_miss_str_complete[complete.cases(data_miss_str_complete), ]
  
  data_miss_txt_complete <-
    data_miss %>% select(x1, x2, x3, x4obstext_sim, x5, A1, A2, R1, R2)
  names(data_miss_txt_complete)[4] = "x4"
  data_miss_txt_complete = data_miss_txt_complete[complete.cases(data_miss_txt_complete), ]
  
  #Introduce error: when observed in text then there are 60% chance of
  #entered in kg in stead of lb.
  x4errprob <- sapply(1:nrow(data_error), function(x)
    rbinom(1, 1, prob = 0.1 * ((data_error$x1[x] < 1) &
                                 is.na(data_miss$x4obs[x])
    )))
  
  data_error$x4obs_sim <-
    ifelse(x4errprob == 0, data_error$x4, data_error$x4 * 100)
  
  data_err_str <-
    data_error %>% select(x1, x2, x3, x4obs_sim, x5, A1, A2, R1, R2)
  names(data_err_str)[4] = "x4"
  
  
  ## FULL DATA case
  
  # conditional mean model using linear regression
  REG2 <- Reg.mu(Y = R2,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2 <-
    DTRtree(
      R2,
      A2,
      H = cbind(X0, A1, R1),
      pis.hat = pis2.hat,
      mus.reg = mus2.reg,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree <- predict.DTR(tree2, newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(R2 ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2))
    mus2.RF[, d] <-
    predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  
  for (m in 1:N) {
    E.R2.tree[m] <- R2[m] + mus2.RF[m, g2.tree[m] + 1] - mus2.RF[m, A2[m] +
                                                                   1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  
  tree1 <-
    DTRtree(
      PO.tree,
      A1,
      H = X0,
      pis.hat = pis1.hat,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  
  
  ## MISSING WO TEXT IMPUTED case
  X0 = data_miss_str %>% select(x2, x3, x4, x5)
  A1 = data_miss_str$A1
  A2 = data_miss_str$A2
  R1 = data_miss_str$R1
  R2 = data_miss_str$R2
  N = length(R1)
  
  pis2.hat <- M.propen(A2, cbind(R1, X0$x2))
  REG2 <- Reg.mu(Y = R2,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2_miss_str <-
    DTRtree(
      R2,
      A2,
      H = cbind(X0, A1, R1),
      pis.hat = pis2.hat,
      mus.reg = mus2.reg,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree_miss_str <-
    predict.DTR(tree2_miss_str, newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(R2 ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2))
    mus2.RF[, d] <-
    predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  
  for (m in 1:N) {
    E.R2.tree[m] <-
      R2[m] + mus2.RF[m, g2.tree_miss_str[m] + 1] - mus2.RF[m, A2[m] + 1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  pis1.hat <- M.propen(A1, X0)
  tree1_miss_str <-
    DTRtree(
      PO.tree,
      A1,
      H = X0,
      pis.hat = pis1.hat,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  ## MISSING WO TEXT COMPLETE case
  X0 = data_miss_str_complete %>% select(x2, x3, x4, x5)
  A1 = data_miss_str_complete$A1
  A2 = data_miss_str_complete$A2
  R1 = data_miss_str_complete$R1
  R2 = data_miss_str_complete$R2
  N = length(R1)
  
  pis2.hat <- M.propen(A2, cbind(R1, X0$x2))
  REG2 <- Reg.mu(Y = R2,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2_miss_str_complete <- DTRtree(
    R2,
    A2,
    H = cbind(X0, A1, R1),
    pis.hat = pis2.hat,
    mus.reg = mus2.reg,
    lambda.pct = 0.02,
    minsplit = max(0.05 * N, 20)
  )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree_miss_str_complete <- predict.DTR(tree2_miss_str_complete,
                                           newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(R2 ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2))
    mus2.RF[, d] <-
    predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  
  for (m in 1:N) {
    E.R2.tree[m] <-
      R2[m] + mus2.RF[m, g2.tree_miss_str[m] + 1] - mus2.RF[m, A2[m] + 1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  pis1.hat <- M.propen(A1, X0)
  tree1_miss_str_complete <-
    DTRtree(
      PO.tree,
      A1,
      H = X0,
      pis.hat = pis1.hat,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  
  ## MISSING With TEXT IMPUTED case
  X0 = data_miss_txt %>% select(x1, x2, x3, x4, x5)
  A1 = data_miss_txt$A1
  A2 = data_miss_txt$A2
  R1 = data_miss_txt$R1
  R2 = data_miss_txt$R2
  N = length(R1)
  
  pis2.hat <- M.propen(A2, cbind(R1, X0$x2))
  REG2 <- Reg.mu(Y = R2,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2_miss_text <-
    DTRtree(
      R2,
      A2,
      H = cbind(X0, A1, R1),
      pis.hat = pis2.hat,
      mus.reg = mus2.reg,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree_miss_text <-
    predict.DTR(tree2_miss_text, newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(R2 ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2))
    mus2.RF[, d] <-
    predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  
  for (m in 1:N) {
    E.R2.tree[m] <-
      R2[m] + mus2.RF[m, g2.tree_miss_text[m] + 1] - mus2.RF[m, A2[m] + 1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  pis1.hat <- M.propen(A1, cbind(X0$x1, X0$x4, X0$x5))
  tree1_miss_text <-
    DTRtree(
      PO.tree,
      A1,
      H = X0,
      pis.hat = pis1.hat,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  ## MISSING With TEXT COMPLETE case
  X0 = data_miss_txt_complete %>% select(x1, x2, x3, x4, x5)
  A1 = data_miss_txt_complete$A1
  A2 = data_miss_txt_complete$A2
  R1 = data_miss_txt_complete$R1
  R2 = data_miss_txt_complete$R2
  N = length(R1)
  
  pis2.hat <- M.propen(A2, cbind(R1, X0$x2))
  REG2 <- Reg.mu(Y = R2,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2_miss_text_complete <-
    DTRtree(
      R2,
      A2,
      H = cbind(X0, A1, R1),
      pis.hat = pis2.hat,
      mus.reg = mus2.reg,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree_miss_text_complete <-
    predict.DTR(tree2_miss_text_complete,
                newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(R2 ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2)) {
    mus2.RF[, d] <-
      predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  }
  for (m in 1:N) {
    E.R2.tree[m] <-
      R2[m] + mus2.RF[m, g2.tree_miss_text_complete[m] + 1] - mus2.RF[m, A2[m] +
                                                                        1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  pis1.hat <- M.propen(A1, cbind(X0$x1, X0$x4, X0$x5))
  tree1_miss_text_complete <-
    DTRtree(
      PO.tree,
      A1,
      H = X0,
      pis.hat = pis1.hat,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  
  ## ERROR WO TEXT case
  X0 = data_err_str %>% select(x2, x3, x4, x5)
  A1 = data_err_str$A1
  A2 = data_err_str$A2
  R1 = data_err_str$R1
  R2 = data_err_str$R2
  N = length(R1)
  
  pis2.hat <- M.propen(A2, cbind(R1, X0$x2))
  REG2 <- Reg.mu(Y = R2,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2_err_str <-
    DTRtree(
      R2,
      A2,
      H = cbind(X0, A1, R1),
      pis.hat = pis2.hat,
      mus.reg = mus2.reg,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree_err_str <-
    predict.DTR(tree2_err_str, newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(R2 ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2))
    mus2.RF[, d] <-
    predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  
  for (m in 1:N) {
    E.R2.tree[m] <-
      R2[m] + mus2.RF[m, g2.tree_err_str[m] + 1] - mus2.RF[m, A2[m] + 1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  pis1.hat <- M.propen(A1, X0)
  tree1_err_str <- DTRtree(
    PO.tree,
    A1,
    H = X0,
    pis.hat = pis1.hat,
    lambda.pct = 0.02,
    minsplit = max(0.05 * N, 20)
  )
  
  
  ## Prediction using new data
  set.seed(10000)
  x1 <- text_data$smk_extract[test_id]#x1=smoke
  x2 <- rnorm(N2)
  x3 <- rnorm(N2)
  x4 <- text_data$wt_full[test_id]#x4=weight
  x5 <- rnorm(N2)
  X0 <- cbind(x1, x2, x3, x4, x5)
  
  
  # true optimal for regime at stage 1
  g1.opt <- (x1 < 1) * ((x2 > -0.5) + (x2 > 0.5))
  
  R1 <- exp(1.5 + 0.003 * x4) + rnorm(N2, 0, 1)
  
  R2 <- exp(1.18 + 0.2 * x2) + rnorm(N2, 0, 1)
  
  ## Stage 1 prediction
  
  # predict selection %
  
  g1.tree <- predict.DTR(tree1, newdata = data.frame(X0))
  g1.tree_test_miss_str <-
    predict.DTR(tree1_miss_str, newdata = data.frame(X0[, -1]))
  g1.tree_test_miss_str_complete <-
    predict.DTR(tree1_miss_str_complete, newdata = data.frame(X0[, -1]))
  g1.tree_test_miss_text <-
    predict.DTR(tree1_miss_text, newdata = data.frame(X0))
  g1.tree_test_miss_text_complete <-
    predict.DTR(tree1_miss_text_complete, newdata = data.frame(X0))
  g1.tree_test_err_str <-
    predict.DTR(tree1_err_str, newdata = data.frame(X0[, -1]))
  
  
  select1 <- mean(g1.tree == g1.opt)
  select1_miss_str <- mean(g1.tree_test_miss_str == g1.opt)
  select1_miss_str_complete <-
    mean(g1.tree_test_miss_str_complete == g1.opt)
  select1_miss_text <- mean(g1.tree_test_miss_text == g1.opt)
  select1_miss_text_complete <-
    mean(g1.tree_test_miss_text_complete == g1.opt)
  select1_err_str <- mean(g1.tree_test_err_str == g1.opt)
  
  
  
  R1.tree <-
    exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) * (g1.tree - g1.opt) ^ 2) + rnorm(N2, 0, 1)
  R1.tree_miss_str <-
    exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) * (g1.tree_test_miss_str - g1.opt)^2) + 
    rnorm(N2, 0, 1)
  R1.tree_miss_str_complete <- exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) *
                                     (g1.tree_test_miss_str_complete - g1.opt) ^
                                     2) +
    rnorm(N2, 0, 1)
  R1.tree_miss_text <- exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) *
                             (g1.tree_test_miss_text - g1.opt) ^ 2) + rnorm(N2, 0, 1)
  R1.tree_miss_text_complete <- exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) *
                                      (g1.tree_test_miss_text_complete - g1.opt) ^
                                      2) +
    rnorm(N2, 0, 1)
  R1.tree_err_str <-
    exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) * (g1.tree_test_err_str - g1.opt) ^
          2) +
    rnorm(N2, 0, 1)
  
  ## Stage 2 prediction
  
  g2.tree <-
    predict.DTR(tree2, newdata = data.frame(X0, A1 = g1.tree, R1 = R1.tree))
  g2.tree_test_miss_str <- predict.DTR(tree2_miss_str,
                                       newdata = data.frame(X0[, -1],
                                                            A1 = g1.tree_test_miss_str,
                                                            R1 = R1.tree_miss_str))
  g2.tree_test_miss_str_complete <-
    predict.DTR(
      tree2_miss_str_complete,
      newdata = data.frame(X0[, -1],
                           A1 = g1.tree_test_miss_str_complete,
                           R1 = R1.tree_miss_str_complete)
    )
  g2.tree_test_miss_text <- predict.DTR(tree2_miss_text,
                                        newdata = data.frame(X0,
                                                             A1 = g1.tree_test_miss_text,
                                                             R1 = R1.tree_miss_text))
  g2.tree_test_miss_text_complete <-
    predict.DTR(
      tree2_miss_text_complete,
      newdata = data.frame(X0,
                           A1 = g1.tree_test_miss_text_complete,
                           R1 = R1.tree_miss_text_complete)
    )
  g2.tree_test_err_str <- predict.DTR(tree2_err_str,
                                      newdata = data.frame(X0[, -1],
                                                           A1 = g1.tree_test_err_str,
                                                           R1 = R1.tree_err_str))
  # true optimal for regime at stage 2
  g2.opt.tree <- (x3 > -1) * ((R1.tree > 0) + (R1.tree > 2))
  
  select2 <- mean(g2.tree == g2.opt.tree)
  select2_miss_str <- mean(g2.tree_test_miss_str == g2.opt.tree)
  select2_miss_str_complete <-
    mean(g2.tree_test_miss_str_complete == g2.opt.tree)
  select2_miss_text <- mean(g2.tree_test_miss_text == g2.opt.tree)
  select2_miss_text_complete <-
    mean(g2.tree_test_miss_text_complete == g2.opt.tree)
  select2_err_str <- mean(g2.tree_test_err_str == g2.opt.tree)
  
  
  selects <- mean(g1.tree == g1.opt & g2.tree == g2.opt.tree)
  selects_miss_str <-
    mean(g1.tree_test_miss_str == g1.opt &
           g2.tree_test_miss_str == g2.opt.tree)
  selects_miss_str_complete <-
    mean(
      g1.tree_test_miss_str_complete == g1.opt &
        g2.tree_test_miss_str_complete == g2.opt.tree
    )
  selects_miss_text <-
    mean(g1.tree_test_miss_text == g1.opt &
           g2.tree_test_miss_text == g2.opt.tree)
  selects_miss_text_complete <-
    mean(
      g1.tree_test_miss_text_complete == g1.opt &
        g2.tree_test_miss_text_complete == g2.opt.tree
    )
  selects_err_str <-
    mean(g1.tree_test_err_str == g1.opt &
           g2.tree_test_err_str == g2.opt.tree)
  
  R2.tree <-
    exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) * (g2.tree - g2.opt.tree) ^ 2) +
    rnorm(N2, 0, 1)
  R2.tree_miss_str <-
    exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) * (g2.tree_test_miss_str - g2.opt.tree) ^
          2) +
    rnorm(N2, 0, 1)
  R2.tree_miss_str_complete <- exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) *
                                     (g2.tree_test_miss_str_complete - g2.opt.tree) ^
                                     2) +
    rnorm(N2, 0, 1)
  R2.tree_miss_text <- exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) *
                             (g2.tree_test_miss_text - g2.opt.tree) ^ 2) +
    rnorm(N2, 0, 1)
  R2.tree_miss_text_complete <- exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) *
                                      (g2.tree_test_miss_text_complete - g2.opt.tree) ^
                                      2) +
    rnorm(N2, 0, 1)
  
  R2.tree_err_str <- exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) *
                           (g2.tree_test_err_str - g2.opt.tree) ^ 2) + rnorm(N2, 0, 1)
  
  EY_full <- mean(R1.tree + R2.tree)
  EY_miss_str <- mean(R1.tree_miss_str + R2.tree_miss_str)
  EY_miss_str_complete <-
    mean(R1.tree_miss_str_complete + R2.tree_miss_str_complete)
  EY_miss_text <- mean(R1.tree_miss_text + R2.tree_miss_text)
  EY_miss_text_complete <-
    mean(R1.tree_miss_text_complete + R2.tree_miss_text_complete)
  EY_err_str <- mean(R1.tree_err_str + R2.tree_err_str)
  return(
    c(
      selects * 100,
      selects_err_str * 100,
      selects_miss_text_complete * 100,
      selects_miss_text * 100,
      selects_miss_str_complete * 100,
      selects_miss_str * 100,
      EY_full,
      EY_err_str,
      EY_miss_text_complete,
      EY_miss_text,
      EY_miss_str_complete,
      EY_miss_str
    )
  )
}

#' Simulation for three stages
#' Returns the opt% and expected mean counterfactual outcome for 
#' Case 1 with text, Case 1 without text,
#' Case 2 with text imputed, and
#' Case 2 without text imputed, respectively.
#' 
#' @param i seed.
#' @param N number of patients.
#' @param text_data input data containing text information.
sim_text_3stage <- function(i, N, text_data) {
  set.seed(i)
  train_id <- sample(1:nrow(text_data), N)
  N2 <- 1000 # sample size of test data
  rem_id <- c(1:nrow(text_data))[-train_id]
  test_id <- sample(rem_id, N2)
  #iter <- 5 # replication
  
  
  #select1 <- select2 <- selects <- rep(NA, iter) # percent of optimality
  #EYs <- rep(NA, iter) # estimated mean counterfactual outcome
  
  set.seed(i * 50)
  x1 <- text_data$smk_extract[train_id]#x1=smoke
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  x4 <- text_data$wt_full[train_id]#x4=weight
  x5 <- rnorm(N)
  X0 <- cbind(x1, x2, x3, x4, x5)
  
  #Variables are coded as (names the paper = names in the code):
  #X1 = X2,
  #X2 = X3,
  #X3 = X5,
  #X4 = X4,
  #X5 = X1
  
  ## stage 1 data simulation
  # simulate A1, stage 1 treatment with K1=3
  pi10 <-
    rep(1, N)
  pi11 <- exp(0.005 * x4 + 0.5 * x1)
  pi12 <- exp(0.5 * x5 - 0.5 * x1)
  # weights matrix
  matrix.pi1 <- cbind(pi10, pi11, pi12)
  
  A1 <- A.sim(matrix.pi1)
  class.A1 <- sort(unique(A1))
  
  # propensity stage 1
  pis1.hat <- M.propen(A1, cbind(x1, x4, x5))
  
  # simulate stage 1 optimal g1.opt
  g1.opt <- (x1 < 1) * ((x2 > -0.5) + (x2 > 0.5))
  # stage 1 outcome
  R1 <-
    exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) * (A1 - g1.opt) ^ 2) + rnorm(N, 0, 1)
  
  ## stage 2 data simulation
  # A2, stage 2 treatment with K2=3
  pi20 <- rep(1, N)
  pi21 <- exp(0.2 * R1 - 0.5)
  pi22 <- exp(0.5 * x2)
  matrix.pi2 <- cbind(pi20, pi21, pi22)
  A2 <- A.sim(matrix.pi2)
  class.A2 <- sort(unique(A2))
  
  # propensity stage 2
  pis2.hat <- M.propen(A2, cbind(R1, x2))
  
  # optimal g2.opt
  g2.opt <- (x3 > -1) * ((R1 > 0) + (R1 > 2))
  # stage 2 outcome R2
  R2 <-
    exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) * (A2 - g2.opt) ^ 2) + rnorm(N, 0, 1)
  
  ## stage 3 data simulation
  # A3, stage 3 treatment with K3=2 (0 or 1)
  pi30 <- rep(1, N)
  pi31 <- exp(0.1 * R2 + 0.3 * x3) # prob A3 being 1.
  
  matrix.pi3 <- cbind(pi30, pi31)
  A3 <- A.sim(matrix.pi3)
  class.A3 <- sort(unique(A3))
  
  # propensity stage 3
  pis3.hat <- M.propen(A3, cbind(R2, x3))
  
  # optimal g3.opt
  g3.opt <- ifelse((R2 > 0), 1, 0)
  # stage 3 outcome R3
  R3 <- exp(0.1 * x2 - (A3 - g3.opt) ^ 2) + rnorm(N, 0, 1)
  
  #Different versions of the data
  data_full <- data.frame(X0, A1, A2, A3, R1, R2, R3)
  #data with missing
  text_data_train <- text_data[train_id, ]
  data_miss <- data_full
  data_error <- data_full
  data_miss$x4obs <-
    ifelse(is.na(text_data_train$weightlb), data_miss$x4, NA)
  #- this is the structured only wt(excluding all information in text)
  #Introduce missing: missing is 90% in cells with x2>-0.5 and 10% randomly in other cells
  x4misprob <- sapply(1:nrow(data_full), function(x)
    rbinom(1, 1, prob = 0.1 + 0.5 * (
      data_full$x1[x] < 1 & is.na(data_miss$x4obs[x])
    )))
  
  data_miss$x4obs_sim <- ifelse(x4misprob == 0, data_miss$x4, NA)
  data_miss$x4obstext_sim <-
    ifelse(!is.na(text_data_train$weightlb),
           data_miss$x4,
           data_miss$x4obs_sim)
  
  # Do imputations missForest
  invisible(capture.output(
    data_miss$x4obs_sim_imp <- missForest(data_miss %>%
                                            select(x1, x2, x3, x4obs_sim, A1, A2, A3, R1, R2, R3))$ximp$x4obs_sim
  ))
  invisible(capture.output(
    data_miss$x4obstext_sim_imp <- missForest(
      data_miss %>%
        select(x1, x2, x3, x4obstext_sim, A1, A2, A3, R1, R2, R3)
    )$ximp$x4obstext_sim
  ))
  
  data_miss_str <-
    data_miss %>% select(x1, x2, x3, x4obs_sim_imp, x5, A1, A2, A3, R1, R2, R3)
  names(data_miss_str)[4] = "x4"
  
  data_miss_txt <-
    data_miss %>% select(x1, x2, x3, x4obstext_sim_imp, x5, A1, A2, A3, R1, R2, R3)
  names(data_miss_txt)[4] = "x4"
  
  #Introduce error: when observed in text then there are 60% chance of
  #entered in kg in stead of lb.
  x4errprob <- sapply(1:nrow(data_error), function(x)
    rbinom(1, 1, prob = 0.1 * ((data_error$x1[x] < 1) &
                                 is.na(data_miss$x4obs[x])
    )))
  
  data_error$x4obs_sim <-
    ifelse(x4errprob == 0, data_error$x4, data_error$x4 * 100)
  
  data_err_str <-
    data_error %>% select(x1, x2, x3, x4obs_sim, x5, A1, A2, A3, R1, R2, R3)
  names(data_err_str)[4] = "x4"
  
  
  ## FULL DATA case
  # stage 3 conditional mean model using linear regression
  REG3 <- Reg.mu(Y = R3,
                 As = cbind(A1, A2, A3),
                 H = cbind(X0, R1, R2))
  mus3.reg <- REG3$mus.reg
  
  tree3 <- DTRtree(
    R3,
    A3,
    H = cbind(X0, A1, A2, R1, R2),
    pis.hat = pis3.hat,
    mus.reg = mus3.reg,
    m.method = "AIPW",
    lambda.pct = 0.02,
    minsplit = max(0.05 * N, 20)
  )
  
  # expected optimal stage 2 outcome
  E.R3.tree <- rep(NA, N)
  
  # estimated optimal regime
  g3.tree <- predict.DTR(tree3, newdata = data.frame(X0, A1, A2, R1, R2))
  
  # random forest for the estimated mean
  RF3 <- randomForest(R3 ~ ., data = data.frame(A3, X0, A1, A2, R1, R2))
  mus3.RF <- matrix(NA, N, length(class.A3))
  for (d in 1L:length(class.A3)) {
    mus3.RF[, d] <-
      predict(RF3, newdata = data.frame(A3 = rep(class.A3[d], N), X0, A1, A2, R1, R2))
  }
  for (m in 1:N) {
    E.R3.tree[m] <-
      R3[m] + mus3.RF[m, g3.tree[m] + 1] - mus3.RF[m, A3[m] + 1] #starts from 0 vs 1
  }
  
  # pseudo outcomes stage 3 optimized.
  PO3.tree <- R2 + E.R3.tree
  
  #Stage 2
  # conditional mean model using linear regression
  REG2 <- Reg.mu(Y = PO3.tree,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2 <-
    DTRtree(
      PO3.tree,
      A2,
      H = cbind(X0, A1, R1),
      pis.hat = pis2.hat,
      mus.reg = mus2.reg,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree <- predict.DTR(tree2, newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(PO3.tree ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2))
    mus2.RF[, d] <-
    predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  
  for (m in 1:N) {
    E.R2.tree[m] <-
      PO3.tree[m] + mus2.RF[m, g2.tree[m] + 1] - mus2.RF[m, A2[m] + 1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  
  tree1 <-
    DTRtree(
      PO.tree,
      A1,
      H = X0,
      pis.hat = pis1.hat,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  
  
  ## MISSING WO TEXT IMPUTED case
  X0 = data_miss_str %>% select(x2, x3, x4, x5)
  A1 = data_miss_str$A1
  A2 = data_miss_str$A2
  A3 = data_miss_str$A3
  R1 = data_miss_str$R1
  R2 = data_miss_str$R2
  R3 = data_miss_str$R3
  N = length(R1)
  
  #Stage 3
  pis3.hat <- M.propen(A3, cbind(R2, X0$x3))
  REG3 <- Reg.mu(Y = R3,
                 As = cbind(A1, A2, A3),
                 H = cbind(X0, R1, R2))
  mus3.reg <- REG3$mus.reg
  
  tree3_miss_str <-
    DTRtree(
      R3,
      A3,
      H = cbind(X0, A1, A2, R1, R2),
      pis.hat = pis3.hat,
      mus.reg = mus3.reg,
      m.method = "AIPW",
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R3.tree <- rep(NA, N)
  
  # estimated optimal regime
  g3.tree_miss_str <-
    predict.DTR(tree3_miss_str, newdata = data.frame(X0, A1, A2, R1, R2))
  
  # random forest for the estimated mean
  RF3 <- randomForest(R3 ~ ., data = data.frame(A3, X0, A1, A2, R1, R2))
  mus3.RF <- matrix(NA, N, length(class.A3))
  for (d in 1L:length(class.A3))
    mus3.RF[, d] <-
    predict(RF3, newdata = data.frame(A3 = rep(class.A3[d], N), X0, A1, A2, R1, R2))
  
  for (m in 1:N) {
    E.R3.tree[m] <-
      R3[m] + mus3.RF[m, g3.tree_miss_str[m] + 1] - mus3.RF[m, A3[m] + 1]
  }
  
  # pseudo outcomes stage 3 optimized
  PO3.tree <- R2 + E.R3.tree
  
  # STAGE 2
  pis2.hat <- M.propen(A2, cbind(R1, X0$x2))
  REG2 <- Reg.mu(Y = PO3.tree,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2_miss_str <- DTRtree(
    PO3.tree,
    A2,
    H = cbind(X0, A1, R1),
    pis.hat = pis2.hat,
    mus.reg = mus2.reg,
    lambda.pct = 0.02,
    minsplit = max(0.05 * N, 20)
  )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree_miss_str <-
    predict.DTR(tree2_miss_str, newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(PO3.tree ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2))
    mus2.RF[, d] <-
    predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  
  for (m in 1:N) {
    E.R2.tree[m] <-
      PO3.tree[m] + mus2.RF[m, g2.tree_miss_str[m] + 1] - mus2.RF[m, A2[m] + 1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  pis1.hat <- M.propen(A1, X0)
  tree1_miss_str <-
    DTRtree(
      PO.tree,
      A1,
      H = X0,
      pis.hat = pis1.hat,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  ## MISSING With TEXT IMPUTED case
  X0 = data_miss_txt %>% select(x1, x2, x3, x4, x5)
  A1 = data_miss_txt$A1
  A2 = data_miss_txt$A2
  A3 = data_miss_txt$A3
  R1 = data_miss_txt$R1
  R2 = data_miss_txt$R2
  R3 = data_miss_txt$R3
  N = length(R1)
  
  
  #Stage 3
  pis3.hat <- M.propen(A3, cbind(R2, X0$x3))
  REG3 <- Reg.mu(Y = R3,
                 As = cbind(A1, A2, A3),
                 H = cbind(X0, R1, R2))
  mus3.reg <- REG3$mus.reg
  
  tree3_miss_text <- DTRtree(
    R3,
    A3,
    H = cbind(X0, A1, A2, R1, R2),
    pis.hat = pis3.hat,
    mus.reg = mus3.reg,
    m.method = "AIPW",
    lambda.pct = 0.02,
    minsplit = max(0.05 * N, 20)
  )
  
  # expected optimal stage 2 outcome
  E.R3.tree <- rep(NA, N)
  
  # estimated optimal regime
  g3.tree_miss_text <-
    predict.DTR(tree3_miss_text, newdata = data.frame(X0, A1, A2, R1, R2))
  
  # random forest for the estimated mean
  RF3 <- randomForest(R3 ~ ., data = data.frame(A3, X0, A1, A2, R1, R2))
  mus3.RF <- matrix(NA, N, length(class.A3))
  for (d in 1L:length(class.A3))
    mus3.RF[, d] <-
    predict(RF3, newdata = data.frame(A3 = rep(class.A3[d], N), X0, A1, A2, R1, R2))
  
  for (m in 1:N) {
    E.R3.tree[m] <-
      R3[m] + mus3.RF[m, g3.tree_miss_str[m] + 1] - mus3.RF[m, A3[m] + 1]
  }
  
  # pseudo outcomes stage 3 optimized
  PO3.tree <- R2 + E.R3.tree
  
  #Stage 2
  pis2.hat <- M.propen(A2, cbind(R1, X0$x2))
  REG2 <- Reg.mu(Y = PO3.tree,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2_miss_text <- DTRtree(
    PO3.tree,
    A2,
    H = cbind(X0, A1, R1),
    pis.hat = pis2.hat,
    mus.reg = mus2.reg,
    lambda.pct = 0.02,
    minsplit = max(0.05 * N, 20)
  )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree_miss_text <-
    predict.DTR(tree2_miss_text, newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(PO3.tree ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2)) {
    mus2.RF[, d] <-
      predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  }
  for (m in 1:N) {
    E.R2.tree[m] <-
      PO3.tree[m] + mus2.RF[m, g2.tree_miss_text[m] + 1] - mus2.RF[m, A2[m] +
                                                                     1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  pis1.hat <- M.propen(A1, cbind(X0$x1, X0$x4, X0$x5))
  tree1_miss_text <-
    DTRtree(
      PO.tree,
      A1,
      H = X0,
      pis.hat = pis1.hat,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  
  ## ERROR WO TEXT case
  X0 = data_err_str %>% select(x2, x3, x4, x5)
  A1 = data_err_str$A1
  A2 = data_err_str$A2
  A3 = data_err_str$A3
  R1 = data_err_str$R1
  R2 = data_err_str$R2
  R3 = data_err_str$R3
  N = length(R1)
  
  #Stage 3
  pis3.hat <- M.propen(A3, cbind(R2, X0$x3))
  REG3 <- Reg.mu(Y = R3,
                 As = cbind(A1, A2, A3),
                 H = cbind(X0, R1, R2))
  mus3.reg <- REG3$mus.reg
  
  tree3_err_str <-
    DTRtree(
      R3,
      A3,
      H = cbind(X0, A1, A2, R1, R2),
      pis.hat = pis3.hat,
      mus.reg = mus3.reg,
      m.method = "AIPW",
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R3.tree <- rep(NA, N)
  
  # estimated optimal regime
  g3.tree_err_str <-
    predict.DTR(tree3_err_str, newdata = data.frame(X0, A1, A2, R1, R2))
  
  # random forest for the estimated mean
  RF3 <- randomForest(R3 ~ ., data = data.frame(A3, X0, A1, A2, R1, R2))
  mus3.RF <- matrix(NA, N, length(class.A3))
  for (d in 1L:length(class.A3))
    mus3.RF[, d] <-
    predict(RF3, newdata = data.frame(A3 = rep(class.A3[d], N), X0, A1, A2, R1, R2))
  
  for (m in 1:N) {
    E.R3.tree[m] <-
      R3[m] + mus3.RF[m, g3.tree_err_str[m] + 1] - mus3.RF[m, A3[m] + 1]
  }
  
  # pseudo outcomes stage 3 optimized
  PO3.tree <- R2 + E.R3.tree
  
  # Stage 2
  pis2.hat <- M.propen(A2, cbind(R1, X0$x2))
  REG2 <- Reg.mu(Y = PO3.tree,
                 As = cbind(A1, A2),
                 H = cbind(X0, R1))
  mus2.reg <- REG2$mus.reg
  
  tree2_err_str <-
    DTRtree(
      PO3.tree,
      A2,
      H = cbind(X0, A1, R1),
      pis.hat = pis2.hat,
      mus.reg = mus2.reg,
      lambda.pct = 0.02,
      minsplit = max(0.05 * N, 20)
    )
  
  # expected optimal stage 2 outcome
  E.R2.tree <- rep(NA, N)
  
  # estimated optimal regime
  g2.tree_err_str <-
    predict.DTR(tree2_err_str, newdata = data.frame(X0, A1, R1))
  
  # random forest for the estimated mean
  RF2 <- randomForest(PO3.tree ~ ., data = data.frame(A2, X0, A1, R1))
  mus2.RF <- matrix(NA, N, length(class.A2))
  for (d in 1L:length(class.A2))
    mus2.RF[, d] <-
    predict(RF2, newdata = data.frame(A2 = rep(class.A2[d], N), X0, A1, R1))
  
  for (m in 1:N) {
    E.R2.tree[m] <-
      PO3.tree[m] + mus2.RF[m, g2.tree_err_str[m] + 1] - mus2.RF[m, A2[m] + 1]
  }
  
  # pseudo outcomes
  PO.tree <- R1 + E.R2.tree
  pis1.hat <- M.propen(A1, X0)
  tree1_err_str <- DTRtree(
    PO.tree,
    A1,
    H = X0,
    pis.hat = pis1.hat,
    lambda.pct = 0.02,
    minsplit = max(0.05 * N, 20)
  )
  
  
  ## Predictions using test data
  set.seed(10000)
  x1 <- text_data$smk_extract[test_id]#x1=smoke
  x2 <- rnorm(N2)
  x3 <- rnorm(N2)
  x4 <- text_data$wt_full[test_id]#x4=weight
  x5 <- rnorm(N2)
  X0 <- cbind(x1, x2, x3, x4, x5)
  
  
  # true optimal for regime at stage 1
  g1.opt <- (x1 < 1) * ((x2 > -0.5) + (x2 > 0.5))
  
  R1 <- exp(1.5 + 0.003 * x4) + rnorm(N2, 0, 1)
  
  R2 <- exp(1.18 + 0.2 * x2) + rnorm(N2, 0, 1)
  
  R3 <- exp(0.1 * x2) + rnorm(N2, 0, 1)
  
  ## Stage 1 prediction
  
  # predict selection %
  
  g1.tree <- predict.DTR(tree1, newdata = data.frame(X0))
  g1.tree_test_miss_str <-
    predict.DTR(tree1_miss_str, newdata = data.frame(X0[, -1]))
  g1.tree_test_miss_text <-
    predict.DTR(tree1_miss_text, newdata = data.frame(X0))
  g1.tree_test_err_str <-
    predict.DTR(tree1_err_str, newdata = data.frame(X0[, -1]))
  
  
  select1 <- mean(g1.tree == g1.opt)
  select1_miss_str <- mean(g1.tree_test_miss_str == g1.opt)
  select1_miss_text <- mean(g1.tree_test_miss_text == g1.opt)
  select1_err_str <- mean(g1.tree_test_err_str == g1.opt)
  
  
  
  R1.tree <-
    exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) * (g1.tree - g1.opt) ^ 2) + rnorm(N2, 0, 1)
  R1.tree_miss_str <-
    exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) * (g1.tree_test_miss_str - g1.opt) ^
          2) +
    rnorm(N2, 0, 1)
  R1.tree_miss_text <- exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) *
                             (g1.tree_test_miss_text - g1.opt) ^ 2) + rnorm(N2, 0, 1)
  R1.tree_err_str <-
    exp(1.5 + 0.003 * x4 - abs(1.5 * x1 - 2) * (g1.tree_test_err_str - g1.opt) ^
          2) +
    rnorm(N2, 0, 1)
  
  ## Stage 2 prediction
  
  g2.tree <-
    predict.DTR(tree2, newdata = data.frame(X0, A1 = g1.tree, R1 = R1.tree))
  g2.tree_test_miss_str <- predict.DTR(tree2_miss_str,
                                       newdata = data.frame(X0[, -1],
                                                            A1 = g1.tree_test_miss_str,
                                                            R1 = R1.tree_miss_str))
  g2.tree_test_miss_text <- predict.DTR(tree2_miss_text,
                                        newdata = data.frame(X0,
                                                             A1 = g1.tree_test_miss_text,
                                                             R1 = R1.tree_miss_text))
  g2.tree_test_err_str <- predict.DTR(tree2_err_str,
                                      newdata = data.frame(X0[, -1],
                                                           A1 = g1.tree_test_err_str,
                                                           R1 = R1.tree_err_str))
  
  # true optimal for regime at stage 2
  g2.opt.tree <- (x3 > -1) * ((R1.tree > 0) + (R1.tree > 2))
  
  select2 <- mean(g2.tree == g2.opt.tree)
  select2_miss_str <- mean(g2.tree_test_miss_str == g2.opt.tree)
  select2_miss_text <- mean(g2.tree_test_miss_text == g2.opt.tree)
  select2_err_str <- mean(g2.tree_test_err_str == g2.opt.tree)
  
  R2.tree <-
    exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) * (g2.tree - g2.opt.tree) ^ 2) +
    rnorm(N2, 0, 1)
  R2.tree_miss_str <-
    exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) * (g2.tree_test_miss_str - g2.opt.tree) ^
          2) +
    rnorm(N2, 0, 1)
  R2.tree_miss_text <- exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) *
                             (g2.tree_test_miss_text - g2.opt.tree) ^ 2) +
    rnorm(N2, 0, 1)
  R2.tree_err_str <- exp(1.18 + 0.2 * x2 - abs(1.5 * x3 + 2) *
                           (g2.tree_test_err_str - g2.opt.tree) ^ 2) + rnorm(N2, 0, 1)
  
  
  ## Stage 3 prediction
  
  g3.tree <- predict.DTR(tree3,
                         newdata = data.frame(
                           X0,
                           A1 = g1.tree,
                           A2 = g2.tree,
                           R1 = R1.tree,
                           R2 = R2.tree
                         ))
  g3.tree_test_miss_str <- predict.DTR(
    tree3_miss_str,
    newdata = data.frame(
      X0[, -1],
      A1 = g1.tree_test_miss_str,
      A2 = g2.tree_test_miss_str,
      R1 = R1.tree_miss_str,
      R2 = R2.tree_miss_str
    )
  )
  
  g3.tree_test_miss_text <- predict.DTR(
    tree3_miss_text,
    newdata = data.frame(
      X0,
      A1 = g1.tree_test_miss_text,
      A2 = g2.tree_test_miss_text,
      R1 = R1.tree_miss_text,
      R2 = R2.tree_miss_text
    )
  )
  
  g3.tree_test_err_str <- predict.DTR(
    tree3_err_str,
    newdata = data.frame(
      X0[, -1],
      A1 = g1.tree_test_err_str,
      A2 = g2.tree_test_err_str,
      R1 = R1.tree_err_str,
      R2 = R2.tree_err_str
    )
  )
  # true optimal for regime at stage 3
  g3.opt.tree <- ifelse((R2 > 0), 1, 0)
  
  select3 <- mean(g3.tree == g3.opt.tree)
  select3_miss_str <- mean(g3.tree_test_miss_str == g3.opt.tree)
  select3_miss_text <- mean(g3.tree_test_miss_text == g3.opt.tree)
  select3_err_str <- mean(g3.tree_test_err_str == g3.opt.tree)
  
  selects <-
    mean(g1.tree == g1.opt &
           g2.tree == g2.opt.tree & g3.tree == g3.opt.tree)
  selects_miss_str <- mean(
    g1.tree_test_miss_str == g1.opt &
      g2.tree_test_miss_str == g2.opt.tree &
      g3.tree_test_miss_str == g3.opt.tree
  )
  
  selects_miss_text <- mean(
    g1.tree_test_miss_text == g1.opt &
      g2.tree_test_miss_text == g2.opt.tree &
      g3.tree_test_miss_text == g3.opt.tree
  )
  
  selects_err_str <- mean(
    g1.tree_test_err_str == g1.opt &
      g2.tree_test_err_str == g2.opt.tree &
      g3.tree_test_err_str == g3.opt.tree
  )
  
  R3.tree <- exp(0.1 * x2 - (g3.tree - g3.opt.tree) ^ 2) + rnorm(N2, 0, 1)
  R3.tree_miss_str <-
    exp(0.1 * x2 - (g3.tree_test_miss_str - g3.opt.tree) ^ 2) +
    rnorm(N2, 0, 1)
  
  R3.tree_miss_text <-
    exp(0.1 * x2 - (g3.tree_test_miss_text - g3.opt.tree) ^ 2) +
    rnorm(N2, 0, 1)
  
  R3.tree_err_str <-
    exp(0.1 * x2 - (g3.tree_test_err_str - g3.opt.tree) ^ 2) + rnorm(N2, 0, 1)
  
  
  EY_full <- mean(R1.tree + R2.tree + R3.tree)
  EY_miss_str <-
    mean(R1.tree_miss_str + R2.tree_miss_str + R3.tree_miss_str)
  
  EY_miss_text <-
    mean(R1.tree_miss_text + R2.tree_miss_text + R3.tree_miss_text)
  
  EY_err_str <-
    mean(R1.tree_err_str + R2.tree_err_str + R3.tree_err_str)
  EY_opt <- mean(R1 + R2 + R3)
  return(
    c(
      selects * 100,
      selects_err_str * 100,
      selects_miss_text * 100,
      selects_miss_str * 100,
      EY_full,
      EY_err_str,
      EY_miss_text,
      EY_miss_str
    )
  )
}

####Simulation true reward generation####
#' Function for generating true optimal rewards for 2 stage simulation.
#' Returns the mean reward and variance of reward.
#' 
#' @param seed random seed.
#' @param N2 number of records.
mean_var_calc2stage = function(seed, N2) {
  set.seed(seed)
  test_id <- sample(1:nrow(text_data), N2)
  x1 <- text_data$smk_extract[test_id]#x1=smk
  x2 <- rnorm(N2)
  x3 <- rnorm(N2)
  x4 <- text_data$wt_full[test_id]#x4=wt
  x5 <- rnorm(N2)
  X0 <- cbind(x1, x2, x3, x4, x5)
  
  R1 <- exp(1.5 + 0.003 * x4) + rnorm(N2)
  
  R2 <- exp(1.18 + 0.2 * x2) + rnorm(N2)#has noting to do with a1
  
  return(c(mean(R1 + R2), var(R1 + R2)))
}

#' Function for generating true optimal rewards for 3 stage simulation.
#' Returns the mean reward and variance of reward.
#' 
#' @param seed random seed.
#' @param N2 number of records.
mean_var_calc3stage = function(seed, N2) {
  set.seed(seed)
  test_id <- sample(1:nrow(text_data), N2)
  x1 <- text_data$smk_extract[test_id]#x1=smk
  x2 <- rnorm(N2)
  x3 <- rnorm(N2)
  x4 <- text_data$wt_full[test_id]#x4=wt
  x5 <- rnorm(N2)
  X0 <- cbind(x1, x2, x3, x4, x5)
  
  R1 <- exp(1.5 + 0.003 * x4) + rnorm(N2, 0, 1)
  
  R2 <- exp(1.18 + 0.2 * x2) + rnorm(N2, 0, 1)
  
  R3 <- exp(0.1 * x2) + rnorm(N2, 0, 1)
  
  return(c(mean(R1 + R2 + R3), var(R1 + R2 + R3)))
}





