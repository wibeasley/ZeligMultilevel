# Copied & adapted from `model-mixed.R`.
znlsmixed <- setRefClass("Zelig-nls-mixed",
                      fields = list(x.beta = "ANY",
                                    z.b = "ANY",
                                    offset = "ANY",
                                    error = "ANY",
                                    formula.full = "ANY",# Zelig formula
                                    mm.RE = "ANY", # group membership
                                    simtype = "ANY"), # linear or probability 
                      contains = "Zelig")

# Copied & adapted from both `model-mixed.R` and `model-ls-mixed.R`.
znlsmixed$methods(
  initialize = function() {
    # Copied & adapted from `model-mixed.R`
    callSuper()
    .self$fn <- quote(lme4::nlmer)
    .self$packageauthors <- "TBD"
    .self$modelauthors <- "TBD" 
    .self$year <- 2016
    .self$mm.RE <- NULL
    .self$acceptweights <- TRUE
    
    # Copied & adapted from `model-ls-mixed.R`.
    .self$name <- "nls.mixed"
    .self$fn <- quote(lme4::nlmer)
    .self$category <- "continuous"
    .self$simtype <- "linear"
    .self$wrapper <- "nls.mixed"
  }
)

# Copied & adapted from `model-mixed.R`.
znlsmixed$methods(
  set = function(...) {
    .self$mm.RE <- NULL # reset group membership
    s <- list(...)
    group <- names(ranef(.self$zelig.out$z.out[[1]]))
    print(group)
    set.RE <- intersect(names(s), group)
    if (length(set.RE) > 0)
      .self$mm.RE <- as.data.frame(s[set.RE])
    callSuper(...)
  }
)

# Copied & adapted from `model-mixed.R`.
znlsmixed$methods(
  param = function(z.out) {
    #changing function from arm::sim to myArmSim from myArmSim.R file
    return(list(simparam = myArmSim(z.out, .self$num), simalpha = z.out))
  }
)

# Copied & adapted from `model-ls-mixed.R`.
znlsmixed$methods(
  zelig = function(formula, start, data, ..., weights = NULL, by = NULL) {
    .self$zelig.call <- match.call(expand.dots = TRUE)
    .self$model.call <- match.call(expand.dots = TRUE)
    callSuper(formula = formula, data = data, ..., weights = weights, by = by)
    .self$formula.full <- .self$formula # fixed and random effects
    .self$formula <- as.formula(paste0(.self$formula.full[[2]][[2]],"~",setdiff(all.vars(.self$formula.full[[2]][[3]]), names(.self$zelig.call$start)),"-1"))
    #TURNOFF  
    #.self$formula <- formula(.self$zelig.out$z.out[[1]], fixed.only = TRUE) # fixed effects only
  }
)

# Copied & adapted from `model-ls-mixed.R`.
znlsmixed$methods(
  qi = function(simparam, mm) {
    regression <- simparam$simalpha;
    sims <- simparam$simparam;
    
    sims.tmp <- sims
    ## Check if group memberships are specified
    if (!is.null(.self$mm.RE)){
      for (group in names(.self$mm.RE)) {
        sims@ranef[[group]][, , ] <- 0
        sims@ranef[[group]][, unlist(.self$mm.RE[group]), ] <- sims.tmp@ranef[[group]][, unlist(.self$mm.RE[group]), ]
      }
    }
    
    numSimulations <- dim(sims@fixef)[1];
    devcomp <- getME(regression, "devcomp");
    dims <- devcomp$dims;#list(n=regression$dims$N,q=regression$dims$ngrps[1]*regression$dims$qvec[1],reTrms=regression$dims)
    
    numRanef  <- dims[["q"]];
    numLevels <- dims[["reTrms"]];
    
    simulatedRanef <- matrix(0, numRanef, numSimulations);
    
    index <- 0;
    for (i in 1:length(sims@ranef)) {
      #i <- 1
      levelSims <- sims@ranef[[i]];
      numCoefficientsPerLevel <- dim(levelSims)[2];
      numGroupsPerLevel <- dim(levelSims)[3];
      for (j in 1:numCoefficientsPerLevel) {
        #j <- 1
        ranefRange <- index + 1:numGroupsPerLevel;
        index <- index + numGroupsPerLevel;
        
        simulatedRanef[ranefRange,] <- t(levelSims[,j,]);
      }
    }
    
    X <- getME(regression, "X");
    #TURNOFF
    #X <- matrix(rep(mm, length(X)), nrow(X), ncol(X), byrow = TRUE)

    Zt <- getME(regression, "Zt");

    ## Linear predictor
    x.beta <- as.matrix(tcrossprod(as.matrix(X), sims@fixef))
    z.b <- crossprod(as.matrix(Zt), simulatedRanef)

    mylp0 <- x.beta + z.b
    
    #startTime <- proc.time()
    mylpFe1 <- lapply(as.data.frame(x.beta), function(x)  {matrix(x,nrow=dims[["n"]],dimnames = list(NULL,colnames(X))) })
    
    
    mylp1 <- lapply(as.data.frame(mylp0), function(x)  {matrix(x,nrow=dims[["n"]],dimnames = list(NULL,colnames(X))) })
    

    lpFeList <- eval(regression@resp$nlmod[[2]],merge(do.call("rbind",mylpFe1),as.data.frame(mm)))

    #lpFeList <- lapply(mylpFe1,function(y) { apply(y,1,function(x) {eval(regression@resp$nlmod[[2]],list(c(x,as.data.frame(mm)))[[1]])} ) } )
    
    lpList <- eval(regression@resp$nlmod[[2]],merge(do.call("rbind",mylp1),as.data.frame(mm))) 
    
    #lpList <- lapply(mylp1,function(y) { apply(y,1,function(x) {eval(regression@resp$nlmod[[2]],list(c(x,as.data.frame(mm)))[[1]])} ) } )
    
    
    # lpFe <- do.call("cbind",lpFeList)
    # lp <- do.call("cbind",lpList)
    
    lpFe <- matrix(lpFeList,nrow=dims[["n"]])
    lp <- matrix(lpList,nrow=dims[["n"]])
    #proc.time() - startTime
    
    #TURNOFF
    #lp <- x.beta + z.b + matrix(getME(regression, "offset"), dims[["n"]], numSimulations);

    ## FE
    #TURNOFF #ev <- as.matrix(colMeans(x.beta, 1))
    ev <- as.matrix(colMeans(lpFe, 1))
    ## FE + RE
    lpm <- as.matrix(colMeans(lp, 1))
    pv <- as.matrix(rnorm(n = length(lpm), mean = lpm, sd = sims@sigma), nrow = length(lpm), ncol = 1)
    
    return(list(ev = ev, pv = pv))
  }
)
