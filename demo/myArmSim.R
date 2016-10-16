#I'm changing the original arm::sim function to allow NLMM posterior sampling
#Not sure why they didn't allow this originally. Seems to be a natural extension of LMM
#Changes made to original function are indicated by #TURNOFF
myArmSim <- 
  
  function (object, ...) 
  {
    .local <- function (object, n.sims = 100) 
    {
      applyLeftFactor <- function(decomp, rhs) {
        c(as.vector(decomp$ul %*% rhs[ranefRange] + decomp$ur %*% 
                      rhs[fixefRange]), as.vector(decomp$lr %*% rhs[fixefRange]))
      }
      getInverseInformationLeftFactor <- function(regression) {
        Lz <- getME(regression, "L")
        Rzx <- getME(regression, "RZX")
        Rx <- getME(regression, "RX")
        solveFunc <- getMethod("solve", signature(a = "CHMfactor", 
                                                  b = "diagonalMatrix"))
        Rz.inv <- t(solveFunc(Lz, Diagonal(Lz@Dim[1]), "L"))
        Rx.inv <- solve(Rx)
        Rzx.inv <- -Rz.inv %*% Rzx %*% Rx.inv
        Lambda <- t(getME(regression, "Lambda"))
        return(list(ul = Lambda %*% Rz.inv, ur = Lambda %*% 
                      Rzx.inv, lr = Rx.inv))
      }
      sampleCommonScale <- function(ignored) {
        return(sqrt(1/rgamma(1, 0.5 * numDoF, 0.5 * devcomp$cmp[["pwrss"]])))
      }
      regression <- object
      devcomp <- getME(regression, "devcomp")
      dims <- devcomp$dims
      #TURNOFF
      #if (dims[["NLMM"]] != 0L) 
        #stop("sim not yet implemented for nlmms")
      numObs <- dims[["n"]]
      numRanef <- dims[["q"]]
      numFixef <- dims[["p"]]
      numLevels <- dims[["reTrms"]]
      isLinearMixedModel <- dims[["GLMM"]] == 0L #TURNOFF #&& dims[["NLMM"]] == #0L
      numEffects <- numRanef + numFixef
      numDoF <- numObs - numFixef
      ranefRange <- 1:numRanef
      fixefRange <- numRanef + 1:numFixef
      groupsPerUniqueFactor <- lapply(regression@flist, levels)
      factorPerLevel <- attr(regression@flist, "assign")
      coefficientNamesPerLevel <- regression@cnms
      numCoefficientsPerLevel <- as.numeric(sapply(coefficientNamesPerLevel, 
                                                   length))
      numGroupsPerLevel <- as.numeric(sapply(groupsPerUniqueFactor[factorPerLevel], 
                                             length))
      numRanefsPerLevel <- numCoefficientsPerLevel * numGroupsPerLevel
      ranefLevelMap <- rep.int(seq_along(numRanefsPerLevel), 
                               numRanefsPerLevel)
      simulatedSD <- if (isLinearMixedModel) {
        rep(NA, n.sims)
      }
      else {
        NA
      }
      simulatedRanef <- vector("list", numLevels)
      names(simulatedRanef) <- names(regression@cnms)
      for (i in 1:numLevels) {
        simulatedRanef[[i]] <- array(NA, c(n.sims, numGroupsPerLevel[i], 
                                           numCoefficientsPerLevel[i]), list(NULL, groupsPerUniqueFactor[[factorPerLevel[i]]], 
                                                                             coefficientNamesPerLevel[[i]]))
      }
      simulatedFixef <- matrix(NA, n.sims, numFixef, dimnames = list(NULL, 
                                                                     names(fixef(regression))))
      effectsMean <- c(getME(regression, "b")@x, getME(regression, 
                                                       "beta"))
      effectsCovLeftFactor <- getInverseInformationLeftFactor(regression)
      for (i in 1:n.sims) {
        if (isLinearMixedModel) {
          simulatedSD[i] <- sampleCommonScale(regression)
          sphericalEffects <- rnorm(numEffects, 0, simulatedSD[i])
        }
        else {
          sphericalEffects <- rnorm(numEffects)
        }
        simulatedEffects <- applyLeftFactor(effectsCovLeftFactor, 
                                            sphericalEffects) + effectsMean
        simulatedFixef[i, ] <- simulatedEffects[fixefRange]
        rawRanef <- simulatedEffects[ranefRange]
        simulatedRanefPerLevel <- split(rawRanef, ranefLevelMap)
        for (k in 1:numLevels) {
          simulatedRanef[[k]][i, , ] <- matrix(simulatedRanefPerLevel[[k]], 
                                               ncol = numCoefficientsPerLevel[k], byrow = TRUE)
        }
      }
      ans <- new("sim.merMod", fixef = simulatedFixef, ranef = simulatedRanef, 
                 sigma = simulatedSD)
      return(ans)
    }
    .local(object, ...)
  }