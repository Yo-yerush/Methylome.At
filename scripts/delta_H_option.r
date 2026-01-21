### op1
# calculate_surp <- function(p) {
#     p * log2(1 / (p + 1e-6))
# }

### op2
calculate_surp <- function(p) {
    -(p * log2(p) + (1 - p) * log2(1 - p))
}

# calculate_surp <- function(x, m_ctrl, n_ctrl, m_trnt, n_trnt) {
#     mutate(x,
#         pi0 = pmin(pmax(m_ctrl / n_ctrl, 1e-6), 1 - 1e-6), # ctrl expectation π0
#         surprisal = -(lchoose(n_trnt, m_trnt) +
#             m_trnt * log(pi0) +
#             (n_trnt - m_trnt) * log(1 - pi0))
#     )
# }

#################################################################

.stopIfNotAll <- function(exprs, errorMsg) {
    for (expr in exprs) {
        if (!expr) {
            stop(errorMsg, call. = FALSE)
        }
    }
}

#################################################################

.movingSum <- function(minPos, maxPos, pos, val, weights = 1, windowSize = 150, normalize = FALSE) {
    .stopIfNotAll(c(length(pos) == length(val)), "pos and val vectors need to have the same length")

    if (length(weights) < length(val)) {
        weights <- rep(weights, length.out = length(val))
    } else if (length(weights) > length(val)) {
        weights <- weights[1:length(val)]
    }

    # Filter out NAs
    keepIndexes <- which(!is.na(pos) & !is.na(val) & !is.na(weights))
    pos <- pos[keepIndexes]
    weights <- weights[keepIndexes]
    val <- val[keepIndexes]


    # set the values
    rawVector <- rep(0, (maxPos - minPos + 1))
    rawVector[pos - minPos + 1] <- weights * val
    rawVector <- c(rawVector, rep(0, (windowSize - 1)))

    # Define the (triangular) kernel.
    kernel <- c(rep(1, (windowSize)))


    smoothedVector <- RcppRoll::roll_sum(rawVector, length(kernel), weights = kernel, normalize = normalize)

    # smoothedVector <- smoothedVector[seq(1,length(smoothedVector), by=windowSize)]

    return(smoothedVector)
}

#################################################################

.analyseReadsInsideBinsReplicates_yo_dH <- function(methylationData, bins, currentRegion,
                                                    condition, pseudocountM, pseudocountN) {
    binSize <- min(unique(width(bins)))
    # Rcpp
    m <- grep("readsM", names(mcols(methylationData)))
    n <- grep("readsN", names(mcols(methylationData)))
    readsM <- matrix(0, ncol = length(m), nrow = length(bins))
    for (i in 1:length(m)) {
        test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[m[i]]], windowSize = binSize)
        readsM[, i] <- test[seq(1, length(test) - binSize, by = binSize)]
    }

    readsN <- matrix(0, ncol = length(n), nrow = length(bins))
    for (i in 1:length(n)) {
        test2 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[n[i]]], windowSize = binSize)
        readsN[, i] <- test2[seq(1, length(test) - binSize, by = binSize)]
    }

    # proportions <- (readsM + pseudocountM)/ (readsN + pseudocountN)
    proportions <- calculate_surp((readsM + pseudocountM) / (readsN + pseudocountN))

    proportions <- as.data.frame(proportions)
    names_prop <- paste0("proportionsR", 1:ncol(proportions))
    colnames(proportions) <- names_prop

    m1 <- m[which(condition == unique(condition)[1])]
    n1 <- n[which(condition == unique(condition)[1])]
    m2 <- m[which(condition == unique(condition)[2])]
    n2 <- n[which(condition == unique(condition)[2])]

    readsM1 <- readsM[, which(condition == unique(condition)[1])]
    # readsM1 <- matrix(0, ncol = length(m1), nrow=length(bins))
    # for(i in 1:length(m1)){
    #   test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[m1[i]]], windowSize = binSize)
    #   readsM1[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
    # }
    sumReadsM1 <- apply(readsM1, 1, sum)



    readsN1 <- readsN[, which(condition == unique(condition)[1])]
    # readsN1 <- matrix(0, ncol = length(n1), nrow=length(bins))
    # for(i in 1:length(n1)){
    #   test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[n1[i]]], windowSize = binSize)
    #   readsN1[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
    # }
    sumReadsN1 <- apply(readsN1, 1, sum)
    proportion1 <- calculate_surp((sumReadsM1 + pseudocountM) / (sumReadsN1 + pseudocountN))
    # proportion1 <- (sumReadsM1 + pseudocountM)/ (sumReadsN1 + pseudocountN)


    readsM2 <- readsM[, which(condition == unique(condition)[2])]
    # readsM2 <- matrix(0, ncol = length(m2), nrow=length(bins))
    # for(i in 1:length(m2)){
    #   test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[m2[i]]], windowSize = binSize)
    #   readsM2[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
    # }
    sumReadsM2 <- apply(readsM2, 1, sum)

    readsN2 <- readsN[, which(condition == unique(condition)[2])]
    # readsN2 <- matrix(0, ncol = length(n2), nrow=length(bins))
    # for(i in 1:length(n2)){
    #   test <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), mcols(methylationData)[[n2[i]]], windowSize = binSize)
    #   readsN2[,i] <- test[seq(1,length(test)-binSize, by=binSize)]
    # }
    sumReadsN2 <- apply(readsN2, 1, sum)

    #
    proportion2 <- calculate_surp((sumReadsM2 + pseudocountM) / (sumReadsN2 + pseudocountN))
    # proportion2 <- (sumReadsM2 + pseudocountM)/ (sumReadsN2 + pseudocountN)

    cytosines <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), rep(1, length(start(methylationData))), windowSize = binSize)
    cytosinesCount <- cytosines[seq(1, length(cytosines) - binSize, by = binSize)]


    bins$sumReadsM1 <- sumReadsM1
    bins$sumReadsN1 <- sumReadsN1
    bins$proportion1 <- proportion1
    bins$sumReadsM2 <- sumReadsM2
    bins$sumReadsN2 <- sumReadsN2
    bins$proportion2 <- proportion2
    bins$cytosinesCount <- cytosinesCount
    mcols(bins) <- cbind(mcols(bins), proportions)
    return(bins)
}

#################################################################

.analyseReadsInsideBins_yo_dH <- function(methylationData, bins, currentRegion) {
    binSize <- min(unique(width(bins)))
    # Rcpp
    readsM1 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsM1, windowSize = binSize)
    sumReadsM1 <- readsM1[seq(1, length(readsM1) - binSize, by = binSize)]

    readsN1 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsN1, windowSize = binSize)
    sumReadsN1 <- readsN1[seq(1, length(readsN1) - binSize, by = binSize)]

    # proportion1 <- sumReadsM1/sumReadsN1
    proportion1 <- calculate_surp(sumReadsM1 / sumReadsN1)

    readsM2 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsM2, windowSize = binSize)
    sumReadsM2 <- readsM2[seq(1, length(readsM2) - binSize, by = binSize)]

    readsN2 <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsN2, windowSize = binSize)
    sumReadsN2 <- readsN2[seq(1, length(readsN2) - binSize, by = binSize)]

    # proportion2 <- sumReadsM2/sumReadsN2
    proportion2 <- calculate_surp(sumReadsM2 / sumReadsN2)


    cytosines <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), rep(1, length(start(methylationData))), windowSize = binSize)
    cytosinesCount <- cytosines[seq(1, length(cytosines) - binSize, by = binSize)]

    bins$sumReadsM1 <- sumReadsM1
    bins$sumReadsN1 <- sumReadsN1
    bins$proportion1 <- proportion1
    bins$sumReadsM2 <- sumReadsM2
    bins$sumReadsN2 <- sumReadsN2
    bins$proportion2 <- proportion2
    bins$cytosinesCount <- cytosinesCount

    return(bins)
}

#################################################################

.analyseReadsInsideBinsOneSample_yo_dH <- function(methylationData, bins, currentRegion) {
    binSize <- min(unique(width(bins)))
    # Rcpp
    readsM <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsM, windowSize = binSize)
    sumReadsM <- readsM[seq(1, length(readsM) - binSize, by = binSize)]

    readsN <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsN, windowSize = binSize)
    sumReadsN <- readsN[seq(1, length(readsN) - binSize, by = binSize)]

    # Proportion <- sumReadsM/sumReadsN
    Proportion <- calculate_surp(sumReadsM / sumReadsN)


    bins$sumReadsM <- sumReadsM
    bins$sumReadsN <- sumReadsN
    bins$Proportion <- Proportion

    return(bins)
}

#################################################################

assignInNamespace(
    x = ".analyseReadsInsideBinsReplicates",
    value = .analyseReadsInsideBinsReplicates_yo_dH,
    ns = "DMRcaller"
)

assignInNamespace(
    x = ".analyseReadsInsideBins",
    value = .analyseReadsInsideBins_yo_dH,
    ns = "DMRcaller"
)

assignInNamespace(
    x = ".analyseReadsInsideBinsOneSample",
    value = .analyseReadsInsideBinsOneSample_yo_dH,
    ns = "DMRcaller"
)
