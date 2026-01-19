calculate_surp <- function(x) {
    x * log2(1 / (x + 1e-6))
}

.stopIfNotAll <- function(exprs, errorMsg) {
    for (expr in exprs) {
        if (!expr) {
            stop(errorMsg, call. = FALSE)
        }
    }
}

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
