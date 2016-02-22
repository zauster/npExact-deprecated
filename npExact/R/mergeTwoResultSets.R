
mergeTwoResultSets <- function(results, resultsGreater, resultsLess,
                               merge.d.alt = FALSE) {
        results[["probrej"]] <- resultsGreater[["probrej"]] + resultsLess[["probrej"]]
        results[["probrej"]] <- ifelse(results[["probrej"]] > 1,
                                       1, results[["probrej"]])

        results[["iterations.taken"]] <- max(resultsGreater[["iterations.taken"]],
                                             resultsLess[["iterations.taken"]])

        if(merge.d.alt == TRUE) {
            results[["theta"]] <- c(resultsGreater[["theta"]],
                                    resultsLess[["theta"]])
            results[["typeIIerror"]] <- c(resultsGreater[["typeIIerror"]],
                                          resultsLess[["typeIIerror"]])
            results[["d.alternative"]] <- c(resultsGreater[["d.alternative"]],
                                            1 - resultsLess[["d.alternative"]])
        }

        return(results)
}
