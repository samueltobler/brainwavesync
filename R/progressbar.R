progressbar <- function(i, n, input = "", lim.output = TRUE) {

  if(lim.output == TRUE){
  if(round(i / n * 100) != round((i+1) / n * 100)) {

  print(paste0(input, ": ", round(i / n * 100), '% completed'))
  }
    if (i == n) {print(paste0(input,': Done'))}
}
}

