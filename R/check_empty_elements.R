
check_empty_elements <- function(...){
  rmlist <- NULL
  iniList <- list(...)

  for(i in 1:length(iniList)){
    if(is.null(iniList[[i]])){
      rmlist <- c(rmlist, i)
    }
  }

  if(!is.null(rmlist)){
    outList <- iniList[-c(rmlist)]
  } else {
    outList <- iniList
  }

  return(outList)
}
