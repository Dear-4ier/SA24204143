#' @import Rcpp
#' @useDynLib SA24204143

#' @title ktimes - Integer Multiplication Function
#' @description This function takes two integers in string format as input and returns their product in string format by using Karatsuba algorithm.
#' @param num1 The first integer in string format to be multiplied.
#' @param num2 The second integer in string format to be multiplied.
#' @return The product of the two input integers in string format.
#' @examples
#' \dontrun{
#' result <- ktimes("-987654876543", "12345672345678")
#' print(result)
#' }
#' @export
ktimes <- function(num1, num2){
  sgn <- sign(as.numeric(num1))*sign(as.numeric(num2))
  if (sgn == 0) {
    return("0")
  }
  result <- karatsubaMultiplyRcpp(chr_abs(num1), chr_abs(num2))
  result <- gsub("^0*", "", result)
  if (sgn == -1){
    return(paste("-", result, sep = "", collapse = NULL))
  }
  return(result)
}

#' @title fktimes - Function for Handling Decimal Multiplication
#' @description This function takes two numbers (which could be integers or decimals in string format) as input and returns their product considering decimal points appropriately.
#' @param num1 The first number in string format.
#' @param num2 The second number in string format.
#' @return The product of the two input numbers with decimal points handled correctly in string format.
#' @examples
#' \dontrun{
#' result <- fktimes("-987654.876543", "123.45672345678")
#' print(result)
#' }
#' @export
fktimes <- function(num1, num2) {
  dot_pos1 <- gregexpr("\\.", num1)[[1]]
  if (dot_pos1[1] == -1) {
    len1 <- 0
  } else {
    len1 <- nchar(num1) - dot_pos1[1]
  }
  dot_pos2 <- gregexpr("\\.", num2)[[1]]
  if (dot_pos2[1] == -1) {
    len2 <- 0
  } else {
    len2 <- nchar(num2) - dot_pos2[1]
  }
  len <- len1 + len2
  if (len == 0){
    return(ktimes(num1, num2))
  }
  num1_no_dot <- gsub("\\.", "", num1)
  num2_no_dot <- gsub("\\.", "", num2)
  result <- ktimes(num1_no_dot, num2_no_dot)
  if (len > 0) {
    result <- paste0(substr(result, 1, nchar(result) - len), ".", substr(result, nchar(result) - len + 1, nchar(result)))
  }
  
  result <- gsub("0*$", "", result)
  if (substr(result, nchar(result), nchar(result)) == ".") {
    result <- substr(result, 1, nchar(result) - 1)
  }
  return(result)
}

chr_abs <- function(num) {
  if (substr(num, 1, 1) == "-") {
    return(substr(num, 2, nchar(num)))
  }
  return(num)
}