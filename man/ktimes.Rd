% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Karatsuba.R
\name{ktimes}
\alias{ktimes}
\title{ktimes - Integer Multiplication Function}
\usage{
ktimes(num1, num2)
}
\arguments{
\item{num1}{The first integer in string format to be multiplied.}

\item{num2}{The second integer in string format to be multiplied.}
}
\value{
The product of the two input integers in string format.
}
\description{
This function takes two integers in string format as input and returns their product in string format by using Karatsuba algorithm.
}
\examples{
\dontrun{
result <- ktimes("-987654876543", "12345672345678")
print(result)
}
}
