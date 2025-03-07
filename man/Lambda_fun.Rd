% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct.R
\name{Lambda_fun}
\alias{Lambda_fun}
\title{Lambda Function Calculation}
\usage{
Lambda_fun(type, par0, t, u)
}
\arguments{
\item{type}{A character string specifying the type of Lambda function (e.g., "pp").}

\item{par0}{A vector of four parameters: alpha_t, beta_t, alpha_u, beta_u.}

\item{t}{Numeric vector representing time values.}

\item{u}{Numeric vector representing another variable (e.g., units or cycles).}
}
\value{
A list containing the Lambda_t, Lambda_u, and time values.
}
\description{
Calculates different Lambda values for the given parameters and type.
}
