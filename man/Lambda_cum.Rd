% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct.R
\name{Lambda_cum}
\alias{Lambda_cum}
\title{Cumulative Lambda Calculation}
\usage{
Lambda_cum(par, type, t, u)
}
\arguments{
\item{par}{A vector containing parameters.}

\item{type}{A string specifying the model type (e.g., "pp").}

\item{t}{Time vector.}

\item{u}{Another variable vector (optional).}
}
\value{
A list containing cumulative Lambda values for both time and another variable.
}
\description{
Computes the cumulative Lambda values based on input parameters and model type.
}
