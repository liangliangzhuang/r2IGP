% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_extension.R
\name{Log.liklihod}
\alias{Log.liklihod}
\title{Compute Log-Likelihood for Various Models}
\usage{
Log.liklihod(
  model = model,
  est_par = em_para_re[[q]],
  data = coat_dat$diff_Y,
  t_list = t_list,
  u_list = u_list,
  types = types
)
}
\arguments{
\item{model}{A string indicating the model type. It can be one of "M0", "M1", "M2", "M3", or "M4".}

\item{est_par}{A vector of estimated parameters. These include parameters like kappa, sigma2, and others
specific to each model.}

\item{data}{A matrix or data frame of size m x n, where m is the number of observations and n is the number of series.
Each column represents a time series.}

\item{t_list}{A list of time points (or other relevant sequence) for each series, passed to `Lambda_fun`.}

\item{u_list}{A list of values corresponding to some parameter for each series, passed to `Lambda_fun`.}

\item{types}{A vector specifying the types or categories for each series, used in the `Lambda_fun` function.}
}
\value{
The computed log-likelihood for the specified model, which is a single numeric value.
}
\description{
This function calculates the log-likelihood for different models, including "M0", "M1", "M2", "M3", and "M4".
The log-likelihood is computed based on the provided data and estimated parameters. The function supports
different likelihood formulas for each model type, where the models depend on different parameters such as
kappa, sigma2, gamma, and others.
}
