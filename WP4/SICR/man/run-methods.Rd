\docType{methods}
\name{run,SICR-method}
\alias{run,SICR-method}
\title{Run the model}
\usage{
\S4method{run}{SICR}(model, threads = NULL, seed = NULL)
}
\arguments{
\item{model}{The model to run.}

\item{threads}{Number of threads. Default is NULL, i.e. to use all
available processors.}

\item{seed}{Random number seed. Default is NULL, i.e. the
simulator uses time to seed the random number generator.}
}
\value{
model with result from simulation.
}
\description{
Run the model
}
