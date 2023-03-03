\name{brainsync}
\alias{brainsync}
\title{Brain wave synchrony clustering}
\usage{
brainsync(
  data_wd = "",

  filepattern = "",
  sep = ",",

  channel.lim = c(),

  freq.lim = c(),
  freq.steps,

  freq.range.min = NA,
  freq.range.max = NA,

  itermax = 10000,
  nstart = 1000,

  alpha_chitest = 0.05,
  pval_adjustment = FALSE,
  method = p.adjust.methods,

  detailed_running_output = FALSE,

  writechisquaretable = TRUE,
  comparisonname = ""
)
}
\description{
This function allows exploring the brain activation patterns over large population samples, estimated through EEG measurements and pre-processed through Welch's analyses, for example.

The data need to be prepared in the format of individual csv files with individual columns for the different EEG channels and different rows for the different frequencies in which the data was assessed.
}
\arguments{
\item{data_wd}{Working directory in which the csv files of interest are located}
\item{filepattern}{Part of the file name to subselect files if needed}
\item{sep}{Separator in the csv files}
\item{channel.lim}{Indicates which EEG channels to include in the analysis}
\item{freq.lim}{Indicates in which intervals the frequencies were originally recorded}
\item{freq.steps}{Indicate how fast the frequency increases (e.g., if the frequency increases by 0.25 steps, then freq.steps = 4)}
\item{freq.range.min}{Indicates the minimal frequency to be included in the analysis}
\item{freq.range.max}{Indicates the maximal frequency to be included in the analysis}
\item{itermax}{Indicates the number of iterations of the k-means clustering algorithm}
\item{nstart}{Indicates the number of starting points in the k-means clustering algorithm}
\item{alpha_chitest}{Determines the alpha value for the chi-test significance testing}
\item{pval_adjustment}{Indicates whether the p-values should be adjusted or not}
\item{method}{Indicates the method for p-value adjustment}
\item{detailed_running_output}{Determines whether the console output should be detailed or not}
\item{writechisquaretable}{Indicates whether the results of the chi square test should be saved as table}
\item{comparisonname}{Indicates the name under which the results will be saved}
}
\examples{}
