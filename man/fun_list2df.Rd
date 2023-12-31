% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun.R
\name{fun_list2df}
\alias{fun_list2df}
\title{change a list into dataframe format}
\usage{
fun_list2df(input_list, multiple_items_solver = "keep", colname = "pathway")
}
\arguments{
\item{input_list}{the list, each element is is_atomic format}

\item{multiple_items_solver}{how to deal with items appear multiple time, accept drop and keep}

\item{colname}{what the colname to show, default is pathway}
}
\value{
df
}
\description{
change a list into dataframe format
}
\examples{
fun_list2df(sig_vec_list)
}
