% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{org_data_space}
\alias{org_data_space}
\title{Process and save data}
\usage{
org_data_space(
  meta_data_file,
  gene_count_file,
  output_csv_file,
  selected_clusters
)
}
\arguments{
\item{output_csv_file}{Path to save the processed data as a CSV file.}

\item{selected_clusters}{A character vector specifying the clusters to select from the data.}

\item{input_rds_file}{Path to the input RDS file containing the original data.}

\item{num_rows_per_cluster}{Number of rows to select per cluster.}
}
\value{
None
}
\description{
This function processes the input data from an RDS file, selects specific clusters, and saves the processed data as a CSV file.
}
