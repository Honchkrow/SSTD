#' generate_type4_data
#'
#' This function generates data and plots for selected clusters.
#'
#' @param meta_data_file Path to the metadata file.
#' @param gene_count_file Path to the gene count file.
#' @param selected_clusters A named vector specifying the clusters to include.
#' @param Prefix Prefix for the output CSV files.
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' Prefix = "type_four"
#' selected_clusters <- c("eL4", "eL5", eL6", "Endo")
#' meta_data_file = "inst/extData/meta_data.csv"
#' gene_count_file = "inst/extData/gene_count.csv"
#' source("R/utilities.R")
#' generate_type4_data(meta_data_file, gene_count_file, selected_clusters, Prefix)
#' }
#'
generate_type4_data <- function(meta_data_file, gene_count_file, selected_clusters, Prefix) {
  org_data <- org_data_space(meta_data_file, gene_count_file, output_csv_file = paste0(Prefix, "_org_data_space.csv"), selected_clusters)

  sim_spatial_spot_true_prop_data <- spatial_spot_true_prop(selected_clusters, filename = paste0(Prefix, "_spatial_spot_true_prop.csv"))

  spatial_spot_loc(filename = paste0(Prefix, "_spatial_spot_loc.csv"))

  # Generate gradation arrays for each cell type
  # celltype1
  array <- get_gradation_2d(1, 1,  14, 20, FALSE)
  array_1 <- get_gradation_3d(12, 20, c(1, 0, 0), c(0, 1, 0), c(TRUE, TRUE, FALSE))
  array_2 <- get_gradation_2d(0, 0,  14, 20, FALSE)
  array_3 <- rbind(array, t(array_1[, , 1]))
  get_gradation_3d_celltype1 <- rbind(array_3, array_2)

  # celltype2
  array <- get_gradation_2d(0, 0,  16, 20, TRUE)
  array_0 <- get_gradation_2d(1, 1,  14, 4, TRUE)
  array_1 <- get_gradation_3d(10, 4, c(0.2, 1, 0), c(1, 0.2, 0), c(TRUE, TRUE, FALSE))
  array_2 <- get_gradation_2d(0, 0,  24, 6, TRUE)
  array_3 <- get_gradation_2d(1/3, 1/3, 10, 10, TRUE)
  array_6 <- get_gradation_3d(10, 14, c(1, 0.2, 0), c(0.2, 1, 0), c(TRUE, TRUE, FALSE))
  array_4 <- rbind(t(rbind(t(array_1[, , 1]), array_0)), t(rbind(array_3, array_6[, , 1])))
  array_5 <- rbind(array_4, t(array_2))
  get_gradation_3d_celltype2 <- rbind(array, t(array_5))

  # celltype3
  array <- get_gradation_2d(0, 0,  16, 20, TRUE)
  array_0 <- get_gradation_2d(1, 1,  2, 6, TRUE)
  array_1 <- get_gradation_3d(8, 6, c(0.2, 1, 0), c(1, 0.2, 0), c(TRUE, TRUE, FALSE))
  array_2 <- get_gradation_2d(0, 0,  6, 20, TRUE)
  array_3 <- get_gradation_3d(8, 6, c(1, 0.2, 0), c(0.2, 1, 0), c(TRUE, TRUE, FALSE))
  array_4 <- get_gradation_2d(0, 0,  18, 6, TRUE)
  array_5 <- get_gradation_3d(8, 6, c(1, 0.2, 0), c(0.2, 1, 0), c(TRUE, TRUE, FALSE))
  array_6 <- get_gradation_2d(0.3, 0.3,  8, 6, TRUE)
  array_9 <- get_gradation_3d(8, 6, c(0.2, 1, 0), c(1, 0.2, 0), c(TRUE, TRUE, FALSE))

  array_8 <- rbind(t(array_6), rbind(array_9[, , 1], t(array_6)))
  array_10 <- rbind(t(array_1[, , 1]), rbind(array_0, t(array_5[, , 1])))
  array_7 <- rbind(t(array_4), rbind(t(array_8), t(array_10)))
  get_gradation_3d_celltype3 <- rbind(rbind(array, t(array_7)), array_2)

  # celltype4
  array <- get_gradation_2d(1, 1,  6, 6, FALSE)
  array_1 <- get_gradation_3d(8, 6, c(0.2, 1, 0), c(1, 0.2, 0), c(TRUE, TRUE, FALSE))
  array_2 <- get_gradation_2d(0, 0,  26, 20, FALSE)
  array_3 <- get_gradation_2d(0, 0,  14, 6, FALSE)
  array_4 <- get_gradation_3d(8, 8, c(0.2, 1, 0), c(1, 0.2, 0), c(TRUE, TRUE, FALSE))
  array_5 <- get_gradation_2d(0.3, 0.3,  8, 6, TRUE)

  array_6 <- rbind(t(array_1[, , 1]), array)
  array_7 <- rbind(t(array_5), array_4[, , 1])
  array_8 <- rbind(t(array_7), t(array_6))
  get_gradation_3d_celltype4 <- rbind(array_2, t(rbind(t(array_3), array_8)))

  result_matrices <- list(result_matrix(get_gradation_3d_celltype1),
                          result_matrix(get_gradation_3d_celltype2),
                          result_matrix(get_gradation_3d_celltype3),
                          result_matrix(get_gradation_3d_celltype4))

  new_data <- spatial_spot_nUMI(sim_spatial_spot_true_prop_data, result_matrices, selected_clusters, org_data, Prefix)

  combined_arrays <- combined_arrays(selected_clusters, get_gradation_3d_celltype1, get_gradation_3d_celltype2, get_gradation_3d_celltype3, get_gradation_3d_celltype4)

  generate_plots_and_save(combined_arrays, selected_clusters, new_data, sim_spatial_spot_true_prop_data, Prefix)

  return(NULL)
}


