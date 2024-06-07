#' generate_type3_data
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
#' Prefix = "type_three"
#' meta_data_file <- system.file(package="SSTD", "extdata", "meta_data.csv")
#' gene_count_file <- system.file(package="SSTD", "extdata", "gene_count.csv")
#' selected_clusters <- c("eL4", "eL5", "eL6", "Endo")
#' generate_type3_data(meta_data_file, gene_count_file, selected_clusters, Prefix)
#' }
#'
generate_type3_data <- function(meta_data_file, gene_count_file, selected_clusters, Prefix) {
  org_data <- org_data_space(meta_data_file, gene_count_file, output_csv_file = paste0(Prefix, "_org_data_space.csv"), selected_clusters)

  sim_spatial_spot_true_prop_data <- spatial_spot_true_prop(selected_clusters, filename = paste0(Prefix, "_spatial_spot_true_prop.csv"))

  spatial_spot_loc(filename = paste0(Prefix, "_spatial_spot_loc.csv"))

  # Generate gradation arrays for each cell type
  # celltype1
  array <- get_gradation_2d(1, 1,  6, 20, FALSE)
  array_1 <- get_gradation_3d(8, 20, c(1, 0, 0), c(0, 1, 0), c(TRUE, TRUE, FALSE))
  array_2 <- get_gradation_2d(0, 0,  26, 20, FALSE)
  array_3 <- rbind(array, t(array_1[, , 1]))
  get_gradation_3d_celltype1 <- rbind(array_3, array_2)

  # celltype2
  array <- get_gradation_2d(0, 0,  6, 20, TRUE)
  array_0 <- get_gradation_2d(1, 1,  2, 20, TRUE)
  array_1 <- get_gradation_3d(8, 20, c(0, 1, 0), c(1, 0, 0), c(TRUE, TRUE, FALSE))
  array_2 <- get_gradation_2d(0, 0,  16, 20, TRUE)
  array_3 <- rbind(array, t(array_1[, , 1]))
  array_4 <- rbind(array_3, array_0)
  array_5 <- get_gradation_3d(8, 20, c(1, 0, 0), c(0, 1, 0), c(TRUE, TRUE, FALSE))
  array_5 <- rbind(array_4, t(array_5[, , 1]))
  get_gradation_3d_celltype2 <- rbind(array_5, array_2)

  # celltype3
  array <- get_gradation_2d(0, 0,  16, 20, TRUE)
  array_0 <- get_gradation_2d(1, 1,  2, 20, TRUE)
  array_1 <- get_gradation_3d(8, 20, c(0, 1, 0), c(1, 0, 0), c(TRUE, TRUE, FALSE))
  array_2 <- get_gradation_2d(0, 0,  6, 20, TRUE)
  array_3 <- rbind(array, t(array_1[, , 1]))
  array_4 <- rbind(array_3, array_0)
  array_5 <- get_gradation_3d(8, 20, c(1, 0, 0), c(0, 1, 0), c(TRUE, TRUE, FALSE))
  array_5 <- rbind(array_4, t(array_5[, , 1]))
  get_gradation_3d_celltype3 <- rbind(array_5, array_2)

  # celltype4
  array <- get_gradation_2d(1, 1,  6, 20, FALSE)
  array_1 <- get_gradation_3d(8, 20, c(0, 1, 0), c(1, 0, 0), c(TRUE, TRUE, FALSE))
  array_2 <- get_gradation_2d(0, 0,  26, 20, FALSE)
  array_3 <- rbind(t(array_1[, , 1]), array)
  get_gradation_3d_celltype4 <- rbind(array_2, array_3)

  result_matrices <- list(result_matrix(get_gradation_3d_celltype1),
                          result_matrix(get_gradation_3d_celltype2),
                          result_matrix(get_gradation_3d_celltype3),
                          result_matrix(get_gradation_3d_celltype4))

  new_data <- spatial_spot_nUMI(sim_spatial_spot_true_prop_data, result_matrices, selected_clusters, org_data, Prefix)

  combined_arrays <- combined_arrays(selected_clusters, get_gradation_3d_celltype1, get_gradation_3d_celltype2, get_gradation_3d_celltype3, get_gradation_3d_celltype4)

  generate_plots_and_save(combined_arrays, selected_clusters, new_data, sim_spatial_spot_true_prop_data, Prefix)

  return(NULL)
  }


