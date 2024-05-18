#' generate_type1_data
#'
#' This function generates data and plots for selected clusters.
#'
#' @param meta_data_file Path to the metadata file.
#' @param gene_count_file Path to the gene count file.
#' @param selected_clusters A named vector specifying the clusters to include.
#' @param type Prefix for the output CSV files.
#' @return NULL
#' @export
#'
#' @examples
#' \dontrun{
#' type = "type_one"
#' selected_clusters <- c(celltype1="eL4", celltype2="eL5", celltype3="eL6", celltype4="Endo")
#' generate_type1_data("meta_data.csv", "gene_count.csv", selected_clusters, "type_one")
#' }
#'
generate_type1_data <- function(meta_data_file, gene_count_file, selected_clusters, type) {
  org_data <- org_data_space("meta_data.csv", "gene_count.csv", output_csv_file = paste0(type, "_org_data_space.csv"), selected_clusters)

  sim_spatial_spot_true_prop_data <- spatial_spot_true_prop(selected_clusters, filename = paste0(type, "_spatial_spot_true_prop.csv"))

  spatial_spot_loc(filename = paste0(type, "_spatial_spot_loc.csv"))

  # Generate gradation arrays for each cell type
  # celltype1
  array <- get_gradation_2d(1, 1,  10, 20, TRUE)
  array_1 <- get_gradation_2d(0, 0,  30, 20, TRUE)
  get_gradation_2d_celltype1 <- rbind(array, array_1)

  # celltype2
  array <- get_gradation_2d(0, 0,  10, 20, TRUE)
  array_0 <- get_gradation_2d(1, 1,  10, 20, TRUE)
  array_1 <- get_gradation_2d(0, 0,  20, 20, TRUE)
  get_gradation_2d_celltype2 <- rbind(array, array_0)
  get_gradation_2d_celltype2 <- rbind(get_gradation_2d_celltype2, array_1)

  # celltype3
  array <- get_gradation_2d(0, 0,  10, 20, TRUE)
  array_0 <- get_gradation_2d(1, 1,  10, 20, TRUE)
  array_1 <- get_gradation_2d(0, 0,  20, 20, TRUE)
  get_gradation_2d_celltype3 <- rbind(array_1, array_0)
  get_gradation_2d_celltype3 <- rbind(get_gradation_2d_celltype3, array)

  # celltype4
  array <- get_gradation_2d(1, 1,  10, 20, TRUE)
  array_1 <- get_gradation_2d(0, 0,  30, 20, TRUE)
  get_gradation_2d_celltype4 <- rbind(array_1, array)

  result_matrices <- list(result_matrix(get_gradation_2d_celltype1),
                          result_matrix(get_gradation_2d_celltype2),
                          result_matrix(get_gradation_2d_celltype3),
                          result_matrix(get_gradation_2d_celltype4))

  new_data <- spatial_spot_nUMI(sim_spatial_spot_true_prop_data, result_matrices, selected_clusters, org_data, type)

  combined_arrays <- list(selected_clusters, get_gradation_2d_celltype1,
                          get_gradation_2d_celltype2,
                          get_gradation_2d_celltype3,
                          get_gradation_2d_celltype4)

  generate_plots_and_save(combined_arrays, selected_clusters, new_data, sim_spatial_spot_true_prop_data, type)

  return(NULL)
}


