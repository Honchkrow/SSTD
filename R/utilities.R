#' Process and save data
#'
#' This function processes the input data from an RDS file, selects specific clusters, and saves the processed data as a CSV file.
#'
#' @param input_rds_file Path to the input RDS file containing the original data.
#' @param output_csv_file Path to save the processed data as a CSV file.
#' @param selected_clusters A character vector specifying the clusters to select from the data.
#' @param num_rows_per_cluster Number of rows to select per cluster.
#' @return None
#' @export
#'
#' @importFrom utils read.csv
#'
org_data_space <- function(meta_data_file, gene_count_file, output_csv_file, selected_clusters) {
    # Read the original CSV data
    meta_data <- read.csv(meta_data_file)
    gene_count <- read.csv(gene_count_file)

    # Remove the first column from gene_count (assuming it's the row names)
    gene_count <- gene_count[, -1]

    # Initialize empty data frames to store filtered meta and count data
    filtered_meta <- data.frame()
    filtered_count <- data.frame()

    # Iterate over selected clusters
    for (cluster in selected_clusters) {
        selected_rows <- which(meta_data$cluster_name == cluster)

        selected_rows <- selected_rows[1:100]

        filtered_meta <- rbind(filtered_meta, meta_data[selected_rows, ])
        filtered_count <- rbind(filtered_count, gene_count[selected_rows, ])
    }

    # Rename gene count column names to exclude the first character 'X' for the first 12 genes
    for (i in 1:12) {
        colnames(filtered_count)[i] <- substr(colnames(filtered_count)[i], 2, nchar(colnames(filtered_count)[i]))
    }

    # Combine meta data and gene count data
    org_data_space <- cbind(filtered_meta, filtered_count)

    # Write processed data to CSV file
    write.csv(org_data_space, file = output_csv_file, row.names = FALSE)

    return(org_data_space)
}


#' Generate spatial spot true proportions data
#'
#' This function generates synthetic spatial spot true proportions data, where each spot
#' is represented by a square grid cell. It creates a data frame containing the proportions
#' of different cell types within each spot, along with the coordinates of the spot locations.
#'
#' @param rows The number of rows in the grid (default is 40)
#' @param cols The number of columns in the grid (default is 20)
#' @param cell_size The size of each grid cell in pixels (default is 500)
#' @param cell_types A character vector specifying the cell types (default is c("eL4", "eL5", "eL6", "Endo"))
#' @param filename The filename to save the generated data to (default is "spatial_spot_true_prop.csv")
#' @return A data frame containing the generated spatial spot true proportions data
#' @export
spatial_spot_true_prop <- function(cell_types, filename) {
    sim_spatial_spot_true_prop_data <- data.frame(matrix(0, ncol = length(cell_types), nrow = 40 * 20))
    colnames(sim_spatial_spot_true_prop_data) <- c(cell_types)

    k = 0
    for (j in 1:40) {
        for (i in 1:20) {
            col <- i * 500
            row <- j * 500
            k <- k + 1
            rownames(sim_spatial_spot_true_prop_data)[k] <- paste0("x", col, "y", row)
            sim_spatial_spot_true_prop_data$x[k] <- col
            sim_spatial_spot_true_prop_data$y[k] <- row
            for (cell_type in cell_types) {
                sim_spatial_spot_true_prop_data[k, cell_type] <- 0
            }
        }
    }

    write.csv(sim_spatial_spot_true_prop_data[, 1:4], filename, row.names = TRUE)
    return(sim_spatial_spot_true_prop_data)
}


#' Generate spatial spot location data and save it to a CSV file
#'
#' This function generates spatial spot location data for a 20x40 grid of spots,
#' where each spot is represented by a unique combination of x and y coordinates.
#' The generated data is saved to a CSV file with the specified filename.
#'
#' @param filename The filename to save the generated data to (default: "spatial_spot_loc.csv")
#' @return A data frame containing the generated spatial spot location data
#' @export
#'
#' @importFrom utils read.csv
#'
spatial_spot_loc <- function(filename) {
    # Generate grid data frame to store x, y coordinates, and labels
    grid <- data.frame(x = numeric(0),
                       y = numeric(0),
                       label = numeric(0))

    # Generate data frame to store imagecol, imagerow, x, and y coordinates
    data <- data.frame(imagecol = numeric(800),
                       imagerow = numeric(800),
                       x = numeric(800),
                       y = numeric(800))

    # Initialize counter
    k <- 0

    # Loop to generate data
    for(j in 1:40) {
        for(i in 1:20) {
            col <- i * 500
            row <- j * 500
            k <- k + 1
            rownames(data)[k] <- paste0("x", col, "y", row)
            data$imagecol[k] <- col
            data$imagerow[k] <- row
            data$x[k] <- (i - 1) %% 20 + 1
            data$y[k] <- j
        }
    }

    # Write data to CSV file
    write.csv(data, filename, row.names = TRUE)

    # Read data from CSV file
    spatial_spot_loc_data <- read.csv(filename, header = TRUE)

    # Generate grid data
    for(j in 1:40){
        for(i in 1:20){
            x <- (i-1)*500 + 500
            y <- (j-1)*500 + 500
            grid <- rbind(grid, data.frame(x, y, label=i+(j-1)*10))
        }
    }

    # Return generated data frame
    return(data)
}


#' Calculate total nUMI for all cell types
#'
#' This function calculates the total nUMI for all specified cell types and writes the result to CSV files.
#'
#' @param result_matrices List of matrices representing the cell counts for each cell type
#' @param cluster_names List of names of the cell types
#' @param org_data Original data frame containing cluster information
#' @param spot_nUMI_filename Name of the CSV file to write the spatial spot nUMI counts
#' @param new_data_filename Name of the CSV file to write the new data
#'
#' @return A data frame containing the spatial spot nUMI counts for all specified cell types
spatial_spot_nUMI <- function(sim_spatial_spot_true_prop_data, result_matrices, cluster_names, org_data, prefix = "") {
    num_cell_types <- length(result_matrices)
    new_data <- data.frame(matrix(ncol = ncol(org_data), nrow = 0))
    new_spot_data <- data.frame(matrix(0, ncol = ncol(org_data) - 9, nrow = 0))

    for (cell_type_index in 1:num_cell_types) {
        result_matrix <- result_matrices[[cell_type_index]]
        cluster_name <- cluster_names[[cell_type_index]]

        # Temporary new spot data for the current cell type
        temp_new_spot_data <- data.frame(matrix(0, ncol = ncol(org_data) - 9, nrow = 0))

        for(i in 1:40){
            for(j in 1:20){
                n <- result_matrix[i,j]
                vector1 <- rep(0, 1020)
                for(k in 1:n){
                    if(n == 0) {
                        next
                    }
                    row <- sample(which(org_data$cluster_name == cluster_name),1)
                    new_row <- org_data[row,]
                    x <- round(sample(1:500,1) + 500*(j-1), digits = 10)
                    y <- round(sample(1:500,1) + 500*(i-1), digits = 10)
                    new_row$X <- x
                    new_row$Y <- y
                    new_data <- rbind(new_data, new_row)
                    vector1 <- vector1 + new_row[10:ncol(new_data)]
                }
                colnames(temp_new_spot_data) <- colnames(org_data)[10:ncol(org_data)]
                temp_new_spot_data <- rbind(temp_new_spot_data, vector1)
            }
        }
        # Add current cell type's new spot data to the existing new_spot_data
        if (nrow(new_spot_data) == 0) {
            new_spot_data <- temp_new_spot_data
        } else {
            new_spot_data <- new_spot_data + temp_new_spot_data
        }
    }

    # Calculate total nUMI for all cell types
    for (i in 1:12) {
        colnames(new_spot_data)[i] <- substr(colnames(new_spot_data)[i], 2, nchar(colnames(new_spot_data)[i]))
    }

    rownames(new_spot_data) <- rownames(sim_spatial_spot_true_prop_data)

    # Write new data to CSV files
    spot_nUMI_filename <- paste0(prefix, "_spatial_spot_nUMI.csv")
    new_data_filename <- paste0(prefix, "_new_data.csv")
    write.csv(new_spot_data, spot_nUMI_filename)
    write.csv(new_data, new_data_filename)

    return(new_data)
}


#' Calculate cell counts within each spot based on combined array
#'
#' This function takes a combined array representing the spatial distribution of cells
#' and calculates the cell counts within each spot based on randomly sampled values.
#' The input combined array should have dimensions of 40x20, representing a grid of 40x20 spots.
#' Each element of the combined array indicates the presence (1) or absence (0) of cells in the corresponding spot.
#' The function randomly samples a cell count between 5 and 20 for each spot and multiplies it with the combined array
#' to obtain the result matrix representing cell counts for each spot.
#'
#' @param combined_array A matrix representing the spatial distribution of cells, with dimensions 40x20
#' @return A matrix containing cell counts within each spot, with dimensions 40x20
#' @export
result_matrix <- function(combined_array) {
    # Initialize a matrix to store the result
    matrix_data <- matrix(NA, nrow = 40, ncol = 20)

    # Loop through each element of the combined array
    for (i in 1:800) {
        # Randomly sample a cell count between 5 and 20 for each spot
        sample_size <- sample(5:20, 1)
        # Assign the sampled cell count to the corresponding position in the matrix
        matrix_data[i] <- sample_size
    }

    # Reshape the vector into a matrix with dimensions 40x20
    matrix_data <- matrix(matrix_data, nrow = 40, ncol = 20, byrow = TRUE)

    # Multiply the sampled cell counts with the combined array to get the result matrix
    result_matrix <- matrix_data * combined_array

    return(result_matrix)
}


#' Generate plots and save
#'
#' This function generates plots of cell locations with a grid and a pie chart of cell type proportions.
#' The plots are saved as PDF files.
#'
#' @param combined_arrays A list containing matrices of data for each cell type
#' @param include_celltypes A character vector indicating which cell types to include
#' @return None
#' @export
#'
#' @importFrom ggplot2 scale_color_discrete scale_fill_discrete
#' @importFrom scatterpie geom_scatterpie
#'
#' @examples
#' \dontrun{
#' generate_plots_and_save(combined_arrays, selected_clusters,
#' new_data, sim_spatial_spot_true_prop_data, type)
#' }
#'
generate_plots_and_save <- function(combined_arrays, include_celltypes, new_data, sim_spatial_spot_true_prop_data, prefix = "") {
    # Filter cell types based on parameters
    cell_types <- include_celltypes

    # Plot locations of cells with grid
    plot_Locations_of_cell_with_grid <- ggplot(new_data[new_data$cluster_name %in% cell_types, 1:9], aes(x=X, y=Y, color=cluster_name)) +
        geom_point() +
        scale_color_discrete() +  # 使用 scale_color_discrete 自动生成颜色映射
        scale_x_continuous(breaks = seq(0, 10000, 1000)) +
        scale_y_continuous(breaks = seq(0, 20000, 1000)) +
        theme_gray() +
        theme(legend.title=element_blank(), legend.text=element_text(size=12),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid'),
              panel.grid.minor = element_line(size = 0.5, linetype = 'solid'))

    # 将 combined_arrays 的所有元素转换为矩阵
    for (name in names(combined_arrays)) {
        combined_arrays[[name]] <- as.matrix(combined_arrays[[name]])
    }

    # 更新 sim_spatial_spot_true_prop_data 中相应的列
    for (ctype in cell_types) {
        if (ctype %in% names(combined_arrays)) {  # 确保 ctype 存在于 combined_arrays 中
            sim_spatial_spot_true_prop_data[, ctype] <- c(t(combined_arrays[[ctype]]))
        } else {
            warning(paste("Cell type not found in combined_arrays:", ctype))
        }
    }

    write.csv(sim_spatial_spot_true_prop_data[, cell_types], paste0(prefix, "_spatial_spot_true_prop.csv"))

    # Pie chart of cell type proportions
    Pie_chart_celltype_proportions <- ggplot() +
        geom_scatterpie(aes(x=x, y=y), data=sim_spatial_spot_true_prop_data, cols=cell_types) +
        scale_fill_discrete() +  # 使用 scale_fill_discrete 自动生成颜色映射
        theme_classic() +
        theme(legend.title=element_blank(), legend.text=element_text(size=12))

    # Save plot objects as PDF files
    ggsave(paste0(prefix, "_plot_Locations_of_cell_with_grid.pdf"), plot = plot_Locations_of_cell_with_grid, device = "pdf", width = 5, height = 8)
    ggsave(paste0(prefix, "_Pie_chart_celltype_proportions.pdf"), plot = Pie_chart_celltype_proportions, device = "pdf", width = 5, height = 8)
}




#' Generate a combined 2D gradation array matrix
#'
#' This function generates a combined array matrix based on the specified parameters.
#'
#' @param start The starting value for the sequence
#' @param stop The ending value for the sequence
#' @param width The number of columns in the matrix
#' @param height The number of rows in the matrix
#' @param is_horizontal Logical value indicating whether the matrix should be horizontal (TRUE) or vertical (FALSE)
#' @return A matrix representing the combined array
#' @export
get_gradation_2d <- function(start, stop, width, height, is_horizontal) {
    if (is_horizontal) {
        return(matrix(rep(seq(start, stop, length.out = width), each = height), nrow = width))
    } else {
        return(t(matrix(rep(seq(start, stop, length.out = height), each = width), nrow = height)))
    }
}


#' Generate a 3D Gradient Matrix
#'
#' Generate an array containing 3D gradient matrices, where each matrix represents a gradient in one dimension.
#'
#' @param width The width of the matrices.
#' @param height The height of the matrices.
#' @param start_list A list of starting values for the gradients.
#' @param stop_list A list of ending values for the gradients.
#' @param is_horizontal_list A list of logical values indicating whether the gradient direction for each dimension is horizontal.
#' @return A 3D array where each matrix represents a gradient in one dimension.
#' @export
get_gradation_3d <- function(width, height, start_list, stop_list, is_horizontal_list) {
    result <- array(0, dim = c(height, width, length(start_list)))
    for (i in seq_along(start_list)) {
        start <- start_list[i]
        stop <- stop_list[i]
        is_horizontal <- is_horizontal_list[i]
        result[,,i] <- get_gradation_2d(start, stop, width, height, is_horizontal)
    }
    return(result)
}


#' Create combined arrays based on selected clusters
#'
#' This function creates a list of combined arrays based on the selected clusters and corresponding data matrices.
#'
#' @param selected_clusters A character vector containing the names of selected clusters.
#' @param celltype1 A matrix containing data for celltype1.
#' @param celltype2 A matrix containing data for celltype2.
#' @param celltype3 A matrix containing data for celltype3.
#' @param celltype4 A matrix containing data for celltype4.
#' @return A list of combined arrays.
#' @export
combined_arrays <- function(selected_clusters, celltype1, celltype2, celltype3, celltype4) {
    combined_arrays <- list(
        celltype1,
        celltype2,
        celltype3,
        celltype4
    )
    names(combined_arrays) <- selected_clusters
    return(combined_arrays)
}


