#' Generating_simulation_data
#'
#' A function that generates simulated data. This function is used to generate simulated experimental data and
#' save the results to a file.And the directory path for the output result is the current folder
#'
#' @param spot_spatial_count A vector representing the UMI count of each cell.
#'
#' spot_spatial_count[1:5, 1:5]
#'   A730017C20Rik Aamp Abat Abhd3 Acat2
#' 0             0    0    1     0     0
#' 1             0    0    1     0     2
#' 2             0    0    0     0     0
#' 3             0    0    0     0     1
#' 4             0    0    0     0     1
#'
#' @param sim_spatial_spot_loc A table representing the simulated spatial location.
#' Each row of the table represents the coordinates of a spatial point, and each column
#' represents a dimension (such as x, y coordinates).
#'
#' You must ensure that the column names "X" and "Y" exist in this file.
#' sim_spatial_spot_loc[1:5, ]
#'          X         Y
#' 1 24.08588  8867.800
#' 2 26.38492 10853.262
#' 3 49.77125  4980.778
#' 4 53.47311  8729.177
#' 5 38.32856  8978.149
#'
#' @param cell_celltype A table representing the type of each cell.Each element represents
#' the type of a cell, which can be a character vector or factor.
#'
#' You must ensure that the column names "Cluster" and "celltype" exist in this file.
#' meta_data[1:5, ]
#'   Cluster         celltype
#' 0       0     ExcitatoryL4
#' 1       1 ExcitatoryL2and3
#' 2       2     ExcitatoryL6
#' 3       0     ExcitatoryL4
#' 4       0     ExcitatoryL4
#'
#' @param grid_interval Spacing size of spatial grids.Used to map the simulated spatial
#' positions onto the grid for subsequent spatial analysis.
#'
#' @examples
#' Generating_simulation_data(spot_spatial_count = "Spatial_count.txt",
#'            sim_spatial_spot_loc = Locations.txt",
#'            cell_celltype = "Spatial_annotate.txt",
#'            grid_interval = 500
#'            prefix = “prefix”)
#'
Generating_simulation_data <- function(spot_spatial_count, sim_spatial_spot_loc, cell_celltype, grid_interval, prefix) {
  library(ggplot2)
  library(scatterpie)
  library(randomcoloR)
  library(igraph)

  ## get data ##
  data <- read.table(sim_spatial_spot_loc, header = TRUE)
  meta_data <- read.table(cell_celltype, header = TRUE)
  gene_count <- read.table(spot_spatial_count, header = TRUE)

  columns_to_load <- c("X", "Y")
  meta_data <- cbind(meta_data, data[columns_to_load])
  colnames(meta_data) <- gsub("Cluster", "cluster_id", colnames(meta_data))
  colnames(meta_data) <- gsub("celltype", "cluster_name", colnames(meta_data))
  num_columns <- dim(meta_data)[1] - 1
  index <- sprintf("Dat1_Cell_%05d", 0:num_columns)
  meta_data <- cbind(meta_data, index)
  rownames(gene_count) <- index

  ## first 12 gene names need to exclude the first Character 'X'
  #for (i in 1:ncol(gene_count)) {
    #if (substr(colnames(gene_count)[i], 1, 1) == "X") {
      #colnames(gene_count)[i] <- substr(colnames(gene_count)[i], 2, nchar(colnames(gene_count)[i]))
    #}
  #}

  # check index
  stopifnot(all(meta_data$index == row.names(gene_count)))

  ## create pseudo spot ##
  dat <- meta_data
  count <- gene_count  # cell x gene : 2002 x 1020
  x <- "X"  # x axis
  y <- "Y"  # y axis
  index <- "index"  # cell ID
  cluster <- "cluster_name"  # cell type
  scaler <- 2 * grid_interval  # scale factor
  breaks <- grid_interval  # grid interval

  # get sequencing depth for each cell
  seq_depth <- rowSums(count)

  # convert cell-type name into a factor
  if (is.factor(dat[, cluster]) == 0) {
    dat[, cluster] <- factor(dat[, cluster])
  }

  ## generate grid ##
  range_x <- range(dat[, x])  # get x axis range
  x_start <- floor(range_x[1])
  x_end <- ceiling(range_x[2])
  grid_x <- seq(x_start, x_end, breaks)

  range_y <- range(dat[, y])
  y_start <- floor(range_y[1])
  y_end <- ceiling(range_y[2])
  grid_y <- seq(y_start, y_end, breaks)

  # make sure the grid is identical
  if (max(grid_x) < x_end) {
    grid_x <- c(grid_x, max(grid_x) + breaks)
  }
  if (max(grid_y) < y_end) {
    grid_y <- c(grid_y, max(grid_y) + breaks)
  }

  # init grad table: number of cells, x, y
  grid_table <- data.frame(ncell = numeric(), x = numeric(), y = numeric())
  # init cell-type table
  ct_table <- grid_table

  spot_exp <- list()

  # number of cells
  ncell_fov <- nrow(dat)

  ## find the cells located in each spot and aggregate
  # length(grid_y): 36, length(grid_x):23
  for (i in 1:(length(grid_y) - 1)) {
    for (j in 1:(length(grid_x) - 1)) {  # grid_x[j], grid_x[j + 1], grid_y[i] and grid_y[i + 1] is the target spot
      # find the cells located in target spot
      fall_id <- which(dat[, x] >= grid_x[j] & dat[, x] < grid_x[j + 1] & dat[, y] >= grid_y[i] & dat[, y] < grid_y[i + 1])

      # add spot information: number of cells, x, y
      grid_table <- rbind(grid_table, c(ncell = length(fall_id), x = grid_x[j], y = grid_y[i]))

      # add spot information: number of cell types, x, y
      ct_table <- rbind(ct_table, c(nct = length(unique(dat[fall_id, cluster])), x = grid_x[j], y = grid_y[i]))

      # generate spot expression data
      spot_exp <- c(spot_exp, list(ceiling(colSums(count[dat[fall_id, index], ]/(seq_depth[dat[fall_id, index]] * length(fall_id)) * scaler))))

      # set names to the spot using coordinate
      names(spot_exp)[length(spot_exp)] <- paste0("x", grid_x[j], "y", grid_y[i])

      # compute proportions for each cell types
      tmp_prop <- data.frame(cbind(t(as.matrix(table(dat[fall_id, cluster])/sum(table(dat[fall_id, cluster])), nrow = 1, drop = F)), x = grid_x[j], y = grid_y[i]))

      if (i == 1 & j == 1) {
        grid_prop <- tmp_prop
      } else {
        grid_prop <- rbind(grid_prop, tmp_prop)
      }

    }
  }

  # rename tables
  colnames(grid_table) <- c("ncell", "x", "y")
  colnames(ct_table) <- c("nct", "x", "y")

  # convert expression list to matrix
  spot_exp_mat <- Reduce(cbind, spot_exp)
  colnames(spot_exp_mat) <- names(spot_exp)

  # filter empty spots
  nonempty_cell <- which(grid_table$ncell != 0)
  grid_table <- grid_table[nonempty_cell, ]
  ct_table <- ct_table[nonempty_cell, ]
  grid_prop <- grid_prop[nonempty_cell, ]
  spot_exp_mat <- spot_exp_mat[, nonempty_cell]

  # aggregate results
  starmap_res <- list(ncell_grid = grid_table, nct_grid = ct_table, prop_grid = grid_prop, spatial_exp = spot_exp_mat)

  my.distinct.colors20 <- randomColor(count = 20)
  plot_Locations_of_cell_with_grid <- ggplot(meta_data, aes(x=X, y=Y, color=cluster_name)) +
    geom_point() +
    scale_x_continuous(minor_breaks=seq(floor(min(meta_data$X) / 10) * 10, max(meta_data$X) + breaks, breaks),
                       breaks = seq(floor(min(meta_data$X) / 10) * 10, max(meta_data$X), grid_interval)) +
    scale_y_continuous(minor_breaks=seq(floor(min(meta_data$Y) / 10) * 10, max(meta_data$Y) + breaks, breaks),
                       breaks = seq(floor(min(meta_data$Y) / 10) * 10, max(meta_data$Y), grid_interval)) +
    theme_gray() +
    theme(legend.title=element_blank(), legend.text=element_text(size=12, ),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid'),
          panel.grid.minor = element_line(size = 0.5, linetype = 'solid'))

  # Save plot object as a PDF file
  ggsave(paste0(file_prefix, "_plot_Locations_of_cell_with_grid.pdf"), plot = plot_Locations_of_cell_with_grid, device = "pdf")

  ## rename cell type eL2.3 to eL2/3
  #colnames(starmap_res$prop_grid)[2] = "eL2/3"

  Pie_chart_of_cell_type_proportions <- ggplot() +
    geom_scatterpie(aes(x=x, y=y), data=starmap_res$prop_grid, cols=colnames(starmap_res$prop_grid)[1:12]) +
    scale_fill_manual(values=my.distinct.colors20) +
    theme_classic() +
    theme(legend.title=element_blank(), legend.text=element_text(size=12))

  # Save plot object as a PDF file
  ggsave(paste0(file_prefix, "_Pie_chart_of_cell_type_proportions.pdf"), plot = Pie_chart_of_cell_type_proportions, device = "pdf")

  # generate spot location dataframe
  local_df = starmap_res$prop_grid[, c('x', 'y')]
  row.names(local_df) = paste0('x', local_df$x, 'y', local_df$y)
  colnames(local_df) = c('imagecol', 'imagerow')


  #getNeighbour = function(array_row, array_col) {
    ## based on the (row, col) of one spot, return the (row, col) of all 6 neighbours
    #return(list(c(array_row-500, array_col),
                #c(array_row+500, array_col),
                #c(array_row+0, array_col-500),
                #c(array_row+0, array_col+500)))
    #}


  # adjacency matrix
  A = matrix(0, nrow = nrow(local_df), ncol = nrow(local_df))
  row.names(A) = rownames(local_df)
  colnames(A) = rownames(local_df)

  #for (i in 1:nrow(local_df)) {
    #barcode = rownames(local_df)[i]
    #array_row = local_df[i, 'imagerow']
    #array_col = local_df[i, 'imagecol']

    ## get neighbors
    #neighbours = getNeighbour(array_row, array_col)

    ## fill the adjacency matrix
    #for (this.vec in neighbours) {
      #tmp.p = rownames(local_df[local_df$imagerow==this.vec[1] & local_df$imagecol==this.vec[2], ])

      #if (length(tmp.p) >= 1) {
        ## target spots have neighbors in selected spots
        #for (neigh.barcode in tmp.p) {
          #A[barcode, neigh.barcode] = 1
        #}
      #}
    #}
  #}

  #g = graph_from_adjacency_matrix(A, 'undirected', add.colnames = NA, add.rownames = NA)
  # manually set nodes x and y coordinates
  #vertex_attr(g, name = 'x') = local_df$imagecol
  #vertex_attr(g, name = 'y') = local_df$imagerow
  #plot(g, vertex.size=5, edge.width=4, margin=-0.05)

  to_save = t(starmap_res$spatial_exp)

  write.csv(to_save, paste0(file_prefix, '_sim_spatial_spot_nUMI.csv'))
  print(sprintf('save %d gene nUMIs of %d simulated spatial pseudo-spots into file %s', ncol(to_save), nrow(to_save), paste0(file_prefix, '_sim_spatial_spot_nUMI.csv')))

  local_df$x = floor(local_df$imagecol/breaks) + 1
  local_df$y = floor(local_df$imagerow/breaks) + 1

  write.csv(local_df, paste0(file_prefix, '_sim_spatial_spot_loc.csv'))
  print(sprintf('save Physical Locations of simulated spatial pseudo-spots into file %s', paste0(file_prefix, '_sim_spatial_spot_loc.csv')))

  #write.csv(A, 'sim_spatial_spot_adjacency_matrix.csv')
  #print(sprintf('save Adjacency Matrix of simulated spatial pseudo-spots into file %s', 'sim_spatial_spot_adjacency_matrix.csv'))

  to_save = starmap_res$prop_grid
  row.names(to_save) = paste0('x', to_save$x, 'y', to_save$y)
  to_save = to_save[, !(colnames(to_save) %in% c('x', 'y'))]

  write.csv(to_save, paste0(file_prefix, '_sim_spatial_spot_true_prop.csv'))
  print(sprintf('save True cell type proportions of simulated spatial pseudo-spots into file %s', paste0(file_prefix, '_sim_spatial_spot_true_prop.csv')))

  write.csv(gene_count, paste0(file_prefix, '_cells_nUMI.csv'))
  print(sprintf('save %d gene nUMIs of %d spot cells into file %s', ncol(gene_count), nrow(gene_count), paste0(file_prefix, '_cells_nUMI.csv')))

  to_save = meta_data['cluster_name']
  row.names(to_save) = meta_data$index
  colnames(to_save) = 'celltype'
  to_save$celltype = as.character(to_save$celltype)
  to_save[1:5, , drop=F]

  write.csv(to_save, paste0(file_prefix, '_cells_celltype.csv'))
  print(sprintf('save cell type annotation of spot cells into file %s', paste0(file_prefix, '_cells_celltype.csv')))
}


