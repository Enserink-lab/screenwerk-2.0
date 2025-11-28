#' Plotting drug-dose response based on the AUC as heatmap
#'
#' @description
#' An complementary component of the modular library \pkg{screenwerk} providing a set of plots for the visualization of drug-dose responses.
#' \emph{\code{heatmapPlots}} is a function that generates a heatmap plot as an overview of the dose-response between drugs and samples based on the area under the curve (AUC).
#' 
#' @param data an object of class 'processedData', 'drm' or 'synMetrics'.
#' @param .export \code{vector};  a vector with predefined options for exporting heatmap annotations .
#' @param .saveto string; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{heatmapPlots} is used to provide an overview of the responses of all treatments between individual drugs and samples as heatmaps. This is of particularly benefit for large drug screens
#' in which a large number of drugs and samples have been screened.
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/graphs/heatmaps'
#' 
#' @examples
#' \donttest{\dontrun{
#' # Plot kinetic dose responses from processed data
#' heatmapPlots(processedData, .export = c("dendrogram", "rownames"), .saveto = "path/to/folder/")
#' 
#' # Plot kinetic dose responses from the dose response model
#' heatmapPlots(doseRespModel, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen analysis dose response curve heatmap
#' 
#' @importFrom utils tail
#' @importFrom stats ave
#' @importFrom circlize colorRamp2
#' @importFrom grDevices png dev.off
#' @importFrom ComplexHeatmap Heatmap draw ht_opt row_order row_order column_order
#' 
#' @export

heatmapPlots <- function(data, .export=c("dendrogram"), .saveto){
  
  # Check, if the data has been provided as an object of class S3:processedData
  if(missing(data)){stop("Data missing! Please provide either a data set as an object of 'processedData' or 'drm'.", call. = TRUE)}
  
  # Extract the required data from the class S3 object
  if(class(data) == "processedData"){
    synDataset <- data[["splitDataset"]]
    heatmapMethod = "singledrugresponse"
  } else if(class(data) == "drm") {
    efs <- data
    heatmapMethod = "singledrugresponse"
  } else if(class(data) == "synMetrics") {
    synMetrics <- data
    heatmapMethod = "synergyscores"
  } else {
    stop("in 'data'. Argument needs to be an object of class 'processedData' or 'drm' for plotting heatmaps of single drug responses or an object of class 'synMetrics' for plotting heatmaps of drug interaction scores.", call. = TRUE)
  }
  
  # Check, if elements of the heatmap should be exported
  if(missing(.export)){
    .export = FALSE
  } else if (!.export %in% c("dendrogram")){
    message("Provided '.export' argument not supported. No heatmap annotations will be exported!", '\n')
    .export = FALSE
  }
  
  # Check, if a folder location has been provided
  if(missing(.saveto)){ .saveto <- getwd() }
  # Create folder, if it does not exist
  if(!file.exists(file.path(.saveto))){ dir.create(file.path(.saveto), showWarnings = FALSE, recursive = TRUE) }
  if(!utils::file_test("-d", .saveto)){stop("in '.saveto'. Argument needs to be a valid folder location.", call. = TRUE)}
  
  
  
  if(class(data) == "drm"){
    
    # Extract the original data from the dose response model
    synDataset <- setNames(lapply(names(efs), function(samplename){ 
      synDataset[[samplename]] <- list(singleDrugResponseData = do.call(rbind, setNames(lapply(names(efs[[samplename]]), function(drugname) subset(efs[[samplename]][[drugname]][["drm"]][["origData"]], Drug == drugname)), NULL))); synDataset[[samplename]]
    }), names(efs))
    
  }
  
  
  
  # Create subfolder if it does not exist
  if(!file.exists(file.path(.saveto, "graphs/heatmaps"))){ dir.create(file.path(.saveto, "graphs/heatmaps"), showWarnings = FALSE, recursive = TRUE) }
  
  
  switch (heatmapMethod,
          singledrugresponse = {
  
            listofSamples <- names(synDataset)
            listofDrugs <- unique(unlist(lapply(synDataset, function(x) unique(x[["singleDrugResponseData"]]$Drug))))
            
            
            # Function calculating the area under the curve (AUC) using a simple trapezoidal summation rule
            AUC <- function(x, y){
              auc <-  0
              for (i in 2:(length(x))){
                auc <- auc + (x[i] - x[i-1]) * (y[i] + y[i-1]) * 0.5
              }
              return(auc)
            }
            
            
            heatmapData = matrix(nrow = length(listofDrugs), ncol = length(synDataset), dimnames = list(listofDrugs, names(synDataset)))
            
            for(samplename in names(synDataset)){
              for(drugname in unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug)){
                
                cat('\r', "Calculating the area under the curve (AUC) for ", samplename, ": ", drugname, ".", strrep(" ", 100), sep = "")
                
                
                .data <- subset(synDataset[[samplename]][["singleDrugResponseData"]], Drug == drugname)
                
                .auc = AUC(1:length(.data$Drug.Concentration), .data$Inhibition)
                
                
                heatmapData[drugname, samplename] <- .auc
                
              }
            }
            
            cat('\r', "Finished calculating the area under the curve (AUC).", strrep(" ", 100), '\n', sep = "")
            
            
            # Create subfolder if it does not exist
            if(!file.exists(file.path(.saveto, "graphs/heatmaps/single drug responses"))){ dir.create(file.path(.saveto, "graphs/heatmaps/single drug responses"), showWarnings = FALSE, recursive = TRUE) }
            
            cat('\r', " > Plotting heatmap based on the area under the curve (AUC).", sep = "")
            
            # Plot heatmap with single drug responses based on the AUC in portrait
            grDevices::png(filename = file.path(resultDirectory, "graphs/heatmaps/single drug responses", "heatmap (single drug responses, auc) portrait.png"), 
                           width = 210, height = 297, units = "mm", res = 300, pointsize = 12)
            
            h1 <- ComplexHeatmap::Heatmap(heatmapData, name = "cmb", col = circlize::colorRamp2(c(-10, -1, 0, 1, 10), c("#2066AC", "#BCDAEA", "#F7F7F7", "#FBC9AF", "#B2182B"), space = "LAB"), na_col = "grey95",
                                          width = ncol(heatmapData)*unit(5, "mm"), 
                                          # height = nrow(heatmapData)*unit(0.5, "mm"),
                                          border = FALSE, rect_gp = gpar(col = "white", lwd = 0.1),
                                          row_title = NULL, column_title = NULL, column_title_side = "bottom",
                                          cluster_rows = TRUE, cluster_columns = TRUE,
                                          show_row_dend = TRUE, show_column_dend = TRUE,
                                          row_dend_side = "left", column_dend_side = "top",
                                          row_dend_width = unit(21, "mm"), column_dend_height = unit(21, "mm"),
                                          
                                          clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                                          row_dend_reorder = TRUE, column_dend_reorder = TRUE,
                                          column_km = 1, column_gap = unit(2, "mm"),
                                          row_km = 1, row_gap = unit(2, "mm"),
                                          
                                          show_row_names = TRUE, show_column_names = TRUE,
                                          row_names_side = "right", column_names_side = "bottom",
                                          row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
                                          column_names_rot = 45,
                                          
                                          heatmap_legend_param = list(title = "AUC\n", at = c(-10, 10), labels = c("Excessive\nProliferation", "Inhibition"), 
                                                                      labels_gp = gpar(fontsize = 8), legend_height = unit(16, "mm"), direction = "vertical"))
            
            
            
            
            ComplexHeatmap::draw(h1, newpage = FALSE, column_title = '\nHeatmap (euclidean clustering) of single drug responses based on the area under the curve (AUC)', 
                                 column_title_gp = gpar(fontsize = 8, fontface = "plain"), column_title_side = "bottom", heatmap_legend_side = "right")
            
            cat('\r', " > Saving (AUC) heampap.", strrep(" ", 100), sep = "")
            
            
            grDevices::dev.off()
            
            cat('\r', "Finished plotting (AUC) heatmap.", strrep(" ", 100), '\n', sep = "")
            
            
            message('\n', "Heatmap saved to: ", '\n', file.path(.saveto, "graphs/heatmaps/single drug responses"), strrep(" ", 100), sep = "")
            
            
            
            if(!.export == FALSE){
              
              for (i in .export){
                switch (i,
                        dendrogram = {
                          
                          cat('\n', " > Exporting dendrogram of the heatmap.", sep = "")
                          
                          
                          # Plot just the dendograms and annotations
                          grDevices::png(filename = file.path(resultDirectory, "graphs/heatmaps/single drug responses", "heatmap (column dendrogram).png"),
                                         width = 800, height = 600, units = "mm", res = 300, pointsize = 12)
                          
                          # grDevices::pdf(file = file.path("graphs/heatmaps/single drug responses", paste("heatmap (column dendrogram).pdf", sep = " ")), onefile = FALSE,
                          #     width = (297/2.54)/sqrt(8), height = (210/2.54)/sqrt(8), paper = "special", pointsize = 12, colormodel = "srgb", useKerning = TRUE, compress = FALSE)
                          
                          .textscalefactor = 2
                          
                          h6 <- ComplexHeatmap::Heatmap(heatmapData, name = "h6", border = FALSE, rect_gp = gpar(type = "none"),
                                                        width = ncol(heatmapData)*unit(25, "mm"), height = nrow(heatmapData)*unit(0, "mm"),
                                                        row_title = NULL, column_title = NULL, column_title_side = "bottom",
                                                        row_title_gp = gpar(fontsize = 18*.textscalefactor), column_title_gp = gpar(fontsize = 10*.textscalefactor),
                                                        cluster_rows = FALSE, cluster_columns = TRUE,
                                                        show_row_dend = FALSE, show_column_dend = TRUE,
                                                        row_dend_side = "left", column_dend_side = "top",
                                                        row_dend_width = unit(30*.textscalefactor, "mm"), column_dend_height = unit(20*.textscalefactor, "mm"),
                                                        row_dend_gp = gpar(lwd = 4), column_dend_gp = gpar(lwd = 4),
                                                        
                                                        clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                                                        row_dend_reorder = FALSE, column_dend_reorder = TRUE,
                                                        # row_km = 3, row_gap = unit(2*.textscalefactor, "mm"),
                                                        # row_split = 3, row_gap = unit(2*.textscalefactor, "mm"),
                                                        # column_split = 2, column_gap = unit(2*.textscalefactor, "mm"),
                                                        # cluster_row_slices = FALSE,
                                                        
                                                        show_row_names = FALSE, show_column_names = TRUE,
                                                        row_names_side = "right", column_names_side = "bottom",
                                                        row_names_gp = gpar(fontsize = 16*.textscalefactor), column_names_gp = gpar(fontsize = 16*.textscalefactor),
                                                        column_names_rot = 45, show_heatmap_legend = FALSE
                                                        
                          )
                          
                          ComplexHeatmap::draw(h6, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
                          
                          # decorate_column_dend("h6", {
                          #   grid.yaxis(gp = gpar(fontsize = 16*.textscalefactor))
                          # })
                          
                          grDevices::dev.off()
                          
                          cat('\r', "Finished exporting dendogram.", strrep(" ", 100), '\n', sep = "")
                          
                          
                          message('\n', "Dendogram saved to: ", '\n', file.path(.saveto, "graphs/heatmaps/single drug responses"), strrep(" ", 100), sep = "")
                          
                        })
              }
              
            }
            
          },
          synergyscores = {
            
            # Aggregate data into a single matrix
            heatmapData = list()
            
            # Extract only relevant columns: drug names and synergy score
            for(samplename in names(synMetrics[["bayesynergy"]])){
              cat('\r', " > Extracting synergy and atagonism scores for ", samplename, ".", strrep(" ", 100), sep = "")            
              
              heatmapData[[samplename]] <- synMetrics[["bayesynergy"]][[samplename]][c("Drug A", "Drug B", "Median (syn)", "Median (ant)", "bayesfactor", "divergent", "max_treedepth")]
              names(heatmapData[[samplename]])[names(heatmapData[[samplename]]) == "Drug A"] <- "Drug.A"
              names(heatmapData[[samplename]])[names(heatmapData[[samplename]]) == "Drug B"] <- "Drug.B"
              names(heatmapData[[samplename]])[names(heatmapData[[samplename]]) == "Median (syn)"] <- "Synergy.Score"
              names(heatmapData[[samplename]])[names(heatmapData[[samplename]]) == "Median (ant)"] <- "Antagonism.Score"
              
              heatmapData[[samplename]] <- transform(heatmapData[[samplename]], Drug.A = pmin(Drug.A, Drug.B), Drug.B = pmax(Drug.A, Drug.B), check.names = FALSE)
              heatmapData[[samplename]] <- transform(heatmapData[[samplename]], Synergy.Score = ifelse(Antagonism.Score > Synergy.Score, Antagonism.Score*-1, Synergy.Score), check.names = FALSE)
              heatmapData[[samplename]] <- transform(heatmapData[[samplename]], Synergy.Score = ifelse(divergent > 100, NA, Synergy.Score), check.names = FALSE)
              heatmapData[[samplename]] <- transform(heatmapData[[samplename]], Synergy.Score = ifelse(max_treedepth > 20, NA, Synergy.Score), check.names = FALSE)
              heatmapData[[samplename]] <- transform(heatmapData[[samplename]], Synergy.Score = ifelse(is.na(Synergy.Score), 0, Synergy.Score), check.names = FALSE)
              
              heatmapData[[samplename]] <- heatmapData[[samplename]][, -which(names(heatmapData[[samplename]]) %in% c("Antagonism.Score", "bayesfactor", "divergent", "max_treedepth"))]
              
              names(heatmapData[[samplename]])[names(heatmapData[[samplename]]) == "Synergy.Score"] <- samplename
    
            }
            
            
            # Reduce list of data frames into a single data frame in long format
            heatmapData <- Reduce(function(x, y) merge(x, y, by = c("Drug.A", "Drug.B"), all = TRUE), heatmapData)
            # Combine drug names by drug pairs
            # Convert column with drug names to row names
            heatmapData <- transform(heatmapData, Drug = do.call(paste, c(heatmapData[,grepl("Drug", names(heatmapData))], sep=" + ")))
            rownames(heatmapData) <- heatmapData[["Drug"]]
            # Drop columns with drug names
            heatmapData <- heatmapData[, -which(grepl("Drug", names(heatmapData)))]
            
            # Convert data frame into a matrix
            heatmapData <- as.matrix(heatmapData)
            
            cat('\r', "Finished extracting synergy and atagonism scores.", strrep(" ", 100), sep = "")            
            
            
            # Create subfolder if it does not exist
            if(!file.exists(file.path(.saveto, "graphs/heatmaps/synergy scores"))){ dir.create(file.path(.saveto, "graphs/heatmaps/synergy scores"), showWarnings = FALSE, recursive = TRUE) }
            
            cat('\r', "Plotting heatmap from synergy and antagonism scores.", sep = "")
            
            .textscalefactor = 2
            
            # Full heatmap with all the synergy scores clustering by synergy and antagonism
            grDevices::png(filename = file.path(.saveto, "graphs/heatmaps/synergy scores", "heatmap (synergy scores, median VUS) portrait.png"),
                width = 800*0.45, height = 1085*0.45, units = "mm", res = 300, pointsize = 12)
            
            # grDevices::pdf(file = file.path(.saveto, "graphs/heatmaps/synergy scores", "heatmap (synergy scores, median VUS) portrait.pdf"), onefile = FALSE,
            #     width = (210/2.54)/sqrt(8), height = (297/2.54)/sqrt(8), paper = "special", pointsize = 12, colormodel = "srgb", useKerning = TRUE, compress = FALSE)
            
            ComplexHeatmap::ht_opt(DENDROGRAM_PADDING = unit(2, "mm"))
            h4 <- ComplexHeatmap::Heatmap(heatmapData, name = "cmb", col = colorRamp2(c(-20, -15, -10, 0, 10, 15, 20), c("#2066AC", "#6EA0CB", "#BCDAEA", "#F7F7F7", "#FBC9AF", "#D7716D", "#B2182B"), space = "LAB"), na_col = "#FFFFFF",
                          width = ncol(heatmapData)*unit(20, "mm"), 
                          height = nrow(heatmapData)*unit(1.25, "mm"),
                          border = FALSE, rect_gp = gpar(col = "white", lwd = 0, lty = 0),
                          row_title = "Drug combinations", column_title = NULL, column_title_side = "bottom",
                          row_title_gp = gpar(fontsize = 22*.textscalefactor), column_title_gp = gpar(fontsize = 22*.textscalefactor),
                          cluster_rows = TRUE, cluster_columns = TRUE,
                          show_row_dend = TRUE, show_column_dend = TRUE,
                          row_dend_side = "left", column_dend_side = "top",
                          row_dend_width = unit(30*.textscalefactor, "mm"), column_dend_height = unit(20*.textscalefactor, "mm"),
                          row_dend_gp = gpar(lwd = 2), column_dend_gp = gpar(lwd = 2),
                          
                          clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                          row_dend_reorder = TRUE, column_dend_reorder = TRUE,
                          # row_km = 3, row_gap = unit(2*.textscalefactor, "mm"),
                          # row_split = 2, row_gap = unit(2*.textscalefactor, "mm"),
                          # column_split = 2, column_gap = unit(2*.textscalefactor, "mm"),
                          # cluster_row_slices = FALSE,
                          
                          show_row_names = FALSE, show_column_names = TRUE,
                          row_names_side = "right", column_names_side = "bottom",
                          row_names_gp = gpar(fontsize = 14*.textscalefactor), column_names_gp = gpar(fontsize = 14*.textscalefactor),
                          column_names_rot = 45,
                          
                          # right_annotation = ha,
                          
                          heatmap_legend_param = list(title = NULL, at = c(-20, -15, -10, 0, 10, 15, 20), title_gp = gpar(fontsize = 16*.textscalefactor, fontface = "plain"),
                                                      labels = c("\n   20\nantagonism", "   15", "   10", "    0", " -10", " -15", "synergy\n -20\n"),
                                                      # labels = c("-5", "-1", "0", "1", "5\n\n"),
                                                      # labels = c("-5", "-1", "0", "1", "           synergy >1\n     antagonism <-1\n5\n\n"), 
                                                      labels_gp = gpar(fontsize = 16*.textscalefactor), grid_width = unit(6*.textscalefactor, "mm"), legend_height = unit(48*.textscalefactor, "mm"), direction = "vertical"))
            
            
            ht <- ComplexHeatmap::draw(h4, heatmap_legend_side = "right", annotation_legend_side = "right")
            
            grDevices::dev.off()
            
            cat('\r', "Finished plotting heatmap.", strrep(" ", 100), '\n', sep = "")
            
            
            #### TOP- AND BOTTOM-HITS CLUSTERING ##########################################################################
            
            rownames(heatmapData)[ComplexHeatmap::row_order(ht)]
            
            
            # reorder the matrix according to the clustering
            heatmapData = heatmapData[ComplexHeatmap::row_order(ht), ComplexHeatmap::column_order(ht)]
            
            heatmapDataTop <- head(heatmapData, 25)
            heatmapDataBottom <- tail(heatmapData, 25)
            
            
            cat('\r', "Plotting top-hits.", sep = "")
            
            # Partial heatmap with synergy scores above 1.0 and the top-hits
            grDevices::png(filename = file.path(.saveto, "graphs/heatmaps/synergy scores", "heatmap (top-25 hits, median VUS) portrait.png"),
                width = 800, height = 600, units = "mm", res = 300, pointsize = 12)
            
            # grDevices::pdf(file = file.path(.saveto, "graphs/heatmaps/synergy scores", "heatmap (top-25 hits, median VUS) portrait.pdf"), onefile = FALSE,
            #     width = (297/2.54)/sqrt(8), height = (210/2.54)/sqrt(8), paper = "special", pointsize = 12, colormodel = "srgb", useKerning = TRUE, compress = FALSE)
            
            ComplexHeatmap::ht_opt(DENDROGRAM_PADDING = unit(2, "mm"))
            ComplexHeatmap::ht_opt(ROW_ANNO_PADDING = unit(2, "mm"))
            
            h1 <- ComplexHeatmap::Heatmap(heatmapDataTop, name = "cmb", col = colorRamp2(c(-20, -15, -10, 0, 10, 15, 20), c("#2066AC", "#6EA0CB", "#BCDAEA", "#F7F7F7", "#FBC9AF", "#D7716D", "#B2182B"), space = "LAB"), na_col = "#FFFFFF",
                          width = ncol(heatmapDataTop)*unit(20, "mm"), height = nrow(heatmapDataTop)*unit(10, "mm"),
                          border = FALSE, rect_gp = gpar(col = "white", lwd = 0.1, lty = 0),
                          row_title = "Drug combinations (Top-25)", column_title = NULL, column_title_side = "bottom",
                          row_title_gp = gpar(fontsize = 22*.textscalefactor), column_title_gp = gpar(fontsize = 22*.textscalefactor),
                          cluster_rows = FALSE, cluster_columns = FALSE,
                          show_row_dend = TRUE, show_column_dend = TRUE,
                          row_dend_side = "left", column_dend_side = "top",
                          row_dend_width = unit(30*.textscalefactor, "mm"), column_dend_height = unit(20*.textscalefactor, "mm"),
                          row_dend_gp = gpar(lwd = 3), column_dend_gp = gpar(lwd = 3),
                          
                          clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                          row_dend_reorder = FALSE, column_dend_reorder = FALSE,
                          # row_km = 3, row_gap = unit(2*.textscalefactor, "mm"),
                          # row_split = 3, row_gap = unit(2*.textscalefactor, "mm"),
                          # column_split = 8, column_gap = unit(2*.textscalefactor, "mm"),
                          # cluster_row_slices = FALSE,
                          
                          show_row_names = TRUE, show_column_names = TRUE,
                          row_names_side = "right", column_names_side = "bottom",
                          row_names_gp = gpar(fontsize = 14*.textscalefactor), column_names_gp = gpar(fontsize = 12*.textscalefactor),
                          column_names_rot = 45,
                          
                          # right_annotation = ha,
                          
                          heatmap_legend_param = list(title = "antagonism synergy ", at = c(-20, -15, -10, 0, 10, 15, 20), title_gp = gpar(fontsize = 12*.textscalefactor, fontface = "plain"),
                                                      # labels = c("-5", "-1 > antagonism", "0", "1 < synergy", "5"),
                                                      labels = c("20 ", "15", "10 ", " 0 ", "-10", "-15", "-20 \n\n"),
                                                      # labels = c("-5", "-1", "0", "1", "           synergy >1\n     antagonism <-1\n5\n\n"), 
                                                      labels_gp = gpar(fontsize = 14*.textscalefactor), grid_width = unit(6*.textscalefactor, "mm"), legend_width = unit(38*.textscalefactor, "mm"), direction = "horizontal"))
            
            ComplexHeatmap::draw(h1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", legend_gap = unit(25, "mm"))
            
            grDevices::dev.off()
            
            cat('\r', "Finished plotting top-hits.", strrep(" ", 100), '\n', sep = "")
            
            
            
            cat('\r', "Plotting bottom-hits.", sep = "")
            
            # Partial heatmap with synergy scores above 1.0 and the top-hits
            grDevices::png(filename = file.path(.saveto, "graphs/heatmaps/synergy scores", "heatmap (bottom-25 hits, median VUS) portrait.png"),
                           width = 800, height = 600, units = "mm", res = 300, pointsize = 12)
            
            # grDevices::pdf(file = file.path(.saveto, "graphs/heatmaps/synergy scores", "heatmap (bottom-25 hits, median VUS) portrait.pdf"), onefile = FALSE,
            #     width = (297/2.54)/sqrt(8), height = (210/2.54)/sqrt(8), paper = "special", pointsize = 12, colormodel = "srgb", useKerning = TRUE, compress = FALSE)
            
            ComplexHeatmap::ht_opt(DENDROGRAM_PADDING = unit(2, "mm"))
            ComplexHeatmap::ht_opt(ROW_ANNO_PADDING = unit(2, "mm"))
            
            h2 <- ComplexHeatmap::Heatmap(heatmapDataBottom, name = "cmb", col = colorRamp2(c(-20, -15, -10, 0, 10, 15, 20), c("#2066AC", "#6EA0CB", "#BCDAEA", "#F7F7F7", "#FBC9AF", "#D7716D", "#B2182B"), space = "LAB"), na_col = "#FFFFFF",
                          width = ncol(heatmapDataBottom)*unit(20, "mm"), height = nrow(heatmapDataBottom)*unit(10, "mm"),
                          border = FALSE, rect_gp = gpar(col = "white", lwd = 0.1, lty = 0),
                          row_title = "Drug combinations (Bottom-25)", column_title = NULL, column_title_side = "bottom",
                          row_title_gp = gpar(fontsize = 22*.textscalefactor), column_title_gp = gpar(fontsize = 22*.textscalefactor),
                          cluster_rows = FALSE, cluster_columns = FALSE,
                          show_row_dend = TRUE, show_column_dend = TRUE,
                          row_dend_side = "left", column_dend_side = "top",
                          row_dend_width = unit(30*.textscalefactor, "mm"), column_dend_height = unit(20*.textscalefactor, "mm"),
                          row_dend_gp = gpar(lwd = 3), column_dend_gp = gpar(lwd = 3),
                          
                          clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                          row_dend_reorder = FALSE, column_dend_reorder = FALSE,
                          # row_km = 3, row_gap = unit(2*.textscalefactor, "mm"),
                          # row_split = 3, row_gap = unit(2*.textscalefactor, "mm"),
                          # column_split = 8, column_gap = unit(2*.textscalefactor, "mm"),
                          # cluster_row_slices = FALSE,
                          
                          show_row_names = TRUE, show_column_names = TRUE,
                          row_names_side = "right", column_names_side = "bottom",
                          row_names_gp = gpar(fontsize = 14*.textscalefactor), column_names_gp = gpar(fontsize = 12*.textscalefactor),
                          column_names_rot = 45,
                          
                          # right_annotation = ha,
                          
                          heatmap_legend_param = list(title = "antagonism synergy ", at = c(-20, -15, -10, 0, 10, 15, 20), title_gp = gpar(fontsize = 12*.textscalefactor, fontface = "plain"),
                                                      # labels = c("-5", "-1 > antagonism", "0", "1 < synergy", "5"),
                                                      labels = c("20 ", "15", "10 ", " 0 ", "-10", "-15", "-20 \n\n"),
                                                      # labels = c("-5", "-1", "0", "1", "           synergy >1\n     antagonism <-1\n5\n\n"), 
                                                      labels_gp = gpar(fontsize = 14*.textscalefactor), grid_width = unit(6*.textscalefactor, "mm"), legend_width = unit(38*.textscalefactor, "mm"), direction = "horizontal"))
            
            ComplexHeatmap::draw(h2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", legend_gap = unit(25, "mm"))
            
            grDevices::dev.off()
            
            cat('\r', "Finished plotting bottom-hits.", strrep(" ", 100), '\n', sep = "")
            
            
            message('\n', "Heatmaps saved to: ", '\n', file.path(.saveto, "graphs/heatmaps/synergy scores"), strrep(" ", 100), sep = "")
            
            
            if(!.export == FALSE){
              
              for (i in .export){
                switch (i,
                        dendrogram = {
                          
                          cat('\n', " > Exporting dendrogram of the heatmap.", sep = "")
                          
                          
                          # Plot just the dendograms and annotations
                          grDevices::png(filename = file.path(resultDirectory, "graphs/heatmaps/synergy scores", "heatmap (column dendrogram).png"),
                                         width = 800, height = 600, units = "mm", res = 300, pointsize = 12)
                          
                          # grDevices::pdf(file = file.path("graphs/heatmaps/synergy scores", paste("heatmap (column dendrogram).pdf", sep = " ")), onefile = FALSE,
                          #     width = (297/2.54)/sqrt(8), height = (210/2.54)/sqrt(8), paper = "special", pointsize = 12, colormodel = "srgb", useKerning = TRUE, compress = FALSE)
                          
                          .textscalefactor = 2
                          
                          h6 <- ComplexHeatmap::Heatmap(heatmapData, name = "h6", border = FALSE, rect_gp = gpar(type = "none"),
                                                        width = ncol(heatmapData)*unit(25, "mm"), height = nrow(heatmapData)*unit(0, "mm"),
                                                        row_title = NULL, column_title = NULL, column_title_side = "bottom",
                                                        row_title_gp = gpar(fontsize = 18*.textscalefactor), column_title_gp = gpar(fontsize = 10*.textscalefactor),
                                                        cluster_rows = FALSE, cluster_columns = TRUE,
                                                        show_row_dend = FALSE, show_column_dend = TRUE,
                                                        row_dend_side = "left", column_dend_side = "top",
                                                        row_dend_width = unit(30*.textscalefactor, "mm"), column_dend_height = unit(20*.textscalefactor, "mm"),
                                                        row_dend_gp = gpar(lwd = 4), column_dend_gp = gpar(lwd = 4),
                                                        
                                                        clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                                                        row_dend_reorder = FALSE, column_dend_reorder = TRUE,
                                                        # row_km = 3, row_gap = unit(2*.textscalefactor, "mm"),
                                                        # row_split = 3, row_gap = unit(2*.textscalefactor, "mm"),
                                                        # column_split = 2, column_gap = unit(2*.textscalefactor, "mm"),
                                                        # cluster_row_slices = FALSE,
                                                        
                                                        show_row_names = FALSE, show_column_names = TRUE,
                                                        row_names_side = "right", column_names_side = "bottom",
                                                        row_names_gp = gpar(fontsize = 16*.textscalefactor), column_names_gp = gpar(fontsize = 16*.textscalefactor),
                                                        column_names_rot = 45, show_heatmap_legend = FALSE
                                                        
                          )
                          
                          ComplexHeatmap::draw(h6, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
                          
                          # decorate_column_dend("h6", {
                          #   grid.yaxis(gp = gpar(fontsize = 16*.textscalefactor))
                          # })
                          
                          grDevices::dev.off()
                          
                          cat('\r', "Finished exporting dendogram.", strrep(" ", 100), '\n', sep = "")
                          
                          
                          message('\n', "Dendogram saved to: ", '\n', file.path(.saveto, "graphs/heatmaps/synergy scores"), strrep(" ", 100), sep = "")
                          
                        },
                        rownames = {
                          
                          cat('\n', " > Exporting rownames of the heatmap.", sep = "")
                          
                          
                          cat(readLines(file.path(.saveto, "graphs/heatmaps/synergy scores", "Heatmap Clustering Drug Row Names.txt")), sep = "\n")
                          
                          message('\n', "Rownames saved to: ", '\n', file.path(.saveto, "graphs/heatmaps/synergy scores"), strrep(" ", 100), sep = "")
                          
                          
                        })
              }
              
            }
            
            
          }
          
  
  )
}