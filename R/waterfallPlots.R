#' Plotting drug-dose response based on the AUC as waterfall plots
#'
#' @description
#' An complementary component of the modular library \pkg{screenwerk} providing a set of kinetic plots for the visualization of drug-dose responses.
#' \emph{\code{waterfallPlots}} is a function that generates a set of waterfall plots that provide an overview of the dose-response between drugs and samples based on the area under the curve (AUC).
#' 
#' @param data an object of class 'processedData' or 'drm'.
#' @param .saveto string; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{waterfallPlots} is used to provide an overview of the responses of all treatments between individual drugs and samples as waterfall plots. This is of particularly benefit for large drug screens
#' in which a large number of drugs and samples have been screened.
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/graphs/waterfall plot'
#' 
#' @examples
#' \donttest{\dontrun{
#' # Plot kinetic dose responses from processed data
#' waterfallPlots(processedData, .saveto = "path/to/folder/")
#' 
#' # Plot kinetic dose responses from the dose response model
#' waterfallPlots(doseRespModel, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen analysis dose response curve waterfall
#' 
#' @importFrom utils tail
#' @importFrom stats ave
#' @importFrom ggplot2 ggplot ggsave labs
#' @importFrom grid gpar textGrob unit viewport
#' @importFrom gridExtra grid.arrange
#' 
#' @export

waterfallPlots <- function(data, .saveto){
  
  # Check, if the data has been provided as an object of class S3:processedData
  if(missing(data)){stop("Data missing! Please provide either a data set as an object of 'processedData' or 'drm'.", call. = TRUE)}
  
  # Extract the required data from the class S3 object
  if(class(data) == "processedData"){
    synDataset <- data[["splitDataset"]]
  } else if(class(data) == "drm") {
    efs <- data
  } else {
    stop("in 'data'. Argument needs to be an object of 'processedData' or 'drm'.", call. = TRUE)
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
  if(!file.exists(file.path(.saveto, "graphs/waterfall plots"))){ dir.create(file.path(.saveto, "graphs/waterfall plots"), showWarnings = FALSE, recursive = TRUE) }
  
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
  
  AUCList <- list(viability = list(), inhibition = list())
  
  for (samplename in listofSamples){
    
    for (drugname in sort(listofDrugs)){
      
      cat('\r', "Calculating the area under the curve (AUC) for ", samplename, ": ", drugname, ".", strrep(" ", 100), sep = "")
      
      .data <- subset(synDataset[[samplename]][["singleDrugResponseData"]], Drug == drugname)
      
      AUCList[["viability"]][[samplename]][[drugname]]   = AUC(1:length(unique(.data$Drug.Concentration)), tapply(.data$Viability, list(.data$Drug, .data$Drug.Concentration), median))
      AUCList[["inhibition"]][[samplename]][[drugname]]  = AUC(1:length(unique(.data$Drug.Concentration)), tapply(.data$Inhibition, list(.data$Drug, .data$Drug.Concentration), median))
      
    }
    
    AUCList[["viability"]][[samplename]] <- AUCList[["viability"]][[samplename]][order(unlist(AUCList[["viability"]][[samplename]]), decreasing = FALSE)]
    AUCList[["inhibition"]][[samplename]] <- AUCList[["inhibition"]][[samplename]][order(unlist(AUCList[["inhibition"]][[samplename]]), decreasing = TRUE)]
    
    cat('\r', strrep(" ", 100), sep = "")
    
  }
  
  cat('\r', "Finished calculating the area under the curve (AUC).", strrep(" ", 100), '\n', sep = "")
  
  
  
  # Plotting dose response as waterfall plot based on AUC
  for(samplename in names(AUCList[["inhibition"]])){
    
    cat('\r', " > Plotting waterfall plots for ", samplename, ".", sep = "")
    
    
    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "graphs/waterfall plots", samplename))){ dir.create(file.path(.saveto, "graphs/waterfall plots", samplename), showWarnings = FALSE, recursive = TRUE) }
    
    
    data <- data.frame("Sample" = samplename, "Drug" = names(AUCList[["inhibition"]][[samplename]]), "AUC" =  unlist(AUCList[["inhibition"]][[samplename]]), row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    data$Drug <- factor(data$Drug, levels = data$Drug)
    data$category <- ifelse(data$AUC < 0, "negative", "positive")
    data$category <- factor(data$category, levels = c("positive", "negative"), labels = c("inhibition", "proliferation"))
    
    # Get the number of unique doses
    noDoses <- unique(as.numeric(tapply(synDataset$`OVCAR-8 R1`$singleDrugResponseData$Drug.Concentration, list(synDataset$`OVCAR-8 R1`$singleDrugResponseData$Sample, synDataset$`OVCAR-8 R1`$singleDrugResponseData$Drug), FUN = function(x) length(unique(x)))))
    
    # Plotting waterfall plot based on full drug list
    ggplot2::ggplot(data, aes(x=Drug, y=AUC, fill = category)) +
      geom_bar(stat="identity", width=0.98, position = position_dodge(width=0)) +
      # scale_fill_gradient2(low = rainbow(2, start = 0.6, end = 0.9)[2], high = rainbow(2, start = 0.6, end = 0.9)[1], midpoint = 0, guide="none") + 
      scale_fill_manual(name="Inhibition", breaks = levels(data$category), values = rainbow(2, start = 0.9, end = 0.6)) + 
      # scale_fill_discrete(name="Inhibition") + 
      scale_color_discrete(guide="none") +
      # coord_cartesian(ylim = c(floor(min(data$AUC)), ceiling(max(data$AUC)))) +
      scale_y_continuous(name = "AUC", breaks = seq(floor(min(data$AUC)), noDoses), limits = c(floor(min(data$AUC)), noDoses)) +
      geom_vline(xintercept=seq(1.5, nrow(data)-0.5, 1), lwd=0.6, colour="white") +
      geom_hline(aes(yintercept=noDoses), linetype = "dotted", lwd=0.6, colour="black") +
      geom_hline(aes(yintercept=ifelse(min(AUC) < 1-noDoses, floor(min(AUC)), 0)), linetype = "dotted", lwd=ifelse(min(data$AUC) < 1-noDoses, 0.6, 0), colour="black") +
      geom_text(aes(nrow(data), noDoses, label = "max.", hjust = 0.5, vjust = -1)) +
      facet_grid(Sample~.) +
      xlab(NULL) +
      ylab("AUC") +
      ggtitle("Drug Response (AUC, Inhibition)") +
      theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 16),
            axis.title.y = element_text(margin = margin(r = 45)),
            axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(25, 75, 0, 110),
            legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(20,0,0,0),
            text = element_text(size = 20), plot.caption = element_text(size = 12), plot.caption.position = "plot",
            plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
            axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1),
            # panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    
    cat('\r', " > Saving waterfall plots for ", samplename, ".", strrep(" ", 100), sep = "")
    
    suppressWarnings(
    ggplot2::ggsave(filename = file.path(.saveto, "graphs/waterfall plots", samplename, paste(samplename, "drug-response (AUC).png", sep = " ")), device = "png", type="cairo-png", width = 819, height = 512, units = "mm", dpi = 300, limitsize = FALSE, bg = "white")
    )
    
    
    
    cat('\r', " > Plotting top-hits as waterfall plots for ", samplename, ".", sep = "")
    
    # Plotting the top-hits, based on N out of full drug list
    ggplot2::ggplot(head(data, floor(nrow(data)*0.3)), aes(x=Drug, y=AUC, fill = category)) +
      geom_bar(stat="identity", width=0.96, position = position_dodge(width=0)) +
      # scale_fill_gradient2(low = rainbow(2, start = 0.6, end = 0.9)[2], high = rainbow(2, start = 0.6, end = 0.9)[1], midpoint = 0, guide="none") + 
      scale_fill_manual(name="Inhibition", breaks = levels(data$category), values = rainbow(2, start = 0.9, end = 0.6)) + 
      # scale_fill_discrete(name="Inhibition") + 
      scale_color_discrete(guide="none") +
      # coord_cartesian(ylim = c(floor(min(data$AUC)), ceiling(max(data$AUC)))) +
      scale_y_continuous(name = "AUC", breaks = seq(0, ceiling(max(data$AUC)), 0.5), limits = c(0, ceiling(max(data$AUC)))) +
      geom_vline(xintercept=seq(1.5, floor(nrow(data)*0.3)-0.5, 1), lwd=0.6, colour="white") +
      geom_hline(aes(yintercept=ifelse(max(head(data$AUC, floor(nrow(data)*0.3))) > noDoses, ceiling(max(AUC)), 0)), linetype = "dotted", lwd=ifelse(max(head(data$AUC, floor(nrow(data)*0.3))) > noDoses, 0.6, 0), colour="black") +
      geom_hline(aes(yintercept=ifelse(min(head(data$AUC, floor(nrow(data)*0.3))) < 1-noDoses, floor(min(AUC)), 0)), linetype = "dotted", lwd=ifelse(min(head(data$AUC, floor(nrow(data)*0.3))) < 1-noDoses, 0.6, 0), colour="black") +
      geom_text(aes(floor(nrow(data)*0.3), ceiling(max(AUC)), label = ifelse(max(head(data$AUC, floor(nrow(data)*0.3))) > noDoses, "max.", ""), hjust = -0.25, vjust = -1)) +
      facet_grid(Sample~.) +
      xlab(NULL) +
      ylab("AUC") +
      ggtitle("Drug Response (AUC, Inhibition)") +
      theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 16),
            axis.title.y = element_text(margin = margin(r = 45)),
            axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(25, 75, 0, 75),
            legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(20,0,0,0),
            text = element_text(size = 20), plot.caption = element_text(size = 12), plot.caption.position = "plot",
            plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
            axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1),
            # panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    
    cat('\r', " > Saving top-hits as waterfall plots for ", samplename, ".", strrep(" ", 100), sep = "")
    
    suppressWarnings(
      ggplot2::ggsave(filename = file.path(.saveto, "graphs/waterfall plots", samplename, paste(samplename, "drug-response (AUC, top-hits).png", sep = " ")), device = "png", type="cairo-png", width = 512, height = 512, units = "mm", dpi = 300, limitsize = FALSE, bg = "white")
    )
    
    
    
    cat('\r', " > Plotting bottom-hits as waterfall plots for ", samplename, ".", sep = "")
    
    # Plotting the bottom-hits, based on N out of full drug list
    ggplot(tail(data, floor(nrow(data)*0.3)), aes(x=Drug, y=AUC, fill = category)) +
      geom_bar(stat="identity", width=0.96, position = position_dodge(width=0)) +
      # scale_fill_gradient2(low = rainbow(2, start = 0.6, end = 0.9)[2], high = rainbow(2, start = 0.6, end = 0.9)[1], midpoint = 0, guide="none") + 
      scale_fill_manual(name="Inhibition", breaks = levels(data$category), values = rainbow(2, start = 0.9, end = 0.6)) + 
      # scale_fill_discrete(name="Inhibition") + 
      scale_color_discrete(guide="none") +
      # coord_cartesian(ylim = c(floor(min(data$AUC)), ceiling(max(data$AUC)))) +
      scale_y_continuous(name = "AUC", breaks = seq(floor(min(tail(data$AUC, floor(nrow(data)*0.3)))), ceiling(max(tail(data$AUC, floor(nrow(data)*0.3)))), 0.5), limits = c(floor(min(tail(data$AUC, floor(nrow(data)*0.3)))), ceiling(max(tail(data$AUC, floor(nrow(data)*0.3)))))) +
      geom_vline(xintercept=seq(1.5, floor(nrow(data)*0.3)-0.5, 1), lwd=0.6, colour="white") +
      geom_hline(aes(yintercept=ifelse(max(tail(data$AUC, floor(nrow(data)*0.3))) > noDoses, ceiling(max(AUC)), 0)), linetype = "dotted", lwd=ifelse(max(tail(data$AUC, floor(nrow(data)*0.3))) > noDoses, 0.6,0), colour="black") +
      geom_hline(aes(yintercept=ifelse(min(tail(data$AUC, floor(nrow(data)*0.3))) < 1-noDoses, floor(min(AUC)), 0)), linetype = "dotted", lwd=ifelse(min(tail(data$AUC, floor(nrow(data)*0.3))) < 1-noDoses, 0.6, 0), colour="black") +
      geom_text(aes(floor(nrow(data)*0.3), ceiling(max(AUC)), label = ifelse(max(tail(data$AUC, floor(nrow(data)*0.3))) > noDoses, "max.", ""), hjust = -0.25, vjust = -1)) +
      facet_grid(Sample~.) +
      xlab(NULL) +
      ylab("AUC") +
      ggtitle("Drug Response (AUC, Inhibition)") +
      theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 16),
            axis.title.y = element_text(margin = margin(r = 45)),
            axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(25, 75, 0, 75),
            legend.position = "bottom", legend.justification="center", legend.direction = "horizontal", legend.box.margin=margin(20,0,0,0),
            text = element_text(size = 20), plot.caption = element_text(size = 12), plot.caption.position = "plot",
            plot.tag = element_text(hjust = 1, size = 18), plot.tag.position = c(1, 0.99),
            axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1),
            # panel.border = element_blank(), panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    
    cat('\r', " > Saving bottom-hits as waterfall plots for ", samplename, ".", strrep(" ", 100), sep = "")
    
    suppressWarnings(
      ggplot2::ggsave(filename = file.path(.saveto, "graphs/waterfall plots", samplename, paste(samplename, "drug-response (AUC, bottom-hits).png", sep = " ")), device = "png", type="cairo-png", width = 512, height = 512, units = "mm", dpi = 300, limitsize = FALSE, bg = "white")
    )
    
    
    
    cat('\r', "Finished plotting waterfall plots for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
    
    
  }
  
  if(samplename == utils::tail(listofSamples, n = 1L)){
    message('\n', "Waterfall plots saved to: ", '\n', file.path(.saveto, "graphs/waterfall plots"), strrep(" ", 100), sep = "")
    rm(samplename, drugname, data)
  }
  
}