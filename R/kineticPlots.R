#' Plotting kinetic drug-dose response
#'
#' @description
#' An complementary component of the modular library \pkg{screenwerk} providing a set of kinetic plots for the visualization of drug-dose responses.
#' \emph{\code{kineticPlots}} is a function that generates a set of plots that provide an overview of the raw and unfitted dose-response between drugs and samples.
#' 
#' @param data an object of class 'processedData' or 'drm'.
#' @param .saveto string; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{kineticPlots} is used to provide an overview of the raw and unfitted responses of all treatments between individual drugs and samples. This is of particularly benefit for large drug screens
#' in which a large number of drugs and samples have been screened.
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/graphs/single drug response'
#' 
#' @examples
#' \donttest{\dontrun{
#' # Plot kinetic dose responses from processed data
#' kineticPlots(processedData, .saveto = "path/to/folder/")
#' 
#' # Plot kinetic dose responses from the dose response model
#' kineticPlots(doseRespModel, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen analysis dose response curve kinetic
#' 
#' @importFrom utils tail
#' @importFrom stats ave
#' @importFrom ggplot2 ggplot ggsave labs
#' @importFrom grid gpar textGrob unit viewport
#' @importFrom gridExtra grid.arrange
#' 
#' @export

kineticPlots <- function(data, .saveto){
  
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
  if(!file.exists(file.path(.saveto, "graphs/single drug response"))){ dir.create(file.path(.saveto, "graphs/single drug response"), showWarnings = FALSE, recursive = TRUE) }
  
  listofSamples <- names(synDataset)
  listofDrugs <- unique(unlist(lapply(synDataset, function(x) unique(x[["singleDrugResponseData"]]$Drug))))
  
  
  # Plotting dose response curves based on viability.
  # Plot dose response curves as a composite.
  # all-replicates, unfitted without points
  # alphabetical order
  
  cat("Plotting kinetic dose response based on viability, all-samples, unfitted w/o points.", '\n', sep = "")
  
  
  for(p in 1:ceiling(length(listofDrugs)/60)){
    
    grDevices::png(filename = file.path(.saveto, "graphs/single drug response", paste("dose-response (viability) (all-samples, unfitted wo points) [", p, "of", ceiling(length(sort(names(efs[[1]])))/60), "].png", sep = "")), 
        width = (7920/10)*10, height = (4080/6)*6, units = "px", pointsize = 24)
    op <- graphics::par(mfrow = c(6, 10), oma = c(4,0,7,0))
    
    for(drugname in na.omit(listofDrugs[((60*p)-60+1):(60*p)])){
      
      cat('\r', " > Plotting kinetic dose-response curve for ", drugname, ".", strrep(" ", 100), sep = "")
      
      
      for(samplename in listofSamples){
        
        drm.data <- do.call(rbind, setNames(lapply(names(synDataset), function(x) synDataset[[x]][["singleDrugResponseData"]]), NULL))
        drm.data <- subset(drm.data, Drug == drugname & Sample == samplename)
        drm.data$`Viability (Median)` <- with(drm.data, ave(Viability, Drug.Concentration, FUN = median))
        drm.data$`Viability (Median)` <- with(drm.data, ifelse(`Viability (Median)` > 1, 1, `Viability (Median)`))
        
        
        if(samplename == head(listofSamples, n = 1)){
          # Draw empty plot by type = n, adding labels and axis
          graphics::plot(drm.data$Drug.Concentration, drm.data$Viability, log = "x", type = "n", main = drugname, xlab = paste("Concentration [", unique(drm.data$Unit) , "]", sep = ""), ylab = "Viability",
               axes=FALSE, pch=1, cex=1, lty=1, lwd=2, ylim=c(0, 1), yaxs="i", xaxs="r", bty="n")
          
          if(samplename == head(names(synDataset), n = 1)){graphics::axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
          if(samplename == head(names(synDataset), n = 1)){graphics::axis(2)}
        }
        
        # graphics::grid(nx = NULL, ny = "", col = "lightgrey", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
        # graphics::abline(v=outer((1:10),(10^(-5:5))), col="grey80", lty = "dotted", lwd = 0.3)
        # graphics::abline(h=seq(0,1,0.1), col="grey80", lty = "dotted", lwd = 0.3)
        
        graphics::lines(drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)], drm.data$`Viability (Median)`[order(drm.data$Drug.Concentration)], lwd = 0.8)
        # graphics::points(drm.data$Drug.Concentration, drm.data$`Viability (Mean)`, type="p", pch=16, cex = 0.8, col="#FF0066")
        
        }
      
    }
    
    cat('\r', "Saving kinetic dose response plots.", strrep(" ", 100), sep = "")
    
    
    grDevices::dev.off()
    
    cat('\r', "Finished plotting kinetic dose response based on viability, all-samples, unfitted w/o points.", strrep(" ", 100), '\n', sep = "")
    
    
    if(drugname == tail(sort(names(efs[[1]])), n = 1)){
      message('\n', "Kinetic dose response (viability) plots (all-samples, unfitted wo points) saved to: ", '\n', file.path(.saveto, "results", "graphs/single drug response"), strrep(" ", 100), sep = "")
      rm(samplename, drugname, drm.data, op)
    }
    
  }
  
  
  
  
  # Plotting dose response curves based on viability.
  # Plot dose response curves as a composite.
  # all-replicates, unfitted with points
  # alphabetical order
  
  cat("Plotting kinetic dose response based on viability, all-samples, unfitted w/ points.", '\n', sep = "")
  
  
  for(p in 1:ceiling(length(listofDrugs)/60)){
    
    grDevices::png(filename = file.path(.saveto, "graphs/single drug response", paste("dose-response (viability) (all-samples, unfitted w points) [", p, "of", ceiling(length(sort(names(efs[[1]])))/60), "].png", sep = "")), 
                   width = (7920/10)*10, height = (4080/6)*6, units = "px", pointsize = 24)
    op <- graphics::par(mfrow = c(6, 10), oma = c(4,0,7,0))
    
    for(drugname in na.omit(listofDrugs[((60*p)-60+1):(60*p)])){
      
      cat('\r', " > Plotting kinetic dose-response curve for ", drugname, ".", strrep(" ", 100), sep = "")
      
      
      for(samplename in listofSamples){
        
        drm.data <- do.call(rbind, setNames(lapply(names(synDataset), function(x) synDataset[[x]][["singleDrugResponseData"]]), NULL))
        drm.data <- subset(drm.data, Drug == drugname & Sample == samplename)
        # drm.data$`Viability (Median)` <- with(drm.data, ave(Viability, Drug.Concentration, FUN = median))
        # drm.data$`Viability (Median)` <- with(drm.data, ifelse(`Viability (Median)` > 1, 1, `Viability (Median)`))
        
        
        if(samplename == head(listofSamples, n = 1)){
          # Draw empty plot by type = n, adding labels and axis
          graphics::plot(drm.data$Drug.Concentration, drm.data$Viability, log = "x", type = "n", main = drugname, xlab = paste("Concentration [", unique(drm.data$Unit) , "]", sep = ""), ylab = "Viability",
                         axes=FALSE, pch=1, cex=1, lty=1, lwd=2, ylim=c(0, 1.2), yaxs="i", xaxs="r", bty="n")
          
          if(samplename == head(names(synDataset), n = 1)){graphics::axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
          if(samplename == head(names(synDataset), n = 1)){graphics::axis(2, at=seq(from = 0, to = 1, by = 0.2), labels=format(seq(from = 0, to = 1, by = 0.2), trim = TRUE))}
        }
        
        # graphics::grid(nx = NULL, ny = "", col = "lightgrey", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
        # graphics::abline(v=outer((1:10),(10^(-5:5))), col="grey80", lty = "dotted", lwd = 0.3)
        # graphics::abline(h=seq(0,1,0.1), col="grey80", lty = "dotted", lwd = 0.3)
        
        # graphics::lines(drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)], drm.data$`Viability (Median)`[order(drm.data$Drug.Concentration)], lwd = 0.6)
        graphics::points(drm.data$Drug.Concentration, drm.data$Viability, type="p", pch=16, cex = 0.8, col="#FF0066")
        
      }
      
    }
    
    cat('\r', "Saving kinetic dose response plots.", strrep(" ", 100), sep = "")
    
    
    grDevices::dev.off()
    
    cat('\r', "Finished plotting kinetic dose response based on viability, all-samples, unfitted w/ points.", strrep(" ", 100), '\n', sep = "")
    
    
    if(drugname == tail(sort(names(efs[[1]])), n = 1)){
      message('\n', "Kinetic dose response (viability) plots (all-samples, unfitted w points) saved to: ", '\n', file.path(.saveto, "results", "graphs/single drug response"), strrep(" ", 100), sep = "")
      rm(samplename, drugname, drm.data, op)
    }
    
  }
  
  
  
  
  # Plotting dose response curves based on viability.
  # Plot dose response curves as a composite.
  # all-replicates, unfitted with median points
  # alphabetical order
  
  cat("Plotting kinetic dose response based on viability, all-samples, unfitted w/ median points.", '\n', sep = "")
  
  
  for(p in 1:ceiling(length(listofDrugs)/60)){
    
    grDevices::png(filename = file.path(.saveto, "graphs/single drug response", paste("dose-response (viability) (all-samples, unfitted w median points) [", p, "of", ceiling(length(sort(names(efs[[1]])))/60), "].png", sep = "")), 
                   width = (7920/10)*10, height = (4080/6)*6, units = "px", pointsize = 24)
    op <- graphics::par(mfrow = c(6, 10), oma = c(4,0,7,0))
    
    for(drugname in na.omit(listofDrugs[((60*p)-60+1):(60*p)])){
      
      cat('\r', " > Plotting kinetic dose-response curve for ", drugname, ".", strrep(" ", 100), sep = "")
      
      
      for(samplename in listofSamples){
        
        drm.data <- do.call(rbind, setNames(lapply(names(synDataset), function(x) synDataset[[x]][["singleDrugResponseData"]]), NULL))
        drm.data <- subset(drm.data, Drug == drugname & Sample == samplename)
        drm.data$`Viability (Median)` <- with(drm.data, ave(Viability, Drug.Concentration, FUN = median))
        drm.data$`Viability (Median)` <- with(drm.data, ifelse(`Viability (Median)` > 1, 1, `Viability (Median)`))
        
        
        if(samplename == head(listofSamples, n = 1)){
          # Draw empty plot by type = n, adding labels and axis
          graphics::plot(drm.data$Drug.Concentration, drm.data$Viability, log = "x", type = "n", main = drugname, xlab = paste("Concentration [", unique(drm.data$Unit) , "]", sep = ""), ylab = "Viability",
                         axes=FALSE, pch=1, cex=1, lty=1, lwd=2, ylim=c(0, 1.2), yaxs="i", xaxs="r", bty="n")
          
          if(samplename == head(names(synDataset), n = 1)){graphics::axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
          if(samplename == head(names(synDataset), n = 1)){graphics::axis(2, at=seq(from = 0, to = 1, by = 0.2), labels=format(seq(from = 0, to = 1, by = 0.2), trim = TRUE))}
        }
        
        # graphics::grid(nx = NULL, ny = "", col = "lightgrey", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
        # graphics::abline(v=outer((1:10),(10^(-5:5))), col="grey80", lty = "dotted", lwd = 0.3)
        # graphics::abline(h=seq(0,1,0.1), col="grey80", lty = "dotted", lwd = 0.3)
        
        # graphics::lines(drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)], drm.data$`Viability (Median)`[order(drm.data$Drug.Concentration)], lwd = 0.6)
        graphics::points(drm.data$Drug.Concentration, drm.data$`Viability (Median)`, type="p", pch=16, cex = 0.8, col="#404040")
        
      }
      
    }
    
    cat('\r', "Saving kinetic dose response plots.", strrep(" ", 100), sep = "")
    
    
    grDevices::dev.off()
    
    cat('\r', "Finished plotting kinetic dose response based on viability, all-samples, unfitted w/ median points.", strrep(" ", 100), '\n', sep = "")
    
    
    if(drugname == tail(sort(names(efs[[1]])), n = 1)){
      message('\n', "Kinetic dose response (viability) plots (all-samples, unfitted w median points) saved to: ", '\n', file.path(.saveto, "results", "graphs/single drug response"), strrep(" ", 100), sep = "")
      rm(samplename, drugname, drm.data, op)
    }
    
  }
  
  
}