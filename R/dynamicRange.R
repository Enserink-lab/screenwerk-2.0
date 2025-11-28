#' Assessing the dynamic drug-activity range
#'
#' @description
#' An complementary component of the modular library \pkg{screenwerk} assessing the dynamic drug-activity range for individual drug-dose responses.
#' \emph{\code{dynamicRange}} is a function that estimates the dynamic drug-activity range (DDAR) across a number of doses for each drug response and generates a set of plots.
#' 
#' @param data an object of class 'processedData' or 'drm'.
#' @param .saveto string; path to a folder location where the results are saved to.
#' 
#' 
#' @details The function \code{dynamicRange} requires either one of the two data sets, either an object of class 'processedData' or an object of class 'drm'.
#' 
#' \code{dynamicRange} is used to provide an assessment of the drug activity range across the selected range of doses. It is an essential indicator, especially for
#' drug combination screens, in which only a selected range of doses are combined to assess synergies. If doses are combined at which drugs enfold their full inhibitory potential, it won't leave enough room for potential combinatory effects. 
#' 
#' The dynamic range is considered the dose-response range between the ED10 and ED90. This function will indicate the expected and the observed drug activity range for each drug and sample. 
#' It will generate plots based on the fitted dose-response models as well as unfitted curves.
#' 
#' Files are saved either to the specified location or the default working environment, with the corresponding folder structure: 'results/Dynamic Drug Activity Range' and 'results/graphs/dynamic range'
#' 
#' @examples
#' \donttest{\dontrun{
#' # Run dose-response model and estimate the dynamic range
#' dynamicRange(processedData, .saveto = "path/to/folder/")
#' 
#' # Estimate the dynamic range based on the provided dose-response models
#' dynamicRange(doseRespModel, .saveto = "path/to/folder/")
#' }}
#'
#' @keywords drug screen analysis dose response curve matrix
#' 
#' @importFrom utils capture.output tail
#' @importFrom stats ave predict
#' @importFrom drc drm LL.4 drmc ED
#' @importFrom grDevices dev.off png rainbow rgb
#' @importFrom graphics abline axis box grid legend lines mtext par plot plot.new points rect title
#' 
#' @export

dynamicRange <- function(data, .saveto){
  
  # Check, if the data has been provided as an object of class S3:processedData
  if(missing(data)){stop("Data missing! Please provide either a data set as an object of 'processedData' or 'drm'.", call. = TRUE)}
  
  # Extract the required data from the class S3 object
  if(class(data) == "processedData"){
    analysisData <- data[["analysisData"]]
    synDataset <- data[["splitDataset"]]
    synDRM <- data[["doserespMatrix"]]
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
  
  
  
  if(class(data) == "processedData"){
    
    # SECTION A: DOSE RESPONSE MODEL AND CURVE FITTING ##########################################################
    
    cat('\n', "Running dose-response (drm) analysis.", strrep(" ", 100), '\n\n', sep = "")
    
    
    efs <- list()
    
    for(samplename in names(synDataset)){
      
      for(drugname in sort(unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug))){
        
        nocalculations <- length(synDataset) * length(unique(synDataset[[samplename]][["singleDrugResponseData"]]$Drug))
        noutput <- 2
        
        if(exists("p")){ p <- p+noutput }else{ p <- 1 }
        
        cat('\r', "[", p, "/", nocalculations*noutput, "] ", "Fitting curve for ", samplename, " and ", drugname, ".", strrep(" ", 50), sep = "")
        
        
        if(median(subset(subset(synDataset[[samplename]][["singleDrugResponseData"]], Drug == drugname), Drug.Concentration == max(Drug.Concentration))$Viability) <= 0.8){
          
          # Run dose response model and curve fitting with R-package drc based on
          # Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response Analysis Using R. PLOS ONE 10 (12), e0146021. DOI: 10.1371/journal.pone.0146021
          # package version: drc v3.0-1
          
          # Note: viability used for viability plotting 
          # and for summary statistics (ED50, Std. Error , p-value)
          utils::capture.output(type="message",
                                efs[[samplename]][[drugname]][["drm"]] <- tryCatch({drc::drm(Viability ~ Drug.Concentration, curveid = Drug, data = synDataset[[samplename]][["singleDrugResponseData"]],
                                                                                             subset = Drug == drugname, fct = drc::LL.4(fixed = c(NA,NA,NA,NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                                                                             lowerl = c(-Inf, 0, 0, 0), upperl = c(Inf, 1, 1, Inf),
                                                                                             separate = TRUE, control = drc::drmc(useD = FALSE))}, error=function(e){
                                                                                               # If model fails with lower and upper limits, same model is run again without limits 
                                                                                               # and with derivatives for estimation, this usually works for most of the problematic dose response curves
                                                                                               drc::drm(Viability ~ Drug.Concentration, curveid = Drug, data = synDataset[[samplename]][["singleDrugResponseData"]],
                                                                                                        subset = Drug == drugname, fct = drc::LL.4(fixed = c(NA,NA,NA,NA), names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                                                                                        separate = TRUE, control = drc::drmc(useD = TRUE))
                                                                                             }) )
          
          
          # Calculate summary and EC50 based on the curve fitting data
          cat('\r', "[", p+1, "/", nocalculations*noutput, "] ", "Estimating the ED50 for ", samplename, " and ", drugname, ".", strrep(" ", 100), sep = "")
          
          efs[[samplename]][[drugname]][["summary"]] <- suppressWarnings(summary(efs[[samplename]][[drugname]][["drm"]]))
          efs[[samplename]][[drugname]][["ED"]] <- suppressWarnings(drc::ED(efs[[samplename]][[drugname]][["drm"]], c(10, 50, 90), interval = "delta", type = "relative", bound = FALSE, display = FALSE))
          
        }else{
          
          efs[[samplename]][[drugname]][["drm"]] <- lm(formula = Viability ~ Drug.Concentration, data = synDataset[[samplename]][["singleDrugResponseData"]], subset = Drug == drugname)
          
          # Calculate summary statistics for linear model
          efs[[samplename]][[drugname]][["drm"]][["origData"]] <- synDataset[[samplename]][["singleDrugResponseData"]]
          efs[[samplename]][[drugname]][["summary"]] <- summary(efs[[samplename]][[drugname]][["drm"]])
          efs[[samplename]][[drugname]][["ED"]] <- NULL
          
        }
        
        if(p == (nocalculations*noutput)-1){cat('\r', "Finished fitting dose response model for all samples.", strrep(" ", 100), '\n\n', sep = "")}
        
      }
      
    }
    
  }
  
  
  
  # Extracting Data from the Analysis. Exporting a list of ECs.
  # Extract the dynamic range for each drug and sample based on the EC10 and EC90.
  
  cat("Extracting Dynmaic Drug Activity Range (DDAR).", '\n', sep = "")
  
  
  ECDRList <- list()
  
  for(samplename in names(efs)){
    
    for(drugname in names(efs[[samplename]])){
      
      cat('\r', "> Extracting dynmaic drug-activity range (DDAR) for ", samplename, ":", drugname, ".", strrep(" ", 100), sep = "")
      
      if(all(!is.null(efs[[samplename]][[drugname]][["ED"]]))){
        
        # Extract the all ED values for each drug and sample
        ECDRList[[samplename]][[drugname]] <- cbind(ED = rownames(efs[[samplename]][[drugname]][["ED"]]), data.frame(efs[[samplename]][[drugname]][["ED"]], row.names = NULL, check.names = FALSE), stringsAsFactors = FALSE)
        ECDRList[[samplename]][[drugname]] <- ECDRList[[samplename]][[drugname]][c("ED", "Estimate")]
        # Clean-up of ED labels
        ECDRList[[samplename]][[drugname]]$ED <- with(ECDRList[[samplename]][[drugname]] , sapply(strsplit(ED, ":"), `[`, 3))
        # Add columns with cell line and drug name
        ECDRList[[samplename]][[drugname]]$Sample <- samplename
        ECDRList[[samplename]][[drugname]]$Drug <- drugname
        # Add the status and model used for curve fitting
        ECDRList[[samplename]][[drugname]]$Status <- ifelse(class(efs[[samplename]][[drugname]]$drm) == "drc", "responder", "failed")
        ECDRList[[samplename]][[drugname]]$Model <- ifelse(class(efs[[samplename]][[drugname]]$drm) == "drc", "LL.4", NA)
        
        # Convert list from long to wide-format
        ECDRList[[samplename]][[drugname]] <- stats::reshape(ECDRList[[samplename]][[drugname]], idvar = c("ED", "Estimate"), v.names = "Estimate", timevar = "ED", direction = "wide", new.row.names = NULL, sep = "")
        names(ECDRList[[samplename]][[drugname]]) <- gsub("Estimate", "ED", names(ECDRList[[samplename]][[drugname]]))
        # Combine rows dropping empty columns
        ECDRList[[samplename]][[drugname]] <- as.data.frame(lapply(ECDRList[[samplename]][[drugname]], function(x) x[!is.na(x)][1]), stringsAsFactors = FALSE)
        
      } else {
        
        ECDRList[[samplename]][[drugname]] <- data.frame(Sample = samplename, Drug = drugname, Status = "failed", Model = NA, ED10 = NA, ED50 = NA, ED90 = NA, stringsAsFactors = FALSE, check.names = FALSE)
        
      }
 
    }
    
    ECDRList[[samplename]] <- do.call(rbind, setNames(ECDRList[[samplename]], NULL)) 
    
    cat('\r', "Finished Extracting dynmaic drug-activity range (DDAR) for ", samplename, ".", strrep(" ", 100), '\n', sep = "")
    

    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "Dynamic Drug Activity Range"))){ dir.create(file.path(.saveto, "Dynamic Drug Activity Range"), showWarnings = FALSE, recursive = TRUE) }
    
    
    if(samplename == utils::tail(names(efs), n = 1)){
      ECDRList <- do.call(rbind, setNames(ECDRList, NULL))
      
      # Calculate the minimum, mean and maximum Efor each ED
      ECDRList$minED10 <- with(ECDRList, stats::ave(ED10, Drug, FUN = function(x) ifelse(all(is.na(x)), NA, min(x, na.rm = TRUE))))
      ECDRList$maxED90 <- with(ECDRList, stats::ave(ED90, Drug, FUN = function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE))))
      ECDRList$avgED10 <- with(ECDRList, stats::ave(ED10, Drug, FUN = function(x) ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))))
      ECDRList$avgED90 <- with(ECDRList, stats::ave(ED90, Drug, FUN = function(x) ifelse(all(is.na(x)), NA, mean(x, na.rm = TRUE))))
      ECDRList$medED10 <- with(ECDRList, stats::ave(ED10, Drug, FUN = function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE))))
      ECDRList$medED90 <- with(ECDRList, stats::ave(ED90, Drug, FUN = function(x) ifelse(all(is.na(x)), NA, median(x, na.rm = TRUE)))) 
      
      ECDRList$Unit <- unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit)
      
      ECDRList <- ECDRList[with(ECDRList, order(Drug)),]

      
      # Export the data set prior EC50 calculation
      write.csv2(ECDRList, file = file.path(.saveto, "Dynamic Drug Activity Range", "Dynmaic Drug Activity Range (DDAR).csv"), row.names = FALSE, quote = FALSE)

    }
    
    if(samplename == utils::tail(names(efs), n = 1L)){
      # cat('\r', "Finished extracting dynamic drug-activity range (DDAR).", strrep(" ", 100), '\n\n', sep = "")
      message('\n', "Dynamic drug-activity range (DDAR) data saved to: ", '\n', file.path(.saveto, "Dynamic Drug Activity Range"), strrep(" ", 100), sep = "")
      rm(samplename, drugname)
    }
    
  }
 
  
  
  # Plotting dose response curves with the dynamic range.
  # Plot curves as composite for each cell line with all the drugs and export as png.
  
  cat("Plotting dynamic drug-activity range (DDAR).", '\n', sep = "")
  
  for(samplename in names(efs)){
    
    cat('\r', " > Plotting dynamic range based on single drug viability curves for ", samplename, ".", strrep(" ", 100), sep = "")
    
    
    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "graphs/dynamic range", "by sample", samplename))){ dir.create(file.path(.saveto, "graphs/dynamic range", "by sample", samplename), showWarnings = FALSE, recursive = TRUE) }
    
    for(p in 1:ceiling(length(sort(names(efs[[samplename]])))/60)){
    

    grDevices::png(filename = file.path(resultDirectory, "graphs/dynamic range/by sample", samplename, paste(samplename, " dynamic range (viability) [", p, "of", ceiling(length(sort(names(efs[[samplename]])))/70), "].png", sep = "")), 
                       width = (7920/10)*10, height = (4080/6)*6, units = "px", pointsize = 24)
      op <- graphics::par(mfrow = c(6, 10), oma = c(4,0,7,0))
    
    for(drugname in na.omit(sort(names(efs[[samplename]]))[((60*p)-60+1):(60*p)])){
      
      cat('\r', " > Plotting dynamic range based on single drug viability curves for ", drugname, strrep(" ", 100), sep = "")
      
      
      drm.data <- subset(efs[[samplename]][[drugname]][["drm"]][["origData"]], Drug == drugname)
      drm.model <- efs[[samplename]][[drugname]][["drm"]]
      
      switch (class(drm.model),
              lm = {
                
                # plot all points (replicates) by type = all, or average of points by type = average
                graphics::par(bg="white", fg="black", col="black", col.axis="black", col.lab="black", col.main="black", col.sub="black") 
                graphics::plot(drm.data$Drug.Concentration, drm.data$Viability, log = "x", type = "p", main = drugname, xlab = paste("Concentration [", unique(drm.data$Unit) , "]", sep = ""), ylab = "Viability",
                               axes=FALSE, col = "grey80", pch=1, cex=1, lty=1, lwd=1, ylim=c(0, max(drm.data$Viability)+0.01), yaxs="i", xaxs="r")
                graphics::clip(min(drm.data$Drug.Concentration), max(drm.data$Drug.Concentration), min(drm.data$Inhibition), 1)
                graphics::abline(lm(formula = Viability ~ Drug.Concentration, data = drm.data), lwd=1, col="grey80")
                # lines(drm.data$Drug.Concentration[order(drm.data$Drug.Concentration)], with(drm.data, ave(Viability, Drug.Concentration, FUN = median))[order(drm.data$Drug.Concentration)], lwd=0.5)
                graphics::axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))
                graphics::axis(2)
                
                # Indicate the median points of all replicates
                graphics::points(drm.data$Drug.Concentration, with(drm.data, ave(Viability, Drug.Concentration, FUN = median)), type="p", pch=16, cex = 1, col="#cccccc")
                
                # Indicating the doses used for a given drug and the theoretical dynamic range
                for(n in 1:length(drm.data$Drug.Concentration)){
                  graphics::abline(v = with(drm.data, sort(unique(Drug.Concentration)))[n], col=rgb(0.8,0.8,0.8,alpha=1), lwd = 1.2, lty = 3)
                }
                
              },
              drc = {
                
                EC10 <- efs[[samplename]][[drugname]][["ED"]][paste("e", drugname, "10", sep = ":"), "Estimate"]
                EC50 <- efs[[samplename]][[drugname]][["ED"]][paste("e", drugname, "50", sep = ":"), "Estimate"]
                EC90 <- efs[[samplename]][[drugname]][["ED"]][paste("e", drugname, "90", sep = ":"), "Estimate"]
                
                # Plot all points (replicates) by type = all, or average of points by type = average
                graphics::par(bg="white", fg="black", col="black", col.axis="black", col.lab="black", col.main="black", col.sub="black") 
                graphics::plot(drm.model, level = drugname, log = "x", type = "average", main = drugname, 
                               xlab = paste("Concentration [", unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit), "]", sep = ""), ylab = "Viability"
                               # set an individual scale for plotting
                               , axes = FALSE, pch=1, cex=1, lty=1, lwd=1, ylim=c(0, max(drm.data$Viability)+0.01), yaxs="i", xaxs="r"
                )
                # graphics::grid(10,10, lwd = 1)
                
                if(samplename == head(names(efs), n = 1)){axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))}
                if(samplename == head(names(efs), n = 1)){axis(2)}
                
                # Indicate the measured EC50
                graphics::points(EC50, efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Upper Limit:(Intercept)", "Estimate"]-((
                  efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Upper Limit:(Intercept)", "Estimate"]-
                    efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Lower Limit:(Intercept)", "Estimate"])/2),
                  type="p", pch=16, cex = 1.0, col="red")
                
                # Indicate the observed dynamic range between the EC10 and EC90
                graphics::abline(v = EC10, col="grey", lwd = 1.2, lty = 3)
                graphics::abline(v = EC90, col="grey", lwd = 1.2, lty = 3)
                rect(EC10,par("usr")[3],EC90,par("usr")[4],col=rgb(0.89,0.89,0.89,alpha=0.3),lty=0)
                
                graphics::legend("bottom", c(paste("EC10: ", format(EC10, digits = 3), " ", unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit), "  EC50: ", format(EC50, digits = 3), " ", unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit), "  EC90: ", format(EC90, digits = 3),  " ", unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit), sep = "")),
                                 inset=c(0,0.99), xpd=TRUE, horiz=TRUE, cex = 0.9, col = "black", bty = "n")
                
                # graphics::legend("bottomleft", c(paste("EC10: ", format(EC10, digits = 3), " ", unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit),  sep = ""), 
                #                        paste("EC50: ", format(EC50, digits = 3), " ", unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit), sep = ""),
                #                        paste("EC90: ", format(EC90, digits = 3), " ", unique(efs[[samplename]][[drugname]][["drm"]][["origData"]]$Unit),  sep = "")),
                #        cex = 0.8, col = "black", bty = "n")
                
                
                # Indicating the doses used for a given drug and the theoretical dynamic range
                for(n in 1:length(unique(drm.data$Drug.Concentration))){
                  graphics::points(sort(unique(drm.data$Drug.Concentration))[n], 
                                   suppressWarnings(drc::PR(efs[[samplename]][[drugname]]$drm, sort(unique(drm.data$Drug.Concentration))[n])), 
                                   type="p", pch=16, cex = 1.2, col=rgb(0.0,0.0,0.0,alpha=1))
                  
                  graphics::abline(v = sort(unique(drm.data$Drug.Concentration))[n], col=rgb(0.0,0.0,0.0,alpha=1), lwd = 1.2, lty = 3)
                  
                }
                
                
                #   # Import a list of maximum possible doses
                #   listofDoses <- read.csv2(file=file.path(libDirectory, "listofDoses.csv"), check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, na.strings="", sep=";", dec=",", skip=0)
                #   # Section checking if the EC50 is above the highest dose
                #   if(EC50 > as.numeric(listofDoses$`6th Dose`[which(listofDoses$Drug == drugname)])){
                #     graphics::box("plot", lwd = 2.4, col="orange")}
                # 
                #   # Section checking if the EC10 is below the lowest dose
                #   if(EC10 < as.numeric(listofDoses$`1st Dose`[which(listofDoses$Drug == drugname)])){
                #     graphics::box("plot", lwd = 2.4, col="orange")}
                
                # Section checking if the Std. Error is larger than the EC50
                if(is.na(efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["ED50:(Intercept)", "Std. Error"]) &
                   EC50 > as.numeric(max(sort(unique(drm.data$Drug.Concentration))))){
                  graphics::box("plot", lwd = 2.4, col="orange")}
                
                # 
                #   # Section checking if there is less than 25% killing
                #   if(efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Upper Limit:(Intercept)", "Estimate"] -
                #      efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Lower Limit:(Intercept)", "Estimate"] < 0.2){
                #     graphics::box("plot", lwd = 2.4, col="red")}
                
                # Section checking if the slope is positive or negative
                if(efs[[samplename]][[drugname]][["summary"]][["coefficients"]]["Slope:(Intercept)", "Estimate"] < 0.0){
                  graphics::box("plot", lwd = 2.4, col="red")}
                
                # Section checking if the EC50 is outside the dose range
                if(EC50 > max(efs[[samplename]][[drugname]][["drm"]][["dataList"]][["dose"]]) | EC50 < min(efs[[samplename]][[drugname]][["drm"]][["dataList"]][["dose"]])){
                  graphics::box("plot", lwd = 2.4, col="orange")}
                
                # Section checking if response is linear
                if(efs[[samplename]][[drugname]][["drm"]][["coefficients"]]["Lower Limit:(Intercept)"] == efs[[samplename]][[drugname]][["drm"]][["coefficients"]]["Upper Limit:(Intercept)"]){
                  graphics::box("plot", lwd = 2.4, col="red")}
                
                
              }
      )
      
      
    }
      
      if(p == ceiling(length(sort(listofDrugs))/60)){
        graphics::plot.new()
        graphics::legend("topleft", c("observed dynamic range (DR)"), fill=c(rgb(0.89,0.89,0.89,alpha=0.3)), xpd = TRUE, horiz = FALSE, inset = c(0, 0.1), bty = "n", cex = 1.5)
        graphics::legend("topleft", c("EC50"), col=rgb(1,0,0,alpha=1), xpd = TRUE, horiz = FALSE, inset = c(0.0125, 0.25), pch=16, bty = "n", cex = 1.5)
      }
      
    graphics::title(paste(samplename, ": Dynamic Drug Activity Range", sep = ""), outer = TRUE, cex.main = 4)
    graphics::mtext(paste("  Package: ", "drc", " (v", packageVersion("drc"), ")   ",  sep = ""), side = 1, line = 2, adj = 1, cex = 1, col = "black", outer = TRUE) 
    grDevices::dev.off()
    
    }
    cat('\r', "Finished plotting dynamic drug-activity range (DDAR) for ", samplename, ".", strrep(" ", 100), '\n', sep = "")

    if(samplename == utils::tail(names(efs), n = 1L)){
      message('\n', "Dynamic drug-activity range (DDAR) plots saved to: ", '\n', file.path(.saveto, "graphs/dynamic range", "by sample"), strrep(" ", 100), sep = "")
      rm(samplename, drugname, EC10, EC50, EC90, n, op)
    }
    
  }
  
  
  
  
  # Plotting the dynamic range for all samples.
  # Plot curves as composite with all cell lines for each of the drugs and export as png.
  
  cat("Plotting dynamic drug-activity range (DDAR) based on single drug viability curves as composite for all samples.", '\n', sep = "")
  
  # Create subfolder if it does not exist
  if(!file.exists(file.path(.saveto, "graphs/dynamic range"))){ dir.create(file.path(.saveto, "graphs/dynamic range"), showWarnings = FALSE, recursive = TRUE) }
  
  
  for(p in 1:ceiling(length(sort(names(efs[[1]])))/60)){
    
    grDevices::png(filename = file.path(.saveto, "graphs/dynamic range", paste("dynamic range (DR) (all cell lines) [", p, "of", ceiling(length(sort(names(efs[[1]])))/70), "] graph.png", sep = "")), 
                   width = (7920/10)*10, height = (4080/6)*6, units = "px", pointsize = 24)
    op <- graphics::par(mfrow = c(6, 10), oma = c(4,0,7,0))
    
    for(drugname in na.omit(sort(names(efs[[1]]))[((60*p)-60+1):(60*p)])){
      
      cat('\r', " > Plotting dynamic range based on single drug viability curves for ", drugname, ".", strrep(" ", 100), sep = "")
      
      
      drm.data <- do.call(rbind, setNames(lapply(names(efs), function(x) efs[[x]][[drugname]][["drm"]][["origData"]]), NULL))
      drm.data <- subset(drm.data, Drug == drugname)
      
      
      # Draw empty plot by type = n, adding labels and axis
      graphics::plot(drm.data$Drug.Concentration, drm.data$Viability, log = "x", type = "n", main = drugname, xlab = paste("Concentration [", unique(drm.data$Unit) , "]", sep = ""), ylab = "Viability",
                     axes=FALSE, pch=1, cex=1, lty=1, lwd=1, 
                     # ylim=c(0, max(drm.data$Viability)+0.1), 
                     ylim=c(0, 1.1), yaxs="i", xaxs="r", bty="n")
      graphics::box(col = "black")
      
      graphics::axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))
      graphics::axis(2)
      
      for(samplename in names(efs)){
        
        drm.model <- efs[[samplename]][[drugname]][["drm"]]
        
        
        switch (class(drm.model),
                lm = {
                  
                  graphics::clip(min(drm.data$Drug.Concentration), max(drm.data$Drug.Concentration), min(drm.data$Inhibition), 1.1)
                  graphics::abline(drm.model, lwd=1, col="black")
                  
                },
                drc = {
                  
                  # Plot all points (replicates) by type = all, or average of points by type = average
                  graphics::plot(drm.model, log = "x", type = "none", main = drugname, xlab = paste("Concentration [", unique(drm.data$Unit), "]", sep = ""), ylab = "Viability", axes = FALSE,
                                 lty = 1, lwd = 1, cex = 1, pch = 16, add = TRUE,
                                 col = FALSE, legend = FALSE)

                }
        )
        
      }
      
      # Indicate the observed dynamic range between the EC10 and EC90
      # Note, if LL2.4 is used, the scale is logarithmic: use log(ED50)
      graphics::abline(v = unique(ECDRList$minED10[which(ECDRList$Drug == drugname)]), col="grey", lwd = 1.2, lty = 3)
      graphics::abline(v = unique(ECDRList$maxED90[which(ECDRList$Drug == drugname)]), col="grey", lwd = 1.2, lty = 3)
      rect(unique(ECDRList$minED10[which(ECDRList$Drug == drugname)]),
           par("usr")[3],
           unique(ECDRList$maxED90[which(ECDRList$Drug == drugname)]),
           par("usr")[4],col=rgb(0.00,0.89,0.00,alpha=0.1),lty=0)
      
      # graphics::legend("bottomleft", adj = c(0.05, 0.5), paste("DR: EC10 (median): ", format(unique(ECDRList$medED10[which(ECDRList$Drug == drugname)]), digits = 3),
      #                      " ", "EC90 (median): ", format(unique(ECDRList$medED90[which(ECDRList$Drug == drugname)]), digits = 3), sep = ""),
      #        cex = 0.8, col = "black", bty = "n")
      
      
      
      # Indicating the doses used for a given drug and the theoretical dynamic range
      for(n in 1:length(unique(drm.data$Drug.Concentration))){
        graphics::abline(v = sort(unique(drm.data$Drug.Concentration))[n], col=rgb(0.8,0.8,0.8,alpha=1), lwd = 1.2, lty = 3)
      }
      
      if(any(unlist(lapply(names(efs), function(x) class(efs[[x]][[drugname]][["drm"]]) == "drc")))){
        graphics::rect(sort(unique(drm.data$Drug.Concentration))[2],
                       par("usr")[3],
                       sort(unique(drm.data$Drug.Concentration))[length(sort(unique(drm.data$Drug.Concentration)))-1],
                       par("usr")[4],col=rgb(1.00,0.00,0.00,alpha=0.05),lty=0)
      }
      
      
      # Indicate the limit by the highest theoretical dose possible based on the stock concentration for a given drug
      graphics::abline(v = with(drm.data, max(unique(Drug.Concentration))), col=rgb(1,0,0,alpha=1), lwd = 1.2, lty = 3)
      
    }
    
    
    if(p == ceiling(length(sort(listofDrugs))/60)){
      plot.new()
      graphics::legend("topleft", c("observed dynamic range (DR)", "expected dynamic range (DR)"), fill=c(rgb(0.00,0.89,0.00,alpha=0.1), rgb(1.00,0.00,0.00,alpha=0.05)), xpd = TRUE, horiz = FALSE, inset = c(0, 0.1), bty = "n", cex = 1.5)
      graphics::legend("topleft", c("max. Dose"), col=c(rgb(1,0,0,alpha=1)), xpd = TRUE, horiz = FALSE, inset = c(0, 0.4), lwd = 1.2, lty = 3, bty = "n", cex = 1.5)
    }
    
    
    graphics::title(paste("Dynamic Drug Activity Range (all cell lines)", sep = " "), outer = TRUE, cex.main = 4)
    graphics::mtext(paste("  Package: ", "drc", " (v", packageVersion("drc"), ")   ",  sep = ""), side = 1, line = 2, adj = 1.0, cex = 1, col = "black", outer = TRUE) 
    grDevices::dev.off()
    
    cat('\r', "Finished plotting dynamic drug-activity range (DDAR) as composite.", strrep(" ", 100), '\n', sep = "")
    
    message('\n', "Dynamic drug-activity range (DDAR) plot saved to: ", '\n', file.path(.saveto, "graphs/dynamic range"), strrep(" ", 100), sep = "")
    rm(drm.data, samplename, drugname, n, op)
    
  }
  
  
  
  
  # Plotting the dynamic range for all drugs
  # Plotting individual dose response curves for all cell line and by each drug
  
  cat("Plotting dynamic drug-activity range (DDAR) based on single drug viability curves by drug.", '\n', sep = "")
  
  for(drugname in sort(names(efs[[1]]))){
    
    cat('\r', " > Plotting dynamic range based on single drug viability curves for ", drugname, ".", strrep(" ", 100), sep = "")
    
    
    # Create subfolder if it does not exist
    if(!file.exists(file.path(.saveto, "graphs/dynamic range", "by drug"))){ dir.create(file.path(.saveto, "graphs/dynamic range", "by drug"), showWarnings = FALSE, recursive = TRUE) }

    drm.data <- do.call(rbind, setNames(lapply(names(efs), function(x) efs[[x]][[drugname]][["drm"]][["origData"]]), NULL))
    drm.data <- subset(drm.data, Drug == drugname)

    grDevices::png(filename = file.path(.saveto, "graphs/dynamic range", "by drug", paste(gsub("/", " ", drugname), " dynamic range (viability)", ".png", sep = "")),
        width = 1920*0.8, height = 1080*0.8, units = "px", pointsize = 24)
    
    graphics::par(mar=c(5.1, 4.1, 3.1, 8.1), xpd=TRUE)
    
    
    # Draw empty plot by type = n, adding labels and axis
    graphics::plot(drm.data$Drug.Concentration, drm.data$Viability, log = "x", type = "n", main = drugname, xlab = paste("Concentration [", unique(drm.data$Unit) , "]", sep = ""), ylab = "Viability",
         axes=FALSE, pch=1, cex=1, lty=1, lwd=1, 
         # ylim=c(0, max(drm.data$Viability)+0.1), 
         ylim=c(0,1.1), yaxs="i", xaxs="r", bty="n")
    
    graphics::axis(1, at=unique(drm.data$Drug.Concentration), labels=format(unique(drm.data$Drug.Concentration), trim = TRUE, digits = 2, scientific = FALSE, drop0trailing = TRUE))
    graphics::axis(2)
    
    graphics::legend("right", inset = c(-0.15, 0), legend = names(efs), xpd = TRUE, 
           horiz = FALSE, col = rainbow(length(names(efs)), start = 0.6, end = 0.9), lty = 1:18, lwd = 1,5, cex = 0.8, bty = "n")
    
    graphics::clip(min(drm.data$Drug.Concentration), max(drm.data$Drug.Concentration), 0, 1)
    
    
    for(samplename in names(efs)){
      
      drm.model <- efs[[samplename]][[drugname]][["drm"]]
      
      switch (class(drm.model),
              lm = {
                
                graphics::clip(min(drm.data$Drug.Concentration), max(drm.data$Drug.Concentration), min(drm.data$Inhibition), 1.1)
                graphics::abline(drm.model, lwd=1, lty = match(samplename, names(efs)), col = rainbow(length(names(efs)), start = 0.6, end = 0.9)[match(samplename, names(efs))])
                
              },
              drc = {

                # Plot all points (replicates) by type = all, or average of points by type = average
                graphics::plot(drm.model, log = "x", type = "none", main = drugname, xlab = paste("Concentration [", unique(drm.data$Unit), "]", sep = ""), ylab = "Viability", axes=FALSE,
                               lty = match(samplename, names(efs)), lwd = 1, cex = 1, pch = 16, add = TRUE,
                               col = rainbow(length(names(efs)), start = 0.6, end = 0.9)[match(samplename, names(efs))], legend = FALSE)
              }
      )
      
    }
  
    
    graphics::legend("right", inset = c(-0.15, 0), legend = unique(names(efs)), xpd = TRUE, 
           horiz = FALSE, col = grDevices::rainbow(18, start = 0.6, end = 0.9), lty = 1:18, lwd = 1,5, cex = 0.8, bty = "n")
    
    # Indicate the observed dynamic range between the EC10 and EC90
    # Note, if LL2.4 is used, the scale is logarithmic: use log(ED50)
    graphics::abline(v = unique(ECDRList$minED10[which(ECDRList$Drug == drugname)]), col="grey", lwd = 1.2, lty = 3)
    graphics::abline(v = unique(ECDRList$maxED90[which(ECDRList$Drug == drugname)]), col="grey", lwd = 1.2, lty = 3)
    graphics::rect(unique(ECDRList$minED10[which(ECDRList$Drug == drugname)]),
         par("usr")[3],
         unique(ECDRList$maxED90[which(ECDRList$Drug == drugname)]),
         par("usr")[4],col=rgb(0.00,0.89,0.00,alpha=0.1),lty=0)
    
    
    
    # Indicating the doses used for a given drug and the theoretical dynamic range
    for(n in 1:length(unique(drm.data$Drug.Concentration))){
      graphics::abline(v = sort(unique(drm.data$Drug.Concentration))[n], col=rgb(0.8,0.8,0.8,alpha=1), lwd = 1.2, lty = 3)
    }
    
    if(any(unlist(lapply(names(efs), function(x) class(efs[[x]][[drugname]][["drm"]]) == "drc")))){
      graphics::rect(with(drm.data, sort(unique(Drug.Concentration)))[2],
           par("usr")[3],
           with(drm.data, sort(unique(Drug.Concentration)))[with(drm.data, length(unique(Drug.Concentration)))-1],
           par("usr")[4],col=rgb(1.00,0.00,0.00,alpha=0.05),lty=0)
    }

    
    grDevices::dev.off()
    
    cat('\r', "Finished plotting dynamic drug-activity range (DDAR) for ", drugname, ".", strrep(" ", 100), sep = "")
    
    
    if(drugname == utils::tail(sort(names(efs[[1]])), n = 1)){
      cat('\r', "Finished plotting dynamic drug-activity range (DDAR) by drug.", strrep(" ", 100), '\n', sep = "")
      message('\n', "Dynamic drug-activity range (DDAR) plots saved to: ", '\n', file.path(.saveto, "graphs/dynamic range", "by drug"), strrep(" ", 100), sep = "")
      rm(drm.data, drm.model, samplename, drugname, n)
    }
  }
  
  
   
  # Return object of class DRM
  class(ECDRList) <- "ECDR"
  return(ECDRList)
  
}