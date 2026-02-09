#' Process consolidated data for downstream analysis 
#'
#' @description
#' An essential component of the modular library \pkg{screenwerk}, and imperative for the downstream analysis of the experimental data from a drug sensitivity screen.
#' \emph{\code{processData}} is a function that normalizes the raw measurements to the positive and negative control, and subsequently splits the data into individual data sets, 
#' one for the controls, single drug treatments and combination treatments, in case of a drug combination screen. Furthermore, a data table and a matrix is assembled for each drug pair, 
#' representing the dose response between two drugs at each dose.
#' 
#' @param consolidatedData an object of class 'consolidatedData'.
#' @param .ctrls \code{list}; a list containing the name of the positive and negative control.
#' 
#' 
#' @details The function \code{processData} is used to prepare the data for various tasks in the downstream analysis of the drug sensitivity screen. It carries out three essential steps in data processing: 
#' normalizing, splitting and assembling (re-formatting) the data.
#' 
#' \emph{.ctrls} is essential for the normalization of the raw measurements to the positive and negative control. The positive control has to be specified, while the negative control will be retrieved from the solvent in which the drug has been dissolved.
#' This requires an individual set of controls for each solvent and combination of solvents. Since drug treatments, such as single drug treatments, are backfilled with one unit of the control corresponding to the solvent in which the drug has been dissolved,
#' each well has the highest concentration of solvent used across all the treatments in the drug screen.

#' @return Returns an object of class S3:processedData
#'
#' @examples
#' \donttest{\dontrun{
#' # Process consolidated data
#' processData(consolidatedData, .ctrls=list(positive="BzCl"))
#' }}
#'
#' @keywords drug screen process normalize split assemble downstream analysis
#' 
#' @importFrom stats median
#' 
#' @export

processData <- function(consolidatedData, .ctrls=list(positive, negative)){
  
  
  # Check, if the data has been provided as an object of class S3:consolidatedData
  if(missing(consolidatedData)){stop("Data missing! Please provide a consolidated data set.", call. = TRUE)}
  if(class(consolidatedData) != "consolidatedData"){
    stop("Provided data not of class 'consolidatedData'!", call. = TRUE)
  } else {
    # Extract all the necessary data from the class S3 object
    analysisData <- consolidatedData[["consolidated"]]
    ctrlList <- consolidatedData[["dispensingData"]][["dataList"]][["ctrlList"]]
  }
  
  
  # Check, if a set of positive and negative controls have been provided
  if(missing(.ctrls)){ 
    stop("in '.ctrls'. Please provide a set of controls to normalize the data to.", call. = TRUE)
  } else if(class(.ctrls) != "list") { stop("in '.ctrls'. Argument needs to be a list! Please provide a list with a name of the positive (and negative) control.", call. = TRUE) 
  } else if(!exists("positive", where=.ctrls)){
    stop("in '.ctrls'. Positive control missing! Please provide a positive control to normalize data to.", call. = TRUE)
  } else {
    
    if(exists("positive", where=.ctrls) & !.ctrls$positive %in% ctrlList$NAME){
    stop("in '.ctrls'. Positive control not found in the list of controls! Please provide a positive control to normalize data to.", call. = TRUE)
    }
    
    if(exists("negative", where=.ctrls)){
    .controls = list(negative = c(unique(sapply(split(subset(analysisData, Drug %in% setdiff(consolidatedData[["dispensingData"]][["dataList"]][["ctrlList"]]$NAME, .ctrls$positive)), 
                                                      subset(analysisData, Drug %in% setdiff(consolidatedData[["dispensingData"]][["dataList"]][["ctrlList"]]$NAME, .ctrls$positive))$Combination.ID), function(x) paste(unique(x$Drug), collapse = "+"), simplify = TRUE, USE.NAMES = FALSE)),
                                  consolidatedData[["dispensingData"]][["origData"]][[".addUntreated"]]$name))
    
    message("Negative control provided with ", gsub(",([^,]*)$"," and\\1", paste(.ctrls$negative, collapse = ", ")), ".")
      
    # Check, if any of the provided negative controls are not found in the data set
      if(any(!.ctrls$negative %in% .controls$negative)){
        message("Negative control(s) ", gsub(",([^,]*)$"," and\\1", paste(setdiff(.ctrls$negative, .controls$negative), collapse = ", ")), " not found in the retrieved list of negative controls in the data set.")
        message("Negative control(s) ", setdiff(gsub(",([^,]*)$"," and\\1", paste(setdiff(.ctrls$negative, .controls$negative), collapse = ", ")), .controls$negative), " will be ignored!")
      }
      # Check, if the provided negative controls are missing some of the negative controls that are found in the data set
      if(!length(setdiff(.controls$negative, .ctrls$negative)) == 0){
        message("In addition to the provided negative controls, the additional control(s) ", gsub(",([^,]*)$"," and\\1", paste(setdiff(.controls$negative, .ctrls$negative), collapse = ", ")), " were found in the data set.")
      }
    }
    
    .ctrls$positive = .ctrls$positive
    .ctrls$negative = c(unique(sapply(split(subset(analysisData, Drug %in% setdiff(consolidatedData[["dispensingData"]][["dataList"]][["ctrlList"]]$NAME, .ctrls$positive)), 
                                            subset(analysisData, Drug %in% setdiff(consolidatedData[["dispensingData"]][["dataList"]][["ctrlList"]]$NAME, .ctrls$positive))$Combination.ID), function(x) paste(unique(x$Drug), collapse = "+"), simplify = TRUE, USE.NAMES = FALSE)),
                        consolidatedData[["dispensingData"]][["origData"]][[".addUntreated"]]$name)
  }
  
  
  
  cat('\n', "Positive control specified as '", .ctrls$positive, "'.", '\n', sep = "")
  cat("Negative controls identified as '", gsub(",([^,]*)$"," and\\1", paste(.ctrls$negative, collapse = ", ")), "'.", '\n', sep = "")
  cat('\n', sep = "")
  cat('\r', "> Normalizing data... This step may take some time.", strrep(" ", 100), sep = "")
  
  # Select only essential columns
  analysisData <- analysisData[,c("Sample", "Destination.Plate.Barcode", "Dispensing.Set", "Plate.Number", "Combination.ID", "Drug", "Drug.Concentration", "Unit", "Destination.Well", "CPS")]
  # Set factors for individual columns to allow proper sorting
  analysisData$Plate.Number <- factor(as.numeric(analysisData$Plate.Number), levels = order(unique(analysisData$Plate.Number)))
  
  # Add solvent to data set
  listofDrugs <- consolidatedData[["dispensingData"]][["origData"]][["listofDrugs"]]
  analysisData <- as.data.frame(append(analysisData, list(Solvent = listofDrugs$SOLVENT[match(analysisData$Drug, listofDrugs$NAME)]), after = grep("^Drug$", colnames(analysisData))))
  analysisData <- transform(analysisData, Solvent = ifelse(Drug %in% ctrlList$NAME, Drug, Solvent))
  analysisData <- transform(analysisData, Solvent = ifelse(Drug == "Untreated", "Untreated", Solvent))
  
  
  # LEGACY: Removal of backfillings at which only DMSO at half the maximum dose have been dispensed
  # Extract the list of controls and remove treatments, in which one of the controls has been added to a single drug treatment
  # for(ctrl in ctrlList$NAME){
  #   analysisData <- subset(analysisData, (Drug != ctrl | Drug == ctrl & Drug.Concentration == ctrlList$DOSE[match(ctrl, ctrlList$NAME)]))
  # }
  
  # Generate a list of indices selecting only rows in which drug treatments are not controls and controls that have been dispensed as such
  # removing treatments, in which one of the controls has been added to a single drug treatment
  # analysisData <- analysisData[Reduce(intersect, sapply(ctrlList$NAME, function(x) which(analysisData$Drug != x | analysisData$Drug == x & analysisData$Drug.Concentration == ctrlList$DOSE[match(x, ctrlList$NAME)]))),]
  
  
  # Removal of backfillings for single drug treatments with different solvents
  # Pairwise comparison at which a drug with a solvent has been dispensed and removing the solvent
  analysisData <- do.call(rbind, setNames(lapply(split(analysisData, list(analysisData$Combination.ID, analysisData$Dispensing.Set)),
                        function(x) if(nrow(x) >= 2 & any(x$Drug %in% ctrlList$NAME) & !all(x$Drug %in% ctrlList$NAME)){
                          subset(x, !(Drug %in% ctrlList$NAME))}else{x}), NULL))
  

    
  # SECTION A: NORMALIZE DATA FOR DOWSTREAM ANALYSIS  #############################################
  # Normalize data to individual controls and calculate response values
  analysisData <- transform(analysisData, NORM = CPS, Viability = NA, `Viability (%)` = NA, Inhibition = NA, `Inhibition (%)` = NA)
  
  
  # Group by plate to normalize the controls to each individual plate
  analysisData <- split(analysisData, analysisData$Destination.Plate.Barcode)
  
  # Calculate the median of the positive and negative controls for each plate, respectively
  # analysisData <- lapply(analysisData, function(x) transform(x, CPS = ifelse(Drug == .ctrls$positive, median(CPS[Drug == .ctrls$positive]), CPS)))
  # analysisData <- lapply(analysisData, function(x) transform(x, CPS = ifelse(Drug == .ctrls$negative, median(CPS[Drug == .ctrls$negative]), CPS)))
  
  
  # Split data set further by separating controls from the drug treatment
  analysisData <- lapply(analysisData, function(x) split(x, x$Drug %in% c(ctrlList$NAME, "Untreated")))
  
  
  # Calculate the median for each negative control, both controls that represent single drug treatments and combinations
  # NB: Controls are grouped single- or pair-wise based on the treatment: unlist(sapply(split(x[["TRUE"]], x[["TRUE"]]$Combination.ID), function(y) cbind(y, grp = paste(y$Drug, collapse=""))$grp))
  analysisData <- lapply(analysisData, function(x){ x[["TRUE"]] <- data.frame(do.call(rbind, lapply(split(x[["TRUE"]], unlist(sapply(split(x[["TRUE"]], x[["TRUE"]]$Combination.ID), function(y) cbind(y, grp = paste(y$Drug, collapse=""))$grp))), 
                                                                                                    function(a) do.call(rbind, lapply(split(a, a$Drug.Concentration), function(b) transform(b, NORM = median(b$CPS)))))), row.names=NULL, check.names=FALSE); x })
  
  
  # Merge data sets with controls and drug treatments together again
  analysisData <- lapply(analysisData, function(x) do.call(rbind, setNames(x, NULL)))

  
  # Normalize the drug treatments to the positive control by
  # subtracting the background level from the measurements
  # analysisData <- lapply(analysisData, function(x) transform(x, CPS = CPS - unique(CPS[Drug == .ctrls$positive])))
  analysisData <- lapply(analysisData, function(x) transform(x, NORM = NORM - unique(NORM[Drug == .ctrls$positive])))
  
  
  
  # :Experimental code:
  # Extract combinations
  # a = analysisData[[1]]
  # a[ave(a$Combination.ID,a$Combination.ID,FUN=length)>1,]
  # match(paste(a$Drug, a$Drug.Concentration, sep = ""), unique(paste(a$Drug, a$Drug.Concentration, sep = "")))
  # Extract pairwise names
  # unlist(sapply(split(a, a$Combination.ID), \(x) cbind(x, paste(x$Drug, collapse=""))$Drug))
  # :\End: 

  # Set a lower limit by replacing negative values with zero and
  # optionally, set an upper limit by replacing values larger than 1 with 1
  # analysisData <- lapply(analysisData, function(x) transform(x, CPS = ifelse(CPS < 0, 0, CPS)))
  # analysisData <- lapply(analysisData, function(x) transform(x, CPS = ifelse(CPS > 1, 1, CPS)))
  
  # Legacy code for normalization, in which the negative control is predefined by the user
  # analysisData <- lapply(analysisData, function(x) transform(x, Viability = CPS / unique(CPS[Drug == .ctrls$negative]), check.names=FALSE))
  
  
  # Normalize the drug treatments to the negative control by
  # dividing the measurements by the reference (negative) control
  # to calculate response as viability and inhibition, and percentages (%)
  analysisData <- lapply(analysisData, function(x){
    
  # Create a list of solvents and drugs based on the combination id for each normalization set
  .solvents <- sapply(split(x, x$Combination.ID), \(x) cbind(x, paste(x$Drug, collapse=""))$Solvent)
  .drugs <- sapply(split(x, x$Combination.ID), \(x) cbind(x, paste(x$Drug, collapse=""))$Drug)

  # Legacy code for normalization, in which all treatments are normalized to the highest concentration of the corresponding control
  # for(y in unique(x$Combination.ID)){
    # Retrieve CPS (norm.) from the corresponding control based on the solvent of the drug treatment
    # and with multiple controls of the same solvent use the highest concentration (double-treatment, backfilled single drug treatments)
    # x[which(x$Combination.ID == y), "Viability"] = x[which(x$Combination.ID == y), "NORM"] / unique(x$NORM[which(x$Combination.ID %in% names(which(sapply(.drugs, setequal, unique(.solvents[[as.character(y)]])))) & 
    #                                                                                                               x$Drug.Concentration == max(x$Drug.Concentration[x$Combination.ID %in% names(which(sapply(.drugs, setequal, unique(.solvents[[as.character(y)]]))))]) )])
  # };x })
  
  for(y in unique(x$Combination.ID)){
    
    cat('\r', "> ", "[", format(sum(sapply(analysisData[1:match(unique(x$Destination.Plate.Barcode), names(analysisData))], function(x) length(unique(x$Combination.ID)))) - sum(sapply(analysisData[match(unique(x$Destination.Plate.Barcode), names(analysisData))], function(x) length(unique(x$Combination.ID)))) + match(y, unique(x$Combination.ID)), nsmall=1, big.mark=","), "/", format(sum(sapply(analysisData, function(x) length(unique(x$Combination.ID)))), nsmall=1, big.mark=","), "]", 
        " Normalizing data set: ", unique(x$Sample), "  Barcode: ",  unique(x$Destination.Plate.Barcode), "  Plate: ", unique(x$Plate.Number), "  ID: ", y, strrep(" ", 100), sep = "")
    
    # Retrieve CPS (norm.) from the corresponding control based on the solvent of the drug treatment
    # and with multiple controls of the same solvent use the highest concentration (combination-treatment, backfilled single drug treatments), or
    # half the highest concentration (single-treatment, non-backfilled single drug treatments).
    x[which(x$Combination.ID == y), "Viability"] = ifelse(x[which(x$Combination.ID == y), "Drug"] %in% ctrlList$NAME, 
                                                          x[which(x$Combination.ID == y), "NORM"] / x[which(x$Combination.ID == y), "NORM"],
                                                          ifelse(consolidatedData[["dispensingData"]][["origData"]][[".backfilling"]] == TRUE,
                                                                 # normalize to the highes concentration of a control
                                                                 x[which(x$Combination.ID == y), "NORM"] / unique(x$NORM[which(x$Combination.ID %in% names(which(sapply(.drugs, setequal, unique(.solvents[[as.character(y)]])))) & 
                                                                                                                                 x$Drug.Concentration == max(x$Drug.Concentration[x$Combination.ID %in% names(which(sapply(.drugs, setequal, unique(.solvents[[as.character(y)]]))))]) )]),
                                                                 
                                                                 # normalize to half the concentration for single drug treatments and the highest concentration for combination treatments
                                                                 ifelse(nrow(x[which(x$Combination.ID == y),]) > 1, 
                                                                        x[which(x$Combination.ID == y), "NORM"] / unique(x$NORM[which(x$Combination.ID %in% names(which(sapply(.drugs, setequal, unique(.solvents[[as.character(y)]])))) & 
                                                                                                                                        x$Drug.Concentration == max(x$Drug.Concentration[x$Combination.ID %in% names(which(sapply(.drugs, setequal, unique(.solvents[[as.character(y)]]))))]) )]),
                                                                        
                                                                        x[which(x$Combination.ID == y), "NORM"] / unique(x$NORM[which(x$Combination.ID %in% names(which(sapply(.drugs, setequal, unique(.solvents[[as.character(y)]])))) & 
                                                                                                                                        x$Drug.Concentration == min(x$Drug.Concentration[x$Combination.ID %in% names(which(sapply(.drugs, setequal, unique(.solvents[[as.character(y)]]))))]) )])
                                                                 )))
  };x })
  
  
  # Normalize the positive control manually
  analysisData <- lapply(analysisData, function(x) transform(x, Viability = ifelse(Drug == .ctrls$positive, 0, Viability), check.names=FALSE))
  
  # Calculate inhibition and percentages (%)
  analysisData <- lapply(analysisData, function(x) transform(x, 'Viability (%)' = Viability * 100, check.names=FALSE))
  analysisData <- lapply(analysisData, function(x) transform(x, Inhibition = 1 - Viability, check.names=FALSE))
  analysisData <- lapply(analysisData, function(x) transform(x, 'Inhibition (%)' = 100 - `Viability (%)`, check.names=FALSE))
  
  
  # Ungroup data by merging each individual plate back together
  analysisData <- do.call(rbind, setNames(analysisData, NULL))
  analysisData <- analysisData[with(analysisData, order(Dispensing.Set, Plate.Number, Destination.Well)),]
  
  cat('\r', "Finished normalizing data.", strrep(" ", 100), "\n\n", sep = "")
  
  
  
  
  # SECTION B: SPLIT DATA INTO INDIVIDUAL DATA SETS FOR DOWSTREAM PROCESSING  ##############################
  # Select only essential data and reformat it to individual data sets
  
  # Create a list of data frames containing individual data sets needed for downstream processing
  # The analysis data is retained and split into a control, single drug response and combination response data set
  
  # Split analysis data by sample and 
  # generate a complete data set for analysis for each sample
  synDataset <- lapply(split(analysisData, analysisData$Sample), list)
  synDataset <- lapply(synDataset, setNames, "analysisData")
  synDataset <- synDataset[unique(analysisData$Sample)]
  
  for (samplename in names(synDataset)){
  
  cat(" > Splitting data set into controls for ", samplename, ".", strrep(" ", 100), sep = "")
  
    
    # Split data set containing only the control treatments ###############################
    synDataset[[samplename]][["controlData"]] <- synDataset[[samplename]][["analysisData"]]
    
    # Select only essential columns
    synDataset[[samplename]][["controlData"]] <- synDataset[[samplename]][["controlData"]][,c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID", "Drug", "Drug.Concentration", "Unit", "Destination.Well", 
                                    "Viability", "Viability (%)", "Inhibition", "Inhibition (%)")]
    
    
    .controls <- c(consolidatedData[["dispensingData"]][["dataList"]][["ctrlList"]]$NAME,
                   consolidatedData[["dispensingData"]][["origData"]][[".addUntreated"]]$name)
    
    # Filter treatments by control
    synDataset[[samplename]][["controlData"]] <- subset(synDataset[[samplename]][["controlData"]], Drug %in% .controls)
    
    # Remove duplicated control treatments to only one replicate
    # synDataset[[samplename]][["controlData"]] <- synDataset[[samplename]][["controlData"]][!duplicated(synDataset[[samplename]][["controlData"]][,c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Drug")]),]
    
    # Select columns and sort rows 
    synDataset[[samplename]][["controlData"]] <- synDataset[[samplename]][["controlData"]][,c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID", "Drug", "Drug.Concentration", "Unit", "Destination.Well", "Viability", "Viability (%)", "Inhibition", "Inhibition (%)")]
    synDataset[[samplename]][["controlData"]] <- synDataset[[samplename]][["controlData"]][with(synDataset[[samplename]][["controlData"]], order(Sample, Destination.Plate.Barcode, Plate.Number, Combination.ID,Drug)),]
    
    
    
  cat('\r', " > Splitting data set into single drug treatments for ", samplename, ".", strrep(" ", 100), sep = "")
    
    
    # Split data set containing only the single drug treatments ######################################
    synDataset[[samplename]][["singleDrugResponseData"]] <- synDataset[[samplename]][["analysisData"]]
    
    # Select only essential columns
    synDataset[[samplename]][["singleDrugResponseData"]] <- synDataset[[samplename]][["singleDrugResponseData"]][,c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID", 
                                                                                                                    "Drug", "Drug.Concentration", "Unit", "Destination.Well", 
                                                                                                                    "Viability", "Viability (%)", "Inhibition", "Inhibition (%)")]
    
    
    # Remove control treatments from actual drug treatments
    synDataset[[samplename]][["singleDrugResponseData"]] <- subset(synDataset[[samplename]][["singleDrugResponseData"]], !Drug %in% .controls)
    # and select only for unique single drug treatments by filtering out drug combinations
    synDataset[[samplename]][["singleDrugResponseData"]] <- synDataset[[samplename]][["singleDrugResponseData"]][!(duplicated(synDataset[[samplename]][["singleDrugResponseData"]][c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID", "Destination.Well")]) | 
                                                                                                                   duplicated(synDataset[[samplename]][["singleDrugResponseData"]][c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID", "Destination.Well")], fromLast = TRUE)),]
    
    
    # Sort final data set
    synDataset[[samplename]][["singleDrugResponseData"]] <- synDataset[[samplename]][["singleDrugResponseData"]][with(synDataset[[samplename]][["singleDrugResponseData"]], order(Sample, Drug, Drug.Concentration)),]
    
    
    
  cat('\r', " > Splitting data set into combination drug treatments for ", samplename, ".", strrep(" ", 100), sep = "")
    
    
    # Split data set containing only the combination drug treatments ###########################
    synDataset[[samplename]][["combinationData"]] <- synDataset[[samplename]][["analysisData"]]
  
    # Select only essential columns
    synDataset[[samplename]][["combinationData"]] <- synDataset[[samplename]][["combinationData"]][,c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID", 
                                                                                                                  "Drug", "Drug.Concentration", "Unit", "Destination.Well", 
                                                                                                                  "Viability", "Viability (%)", "Inhibition", "Inhibition (%)")]


    # Remove control treatments from actual drug treatments
    synDataset[[samplename]][["combinationData"]] <- subset(synDataset[[samplename]][["combinationData"]], !Drug %in% .controls)
    # and select only for drug combination treatments by filtering out single drug treatments
    synDataset[[samplename]][["combinationData"]] <- synDataset[[samplename]][["combinationData"]][duplicated(synDataset[[samplename]][["combinationData"]][c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID", "Destination.Well")]) | 
                                                                                                   duplicated(synDataset[[samplename]][["combinationData"]][c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID", "Destination.Well")], fromLast = TRUE),]
    
    
    # Sort final data set
    synDataset[[samplename]][["combinationData"]] <- synDataset[[samplename]][["combinationData"]][with(synDataset[[samplename]][["combinationData"]], order(Sample, Combination.ID)),]
    

    
  cat('\r', "Finished splitting data set for ", samplename, strrep(" ", 100), "\n", sep = "")
  
  }
  
  
  
  
  cat('\n')
  
  # SECTION C: RE-FORMAT INDIVIDUAL DATA SETS FOR DOWSTREAM PROCESSING  ##############################
  # Reformat the data set into a data frame and a matrix
  
  # Create a list of matrices with the dose response data
  # A dose response matrix is assembled from the combination and single drug response data
  
  # Split combination treatments by individual drug
  synDRM <- lapply(synDataset, function(x) split(x[["combinationData"]], x[["combinationData"]]$Drug))
  # synDRM <- setNames(lapply(vector("list", length(synDataset)), function(x) list()), names(synDataset))

  
  
  for (samplename in names(synDRM)){
    
    for(drugname in names(synDRM[[samplename]])){
     
      # Add to each primary drug the corresponding drug pair
      synDRM[[samplename]][[drugname]] <- merge(synDRM[[samplename]][[drugname]], synDataset[[samplename]][["combinationData"]], 
                                                by = c("Sample", "Destination.Plate.Barcode", "Plate.Number", "Combination.ID",
                                                       "Destination.Well", "Inhibition"), suffixes = c(".1", ".2"))
      
      # Remove same drug combinations
      synDRM[[samplename]][[drugname]] <- subset(synDRM[[samplename]][[drugname]], Drug.1 != Drug.2)
      # Select only essential columns and rename columns by name
      synDRM[[samplename]][[drugname]] <- synDRM[[samplename]][[drugname]][c("Inhibition", grep("Drug|Unit", names(synDRM[[samplename]][[drugname]]), ignore.case = TRUE, value = TRUE))]
      names(synDRM[[samplename]][[drugname]])[match("Inhibition", names(synDRM[[samplename]][[drugname]]))] <- "Response"
      names(synDRM[[samplename]][[drugname]])[match("Drug.2", names(synDRM[[samplename]][[drugname]]))] <- "DrugRow"
      names(synDRM[[samplename]][[drugname]])[match("Drug.1", names(synDRM[[samplename]][[drugname]]))] <- "DrugCol"
      names(synDRM[[samplename]][[drugname]])[match("Drug.Concentration.2", names(synDRM[[samplename]][[drugname]]))] <- "ConcRow"
      names(synDRM[[samplename]][[drugname]])[match("Drug.Concentration.1", names(synDRM[[samplename]][[drugname]]))] <- "ConcCol"
      names(synDRM[[samplename]][[drugname]])[match("Unit.2", names(synDRM[[samplename]][[drugname]]))] <- "ConcRowUnit"
      names(synDRM[[samplename]][[drugname]])[match("Unit.1", names(synDRM[[samplename]][[drugname]]))] <- "ConcColUnit"
      
      
      # Split the data set for each primary drug into individual data sets for each drug pair
      synDRM[[samplename]][[drugname]] <- split(synDRM[[samplename]][[drugname]] , synDRM[[samplename]][[drugname]]$DrugRow)
   
      
      for(drugpair in names(synDRM[[samplename]][[drugname]])){
        
        cat('\r', " > Generating dose response matrix for ", samplename, ": ", drugname, " + ", drugpair, strrep(" ", 50), sep = "")
        
        
        # Assemble data frame to create a dose response matrix for each drug combination
        # Note: Column: primary drug, Row: drug pair
        
        
        # Create a data frame with the single drug response for the primary drug (DrugCol)
        .sdrCol <- subset(synDataset[[samplename]][["singleDrugResponseData"]], Drug == drugname)
        # Add additional columns to build a data frame that matches the existing data frame
        .sdrCol <- transform(.sdrCol, Response = Inhibition, DrugRow = drugpair, ConcRow = 0, ConcRowUnit = Unit, stringsAsFactors = FALSE)
        # Rename individual columns
        names(.sdrCol)[names(.sdrCol) %in%  c("Drug", "Drug.Concentration", "Unit")] <- c("DrugCol", "ConcCol", "ConcColUnit")
        # Calculate the mean response for each concentration
        .sdrCol <- merge(stats::aggregate(Response ~ ConcCol, .sdrCol, median), .sdrCol, by = c("Response", "ConcCol"), all = FALSE)
        .sdrCol <- .sdrCol[c("Response", "DrugRow", "DrugCol", "ConcRow", "ConcCol", "ConcRowUnit", "ConcColUnit")]
        # Remove duplicate rows
        .sdrCol <- .sdrCol[!duplicated(.sdrCol),]
        
        # Create a data frame with the single drug response for the drug pair (DrugRow)
        .sdrRow <- subset(synDataset[[samplename]][["singleDrugResponseData"]], Drug == drugpair)
        # Add additional columns to build a data frame that matches the existing data frame
        .sdrRow <- transform(.sdrRow, Response = Inhibition, DrugCol = drugname, ConcCol = 0, ConcColUnit = Unit, stringsAsFactors = FALSE)
        # Rename individual columns
        names(.sdrRow)[names(.sdrRow) %in%  c("Drug", "Drug.Concentration", "Unit")] <- c("DrugRow", "ConcRow", "ConcRowUnit")
        # Calculate the mean response for each concentration
        .sdrRow <- merge(stats::aggregate(Response ~ ConcRow, .sdrRow, median), .sdrRow, by = c("Response", "ConcRow"), all = FALSE)
        .sdrRow <- .sdrRow[c("Response", "DrugRow", "DrugCol", "ConcRow", "ConcCol", "ConcRowUnit", "ConcColUnit")]
        # Remove duplicate rows
        .sdrRow <- .sdrRow[!duplicated(.sdrRow),]

        
        # Adding single drug response data for both the primary drug and drug pair to the combination treatment
        # and add a dose zero for both drugs
        synDRM[[samplename]][[drugname]][[drugpair]] <- rbind(synDRM[[samplename]][[drugname]][[drugpair]], .sdrCol, .sdrRow, 
                                                              data.frame(Response = 0, DrugRow = drugpair, DrugCol = drugname, ConcRow = 0, ConcCol = 0, 
                                                                         ConcRowUnit = unique(synDRM[[samplename]][[drugname]][[drugpair]]$ConcRowUnit), 
                                                                         ConcColUnit = unique(synDRM[[samplename]][[drugname]][[drugpair]]$ConcColUnit)))

        
        # Set the order of columns and rows
        synDRM[[samplename]][[drugname]][[drugpair]] <- synDRM[[samplename]][[drugname]][[drugpair]][c("Response", "DrugRow", "DrugCol", "ConcRow", "ConcCol", "ConcRowUnit", "ConcColUnit")]
        synDRM[[samplename]][[drugname]][[drugpair]] <- synDRM[[samplename]][[drugname]][[drugpair]][with(synDRM[[samplename]][[drugname]][[drugpair]], order(ConcCol, ConcRow)),]
        rownames(synDRM[[samplename]][[drugname]][[drugpair]]) <- NULL



        
        
        # Create an empty data frame with all possible drug combinations and fill it with the data set from above
        # Note: missing combinations are retained with NAs
        # Create an empty data set with all possible drug dose combinations
        .dfs <- data.frame(Response = as.numeric(NA), DrugRow = drugpair, DrugCol = drugname,
                           expand.grid(ConcRow = c(sort(unique(synDRM[[samplename]][[drugname]][[drugpair]]$ConcRow))),
                                       ConcCol = c(sort(unique(synDRM[[samplename]][[drugname]][[drugpair]]$ConcCol))), stringsAsFactors = FALSE),
                           ConcRowUnit = unique(synDRM[[samplename]][[drugname]][[drugpair]]$ConcRowUnit), 
                           ConcColUnit = unique(synDRM[[samplename]][[drugname]][[drugpair]]$ConcColUnit), 
                           check.names = FALSE, stringsAsFactors = FALSE)
        

        # Fill in response by mapping the measurements between the data frames
        # Populate empty set with the response for each combination
        .dfs <- transform(.dfs, Response = synDRM[[samplename]][[drugname]][[drugpair]]$Response[match(paste(.dfs$DrugRow, .dfs$DrugCol, .dfs$ConcRow, .dfs$ConcCol),
                                                                                                    paste(synDRM[[samplename]][[drugname]][[drugpair]]$DrugRow,
                                                                                                          synDRM[[samplename]][[drugname]][[drugpair]]$DrugCol,
                                                                                                          synDRM[[samplename]][[drugname]][[drugpair]]$ConcRow,
                                                                                                          synDRM[[samplename]][[drugname]][[drugpair]]$ConcCol))])
        
        # Add indices for each of the doses, designating a column and row position in the matrix 
        .dfs$Row <- with(.dfs, ave(seq_along(ConcCol), ConcCol, FUN=seq_along))
        .dfs$Column <- with(.dfs, ave(seq_along(ConcRow), ConcRow, FUN=seq_along))

        # Set the order of rows
        .dfs <- .dfs[with(.dfs, order(ConcCol, ConcRow)),]
        
        
        # Assemble a matrix from the data frame
        .drm <- matrix(.dfs$Response, nrow = length(unique(.dfs$ConcRow)), ncol = length(unique(.dfs$ConcCol)), 
                       dimnames = list(c(sort(unique(.dfs$ConcRow))), c(sort(unique(.dfs$ConcCol)))))
        
        
        # Assign the matrix and data frame to the initial list
        synDRM[[samplename]][[drugname]][[drugpair]] <- list(dfs = .dfs, drm = .drm)

      }
      
    }
    
    cat('\r', "Finished assembling data set for ", samplename, ".", strrep(" ", 100), "\n", sep = "")
    
  }
  
  
  cat('\n', "Successfully normalized, split and assembled data set.", strrep(" ", 100), "\n", sep = "")
  
  
  # Create an object of S3 class with a list from the attributes of all three vectors
  object <- list(analysisData = analysisData, splitDataset = synDataset, doserespMatrix = synDRM, .ctrls=list(positive=.ctrls$positive, negative=.ctrls$negative))
  class(object) <- "processedData"
  
  return(object)
  
  
}
