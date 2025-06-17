#' Reduce treatments to a specific matrix design 
#'
#' @description
#' An complementary component of the modular library \pkg{screenwerk}, invaluable for the initial set-up of a drug sensitivity screen.
#' #' \emph{\code{reduceDesign}} is a function that reduces the experimental design by removing combination treatments following a given pattern. 
#' This results into a reduced full matrix design, with only combination treatments being dispensed following a given matrix design.
#' One reason for utilizing a reduction in treatments is to reduce the experimental load and save on resources, or due to restrictions or limitations in experimental design.
#' 
#' @param listofCombinations an object of class 'listofCombinations'.
#' @param .design \code{character}; a predefined statement specifying to which matrix design the treatments should be reduced.
#' 
#' 
#' @details The function \code{reduceDesign} is used to reduced the experimental design of a drug screen by removing combination treatments following a given matrix pattern.
#' In that case individual concentrations of drug combinations will not be combined with each other.#' 
#' 
#' \emph{.design} specifies what drug combinations should be removed from a full list of combination treatments to a given matrix design. At current, only the \strong{'X'} matrix design is supported.

#' @return Returns an object of class S3:listofCombinations
#'
#' @examples
#' \donttest{\dontrun{
#' # Reduce treatments to a given experimental design
#' reduceDesign(listofCombinations, .design="X"))
#' }}
#'
#' @keywords drug screen reduce design matrix
#' 
#' @export

reduceDesign <- function(listofCombinations, .design="x"){
  
  # Check, if arguments are provided and in the proper format
  if(missing(listofCombinations)){stop("Data missing! Please provide a list of combinations.", call. = TRUE)}
  if(!all(sapply(c("id", "drug.number", "drug", "drug.concentration", "unit"), function(x) any(grepl(x, names(listofCombinations), ignore.case = TRUE))))){stop("in 'listofCombinations'. List of combinations has missing data! Please provide a data frame containing a column with a drug ID ['ID'] and drug number ['Drug.Number'], the drug name ['Drug'], concentration ['Drug.Concentration'] and the unit ['Unit'].", call. = TRUE)
  } else {
    colnames(listofCombinations)[grepl("id", names(listofCombinations), ignore.case = TRUE)] <- "ID"
    colnames(listofCombinations)[grepl("number", names(listofCombinations), ignore.case = TRUE)] <- "Drug.Number"
    colnames(listofCombinations)[grepl("^drug$", names(listofCombinations), ignore.case = TRUE)] <- "Drug"
    colnames(listofCombinations)[grepl("concentration", names(listofCombinations), ignore.case = TRUE)] <- "Drug.Concentration"
    colnames(listofCombinations)[grepl("unit", names(listofCombinations), ignore.case = TRUE)] <- "Unit"
  }
  
  # Check, if a specific matrix design has been provided
  if(any(missing(.design), all(!grepl(.design, c("x"), ignore.case = TRUE)))){
    stop("Matrix design missing! Please provide a valid matrix design: 'X'.", call. = TRUE)
  }else if(grepl(.design, c("x"), ignore.case = TRUE)){ 
    .design <- c("x")
  }else{ 
    .design <- .design[.design %in% c("x")]
  }
  
  
  # Filter out single drug treatments and retain only the combination matrix
  .listofCombinationDoses <- subset(listofCombinations, as.logical(ave(Drug, ID, FUN = \(x) all(table(x) == 1))))
  
  for(.id in unique(listofCombinations$ID)){
    
    cat('\r', "[", .id, "/", max(listofCombinations$ID), "] ", "Reducing design for ",  paste(listofCombinations[which(listofCombinations$ID == .id), "Drug"], paste("(", listofCombinations[which(listofCombinations$ID == .id), "Drug.Concentration"], " ", listofCombinations[which(listofCombinations$ID == .id), "Unit"], ")", sep = ""), collapse = " and "), ".", strrep(" ", 50), sep = "")
    
    # Check, if both drugs are the same, skip / do not reduce single drug treatments  
    if(any(duplicated(listofCombinations[which(listofCombinations$ID == .id), "Drug"]))){next}
    
    
    # Build a matrix of doses based on the two drugs in question
    
    # Extract drug names for a given ID
    # listofCombinations[which(listofCombinations$ID == .id), "Drug"]
    
    # Extract all drug concentrations for a given drug pair
    .doseList <- setNames(lapply(listofCombinations[which(listofCombinations$ID == .id), "Drug"], function(x) sprintf("%#.4f", sort(unique(listofCombinations[which(listofCombinations$Drug == x), "Drug.Concentration"])))[1:length(unique(listofCombinations[which(listofCombinations$Drug == x), "Drug.Concentration"]))]),
                          listofCombinations[which(listofCombinations$ID == .id), "Drug"])
    
    # Retrieve the dose range for the combination matrix based on the indices of the minimum and maximum combination dose in relation to the full range of doses
    .doseList <- lapply(.doseList, function(x) unlist(x)[
      unique(sapply(listofCombinations[which(listofCombinations$ID == .id), "Drug"], function(x)
        match(sprintf("%#.4f", min(.listofCombinationDoses[which(.listofCombinationDoses$Drug %in% x), "Drug.Concentration"])), .doseList[[x]])
      )):
        unique(sapply(listofCombinations[which(listofCombinations$ID == .id), "Drug"], function(x)
          match(sprintf("%#.4f", max(.listofCombinationDoses[which(.listofCombinationDoses$Drug %in% x), "Drug.Concentration"])), .doseList[[x]])
      ))] )
    
    # Get the number of doses for a drug pair in order to build the matrix of doses
    .noDoses <- max(sapply(.doseList, length))
    
    
    # Generate a full matrix, and
    m <- matrix(1:.noDoses^2, ncol = .noDoses, byrow = TRUE, dimnames = .doseList)
    # reduce the matrix to a X design, by setting the upper/lower triangular part of a matrix to NA
    # by generating the first diagonal pattern from top left (m1) to bottom right (m2)
    m1 <- m
    m1[lower.tri(m)] <- NA
    m1[upper.tri(m)] <- NA
    
    m2 <- apply(m, MARGIN = 2, rev)
    m2[lower.tri(m2)] <- NA
    m2[upper.tri(m2)] <- NA
    m2 <- apply(m2, MARGIN = 2, rev)
    
    
    
    # Retrieve key indices for the reduced design
    # sort(na.omit(as.vector(pmax(m1,m2,na.rm=TRUE))))
    
    # Parse key indices into coordinates
    # do.call("[", c(list(m), as.list(as.character(listofCombinations[which(listofCombinations$ID == id), "Drug.Concentration"]))))
    # m[t(as.character(listofCombinations[which(listofCombinations$ID == id), "Drug.Concentration"]))]
    
    # Extract the doses for a given drug pair, and retrieve the index from the matrix
    # by matching the doses: listofCombinations[which(listofCombinations$ID == id), "Drug.Concentration"]
    # and compare if the selected dose combination is part of the reduced design,
    # if not, remove the drug combination from the list of combinations 
    if(!(m[t(sprintf("%#.4f", listofCombinations[which(listofCombinations$ID == .id), "Drug.Concentration"]))] %in% sort(na.omit(as.vector(pmax(m1,m2,na.rm=TRUE)))))){
      listofCombinations <- subset(listofCombinations, ID != .id) 
    }
    
    # Retrieve drug and dose names from indices
    # mapply(`[[`, dimnames(m), arrayInd(36, dim(m)))
    
    rm(.id, m, m1, m2)
    
    
  }
  
  listofCombinations <- utils::type.convert(listofCombinations, as.is = TRUE)
  
  
  
  return(listofCombinations)
  
}