#' Combine drugs for a drug combination screen
#'
#' @description
#' An complementary component of the modular library \pkg{screenwerk}, invaluable for the initial set-up of a drug sensitivity screen.
#' \emph{\code{combineDrugs}} is a function that combines drugs at selected doses, with each other. These drug combinations are used for the generation of a dispensing file, 
#' which is used to dispense drugs accordingly for a high-throughput drug combination screen. Drugs are combined either at all doses from a list or at defined ranges.
#' In addition, replicates can be generated either for all drug treatments, or for single drug treatments only.
#' 
#' @param listofDrugs \code{data frame}; a list of drugs in long-format. The data frame needs to include a single column with drug names, the CAS number and optionally the solvent or group.
#' @param listofDoses \code{data frame}; a list of doses in long-format. The data frame needs to include a single column with drug names, doses and units.
#' @param .combineDoses \code{numeric vector}; an individual selection or a range of doses at which drugs are combined.
#' @param .noReplicates \code{numeric}; a value providing the number of replicates for selected drug combinations, representing technical replicates in a drug sensitivity screen.
#' @param .drugRepAttrib \code{character}; a predefined statement specifying, which drug treatments should be replicated. "all", both single drug treatments and drug combinations to be replicated,
#' "single", only single drug treatments to be replicated.
#' @param .pairBy \code{character}; a predefined statement specifying, if drug treatments should be combined by solvent or by a custom group. "solvent", only drugs belonging to the same solvent group will be combined with each other.
#' "group", only drugs belonging to the same group will be combined with each other.
#' @param .inclusive \code{vector}; a selective list of drugs to be combined. Only drugs in this list will be combined with each other.
#' @param .exclusive \code{vector}; a selective list of drugs not to be combine. Drugs in this list won't be combined with other drugs.
#'
#' @details This function is useful for combining drugs at given doses with each other. The combinations are generated from a list of drugs and doses, such as generated with \code{generateListofDoses()}. 
#' Drugs can be combined at individual doses, or at a certain dose range, in which a low index corresponds with a lower dose. i.e. a value of 1 selects for the first, and consequently, for the lowest dose among the list of doses. 
#' The range of 2:5 selects the second lowest up to the fifth dose. A negative number will selected all doses and excluded the dose at the provided index.
#' The default is to combine drugs at all doses with each other.
#' 
#' Each drug combination or drug treatment can be replicated into technical replicates. Depending on the selection, either all drug treatments are replicated or only single drug treatment.
#' If the replication number is not provided, then drugs will not be replicated and only one occurrence of each unique drug dose combination is retained.
#' If no group for replication is selected, all drug combinations will be replicated by default.
#' 
#' The function then generates a list of combinations containing a column with the name, the doses and its corresponding unit for each of the drug pairs. 
#' This list can then be used for the generation of a dispensing file with \code{generateDispensingData()}.
#' 
#' @return Returns an object of class "data.frame", with an essential column of ['drug'] names, ['doses'] and ['units'] for each drug pair.
#'
#' @references Robert Hanes et al.
#' @author Robert Hanes
#' @note Version: 2021.03
#'
#' @examples
#' \donttest{\dontrun{
#' library(metascreen)
#' combineDrugs(listofDoses, .combineDoses=c(2:5), .noReplicates = 3, .drugRepAttrib = "single", .pairBy = c("solvent"))
#' }}
#'
#' @keywords drug screen list drugs combination
#'
#'
#' @importFrom stats setNames
#' @importFrom utils type.convert
#' 
#' @export

combineDrugs <- function(listofDrugs, listofDoses, .combineDoses, .noReplicates = 1, .drugRepAttrib = c("all", "single"), 
                         .pairBy = c("solvent", "group"), .inclusive, .exclusive){

  # Check, if arguments are provided and in the proper format
  
  if(missing(listofDrugs)){stop("Data missing! Please provide a list of drugs.", call. = TRUE)}
  if(!is.list(listofDrugs)){
    stop("in 'listofDrugs'. Argument needs to be an object of class data.frame! Please provide a data frame containing a column with the drug name ['Name'], a unique drug number ['Number'] and the CAS number ['CAS number'].", call. = TRUE)
  } else if (!all(sapply(c("name", "number", "cas", "solvent"), function(x) any(grepl(x, names(listofDrugs), ignore.case = TRUE))))){stop("in 'listofDrugs'. List of drugs has missing data! Please provide a data frame containing a column with the drug name ['Name'], a unique drug number ['Number'], the CAS number ['CAS number'] and the solvent ['Solvent'] in which the drug is dissolved.", call. = TRUE)
  } else {
    colnames(listofDrugs)[grepl("name", names(listofDrugs), ignore.case = TRUE)] <- "NAME"
    colnames(listofDrugs)[grepl("number", names(listofDrugs), ignore.case = TRUE) & !grepl("cas", names(listofDrugs), ignore.case = TRUE)] <- "ID"
    colnames(listofDrugs)[grepl("cas", names(listofDrugs), ignore.case = TRUE)] <- "CAS_NUMBER"
    colnames(listofDrugs)[grepl("solvent", names(listofDrugs), ignore.case = TRUE)] <- "SOLVENT"
    
  }
  
  if(all(sapply(c("group"), function(x) any(grepl(x, names(listofDrugs), ignore.case = TRUE))))){colnames(listofDrugs)[grepl("group", names(listofDrugs), ignore.case = TRUE)] <- "GROUP"}
  
  
  if(missing(listofDoses)){stop("Data missing! Please provide a list of doses.", call. = TRUE)}
  if(!all(sapply(c("Drug", "Dose", "Unit"), function(x) any(grepl(x, names(listofDoses), ignore.case = TRUE))))){stop("List of doses has missing data! Please provide a data frame containing a column for the drug name ['Drug'], the dose ['Dose'] and the unit ['Unit'].", call. = TRUE)}

  # If arguments are missing, set them to default
  if(missing(.combineDoses)){
    message("Doses to be combined not specified. Combining all doses.")
    .combineDoses = 1:max(unique(table(listofDoses$Drug)))
  } else if (max(.combineDoses) > max(unique(table(listofDoses$Drug)))) {
      message("Doses to be combined outside available dose range. Adjusting dose range to ", min(.combineDoses), ":", max(unique(table(listofDoses$Drug))), ".")
      .combineDoses = min(.combineDoses):max(unique(table(listofDoses$Drug)))
    }
    if(any(missing(.noReplicates), is.null(.noReplicates), .noReplicates == 0)){
    message("Number of drug replicates not specified. Using default: 1")
    .noReplicates = 1
  } else if (class(.noReplicates) != "numeric") {
    stop("in '.noReplicates'. Not a number! Please provide the number of replicates.")
  } else if (.noReplicates < 0){
    stop("in '.noReplicates'. Not a positive number! Please provide a valid number of replicates.")
  } else if(.noReplicates%%1 != 0){
    message("Number of replicates not a whole number. Rounding number to ", ifelse(.noReplicates < 1, 1, round(.noReplicates)), ".")
    .noReplicates <- ifelse(.noReplicates < 1, 1, round(.noReplicates))
    }
  if(missing(.drugRepAttrib)){message("Attribute for drug replicates not set. Applying replicates to all drugs.")}
  
  
  # If arguments are missing, set them to default
  if(missing(.pairBy)){
    .pairBy = "all"
    message("Drug pairing argument not specified. No specific drug pairing will be applied.")
  } else if (grepl("solvent", .pairBy, ignore.case = TRUE)){
    if(!all(sapply("solvent", function(x) any(grepl(x, names(listofDrugs), ignore.case = TRUE))))){stop("Drugs are set to be combined by solvent. However, solvents not found in the 'listofDrugs'. Please provide a data frame containing a column with the drug name ['Solvent'].", call. = TRUE)}
    .pairBy = "solvent"
  } else if (grepl("group", .pairBy, ignore.case = TRUE)){
      if(!all(sapply("group", function(x) any(grepl(x, names(listofDrugs), ignore.case = TRUE))))){stop("Drugs are set to be combined by group. However, solvents not found in the 'listofDrugs'. Please provide a data frame containing a column with the drug name ['Group'].", call. = TRUE)}
    .pairBy = "group"
  }
  
  # Check, if drugs are listed by both arguments
  if(!any(missing(.inclusive), missing(.exclusive))){
  if(length(intersect(.inclusive, .exclusive)) != 0){
    message("The drug(s) ", ifelse(length(intersect(.inclusive, .exclusive)) > 1, gsub(",([^,]*)$"," and\\1", paste(intersect(.inclusive, .exclusive), collapse = ", ")), intersect(.inclusive, .exclusive)), " are listed to be combined as well as excluded from being combined. First argument will be ignored and those drugs will be excluded from being combined with other drugs.")
  }}
  
  if(missing(.inclusive)){.inclusive = NULL}
  if(missing(.exclusive)){.exclusive = NULL}
  
  
  .drugRepAttrib <- match.arg(.drugRepAttrib)

  

  cat("\n > Generating a list of combinations. Combining drugs...")
  
  # Create a list of drugs
  # listofDrugs <- unique(listofDoses[[grep("Drug", names(listofDoses), ignore.case = TRUE)]])
  
  # Create a list of doses for each drug
  doseList <- sapply(listofDrugs$NAME, function(x) list(listofDoses[which(listofDoses$Drug == x), grepl("Dose", names(listofDoses))]))
  # Select only the designated four doses for the combination treatment
  doseList <- lapply(doseList, function(x) sort(as.numeric(x))[.combineDoses])
  


  # Combine each drug~dose combination with each other
  listofCombinations <- expand.grid(Drug.1 = paste(listofDoses$Drug, listofDoses$Dose, listofDoses$Unit, sep = ":"), Drug.2 = paste(listofDoses$Drug, listofDoses$Dose, listofDoses$Unit, sep = ":"), stringsAsFactors = FALSE)
  
  # Separate drug and dose from the drug~dose combinations
  listofCombinations <- cbind(setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug.1),':', fixed=TRUE)), stringsAsFactors = FALSE), paste(c("Drug", "Dose", "Unit"), 1, sep = ".")),
                              setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug.2),':', fixed=TRUE)), stringsAsFactors = FALSE), paste(c("Drug", "Dose", "Unit"), 2, sep = ".")))
  
  # Transform drug doses into numeric values
  listofCombinations <- transform(listofCombinations, Dose.1 = as.numeric(Dose.1), Dose.2 = as.numeric(Dose.2))
  
  # Remove same drug combinations, unless they share the same dose
  # Same drug combinations at the same dose are considered single drug treatments
  listofCombinations <- subset(listofCombinations, (Drug.1 == Drug.2 & Dose.1 == Dose.2 | Drug.1 != Drug.2))

  # Remove drug combinations at doses, that are not meant to be combined
  # Combine only doses that have been selected, while keeping the other doses as a single drug treatment
  # listofCombinations <- listofCombinations[sapply(1:nrow(listofCombinations), function(x) (listofCombinations[x, "Drug.1"] == listofCombinations[x, "Drug.2"] | 
  #                                                                                          listofCombinations[x, "Drug.1"] != listofCombinations[x, "Drug.2"] & 
  #                                                                                            listofCombinations[x,"Dose.1"] %in% doseList[[match(listofCombinations[x,"Drug.1"], names(doseList))]] & 
  #                                                                                            listofCombinations[x,"Dose.2"] %in% doseList[[match(listofCombinations[x,"Drug.2"], names(doseList))]] )),]
  listofCombinations <- listofCombinations[
    listofCombinations[, "Drug.1"] == listofCombinations[, "Drug.2"] | 
      sapply(1:nrow(listofCombinations), function(x) {
        listofCombinations[x,"Dose.1"] %in% doseList[[listofCombinations[x,"Drug.1"]]] &
          listofCombinations[x,"Dose.2"] %in% doseList[[listofCombinations[x,"Drug.2"]]]
      }), ]
  
  
  
  # Remove combination that do not meet combination requirements specified by the user
  # Combine only drugs within its group, solvent or custom group
  # NOTE: BOTH, SINGLE AND MULTILEVEL GROUPS ARE NOW SUPPORTED
  # COMBINING (A,B,C) WITH (A), (B) AND (C), BUT NOT (A) WITH (B) OR WITH (C)
  switch(as.character(.pairBy),
         "solvent"={
           # Combine only drugs that belong to the same solvent group
           listofCombinations <- subset(listofCombinations, listofDrugs$SOLVENT[match(listofCombinations$Drug.1, listofDrugs$NAME)] == listofDrugs$SOLVENT[match(listofCombinations$Drug.2, listofDrugs$NAME)])},
         "group"={
           # Combine only drugs that belong to the same user-defined group
           listofCombinations <- listofCombinations[
             listofCombinations[, "Drug.1"] == listofCombinations[, "Drug.2"] | 
               sapply(1:nrow(listofCombinations), function(x) {
                 all(unlist(strsplit(listofDrugs$GROUP[match(listofCombinations[x,"Drug.1"], listofDrugs$NAME)],',|\\/|\\s', fixed=FALSE)) %in% unlist(strsplit(listofDrugs$GROUP[match(listofCombinations[x,"Drug.2"], listofDrugs$NAME)],',|\\/|\\s', fixed=FALSE))) |
                   all(unlist(strsplit(listofDrugs$GROUP[match(listofCombinations[x,"Drug.2"], listofDrugs$NAME)],',|\\/|\\s', fixed=FALSE)) %in% unlist(strsplit(listofDrugs$GROUP[match(listofCombinations[x,"Drug.1"], listofDrugs$NAME)],',|\\/|\\s', fixed=FALSE)))
               }), ] },
         "all"={}
         )
  
  

  
  # Combine only selective drugs; inclusive combination of drugs
  if(!is.null(.inclusive)){
    # Combine only drugs that are specified to be combined
    listofCombinations <- subset(listofCombinations, listofCombinations$Drug.1 %in% .inclusive & listofCombinations$Drug.2 %in% .inclusive)
    if(dim(listofCombinations)[1] == 0){stop('\r', paste("Inclusive drugs to be combined not found after pairing by ", .pairBy, ".", " No drugs combined!", sep = ""), call. = FALSE)}
  }
  
  # Exclude selective drugs from being combined; exclusive combination of drugs 
  if(!is.null(.exclusive)){
    listofCombinations <- subset(listofCombinations, !(listofCombinations$Drug.1 %in% .exclusive & listofCombinations$Drug.2 %in% .exclusive &
                                                         listofCombinations$Drug.1 != listofCombinations$Drug.2))
    if(dim(listofCombinations)[1] == 0){stop('\r', paste("No drugs left to combine after exlusion of drugs and pairing by ", .pairBy, ".", " No drugs combined!", sep = ""), call. = FALSE)}
  }
  
  cat('\r', "Drugs succesfully combined!", rep(" ", 29), "\n")
  
  
  # Preserve numeric formatting of drug concentrations by converting them to characters
  listofCombinations[grep("Dose", names(listofCombinations), value=TRUE)] <- apply(listofCombinations[grep("Dose", names(listofCombinations), value=TRUE)], MARGIN = 2, FUN = function(i) paste(i))
  # Merge name, dose and unit for each drug pair 
  listofCombinations <- data.frame(Drug.1 = apply(listofCombinations[grep("1", names(listofCombinations), value=TRUE)], MARGIN = 1, FUN = function(i) paste(i, collapse = ":")),
                                   Drug.2 = apply(listofCombinations[grep("2", names(listofCombinations), value=TRUE)], MARGIN = 1, FUN = function(i) paste(i, collapse = ":")))
  # Reset row names to default
  rownames(listofCombinations) <- NULL
  
  
  # Removing duplicate combinations:
  # by sorting the rows and eliminating duplicates from a transposed matrix
  listofCombinations <- listofCombinations[!duplicated(t(apply(listofCombinations, 1, sort))), ]
  
  
  
  # CHECK number of combinations
  # Note: combine all other drugs at all their doses with one drug and all it's doses. do that for all the drugs
  #       divide that by half avoiding transposed combinations
  #       add single drug combinations each drug at each dose combined with itself
  # [ ( no. of drugs - 1 * no. of doses ) * no. of doses * no. of drugs ] / [ 2 ] + [ no. of drugs * no. of single drug doses ]
  # cat('\r', ifelse(
  #   ((length(listofDrugs$NAME)-1) * unique(sapply(doseList, length)) * unique(sapply(doseList, length)) * length(listofDrugs$NAME)) / 2 + (length(listofDrugs$NAME) * unique(table(listofDoses$Drug)))
  #   == nrow(listofCombinations), "Combinations passed check! A list of combinations has been succesfully generated.", "Combinations failed check!"), "\n")
  
  
  
  # Adding replicates for an individual treatment group
  switch(as.character(.drugRepAttrib),
         "all"={
           # Adding replicates to all drug treatments
           listofCombinations <- listofCombinations[rep(row.names(listofCombinations), .noReplicates), ]},
         "single"={
           # Adding triplicates to single drug treatments,
           # only if the occurrence of the single drug treatments is equal 1 to avoid multiple generation of triplicates
           if(table(unlist(listofCombinations[which(listofCombinations$Drug.1==listofCombinations$Drug.2),][1]))[[2]] == 1){
             listofCombinations <- listofCombinations[rep(row.names(listofCombinations), ifelse(listofCombinations$Drug.1==listofCombinations$Drug.2, .noReplicates, 1)), ]
           }},{ message("Please select a valid treatment group to add replicates to. The treatment groups supported are \"all\" or \"single\" drug treatments.") })
  
  
  
  # Separate name, dose and unit for each drug pair
  # listofCombinations <- cbind(setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug.1),':', fixed=TRUE)), stringsAsFactors = FALSE), paste(c("Drug", "Dose", "Unit"), 1, sep = ".")),
  #                             setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug.2),':', fixed=TRUE)), stringsAsFactors = FALSE), paste(c("Drug", "Dose", "Unit"), 2, sep = ".")))
  # listofCombinations <- utils::type.convert(listofCombinations, as.is = TRUE)
  
  
  
  
  # Prepare data for format conversion
  listofCombinations <- as.data.frame(listofCombinations)
  
  # Convert list from wide to long format, adding columns ID and Drug
  listofCombinations <- stats::reshape(listofCombinations, idvar = "ID", varying = list(grepl("Drug", names(listofCombinations), ignore.case = TRUE)), v.names = c("Drug"), timevar = "Drug.Number", times = paste("Drug", LETTERS[grep("Drug", names(listofCombinations))]), direction = "long", new.row.names = NULL)
  # Reorder columns and rows, reset row names
  listofCombinations <- listofCombinations[,c("ID", sort(grep("Drug", names(listofCombinations), value = TRUE), decreasing = TRUE))]
  listofCombinations <- listofCombinations[with(listofCombinations, order(ID, Drug.Number)),]
  rownames(listofCombinations) <- NULL
  
  # Replace initial drug column with individual columns for drug name, dose and unit
  listofCombinations <- cbind(within(listofCombinations, rm("Drug")), setNames(data.frame(do.call('rbind', strsplit(as.character(listofCombinations$Drug),':', fixed=TRUE)), stringsAsFactors = FALSE), c("Drug", "Drug.Concentration", "Unit")))
  listofCombinations <- transform(listofCombinations, Drug.Concentration = as.numeric(Drug.Concentration))
  attr(listofCombinations$Drug.Concentration, "Unit") <- listofCombinations$Unit
  
  listofCombinations <- utils::type.convert(listofCombinations, as.is = TRUE)
  
  
  
  return(listofCombinations)
  
}
