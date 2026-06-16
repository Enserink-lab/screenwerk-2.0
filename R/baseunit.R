#' Converts a value with a given unit to its base unit
#'
#' @description
#' A complementary component of the modular library \pkg{screenwerk}, valuable for a wide range of analytical tasks in a drug sensitivity screen. 
#' \emph{\code{baseunit}} is a function that converts a given number from a given unit to its base unit. 
#' This tool is helpful for converting and matching values with different units to a base unit that can be compared.
#' \emph{baseunit} was designed with the objective in mind to address the incompatibility of different volumes, weights and values.
#'
#' @param value \code{integer}; a numeric value designating the value to be converted.
#' @param unit \code{character}; a character denoting the unit of the given value.
#' @param .simplify \code{logical}; if TRUE, then only the converted value will be returned. The default is FALSE, in which, both the value and the base unit will be returned as a list.
#'
#' @details The value for \emph{unit} is restricted to units with the prefixes specified by the twenty-four prefixes for the International System of Units (SI).
#' Those units are defined with the following symbols as prefixes: Q (quetta), R (ronna), Y (yotta), Z (zetta), E (exa), P (peta), T (tera), G (giga), M (mega), k (kilo), h (hecto), da (deca),
#' d (deci), c (centi), m (milli), μ (micro), n (nano), p (pico), f (femto), a (atto), z (zepto), y (yocto), r (ronto), q (quecto).
#'
#' @return Returns an object with a list of \emph{value} and \emph{unit}.
#' \preformatted{
#' }
#'
#' Each list contains two objects: \emph{value}, an object of type "numeric" containing the converted value and the unit \emph{unit} an object of type "character" containing the base unit.
#'
#' @references Robert Hanes et al.
#' @author Robert Hanes
#' @note Version: 2021.03
#'
#' @examples
#' baseunit(3, "nl")
#' baseunit(3, "nL")
#' baseunit(3, "mg")
#' baseunit(3, "mg", .simplify=TRUE)
#' baseunit(3.48, "dag")
#' baseunit(3.48, "kg")
#' 
#' value = list(value = 3, unit = "mg")
#' baseunit(value, .simplify=FALSE)
#'
#' @keywords drug screen base unit conversion
#'
#'
#' @export

baseunit <- function(value, unit, .simplify=FALSE){
  
  # Check, if arguments are in the accepted format
  if(missing(value)){ stop("Data missing! Please provide a list or a value with a unit.", call. = TRUE) }
  if(!any(is.numeric(value), is.list(value))){ stop("The value needs to be a number or a list with the first element being a number.", call. = TRUE) 
  } else if(is.list(value)){
      if(!is.numeric(value[[1]])){ stop("The first element of the list needs to be a number.", call. = TRUE) }
      else if(!is.character(value[[2]])){ stop("The second element of the list needs to be a character", call. = TRUE) }
      else if(!any(nchar(value[[2]]) == 1, grepl(ifelse(nchar(value[[2]]) == 3, substring(value[[2]], 1, 2), substring(value[[2]], 1, 1)), c("Q", "R", "Y", "Z", "E", "P", "T", "G", "M", "k", "h", "da", "d", "c", "m", "μ", "u", "n", "p", "f", "a", "z", "y", "r", "q"), ignore.case = TRUE))){ stop("The unit is not supported! Please provide a valid unit. Te following prefixes for units are supported: Q (quetta), R (ronna), Y (yotta), Z (zetta), E (exa), P (peta), T (tera), G (giga), M (mega), k (kilo), h (hecto), da (deca), d (deci), c (centi), m (milli), μ (micro), n (nano), p (pico), f (femto), a (atto), z (zepto), y (yocto), r (ronto), q (quecto)", call. = TRUE) } 
      else {
        unit = value[[2]]
        value = value[[1]]
      }
    
    } else {
    if(missing(unit)){ stop("Data missing! Please provide a unit.", call. = TRUE) }
    if(!is.character(unit)){stop("The unit needs to be a string.", call. = TRUE)}
    else if(!any(grepl(ifelse(nchar(unit) == 3, substring(unit, 1, 2), substring(unit, 1, 1)), c("Q", "R", "Y", "Z", "E", "P", "T", "G", "M", "k", "h", "da", "d", "c", "m", "μ", "u", "n", "p", "f", "a", "z", "y", "r", "q"), ignore.case = TRUE))){ stop("The unit is not supported! Please provide a valid unit. Te following prefixes for units are supported: Q (quetta), R (ronna), Y (yotta), Z (zetta), E (exa), P (peta), T (tera), G (giga), M (mega), k (kilo), h (hecto), da (deca), d (deci), c (centi), m (milli), μ (micro), n (nano), p (pico), f (femto), a (atto), z (zepto), y (yocto), r (ronto), q (quecto)", call. = TRUE) }
    }
  
  if(missing(.simplify)){ .simplify = FALSE }
    
  
  
  # Check, whether the given unit has a prefix,
  # otherwise it will be treated as a base unit
  if(!nchar(unit) < 2){ 
    
    # Check, if the unit has the two-letter prefix for deka
    if(all(nchar(unit) == 3, substring(unit, 1, 2) == "da")){
      
      c_value = value * `^`(10,1)
      suffix_u <- substring(unit, nchar(unit))
      
    } else {
    
      # otherwise treat unit with a one-letter prefix
      prefix_u <- substring(unit, 1, 1)
      suffix_u <- substring(unit, nchar(unit))
      
      # Correct prefix for the unit micro
      prefix_u <- ifelse(prefix_u == "u", "μ", prefix_u)
      
      # If unit has a prefix, calculate the converted value base on its prefix
      switch (prefix_u,
              
              "Q" = {c_value = value * `^`(10,30)},
              "R" = {c_value = value * `^`(10,27)},
              "Y" = {c_value = value * `^`(10,24)},
              "Z" = {c_value = value * `^`(10,21)},
              "E" = {c_value = value * `^`(10,18)},
              "P" = {c_value = value * `^`(10,15)},
              "T" = {c_value = value * `^`(10,12)},
              "G" = {c_value = value * `^`(10,9)},
              "M" = {c_value = value * `^`(10,6)},
              "k" = {c_value = value * `^`(10,3)},
              "h" = {c_value = value * `^`(10,2)},
              
              "d" = {c_value = value * `^`(10,-1)},
              "c" = {c_value = value * `^`(10,-2)},
              "m" = {c_value = value* `^`(10,-3)},
              "μ" = {c_value = value* `^`(10,-6)},
              "n" = {c_value = value* `^`(10,-9)},
              "p" = {c_value = value * `^`(10,-12)},
              "f" = {c_value = value * `^`(10,-15)},
              "a" = {c_value = value * `^`(10,-18)},
              "z" = {c_value = value * `^`(10,-21)},
              "y" = {c_value = value * `^`(10,-24)},
              "r" = {c_value = value * `^`(10,-27)},
              "q" = {c_value = value * `^`(10,-30)},
              
              {c_value = value* `^`(10,0)}
      )
    }
    
    
  }else{ 
    
    # If the given unit has no prefix, it will be treated as a base unit
    # message("Unit has no prefix and will be treated as base unit. Value was not converted and remains the same.")
    c_value <- value
    
    }

  
  # Create an list object "baseunit"
  if (.simplify == TRUE){
    .baseunit <- c_value
  } else {
    .baseunit <- list(value = c_value, unit = suffix_u)
  }
  
   
  return(.baseunit)
  
}
