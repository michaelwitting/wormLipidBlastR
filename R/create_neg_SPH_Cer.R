#' Generation of fragmentation spectra for PC [M+FA-H]-
#' 
#' @import MSnbase
#' @import S4Vectors
#'
#' @export
create_neg_Cer <- function(lipid_info, adduct, template = NA, ...) {
  
  ## some sanity checks here ---------------------------------------------------
  # TODO add checks for correct lipid class here
  if(!adduct %in% c("[M-H]-", "[M+FA-H]-", "[M+HAc-H]-")) {
    stop("unsupported adduct")
  }
  
  ## check if template file is supplied, other wise use hard coded template ----
  if(is.na(template) & adduct == "[M-H]-") {
    template <- .template_neg_Cer_MH()
  } else if(is.na(template) & adduct == "[M+FA-H]-") {
    template <- .template_neg_Cer_MFAH()
  } else if(is.na(template) & adduct == "[M+HAc-H]-") {
    template <- .template_neg_Cer_MHAcH()
  }
  
  ## get lipid information------------------------------------------------------
  lipid <- lipid_info$lipid
  id <- lipid_info$id
  chemFormula <- lipid_info$chemFormula
  
  ## get masses for calculation ------------------------------------------------
  water_mass <- lipidomicsUtils:::water_mass
  proton_mass <- lipidomicsUtils:::proton_mass
  
  ## get fatty acids and calculate mass ----------------------------------------
  fattyAcids <- lipidomicsUtils::isolate_fatty_acyls(lipid)
  acyl_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[1])
  
  # get sphingoid base ---------------------------------------------------------
  sphingoidBase <- lipidomicsUtils::isolate_sphingoid_base(lipid)
  sphingoid_mass <- lipidomicsUtils::calc_sphingoid_mass(sphingoidBase)
  
  # calculate mass of intact PC ------------------------------------------------
  lipid_mass <- lipidomicsUtils::calc_lipid_mass(lipid)
  
  # calculate adduct mass ------------------------------------------------------
  adduct_mass <- lipidomicsUtils::calc_adduct_mass(lipid_mass, adduct)
  
  # generate MS2 spectrum ------------------------------------------------------
  mz <- unlist(lapply(names(template), function(x) {eval(parse(text = x))}))
  
  # get intensity values
  int <- unlist(lapply(unname(template), function(x) {eval(parse(text = x))}))
  
  # generate splash ------------------------------------------------------------
  splash <- splashR::getSplash(cbind(mz = mz,
                                     intensity = int))
  
  # check additional arguments
  add_args <- list(...)
  
  # collision energiy
  if(!"collisionEnergy" %in% names(add_args)) {
    collisionEnergy <- 40
  } else {
    collisionEnergy <- add_args[["collisionEnergy"]]
  }
  
  # create MS2 spectrum
  ms2Spec <- new(
    "Spectrum2",
    merged = 0,
    precScanNum = as.integer(1),
    precursorMz = adduct_mass,
    precursorIntensity = 100,
    precursorCharge = as.integer(1),
    mz = mz,
    intensity = int,
    collisionEnergy = collisionEnergy,
    centroided = TRUE
  )
  
  # create spectra object for return
  lipidSpectrum <- Spectra(ms2Spec)
  
  # fill basic metadata
  mcols(lipidSpectrum)$id <- id
  mcols(lipidSpectrum)$name <- lipid
  mcols(lipidSpectrum)$formula <- chemFormula
  mcols(lipidSpectrum)$exactMass <- lipid_mass
  mcols(lipidSpectrum)$smiles <- NA
  mcols(lipidSpectrum)$inchi <- NA
  mcols(lipidSpectrum)$instrument <- "prediction"
  mcols(lipidSpectrum)$instrumentType <- "prediction"
  mcols(lipidSpectrum)$msType <- "MS2"
  mcols(lipidSpectrum)$ionMode <- "NEGATIVE"
  mcols(lipidSpectrum)$precursorMz <- adduct_mass
  mcols(lipidSpectrum)$precursorType <- adduct
  mcols(lipidSpectrum)$splash <- splash
  mcols(lipidSpectrum)$numPeak <- peaksCount(ms2Spec)
  
  # additional metadata supplied via ...
  # RT
  if(!"RT" %in% names(add_args)) {
    mcols(lipidSpectrum)$rt <- 0
  } else {
    mcols(lipidSpectrum)$rt <- add_args[["RT"]]
  }
  
  # CCS
  if(!"CCS" %in% names(add_args)) {
    mcols(lipidSpectrum)$ccs <- 0
  } else {
    mcols(lipidSpectrum)$ccs <- add_args[["CCS"]]
  }
  
  # return spectrum
  return(lipidSpectrum)
  
}

#'
#'
#'
.template_neg_Cer_MH <- function() {
  
  template <- list(
    
  )
  
  # return template
  return(template)
}


#'
#'
#'
.template_neg_Cer_MFAH <- function() {
  
  template <- list(

  )
  
  # return template
  return(template)
}

#'
#'
#'
.template_neg_Cer_MHAcH <- function() {
  
  template <- list(
    
  )
  
  # return template
  return(template)
}


#'
#'
#' @export
buildingblocks_neg_Cer <- function() {
  
  building_blocks <- c("adduct_mass", "sphingoid_mass", "acyl_mass",
                       "proton_mass")
  
  # return values
  return(building_blocks)
}

