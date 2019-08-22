#' Generation of PS [M+H]+ MS2 spectrum
#'
#' @import MSnbase
#' @import S4Vectors
#'
#' @export
create_pos_PS <- function(lipid_info, adduct, template = NA, ...) {
  
  ## some sanity checks here ---------------------------------------------------
  # TODO add checks for correct lipid class here
  if(!adduct %in% c("[M+H]+", "[M+Na]+")) {
    stop("unsupported adduct")
  }
  
  ## check if template file is supplied, other wise use hard coded template ----
  if(is.na(template) & adduct == "[M+H]+") {
    template <- .template_pos_PS_MH()
  } else if(is.na(template) & adduct == "[M+Na]+") {
    template <- .template_pos_PS_MNa()
  }
  
  ## get lipid information------------------------------------------------------
  lipid <- lipid_info$lipid
  id <- lipid_info$id
  chemFormula <- lipid_info$chemFormula
  
  ## get masses for calculation ------------------------------------------------
  gps_mass <- lipidomicsUtils:::gps_mass
  ps_mass <- lipidomicsUtils:::ps_mass
  water_mass <- lipidomicsUtils:::water_mass
  serine_mass <- lipidomicsUtils:::serine_mass
  proton_mass <- lipidomicsUtils:::proton_mass
  sodium_ion_mass <- lipidomicsUtils:::sodium_ion_mass
  
  ## get fatty acids -----------------------------------------------------------
  fattyAcids <- lipidomicsUtils::isolate_fatty_acyls(lipid)
  
  # calculate masses of different fatty acids ----------------------------------
  sn1_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[1])
  sn2_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[2])
  
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
  mcols(lipidSpectrum)$ionMode <- "POSITIVE"
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
.template_pos_PS_MH <- function() {
  
  template <- list(
    "adduct_mass" = 20,
    "adduct_mass - ps_mass" = 999,
    "sn1_mass - rcdk::get.formula('OH', charge = -1)@mass" = 15,
    "sn2_mass - rcdk::get.formula('OH', charge = -1)@mass" = 15
  )
  
  # return template
  return(template)
  
}


#'
#'
#'
.template_pos_PS_MNa <- function() {
  
  template <- list(
    "adduct_mass" = 100,
    "adduct_mass - ps_mass" = 100,
    "adduct_mass - ps_mass - sodium_ion_mass + proton_mass" = 100,
    "ps_mass + sodium_ion_mass" = 999

  )
  
  # return template
  return(template)
  
}

#'
#'
#' @export
buildingblocks_pos_PS <- function() {
  
  building_blocks <- c("adduct_mass", "gps_mass", "ps_mass",
                       "water_mass", "serine_mass", "proton_mass",
                       "sodium_ion_mass","sn1_mass", "sn2_mass")
  
  # return values
  return(building_blocks)
}
