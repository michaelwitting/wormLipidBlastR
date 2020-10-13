#' Generation of PS [M+H]+ MS2 spectrum
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @import Spectra
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
  int <- unlist(lapply(unname(template), function(x) {eval(parse(text = x))}))
  
  spec <- DataFrame(mz = mz,
                    int = int)
  
  spec <- spec[order(mz),]
  
  # generate splash ------------------------------------------------------------
  splash <- splashR::getSplash(cbind(mz = mz,
                                     intensity = int))
  
  # check additional arguments
  add_args <- list(...)
  
  # collision energy
  if(!"collisionEnergy" %in% names(add_args)) {
    collisionEnergy <- 40
  } else {
    collisionEnergy <- add_args[["collisionEnergy"]]
  }
  
  # additional metadata supplied via ...
  # RT
  if(!"RT" %in% names(add_args)) {
    rtime <- 0
  } else {
    rtime <- as.numeric(add_args[["RT"]])
  }
  
  # CCS
  if(!"CCS" %in% names(add_args)) {
    ccs <- 0
  } else {
    ccs <- as.numeric(add_args[["CCS"]])
  }
  
  spd <- DataFrame(
    msLevel = 2L,
    precursorMz = adduct_mass,
    rtime = rtime,
    accession = id,
    name = lipid,
    exactmass = lipid_mass,
    formula = chemFormula,
    smiles = NA,
    inchi = NA,
    instrument = "prediction",
    instrument_type = "prediction",
    ms_ms_type = "MS2",
    polarity = 1L,
    adduct = adduct,
    splash = splash,
    pknum = length(mz)
  )
  
  spd$mz <- IRanges::NumericList(spec$mz)
  spd$intensity <- IRanges::NumericList(spec$int)
  
  Spectra(spd)
  
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
