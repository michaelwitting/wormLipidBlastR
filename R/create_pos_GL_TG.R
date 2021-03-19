#' Generation of TG [M+Na]+ or [M+NH4]+ MS2 spectrum
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @import Spectra
#'
#' @export
create_pos_TG <- function(lipid_info, adduct, template = NA, ...) {
  
  ## some sanity checks here ---------------------------------------------------
  # TODO add checks for correct lipid class here
  if(!adduct %in% c("[M+Na]+", "[M+NH4]+")) {
    stop("unsupported adduct")
  }
  
  ## check if template file is supplied, other wise use hard coded template ----
  if(is.na(template) & adduct == "[M+Na]+") {
    template <- .template_pos_TG_MNa()
  } else if(is.na(template) & adduct == "[M+NH4]+") {
    template <- .template_pos_TG_MNH4()
  }
  
  ## get lipid information------------------------------------------------------
  lipid <- lipid_info$lipid
  id <- lipid_info$id
  chemFormula <- lipid_info$chemFormula
  
  ## get masses for calculation ------------------------------------------------
  glycerol_mass <- lipidomicsUtils:::glycerol_mass
  water_mass <- lipidomicsUtils:::water_mass
  proton_mass <- lipidomicsUtils:::proton_mass
  sodium_ion_mass <- lipidomicsUtils:::sodium_ion_mass
  
  ## get fatty acids -----------------------------------------------------------
  fattyAcids <- unlist(lipidomicsUtils::isolate_radyls(lipid))
  
  # calculate masses of different fatty acids ----------------------------------
  sn1_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[1])
  sn2_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[2])
  sn3_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[3])
  
  # calculate mass of intact PE ------------------------------------------------
  lipid_mass <- unlist(lipidomicsUtils::calc_lipid_mass(lipid))
  
  # calculate adduct mass ------------------------------------------------------
  adduct_mass <- as.numeric(MetaboCoreUtils::mass2mz(lipid_mass, adduct))
  
  # generate MS2 spectrum ------------------------------------------------------
  mz <- unname(unlist(lapply(names(template), function(x) {eval(parse(text = x))})))
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
.template_pos_TG_MNa <- function() {
  
  template <- list(

  )
  
  # return template
  return(template)
  
}

#'
#'
#'
.template_pos_TG_MNH4 <- function() {
  
  template <- list(
    "adduct_mass" = 999,
    "adduct_mass - rcdk::get.formula('NH3')@mass" = 100,
    "adduct_mass - sn1_mass - rcdk::get.formula('NH3')@mass" = 900,
    "adduct_mass - sn2_mass - rcdk::get.formula('NH3')@mass" = 500,
    "adduct_mass - sn3_mass - rcdk::get.formula('NH3')@mass" = 900,
    "adduct_mass - sn1_mass - sn2_mass - rcdk::get.formula('NH3')@mass" = 10,
    "adduct_mass - sn1_mass - sn3_mass - rcdk::get.formula('NH3')@mass" = 10,
    "adduct_mass - sn1_mass - sn2_mass - rcdk::get.formula('NH3')@mass" = 10,
    "adduct_mass - sn2_mass - sn3_mass - rcdk::get.formula('NH3')@mass" = 10,
    "sn1_mass - rcdk::get.formula('OH', charge = -1)@mass" = 10,
    "sn2_mass - rcdk::get.formula('OH', charge = -1)@mass" = 10,
    "sn3_mass - rcdk::get.formula('OH', charge = -1)@mass" = 10
  )
  
  # return template
  return(template)
  
}

#'
#'
#' @export
buildingblocks_pos_TG <- function() {
  
  building_blocks <- c("adduct_mass", "glycerol_mass", "water_mass",
                       "proton_mass", "sodium_ion_mass","sn1_mass", "sn2_mass",
                       "sn3_mass")
  
  # return values
  return(building_blocks)
}
