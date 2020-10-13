#' Generation of LPC [M+H]+ MS2 spectrum
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @import Spectra
#'
#' @export
create_pos_LPC <- function(lipid_info, adduct, template = NA, ...) {
  
  ## some sanity checks here ---------------------------------------------------
  # TODO add checks for correct lipid class here
  if(!adduct %in% c("[M+H]+", "[M+Na]+")) {
    stop("unsupported adduct")
  }
  
  ## check if template file is supplied, other wise use hard coded template ----
  if(is.na(template) & adduct == "[M+H]+") {
    template <- .template_pos_LPC_MH()
  } else if(is.na(template) & adduct == "[M+Na]+") {
    template <- .template_pos_LPC_MNa()
  }
  
  ## get lipid information------------------------------------------------------
  lipid <- lipid_info$lipid
  id <- lipid_info$id
  chemFormula <- lipid_info$chemFormula
  
  ## get masses for calculation ------------------------------------------------
  gpc_mass <- lipidomicsUtils:::gpc_mass
  pc_mass <- lipidomicsUtils:::pc_mass
  water_mass <- lipidomicsUtils:::water_mass
  choline_mass <- lipidomicsUtils:::choline_mass
  proton_mass <- lipidomicsUtils:::proton_mass
  sodium_ion_mass <- lipidomicsUtils:::sodium_ion_mass
  
  ## get fatty acids -----------------------------------------------------------
  fattyAcids <- lipidomicsUtils::isolate_fatty_acyls(lipid)
  
  # calculate masses of different fatty acids ----------------------------------
  if(fattyAcids[1] == "0:0") {
    acyl_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[2])
  } else if (fattyAcids[2] == "0:0") {
    acyl_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[2])
  }

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
.template_pos_LPC_MH <- function() {
  
  template <- list(
    "adduct_mass" = 100,
    "adduct_mass - water_mass" = 100,
    "pc_mass + proton_mass" = 999,
    "choline_mass + proton_mass" = 700,
    "choline_mass + proton_mass - water_mass" = 100

  )
  
  # return template
  return(template)
  
}

#'
#'
#'
.template_pos_LPC_MNa <- function() {
  
  template <- list(

  )
  
  # return template
  return(template)
  
}

#'
#'
#' @export
buildingblocks_pos_LPC <- function() {
  
  building_blocks <- c("adduct_mass", "gpc_mass", "pc_mass",
                       "water_mass", "choline_mass", "proton_mass",
                       "sodium_ion_mass","acyl_mass")
  
  # return values
  return(building_blocks)
}
