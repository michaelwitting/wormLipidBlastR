#' Generation of fragmentation spectra for PC [M+FA-H]-
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @import Spectra
#'
#' @export
create_neg_PCO <- function(lipid_info, adduct, template = NA, ...) {
  
  ## some sanity checks here ---------------------------------------------------
  # TODO add checks for correct lipid class here
  if(!adduct %in% c("[M+CHO2]-", "[M+C2H3O2]-")) {
    stop("unsupported adduct")
  }
  
  ## check if template file is supplied, other wise use hard coded template ----
  if(is.na(template) & adduct == "[M+CHO2]-") {
    template <- .template_neg_PCO_MFAH()
  } else if(is.na(template) & adduct == "[M+C2H3O2]-") {
    template <- .template_neg_PCO_MHAcH()
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
  
  ## get fatty acids -----------------------------------------------------------
  fattyAcids <- unlist(lipidomicsUtils::isolate_radyls(lipid))
  
  # calculate masses of different fatty acids ----------------------------------
  sn1_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[1])
  sn2_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[2])
  
  # calculate masses of different fatty acids ----------------------------------
  if(stringr::str_detect(fattyAcids[1], "^O-")) {
    alkyl_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[1])
    acyl_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[2])
  } else {
    alkyl_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[2])
    acyl_mass <- lipidomicsUtils::calc_intact_acyl_mass(fattyAcids[1])
  }
  
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
    polarity = 0L,
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
.template_neg_PCO_MFAH <- function() {
  
  template <- list(

  )
  
  # return template
  return(template)
}



#'
#'
#'
.template_neg_PCO_MHAcH <- function() {
  
  template <- list(

  )
  
  # return template
  return(template)
}

#'
#'
#' @export
buildingblocks_neg_PCO <- function() {
  
  building_blocks <- c("adduct_mass", "gpc_mass", "pc_mass",
                       "water_mass", "choline_mass", "proton_mass",
                       "alkyl_mass", "acyl_mass")
  
  # return values
  return(building_blocks)
}

