#' Generation of fragmentation spectra for PC [M+FA-H]-
#' 
#' @importFrom S4Vectors DataFrame
#' @importFrom IRanges NumericList
#' @import Spectra
#'
#' @export
create_neg_HexCer <- function(lipid_info, adduct, template = NA, ...) {
  
  ## some sanity checks here ---------------------------------------------------
  # TODO add checks for correct lipid class here
  if(!adduct %in% c("[M-H]-", "[M+FA-H]-")) {
    stop("unsupported adduct")
  }
  
  ## check if template file is supplied, other wise use hard coded template ----
  if(is.na(template) & adduct == "[M-H]-") {
    template <- .template_Hexneg_Cer_MH()
  } else if(is.na(template) & adduct == "[M+FA-H]-") {
    template <- .template_neg_HexCer_MFAH()
  } #else if(is.na(template) & adduct == "[M+HAc-H]-") {
    #template <- .template_neg_HxCer_MHAcH()
  #}
  
  ## get lipid information------------------------------------------------------
  lipid <- lipid_info$lipid
  id <- lipid_info$id
  chemFormula <- lipid_info$chemFormula
  
  ## get masses for calculation ------------------------------------------------
  water_mass <- lipidomicsUtils:::water_mass
  proton_mass <- lipidomicsUtils:::proton_mass
  hexose_mass <- lipidomicsUtils:::hexose_mass
  
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
.template_neg_HexCer_MH <- function() {
  
  template <- list(
    
  )
  
  # return template
  return(template)
}


#'
#'
#'
.template_neg_HexCer_MFAH <- function() {
  
  template <- list(
    "adduct_mass" = 5,
    "adduct_mass - rcdk::get.formula('CH2O2')@mass" = 600,
    "adduct_mass - rcdk::get.formula('CH2O2')@mass - hexose_mass" = 999
  )
  
  # return template
  return(template)
}

#'
#'
#'
.template_neg_HexCer_MHAcH <- function() {
  
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
                       "hexose_mass", "proton_mass")
  
  # return values
  return(building_blocks)
}

