#' Function to write .mgf file
#'
#' @export
writeMsp <- function(spectrum, file) {
  
  ## prepare for writing
  con <- file(file, "w")
  on.exit(close(con))
  
  # custom cat function
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  # write content
  .cat("NAME: ", mcols(spectrum)$id, ", ", paste(mcols(spectrum)$name, "predicted", "MS2", collapse = ";"), "\n")
  .cat("PRECURSORMZ: ", mcols(spectrum)$precursorMz, "\n")
  .cat("PRECURSORTYPE: ", mcols(spectrum)$precursorType, "\n")
  .cat("SMILES: ", mcols(spectrum)$smiles, "\n")
  .cat("INCHIKEY: NA\n")
  .cat("FORMULA: ", mcols(spectrum)$formula, "\n")
  .cat("RETENTIONTIME: NA\n")
  .cat("IONMODE: ", mcols(spectrum)$ionMode, "\n")
  .cat("COMPOUNDCLASS: NA\n")
  .cat("Comments: MS2 spectrum predicted with WormLipidBlastR\n")
  .cat("Num Peaks:", mcols(spectrum)$numPeak, "\n")
  .cat(paste(mz(spectrum[[1]]),
             intensity(spectrum[[1]]), collapse = "\n"))
  .cat("\n")
  
}

#' Function to write .mgf file
#'
#' @export
writeMgf <- function(spectrum, file) {
  
  ## prepare for writing
  con <- file(file, "w")
  on.exit(close(con))
  
  # custom cat function
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  # write content
  
}

#' Function to write .mb file
#'
#' @export
writeMassBank <- function(spectrum, file) {
  
  ## prepare for writing
  con <- file(file, "w")
  on.exit(close(con))
  
  # custom cat function
  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  # write content
  .cat("ACCESSION: ", mcols(spectrum)$id, "\n")
  .cat("RECORD_TITLE: ", paste(mcols(spectrum)$name, "predicted", "MS2", collapse = ";"), "\n")
  .cat("DATE: ", Sys.Date(), "\n")
  .cat("AUTHORS: Michael Witting, Helmholtz Zentrum Muenchen\n")
  .cat("LICENSE: CC BY\n")
  .cat("COPYRIGHT: Copyright (C) 2018\n")
  .cat("COMMENT: MS2 spectrum predict with WormLipidBlastR\n")
  .cat("CH$NAME: ", mcols(spectrum)$name, "\n")
  .cat("CH$FORMULA: ", mcols(spectrum)$formula, "\n")
  .cat("CH$EXACT_MASS: ", mcols(spectrum)$exactMass, "\n")
  .cat("CH$SMILES: ", mcols(spectrum)$smiles, "\n")
  .cat("CH$IUPAC: ", mcols(spectrum)$inchi, "\n")
  .cat("AC$INSTRUMENT: prediction\n")
  .cat("AC$INSTRUMENT_TYPE: prediction\n")
  .cat("AC$MASS_SPECTROMETRY: MS_TYPE MS2\n")
  .cat("AC$MASS_SPECTROMETRY: ION_MODE ", mcols(spectrum)$ionMode,"\n")
  .cat("AC$MASS_SPECTROMETRY: COLLISION_ENERGY ", mcols(spectrum)$collisionEnergy, "\n")
  .cat("MS$FOCUSED_ION: PRECURSOR_M/Z ", mcols(spectrum)$precursorMz, "\n")
  .cat("MS$FOCUSED_ION: PRECURSOR_TYPE ", mcols(spectrum)$precursorType, "\n")
  .cat("PK$SPLASH: ", mcols(spectrum)$splash ,"\n")
  .cat("PK$NUM_PEAK: ", mcols(spectrum)$numPeak, "\n")
  .cat("PK$PEAK: m/z int. rel.int.\n")
  .cat(paste("", mz(spectrum[[1]]), intensity(spectrum[[1]]), intensity(spectrum[[1]]), collapse = "\n"))
  .cat("\n//")
}

