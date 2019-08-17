#' @title Prediction of spectra from input file
#'
#' @export
predict_spectra <- function(file) {
  
  
  
  ## read lipids in file -------------------------------------------------------
  lipids <- data.frame(read.table(file, header = TRUE, sep = "\t"),
                       stringsAsFactors = FALSE)
  
  ## make Spectra objects to store predicted spectra ---------------------------
  wormLipidBlast_Neg_Spectra <- MSnbase::Spectra()
  wormLipidBlast_Pos_Spectra <- MSnbase::Spectra()
  
  # make descision on lipids
  for(i in 1:nrow(lipids)) {
    
    # get sub class for decision
    lipid_subclass <- lipidomicsUtils::get_lipid_subclass(lipids$lipid[i])
    
    # create lipid info
    lipid_info <- list(lipid =lipids$lipid[i],
                      id = lipids$id[i],
                      formula = lipids$formula[i])

    # select based on lipid sub class
    if(lipid_subclass == "FA") {
      
    } else if(lipid_subclass == "NAE") {
      
    } else if(lipid_subclass == "CoA") {
      
    } else if(lipid_subclass == "PA") {
      
    } else if(lipid_subclass == "PA-P") {
      
    } else if(lipid_subclass == "PA-O") {
      
    } else if(lipid_subclass == "PC") {

      # create [M+H]+ spectrum
      pos_spectrum <- create_pos_PC_MH(lipid_info, template_file = NA)
      wormLipidBlast_Pos_Spectra <- append(wormLipidBlast_Pos_Spectra,
                                           pos_spectrum)
      
      # create [M+FA-H]- spectrum
      neg_spectrum <- create_neg_PC_MFAH(lipid_info, template_file = NA)
      wormLipidBlast_Neg_Spectra <- append(wormLipidBlast_Neg_Spectra,
                                           neg_spectrum)
      
      
      
    } else if(lipid_subclass == "PC-P") {
      
    } else if(lipid_subclass == "PC-O") {
      
    } else if(lipid_subclass == "PE") {
      
      # create [M+H]+ spectrum
      pos_spectrum <- create_pos_PE_MH(lipid_info, template_file = NA)
      wormLipidBlast_Pos_Spectra <- append(wormLipidBlast_Pos_Spectra,
                                           pos_spectrum)
      
    } else if(lipid_subclass == "PE-P") {
      
    } else if(lipid_subclass == "PE-O") {
      
    } else if(lipid_subclass == "PI") {
      
      # create [M+H]+ spectrum
      pos_spectrum <- create_pos_PI_MH(lipid_info, template_file = NA)
      wormLipidBlast_Pos_Spectra <- append(wormLipidBlast_Pos_Spectra,
                                           pos_spectrum)
      
    } else if(lipid_subclass == "PS") {
      
      # create [M+H]+ spectrum
      pos_spectrum <- create_pos_PS_MH(lipid_info, template_file = NA)
      wormLipidBlast_Pos_Spectra <- append(wormLipidBlast_Pos_Spectra,
                                           pos_spectrum)
      
    } else if(lipid_subclass == "MG") {
      
    } else if(lipid_subclass == "MG-P") {
      
    } else if(lipid_subclass == "MG-O") {
      
    } else if(lipid_subclass == "DG") {
      
    } else if(lipid_subclass == "DG-P") {
      
    } else if(lipid_subclass == "DG-O") {
      
    } else if(lipid_subclass == "TG") {
      
    } else if(lipid_subclass == "TG-P") {
      
    } else if(lipid_subclass == "TG-O") {
      
    } else if(lipid_subclass == "PG") {
      
    } else if(lipid_subclass == "PGP") {
      
    } else if(lipid_subclass == "CDPDG") {
      
    } else if(lipid_subclass == "NAPE") {
      
    } else if(lipid_subclass == "Cer") {
      
    } else if(lipid_subclass == "SM") {
      
    } else if(lipid_subclass == "HexCer") {
      
    } else if(lipid_subclass == "C1P") {
      
    } else if(lipid_subclass == "CL") {
      
    }
  }
  
  # return results
  return(list(wormLipidBlast_Neg_Spectra,
              wormLipidBlast_Pos_Spectra))
}
