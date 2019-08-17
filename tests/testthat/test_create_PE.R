library(wormLipidBlastR)
library(MSnbase)
context("Generation of PE spectra")

################################################################################
# positive mode
################################################################################

# Std template (hardcoded) -----------------------------------------------------
test_that("Correct prediction for PE(16:0/18:1(9Z)), [M+H]+, std template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PE(16:0/18:1(9Z))",
                    id = "0001",
                    chemFormula = "C39H76NO8P")
  
  ## generate spectrum
  pe_spec <- create_pos_PE_MH(lipid_info)
  
  ## tests
  expect_equal(round(precursorMz(pe_spec[[1]]), 4), 718.5381)
  expect_equal(round(mz(pe_spec[[1]]), 4), c(239.2369, 265.2526, 313.2737,
                                             339.2894, 577.5190, 718.5381))
  
})

test_that("Correct prediction for PE(16:0/18:1(9Z)), [M+Na]+", {

  ## precalculations
  lipid_info <- list(lipid = "PE(16:0/18:1(9Z))",
                    id = "0001",
                    chemFormula = "C39H76NO8P")

  ## generate spectrum
  pe_spec <- create_pos_PE_MNa(lipid_info)

  ## tests
  expect_equal(round(precursorMz(pe_spec[[1]]), 4), 740.5201)
  expect_equal(round(mz(pe_spec[[1]]), 4), c(164.0083, 415.2220, 441.2376,
                                             577.5190, 599.5010, 697.4779,
                                             740.5201))

})

################################################################################
# negative mode
################################################################################

# Std template (hardcoded) -----------------------------------------------------
test_that("Correct prediction for PE(16:0/18:1(9Z)), [M-H]-, std template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PE(16:0/18:1(9Z))",
                    id = "0001",
                    chemFormula = "C39H76NO8P")
  
  ## generate spectrum
  pe_spec <- create_neg_PE_MH(lipid_info)
  
  ## tests
  expect_equal(round(precursorMz(pe_spec[[1]]), 4), 716.5236)
  expect_equal(round(mz(pe_spec[[1]]), 4), c(122.0013, 140.0118, 255.2330,
                                             281.2486, 434.2677, 452.2783,
                                             460.2833, 478.2939, 716.5236))
  
})

# user defined template --------------------------------------------------------
# Std template (hardcoded) -----------------------------------------------------
test_that("Correct prediction for PE(16:0/18:1(9Z)), [M-H]-, user-def template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PE(16:0/18:1(9Z))",
                     id = "0001",
                     chemFormula = "C39H76NO8P")
  
  ## create template with only one entry
  template <- list("adduct_mass" = 999)
  
  ## generate spectrum
  pe_spec <- create_neg_PE_MH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pe_spec[[1]]), 4), 716.5236)
  expect_equal(round(mz(pe_spec[[1]]), 4), c(716.5236))
  expect_equal(intensity(pe_spec[[1]]), c(999))
  
  ## create template with several entries
  template <- list("adduct_mass" = 999,
                   "adduct_mass - sn1_mass + water_mass" = 100,
                   "adduct_mass - sn2_mass + water_mass" = 100)
  
  ## generate spectrum
  pe_spec <- create_neg_PE_MH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pe_spec[[1]]), 4), 716.5236)
  expect_equal(round(mz(pe_spec[[1]]), 4), c(452.2783, 478.2939, 716.5236))
  expect_equal(intensity(pe_spec[[1]]), c(100, 100, 999))
  
  ## create template with several entries
  template <- list("adduct_mass" = 999,
                   "adduct_mass - sn1_mass + water_mass" = 100,
                   "adduct_mass - sn2_mass + water_mass" = 100,
                   "sn1_mass - proton_mass" = 999,
                   "sn2_mass - proton_mass" = 999)
  
  ## generate spectrum
  pe_spec <- create_neg_PE_MH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pe_spec[[1]]), 4), 716.5236)
  expect_equal(round(mz(pe_spec[[1]]), 4), c( 255.2330, 281.2486, 452.2783,
                                              478.2939, 716.5236))
  expect_equal(intensity(pe_spec[[1]]), c(999, 999, 100, 100, 999))
})