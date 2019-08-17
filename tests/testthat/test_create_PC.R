library(wormLipidBlastR)
library(MSnbase)
context("Generation of PC spectra")

################################################################################
# positive mode
################################################################################

# Std template (hardcoded) -----------------------------------------------------
test_that("Correct prediction for PC(16:0/18:1(9Z)), [M+H]+, std template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PC(16:0/18:1(9Z))",
                    id = "0001",
                    chemFormula = "C42H82NO8P")
  
  ## generate spectrum
  pc_spec <- create_pos_PC_MH(lipid_info)
  
  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 760.5851)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(184.0733, 760.5851))
  
})

test_that("Correct prediction for PC(16:0/18:1(9Z)), [M+Na]+, std template", {

  ## precalculations
  lipid_info <- list(lipid = "PC(16:0/18:1(9Z))",
                    id = "0001",
                    chemFormula = "C42H82NO8P")

  ## generate spectrum
  pc_spec <- create_pos_PC_MNa(lipid_info)

  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 782.5670)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(146.9818, 184.0733, 239.2369, 
                                             265.2526, 441.2376, 467.2533,
                                             478.3292, 500.3111, 504.3449,
                                             526.3268, 599.5010, 723.4935,
                                             782.5670))

})

################################################################################
# negative mode
################################################################################

# Std template (hardcoded) -----------------------------------------------------
test_that("Correct prediction for PC(16:0/18:1(9Z)), [M+FA-H]-, std template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PC(16:0/18:1(9Z))",
                    id = "0001",
                    chemFormula = "C42H82NO8P")
  
  ## generate spectrum
  pc_spec <- create_neg_PC_MFAH(lipid_info)
  
  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 804.5760)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(255.2330, 281.2486, 462.2990,
                                             480.3096, 488.3146, 506.3252,
                                             744.5549, 804.5760))
})

test_that("Correct prediction for PC(16:0/18:1(9Z)), [M+HAc-H]-, std template", {

  ## precalculations
  lipid_info <- list(lipid = "PC(16:0/18:1(9Z))",
                    id = "0001",
                    chemFormula = "C42H82NO8P")

  ## generate spectrum
  pc_spec <- create_neg_PC_MHAcH(lipid_info)

  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 818.5917)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(255.2330, 281.2486, 462.2990,
                                             480.3096, 488.3146, 506.3252,
                                             744.5549, 818.5917))
})

# user defined template --------------------------------------------------------
test_that("Correct prediction for PC(16:0/18:1(9Z)), [M+FA-H]-, user-def template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PC(16:0/18:1(9Z))",
                     id = "0001",
                     chemFormula = "C42H82NO8P")
  
  ## create template with only one entry
  template <- list("adduct_mass" = 999)
  
  ## generate spectrum
  pc_spec <- create_neg_PC_MFAH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 804.5760)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(804.5760))
  expect_equal(intensity(pc_spec[[1]]), c(999))
  
  ## create template with several entries
  template <- list("adduct_mass" = 999,
                   "adduct_mass - rcdk::get.formula('C2H4O2')@mass" = 100)
  
  ## generate spectrum
  pc_spec <- create_neg_PC_MFAH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 804.5760)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(744.5549, 804.5760))
  expect_equal(intensity(pc_spec[[1]]), c(100, 999))
  
  ## create template with several entries
  template <- list("adduct_mass" = 999,
                   "adduct_mass - rcdk::get.formula('C2H4O2')@mass" = 100,
                   "sn1_mass - proton_mass" = 999,
                   "sn2_mass - proton_mass" = 999)
  
  ## generate spectrum
  pc_spec <- create_neg_PC_MFAH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 804.5760)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(255.2330, 281.2486, 744.5549,
                                             804.5760))
  expect_equal(intensity(pc_spec[[1]]), c(999, 999, 100, 999))
})

test_that("Correct prediction for PC(16:0/18:1(9Z)), [M+HAc-H]-, std template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PC(16:0/18:1(9Z))",
                     id = "0001",
                     chemFormula = "C42H82NO8P")
  
  ## create template with only one entry
  template <- list("adduct_mass" = 999)
  
  ## generate spectrum
  pc_spec <- create_neg_PC_MHAcH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 818.5917)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(818.5917))
  expect_equal(intensity(pc_spec[[1]]), c(999))
  
  ## create template with several entries
  template <- list("adduct_mass" = 999,
                   "adduct_mass - rcdk::get.formula('C3H6O2')@mass" = 100)
  
  ## generate spectrum
  pc_spec <- create_neg_PC_MHAcH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 818.5917)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(744.5549, 818.5917))
  expect_equal(intensity(pc_spec[[1]]), c(100, 999))
  
  ## create template with several entries
  template <- list("adduct_mass" = 999,
                   "adduct_mass - rcdk::get.formula('C3H6O2')@mass" = 100,
                   "sn1_mass - proton_mass" = 999,
                   "sn2_mass - proton_mass" = 999)
  
  ## generate spectrum
  pc_spec <- create_neg_PC_MHAcH(lipid_info, template = template)
  
  ## tests
  expect_equal(round(precursorMz(pc_spec[[1]]), 4), 818.5917)
  expect_equal(round(mz(pc_spec[[1]]), 4), c(255.2330, 281.2486, 744.5549,
                                             818.5917))
  expect_equal(intensity(pc_spec[[1]]), c(999, 999, 100, 999))
})