library(wormLipidBlastR)
context("Generation of PS spectra")

################################################################################
# positive mode
################################################################################

# Std template (hardcoded) -----------------------------------------------------
test_that("Correct prediction for PS(16:0/18:1(9Z)), [M+H]+, std template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PS(16:0/18:1(9Z))",
                     id = "0001",
                     chemFormula = "C42H82NO8P")
  
  ## generate spectrum
  ps_spec <- create_pos_PS(lipid_info, adduct = "[M+H]+")
  
  ## tests
  expect_equal(round(precursorMz(ps_spec), 4), 762.5280)
  expect_equal(unlist(round(mz(ps_spec), 4)), c(239.2369, 265.2526, 577.5190,
                                                762.5280))
  
})

test_that("Correct prediction for PS(16:0/18:1(9Z)), [M+Na]+, std template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PS(16:0/18:1(9Z))",
                     id = "0001",
                     chemFormula = "C42H82NO8P")
  
  ## generate spectrum
  ps_spec <- create_pos_PS(lipid_info, adduct = "[M+Na]+")
  
  ## tests
  expect_equal(round(precursorMz(ps_spec), 4), 784.5099)
  expect_equal(unlist(round(mz(ps_spec), 4)), c(207.9981, 577.5190, 599.5010,
                                                784.5099))
  
})

################################################################################
# negative mode
################################################################################

# Std template (hardcoded) -----------------------------------------------------
test_that("Correct prediction for PS(16:0/18:1(9Z)), [M-H]-, std template", {
  
  ## precalculations
  lipid_info <- list(lipid = "PS(16:0/18:1(9Z))",
                     id = "0001",
                     chemFormula = "C39H76NO8P")
  
  ## generate spectrum
  ps_spec <- create_neg_PS(lipid_info, adduct = "[M-H]-")
  
  ## tests
  expect_equal(round(precursorMz(ps_spec), 4), 760.5134)
  expect_equal(unlist(round(mz(ps_spec), 4)), c(255.2330, 281.2486, 391.2255,
                                                409.2361, 417.2411, 435.2517,
                                                673.4814, 760.5134))
  
})