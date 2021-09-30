# Test that the output table does not contain any MOCK community composition column


test_that("MOCK community columns are removed", {
  output <- MockCommunity(ASV_table_classified_raw, mock_composition, "ASV_names", choose.first=TRUE)


  expect_false(any(grep("MOCK", colnames(output), TRUE)))
})


# Test that the output table is a data.frame


test_that("Output is a data.frame", {
  output <- MockCommunity(ASV_table_classified_raw, mock_composition, "ASV_names",choose.first=TRUE)


  expect_s3_class(output, "data.frame")
})
