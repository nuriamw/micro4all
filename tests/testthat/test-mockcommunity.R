# Test that the output table does not contain any MOCK community composition column


test_that("MOCK community columns are removed", {
  output <- MockCommunity(ASV_table, MOCK_table, "ASV_names")


  expect_false(any(grep("MOCK", colnames(output), TRUE)))
})


# Test that the output table is a data.frame


test_that("Output is a data.frame", {
  output <- MockCommunity(ASV_table, MOCK_table, "ASV_names")


  expect_s3_class(output, "data.frame")
})
