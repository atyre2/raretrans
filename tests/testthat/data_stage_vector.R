context("data integrity and returning state vector")

data("L_elto")
test_that("Data haven't changed",{
  expect_equal_to_reference(L_elto)
})
