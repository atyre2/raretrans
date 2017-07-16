context("data integrity and returning state vector")

data("L_elto")

test_that("Data haven't changed",{
  expect_equal_to_reference(L_elto)
})

testdata <- L_elto %>%
  mutate(stage = factor(stage, levels=c("p","j","a","m")),
         fate = factor(next_stage, levels=c("p","j","a","m"))) %>%
  filter(year == 1, POPNUM == 926)
test_that("State vector is correct",{
  expect_equal(get_state_vector(testdata))
}
