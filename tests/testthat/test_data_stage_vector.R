context("data integrity and returning state vector")

data("L_elto")

test_that("Data haven't changed", {
  expect_equal_to_reference(L_elto, "one.rds")
})

library(dplyr)

testdata <- L_elto %>%
  mutate(
    stage = factor(stage, levels = c("p", "j", "a", "m")),
    fate = factor(next_stage, levels = c("p", "j", "a", "m")),
    strangestage = stage
  ) %>%
  filter(year == 1, POPNUM == 926)

test_that("State vector is correct", {
  expect_equal(get_state_vector(testdata), expected = c(1, 2, 3, 0))
  expect_equal(get_state_vector(select(testdata, -stage), stage = "strangestage"), expected = c(1, 2, 3, 0))
  expect_equal(get_state_vector(testdata, sort = c("j", "m", "a", "p")), expected = c(2, 0, 3, 1))
})

testdata_df <- as.data.frame(testdata)

test_that("State vector works with dataframes too", {
  expect_equal(get_state_vector(testdata_df), expected = c(1, 2, 3, 0))
  expect_equal(get_state_vector(select(testdata_df, -stage), stage = "strangestage"), expected = c(1, 2, 3, 0))
  expect_equal(get_state_vector(testdata, sort = c("j", "m", "a", "p")), expected = c(2, 0, 3, 1))
})

