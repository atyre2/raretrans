library(dplyr)

context("fill fecundity")

data(L_elto)

test_that("data hasn't changed", {
  expect_equal_to_reference(L_elto, "one.rds")
  })

onepop <- L_elto %>% # Filter out population # 250, period (year) 5
  filter(POPNUM == 250, year == 5) %>% # redefine p for el plantón to s for seedling
  mutate(stage = case_when(stage == "p" ~ "s", TRUE ~ stage), next_stage = case_when(next_stage == "p" ~ "s", TRUE ~ next_stage))
TF <- popbio::projection.matrix(as.data.frame(onepop), stage = stage, fate = next_stage,
                                fertility = "fertility", sort = c("s", "j", "a"), TF = TRUE)

N <- get_state_vector(onepop, stage = stage, sort = c("s", "j", "a"))
test_that("args are correct", {
  expect_length(TF, 2)
  expect_type(TF, "list")
  expect_equal(N, c(11, 47, 34))
})

test_that("fill fertility behaves", {
  expect_length(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05)), 9)
  expect_type(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05)), "double")
  expect_vector(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05)))
  expect_vector(fill_fertility(TF = TF, N = N, alpha = c(NA_real_, NA_real_, 1e-05), beta = c(NA_real_, NA_real_, 1e-05)))
  expect_vector(fill_fertility(N = N, TF = TF, alpha = c(NA_real_, NA_real_, 1e-05), beta = c(NA_real_, NA_real_, 1e-05)))
  expect_vector(fill_fertility(N = N, TF, alpha = c(NA_real_, NA_real_, 1e-05), beta = c(NA_real_, NA_real_, 1e-05)))
  expect_vector(fill_fertility(TF, alpha = c(NA_real_, NA_real_, 1e-05), beta = c(NA_real_, NA_real_, 1e-05), N = N))
})

test_that("priorweight can be changed", {
  expect_vector(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05), priorweight = 10))
  expect_vector(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05), priorweight = -3))
  expect_vector(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05), priorweight = 500))
  expect_vector(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05), priorweight = 0))
})

test_that("returnType can be changed", {
  expect_type(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05), returnType=ab), "list")
  expect_vector(fill_fertility(TF, N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05), returnType=A))
})

test_that("fill fertility throws errors", {
  expect_error(fill_fertility(N, TF, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05)))
  expect_error(fill_fertility(N, TF, c(NA_real_, NA_real_, 1e-05)))
  expect_error(fill_fertility(N, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05)))
  expect_error(fill_fertility(TF, c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05)))
  expect_error(fill_fertility(TF, N))
  expect_error(fill_fertility(c(NA_real_, NA_real_, 1e-05), c(NA_real_, NA_real_, 1e-05)))
  expect_error(fill_fertility())
  expect_error(fill_fertility(c(NA_real_, NA_real_, 1e-05)))
  expect_error(fill_fertility(N = TF, TF = N, alpha = c(NA_real_, NA_real_, 1e-05), beta = c(NA_real_, NA_real_, 1e-05)))
  expect_error(fill_fertility(TF, TF, TF, TF))
  expect_error(fill_fertility(TF, N, TF, N))
  expect_error(fill_fertility(TF, N, 1, 1))
})
