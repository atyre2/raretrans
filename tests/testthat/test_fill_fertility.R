library(dplyr)

context("fill fertility")

data(L_elto)

onepop <- L_elto %>% # Filter out population # 250, period (year) 5
  filter(POPNUM == 250, year == 5) %>% # redefine p for el plantón to s for seedling
  mutate(stage = case_when(stage == "p" ~ "s", TRUE ~ stage),
         next_stage = case_when(next_stage == "p" ~ "s", TRUE ~ next_stage))

TF <- popbio::projection.matrix(as.data.frame(onepop), stage = stage, fate = next_stage,
                                fertility = "fertility", sort = c("s", "j", "a"), TF = TRUE)

N <- get_state_vector(onepop, stage = stage, sort = c("s", "j", "a"))


alpha <- matrix(c(NA_real_, NA_real_, 1,
                  NA_real_, NA_real_, NA_real_,
                  NA_real_, NA_real_, NA_real_), nrow=3, ncol = 3, byrow = TRUE)
beta <- matrix(c(NA_real_, NA_real_, 1,
                 NA_real_, NA_real_, NA_real_,
                 NA_real_, NA_real_, NA_real_), nrow=3, ncol = 3, byrow = TRUE)


test_that("args are correct", {
  expect_length(TF, 2)
  expect_type(TF, "list")
  expect_equal(N, c(11, 47, 34))
})

test_that("fill fertility behaves", {
  expect_length(fill_fertility(TF, N,  alpha,  beta), 9)
  expect_type(fill_fertility(TF, N, alpha,  beta), "double")
  expect_vector(fill_fertility(TF, N, alpha,  beta))
  expect_vector(fill_fertility(TF = TF, N = N, alpha = alpha, beta = beta))
  expect_vector(fill_fertility(N = N, TF = TF, alpha = alpha, beta = beta))
  expect_vector(fill_fertility(N = N, TF, alpha = alpha, beta = beta))
  expect_vector(fill_fertility(TF, alpha = alpha, beta = beta, N = N))
})

test_that("priorweight can be changed", {
  expect_vector(fill_fertility(TF, N, alpha, beta, priorweight = 10))
  expect_vector(fill_fertility(TF, N, alpha, beta, priorweight = -3))
  expect_vector(fill_fertility(TF, N, alpha, beta, priorweight = 500))
  expect_vector(fill_fertility(TF, N, alpha, beta, priorweight = 0))
})

test_that("some N == 0 not a problem", {
  expect_vector(fill_fertility(TF, N = c(0,1,1), alpha, beta))
  expect_vector(fill_fertility(TF, N = c(1,0,1), alpha, beta))
  expect_vector(fill_fertility(TF, N = c(1,1,0), alpha, beta))
})

test_that("returnType can be changed", {
  expect_type(fill_fertility(TF, N, alpha, beta, returnType="ab"), "list")
  expect_vector(fill_fertility(TF, N, alpha, beta, returnType="A"))
  expect_error(fill_fertility(TF, N, alpha, beta, returnType=""))
})

test_that("fill fertility throws errors and warnings with invalid arguments", {
  expect_error(fill_fertility(N, TF, alpha, beta))
  expect_error(fill_fertility(N, TF, 1e-05))
  expect_error(fill_fertility(N, alpha, beta))
  expect_error(fill_fertility(TF, alpha, beta))
  expect_warning(fill_fertility(TF, N))
  expect_error(fill_fertility(alpha, beta))
  expect_error(fill_fertility())
  expect_error(fill_fertility(1e-05))
  expect_error(fill_fertility(N = TF, TF = N, alpha = 1e-05, beta = 1e-05))
  expect_error(fill_fertility(TF, TF, TF, TF))
  expect_error(fill_fertility(TF, N, TF, N))
  expect_warning(fill_fertility(TF, N, 1, 1))
})
