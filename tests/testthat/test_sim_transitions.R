library(dplyr)

context("simulate transitions")

onepop <- L_elto %>% # Filter out population # 250, period (year) 5
  filter(POPNUM == 250, year == 5) %>% # redefine p for el plantÃ³n to s for seedling
  mutate(
    stage = case_when(stage == "p" ~ "s", TRUE ~ stage),
    next_stage = case_when(next_stage == "p" ~ "s", TRUE ~ next_stage)
  )

TF <- popbio::projection.matrix(as.data.frame(onepop),
  stage = stage, fate = next_stage,
  fertility = "fertility", sort = c("s", "j", "a"), TF = TRUE
)

N <- get_state_vector(onepop, stage = stage, sort = c("s", "j", "a"))

RLT_Tprior <- matrix(c(0.25, 0.025, 0, 0.05, 0.9, 0.025, 0.01, 0.025, 0.95, 0.69, 0.05, 0.025), byrow = TRUE, nrow = 4, ncol = 3)

alpha <- matrix(c(NA_real_, NA_real_, 1e-05, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_), nrow = 3, ncol = 3, byrow = TRUE)
beta <- matrix(c(NA_real_, NA_real_, 1e-05, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_), nrow = 3, ncol = 3, byrow = TRUE)
set.seed(8390278)

test_that("sim_transitions behaves", {
  expect_equal_to_reference(sim_transitions(TF, N, P = RLT_Tprior, alpha = alpha, beta = beta, priorweight = 0.5), "simtransitions.rds")
  expect_silent(sim_transitions(TF, N, P = RLT_Tprior, alpha = alpha, beta = beta, priorweight = 0.5))
})

set.seed(8390278)
test_that("# of samples can change", {
  expect_equal_to_reference(sim_transitions(TF, N, P = RLT_Tprior, alpha = alpha, beta = beta, priorweight = 0.5, samples = 25), "simtransitions2.rds")
})

test_that("warnings work", {
  expect_warning(sim_transitions(TF, N, P = RLT_Tprior, alpha = 1e-05, beta = 1e-05, priorweight = 0.5))
  expect_error(sim_transitions(TF, N, P = RLT_Tprior, alpha = N, beta = beta))
})

test_that("errors are thrown", {
  expect_error(sim_transitions())
  expect_error(sim_transitions(TF))
  expect_error(sim_transitions(N))
  expect_error(sim_transitions(TF, N, P = N))
  expect_error(sim_transitions(TF, N,
    P = RLT_Tprior, alpha = c(NA_real_, NA_real_, 1e-05, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
    beta = c(NA_real_, NA_real_, 1e-05, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
  ))
  expect_error(sim_transitions(TF, N, P = RLT_Tprior, alpha = c(NA_real_, NA_real_, 1e-05, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_), beta = beta))
  expect_error(sim_transitions(TF, N, P = RLT_Tprior, alpha = alpha, beta = c(NA_real_, NA_real_, 1e-05, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)))
  expect_error(sim_transitions(TF, N, P = matrix(c(1, 2, 3, 4), ncol = 2, nrow = 2, byrow = TRUE), alpha = alpha, beta = beta))
})
