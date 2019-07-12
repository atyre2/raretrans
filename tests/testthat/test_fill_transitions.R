library(dplyr)

context("fill transitions")

data(L_elto)

test_that("data hasn't changed", {
  expect_equal_to_reference(L_elto, "one.rds")
})

# setting up args
onepop <- L_elto %>% # Filter out population # 250, period (year) 5
  filter(POPNUM == 250, year == 5) %>% # redefine p for el plantón to s for seedling
  mutate(stage = case_when(stage == "p" ~ "s", TRUE ~ stage), next_stage = case_when(next_stage == "p" ~ "s", TRUE ~ next_stage))
TF <- popbio::projection.matrix(as.data.frame(onepop),
  stage = stage, fate = next_stage,
  fertility = "fertility", sort = c("s", "j", "a"), TF = TRUE
)
N <- get_state_vector(onepop, stage = stage, sort = c("s", "j", "a"))
Tprior <- matrix(0.25, byrow = TRUE, ncol = 3, nrow = 4)
RLT_Tprior <- matrix(c(0.25, 0.025, 0, 0.05, 0.9, 0.025, 0.01, 0.025, 0.95, 0.69, 0.05, 0.025), byrow = TRUE, nrow = 4, ncol = 3)

test_that("args are correct", {
  expect_equal_to_reference(TF, "TF.rds")
  expect_equal(N, c(11, 47, 34))
})

test_that("fill_transitions behaves", {
  expect_length(fill_transitions(TF, N, P = Tprior), 9)
  expect_length(fill_transitions(TF, N), 9)
  expect_equal_to_reference(fill_transitions(TF, N, P = Tprior), "transmatrix.rds")
  expect_equal_to_reference(fill_transitions(TF, N, P = RLT_Tprior, priorweight = 0.5), "transmatrix2.rds")
})

test_that("priorweight can be changed", {
  expect_vector(fill_transitions(TF, N, P = Tprior, priorweight = 1))
  expect_vector(fill_transitions(TF, N, P = Tprior, priorweight = 531))
  expect_vector(fill_transitions(TF, N, P = Tprior, priorweight = -531))
})

test_that("returnType can be changed", {
  expect_vector(fill_transitions(TF, N, P = Tprior, returnType = "A"))
  expect_vector(fill_transitions(TF, N, P = Tprior, returnType = "TN"))
})

test_that("fill_transitions throws errors", {
  expect_error(fill_transitions(N = N, P = Tprior))
  expect_error(fill_transitions(N))
  expect_error(fill_transitions())
  expect_error(fill_transitions(N, N))
  expect_error(fill_transitions(TF, N, P = matrix(-.25, byrow = TRUE, ncol = 4, nrow = 4)))
  expect_error(fill_transitions(TF, N, returnType = "W"))
})
