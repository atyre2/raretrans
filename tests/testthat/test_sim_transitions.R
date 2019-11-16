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

RLT_Tprior <- matrix(c(0.9474, 0.0, 0.0, 0,
                       0.0211, 0.9371, 0,   0.0,
                       0.0,     0,      0.658,0.1830,
                       0,    0.0429,0.3364,0.8157,
                       0.0315, 0.02, 0.0056, 0.0013),
                     byrow = TRUE, nrow = 5, ncol = 4)
N=c(39,	317,	178,	170)
Trans=matrix(c(
  0.717948718,	0,	0,	0,
  0.256410256,	0.911671924,	0,	0,
  0,	0.009463722,	0.735955056,	0.241176471,
  0,	0.044164038,	0.230337079,	0.752941176), byrow=TRUE,ncol=4)

Fert=matrix(c(
  0.0000000, 0.0000000, 0,	0.111764706,
  0.0000000, 0.0000000, 0, 0.207865169,
  0.0000000, 0.0000000,0,  0.0000000,
  0,0,0,0), byrow=TRUE,ncol=4)


TF = list(T = Trans, F = Fert)

alpha2 <- matrix(c(NA_real_, NA_real_, NA_real_,.001,
                   NA_real_, NA_real_, NA_real_,.001,
                   NA_real_, NA_real_, NA_real_, NA_real_,
                   NA_real_,NA_real_,NA_real_,NA_real_), nrow=4, ncol = 4, byrow = TRUE)
beta2 <- matrix(c(NA_real_, NA_real_, NA_real_,0.001,
                  NA_real_, NA_real_, NA_real_,.001,
                  NA_real_, NA_real_, NA_real_,NA_real_,
                  NA_real_,NA_real_,NA_real_,NA_real_), nrow=4, ncol = 4, byrow = TRUE)

# generate 13 matrices based on the prior transitions and fertilities plus the data
RLT_0.5 <- sim_transitions(TF, N, P = RLT_Tprior,
                           alpha = alpha2, beta = beta2, priorweight = 0,
                           samples = 13)

test_that("return is correct length and type, with no missing values", {
  expect_length(RLT_0.5, 13)
  expect_type(RLT_0.5, "list")
  expect_true(all(!is.na(RLT_0.5)))
  expect_true(all(!is.null(RLT_0.5)))
})




