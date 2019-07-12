
library(dplyr)
library(raretrans)

data(L_elto)

onepop <- L_elto %>% # Filter out population # 250, period (year) 5
  filter(POPNUM == 250, year == 5) %>% # redefine p for el plantón to s for seedling
  mutate(
    stage = case_when(stage == "p" ~ "s", TRUE ~ stage),
    next_stage = case_when(next_stage == "p" ~ "s", TRUE ~ next_stage)
  )

TF <- popbio::projection.matrix(as.data.frame(onepop),
  stage = stage, fate = next_stage,
  fertility = "fertility", sort = c("s", "j", "a"), TF = TRUE
)

N <- get_state_vector(onepop, stage = stage, sort = c("s", "j", "a"))

fill_fertility(TF, N, 1e-1, 1e-05)
alpha <- matrix(c(
  0, 0, 1,
  0, 0, 0,
  0, 0, 0
), nrow = 3, ncol = 3, byrow = TRUE)
beta <- matrix(c(
  0, 0, 1,
  0, 0, 0,
  0, 0, 0
), nrow = 3, ncol = 3, byrow = TRUE)
alpha <- matrix(c(
  NA_real_, NA_real_, 1,
  NA_real_, NA_real_, NA_real_,
  NA_real_, NA_real_, NA_real_
), nrow = 3, ncol = 3, byrow = TRUE)
beta <- matrix(c(
  NA_real_, NA_real_, 1,
  NA_real_, NA_real_, NA_real_,
  NA_real_, NA_real_, NA_real_
), nrow = 3, ncol = 3, byrow = TRUE)

fill_fertility(TF, c(1, 0, 1), alpha, beta)
fill_fertility(TF, c(1, 0, 0), alpha, beta)

fill_fertility(TF, N,
  alpha = c(NA_real_, NA_real_, 1e-05, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
  beta = beta
)
TF$F
