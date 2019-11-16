
context("sim transitions")

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



