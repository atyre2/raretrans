## ----setup, message=FALSE, echo=-7---------------------------------------
library(tidyverse)
library(popbio)
library(raretrans)

# Raymond's theme modifications
rlt_theme <- theme(axis.title.y = element_text(colour="grey20",size=15,face="bold"),
        axis.text.x = element_text(colour="grey20",size=15, face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,face="bold"),  
        axis.title.x = element_text(colour="grey20",size=15,face="bold"))
data("L_elto")
A <- projection.matrix(as.data.frame(L_elto), 
                       stage="stage", fate="next_stage", 
                       fertility="fertility", sort=c("p","j","a"))
knitr::kable(A, digits=2, 
             caption = "Average transition matrix over all census periods and populations.")

## ----hidethis, include = FALSE-------------------------------------------
knitting <- isTRUE(getOption('knitr.in.progress'))
if(!knitting){
  pathtoimages="images"
#   # login to imgur
#   library(imguR)
#   imgurtkn <- imgur_login()
#   imgurtkn$refresh()
}
knitr::opts_chunk$set(fig.width = 5, fig.height = 5)

## ----pop_by_time, fig.cap = "Population size over time."-----------------
plotdata <- L_elto %>% group_by(POPNUM, year, stage) %>%
  tally() %>%
  group_by(POPNUM, year) %>%
  summarise(N = sum(n))

gg1 <- ggplot(plotdata, aes(x=year, y=N, group=POPNUM)) +
  geom_line() + rlt_theme
plotsumm <- plotdata %>% group_by(year) %>%
  summarize(N=sum(N))
gg2 <- ggplot(plotsumm, aes(x=year, y=N)) +
  geom_line() + rlt_theme
ggboth <- gridExtra::arrangeGrob(gg1, gg2)
grid::grid.draw(ggboth)


## ---- eval=FALSE, include = FALSE----------------------------------------
#  ggsave("pop_by_time.png", plot = ggboth, device="png", path=pathtoimages)

## ------------------------------------------------------------------------
allA <- L_elto %>% 
  mutate(stage = factor(stage, levels=c("p","j","a","m")),
         fate = factor(next_stage, levels=c("p","j","a","m"))) %>%
  group_by(POPNUM, year) %>%
  do(A = projection.matrix(as.data.frame(.), 
                           stage="stage", fate="fate", 
                           fertility="fertility", sort=c("p","j","a")),
     TF = projection.matrix(as.data.frame(.), 
                           stage="stage", fate="fate", 
                           fertility="fertility", sort=c("p","j","a"), TF = TRUE),
     N = get_state_vector(as.data.frame(.), 
                           stage="stage", sort=c("p","j","a"))) %>%
  
  filter(year < 12) # last time period all zero for obvious reason


## ----histLmbdaByYear, message=FALSE, fig.cap="histograms of asymptotic population growth for each year."----

library(popdemo)
allA <- ungroup(allA) %>% 
  mutate(A = map(A, matrix, nrow=3, ncol=3,
                 dimnames=list(c("p","j","a"),c("p","j","a")))) %>%
  mutate(lmbda = map_dbl(A, lambda),
         ergodic = map_dbl(A, isErgodic),
         irreduc = map_dbl(A, isIrreducible),
         primitv = map_dbl(A, isPrimitive))
ggplot(allA, aes(x=lmbda)) + geom_histogram(bins = 25) + facet_wrap(~year)+
 ylab("Frequency")+
  xlab(expression(paste(lambda, ", asymptotic population growth"))) + rlt_theme


## ---- eval=FALSE, include = FALSE----------------------------------------
#  ggsave("histLmbdaByYear.png", path=pathtoimages)

## ----checkErgIrr1, echo=-4-----------------------------------------------
allA <- ungroup(allA) %>% mutate(missing = map(A, ~which(.x==0)),
                                 n_missing = map_dbl(missing, length))
missing_summary <- summary(allA$n_missing)
ergo_irr <- as.data.frame(with(allA, table(ergodic, irreduc)))

knitr::kable(ergo_irr, caption = "Ergodicity and irreducibility of individual transition matrices.")

## ----exampleMatrices, echo=FALSE-----------------------------------------
#144 popn 905 year 9 10 individuals go extinct
badA144 <- allA$A[[144]]
badN144 <- allA$N[[144]]
#205 popn 914 year 6 lambda == 1 
badA205 <- allA$A[[205]]
badN205 <- allA$N[[205]]
badA71 <- allA$A[[71]]
badN71 <- allA$N[[71]]
#71 popn 250 year 5 lambda = 0.93
BaseTable <- 
  tribble(~Population, ~Year, ~Stage,     ~N,        ~seedling,    ~juvenile,   ~adult,
          905,             9, "seedling", badN144[1], badA144[1,1],badA144[1,2],badA144[1,3],
          905,             9, "juvenile", badN144[2], badA144[2,1],badA144[2,2],badA144[2,3],
          905,             9, "adult",    badN144[3], badA144[3,1],badA144[3,2],badA144[3,3],
          905,             9, "dead",     NA,         1-sum(badA144[,1]), 1-sum(badA144[,2]), 1-sum(badA144[,3]), 
          914,             6, "seedling", badN205[1], badA205[1,1],badA205[1,2],badA205[1,3],
          914,             6, "juvenile", badN205[2], badA205[2,1],badA205[2,2],badA205[2,3],
          914,             6, "adult",    badN205[3], badA205[3,1],badA205[3,2],badA205[3,3],
          914,             6, "dead",     NA,         1-sum(badA205[,1]), 1-sum(badA205[,2]), 1-sum(badA205[,3]), 
          250,             5, "seedling", badN71[1], badA71[1,1],badA71[1,2],badA71[1,3],
          250,             5, "juvenile", badN71[2], badA71[2,1],badA71[2,2],badA71[2,3],
          250,             5, "adult",    badN71[3], badA71[3,1],badA71[3,2],badA71[3,3],
          250,             5, "dead",     NA,         1-sum(badA71[,1]), 1-sum(badA71[,2]), 1-sum(badA71[,3]))


# make some space for the prior versions
RLT_Tprior <- matrix(c(0.25, 0.025, 0.0,
                       0.05, 0.9,   0.025,
                       0.01, 0.025, 0.95,
                       0.69, 0.05,  0.025), 
                     byrow = TRUE, nrow = 4, ncol = 3)
RLT_Fprior <- matrix(c(0.0, 0.0, 0.025,
                       0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0), 
                     byrow = TRUE, nrow = 3, ncol = 3)

wUnifTable <- bind_cols(BaseTable, BaseTable[,5:7], BaseTable[,5:7])
wUnifTable[1:3, 8:10] <- fill_transitions(allA$TF[[144]], allA$N[[144]])
wUnifTable[4, 8:10] <- 1-colSums(wUnifTable[1:3, 8:10])
wUnifTable[5:7, 8:10] <- fill_transitions(allA$TF[[205]], allA$N[[205]])
wUnifTable[8, 8:10] <- 1-colSums(wUnifTable[5:7, 8:10])
wUnifTable[9:11, 8:10] <- fill_transitions(allA$TF[[71]], allA$N[[71]])
wUnifTable[12, 8:10] <- 1-colSums(wUnifTable[9:11, 8:10])
wUnifTable[1:3, 11:13] <- fill_transitions(allA$TF[[144]], allA$N[[144]], P = RLT_Tprior, priorweight = -1)
wUnifTable[4, 11:13] <- 1-colSums(wUnifTable[1:3, 11:13])
wUnifTable[5:7, 11:13] <- fill_transitions(allA$TF[[205]], allA$N[[205]], P = RLT_Tprior, priorweight = -1)
wUnifTable[8, 11:13] <- 1-colSums(wUnifTable[5:7, 11:13])
wUnifTable[9:11, 11:13] <- fill_transitions(allA$TF[[71]], allA$N[[71]], P = RLT_Tprior, priorweight = -1)
wUnifTable[12, 11:13] <- 1-colSums(wUnifTable[9:11, 11:13])
knitr::kable(wUnifTable, caption = "Two problematic matrices and the matrix with the largest sample size; $\\lambda$ = 0, 1, and 0.92 for population 905, 914 and 205 respectively. The six columns on the right show the effect of assuming a uniform prior for transitions with a total weight of 1, and an RLT prior with total weight equal to 1.  Where columns do not sum to 1 the remaining probability represents mortality.", digits=3, escape = TRUE)

## ----stochGrowthRateBad, fig.cap="Stochastic population growth rates for 17 out of 23 orchid populations. "----
sgr <- allA %>% split(.$POPNUM) %>%
  map(~stoch.growth.rate(.x$A, verbose = FALSE))

# distribute results of stoch.growth.rate() into a tbl
xtrc <- function(x){
  data.frame(approx=x$approx,
                sim = x$sim,
                lowerci = x$sim.CI[1],
                upperci = x$sim.CI[2])
}
sgr <- bind_rows(map(sgr, xtrc), .id="POPNUM") 
filter(sgr, !is.na(sim)) %>%
ggplot(aes(x=approx, y = sim)) + geom_point() + 
  geom_errorbar(aes(ymin=lowerci, ymax=upperci)) + 
  geom_abline(slope = 1, intercept=0) + 
  rlt_theme +
  xlab(expression(paste("Tuljapurkar's approximate stochastic ", log(lambda)))) +
  ylab(expression(paste("Stochastic ", log(lambda), " by simulation")))

## ---- eval=FALSE, include = FALSE----------------------------------------
#  ggsave("stochGrowthRateBad.png", path=pathtoimages)

## ------------------------------------------------------------------------
# 
tmat <- allA[23,"TF"][[1]][[1]]$T ## wow! hard to get at ... 
fmat <- allA[23,"TF"][[1]][[1]]$F
N <- allA[23,"N"][[1]][[1]]
tmat + fmat

## ------------------------------------------------------------------------
isErgodic((tmat + fmat)[-2,-2])
isIrreducible((tmat + fmat)[-2,-2])
all.equal(lambda((tmat + fmat)[-2,-2]), lambda(tmat + fmat))

## ------------------------------------------------------------------------
# get matrix of observed fates, including a row for mortality
# apply(tmat, 1, function(x)x * N)
(TNmat <- L_elto %>% 
  mutate(stage = factor(stage, levels=c("p","j","a","m")),
         fate = factor(next_stage, levels=c("p","j","a","m"))) %>%
  filter(POPNUM == allA$POPNUM[23], year == allA$year[23]+1) %>%
  as.data.frame() %>% 
  xtabs(~fate + stage, data = .) %>% 
  as.matrix())
(TNmat2 <- (TNmat + 0.25)[,-4]) # re-cycles automatically, and drop last column (no one starts as m)

## ------------------------------------------------------------------------
columnsums <- colSums(TNmat2)
(tmat2 <- sweep(TNmat2, 2, columnsums, FUN = "/")[-4,])
lambda(tmat2 + fmat)
isErgodic(tmat2 + fmat)
isIrreducible(tmat2 + fmat)


## ----fillDirichlet, fig.cap = "Asymptotic growth rate using a uniform prior with a total weight of 1 vs. the asymptotic growth rate for the observed data."----
testing <- allA %>% mutate(Adirch = map2(TF, N, fill_transitions, returnType = "A"),
                           ldirch = map_dbl(Adirch, lambda),                 
         irrdirch = map_dbl(Adirch, isIrreducible),
         ergdirch = map_dbl(Adirch, isErgodic),
         TN = map2(TF, N, fill_transitions, returnType = "TN"))

ggplot(testing, aes(x=lmbda, y=ldirch)) + geom_point() + 
  geom_abline(slope = 1, intercept=0) +
  xlab(expression(paste(lambda, " ignoring zeros"))) + 
  ylab(expression(paste(lambda, " using a uniform prior"))) +
  rlt_theme

## ---- eval=FALSE, include = FALSE----------------------------------------
#  ggsave("fillDirichlet.png", path=pathtoimages)

## ---- eval=FALSE, include=FALSE------------------------------------------
#  testing %>% summarize(m_ldirch = mean(ldirch)) %>% kableExtra::kable()

## ----checkErgIrr---------------------------------------------------------
ergo_irr <- as.data.frame(with(testing, table(ergdirch, irrdirch)))

knitr::kable(ergo_irr, caption = "Ergodicity and irreducibility of individual transition matrices.")


## ----stochDirchVsBad-----------------------------------------------------
sgr_unif <- testing %>% split(.$POPNUM) %>%
  map(~stoch.growth.rate(.x$Adirch, verbose = FALSE))

sgr_unif <- bind_rows(map(sgr_unif, xtrc), .id="POPNUM") 
compare_approx <- tibble(unif = sgr_unif$approx,
                             ignored = sgr$approx)
ggplot(compare_approx, aes(x=ignored, y = unif)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept=0) + 
  xlab(expression(paste("log(",lambda,") ignoring zeros"))) + 
  ylab(expression(paste("log(",lambda,") with uniform prior"))) +
  rlt_theme

## ---- eval = FALSE, include=FALSE----------------------------------------
#  ggsave("stochDirchVsBad.png", path=pathtoimages)

## ------------------------------------------------------------------------
post_alpha = 0.00001 + 2
post_beta = 0.00001 + N[3]

## ------------------------------------------------------------------------
fmat2 <- fmat
fmat2[1,3] <- post_alpha / post_beta
fmat2

lambda(tmat + fmat2)
isErgodic(tmat+fmat2)
isIrreducible(tmat+fmat2)

