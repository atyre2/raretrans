## ----setup, message = FALSE----------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy = TRUE)
library(tidyverse)
library(popbio) # for projection.matrix()
library(raretrans)
# Raymond's theme modifications
rlt_theme <- theme(axis.title.y = element_text(colour="grey20",size=15,face="bold"),
        axis.text.x = element_text(colour="grey20",size=15, face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,face="bold"),  
        axis.title.x = element_text(colour="grey20",size=15,face="bold"))

## ----viewData------------------------------------------------------------
data("L_elto") # load the dataset `L_elto` into memory (from package `raretrans`)
head(L_elto) 

## ----mungeData-----------------------------------------------------------
onepop <- L_elto %>%   
# Filter out population # 250, period (year) 5
  filter(POPNUM == 250, year == 5) %>% 
  # redefine p for el plantón to s for seedling
  mutate(stage = case_when(stage == "p" ~ "s",
                           TRUE ~ stage),
         next_stage = case_when(next_stage == "p"~ "s",
                                TRUE ~ next_stage))


# popbio::projection.matrix doesn't like tibbles
# set TF = TRUE to get the right format.
TF <- popbio::projection.matrix(as.data.frame(onepop), 
                        stage = stage, fate = next_stage, 
                        fertility="fertility", sort=c("s","j","a"), TF = TRUE)
TF # This is the stage-structured population model


## ------------------------------------------------------------------------
N <- get_state_vector(onepop, stage = stage, sort=c("s","j","a")) 
N # A vector of # of starting individuals for each stage

## ---- echo=FALSE---------------------------------------------------------
Tmat <- L_elto %>% 
  mutate(stage = factor(stage, levels=c("p","j","a","m")),
         fate = factor(next_stage, levels=c("p","j","a","m"))) %>%
  filter(POPNUM == 231, year == 2) %>%
  as.data.frame() %>% 
  xtabs(~fate + stage, data = .) %>% 
  as.matrix()
Tmat <- Tmat[,-4] #throw away the column for transitions to death
N2 <- colSums(Tmat) #get the total number ... CHECK should be before 86
Tmat <- sweep(Tmat[-4,], 1, N2, "/") #normalize to 1, discarding transitions FROM death
# figure out how much reproduction happened in year 2 by looking at year 3
# L_elto %>% 
#   mutate(stage = factor(stage, levels=c("p","j","a","m")),
#          fate = factor(next_stage, levels=c("p","j","a","m"))) %>%
#   filter(POPNUM == 231, year == 3, stage == "p") %>%
#   count(stage) 
# 2 offspring from N[3] == 16 adults
Fmat <- matrix(0, nrow=3, ncol=3) # create a matrix full of zeros
Fmat[1,3] <- 2/16 # counted 16 adults in this stage, and 2 seedlings in next year
TF2 <- list(Tmat = Tmat, Fmat = Fmat)

## ------------------------------------------------------------------------
TF2
N2

## ----posterior1----------------------------------------------------------
Tprior <- matrix(0.25, byrow = TRUE, ncol = 3, nrow=4)
fill_transitions(TF, N, P = Tprior) # returns transition matrix by default

## ----manualTransitionPost------------------------------------------------
Tobs <- sweep(TF$T, 2, N, "*") # get the observed transitions
Tobs <- rbind(Tobs, N - colSums(Tobs)) # add the mortality row
Tobs <- Tobs + 0.25 # add the prior
sweep(Tobs, 2, colSums(Tobs), "/")[-4,] # divide by the column sum, and discard mortalityrow

