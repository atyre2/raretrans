# load and process raw data on Leopinthes elto
library(tidyverse)
df <- read_csv("data-raw/L.elto.george-30-08-2004_Drew.csv")
newnames <- toupper(gsub("[^a-zA-Z0-9_]", "", names(df)))
# newnames[5:6] <- c("NUMLEA1","NUMINF1")
newnames[1] <- "POPNUM"
newnames[2] <- "IND_NUM"
names(df) <- newnames
df_stages <- select(df, 1:2, matches("STAGE0[1-9]+")) %>%
  gather(key = "year", value = "stage", matches("STAGE[0-9]+"))
df_stages$year <- as.numeric(regmatches(df_stages$year, regexpr("[0-9]+", df_stages$year)))
df_stages <- group_by(df_stages, POPNUM, IND_NUM) %>%
  mutate(next_stage = lead(stage)) %>%
  ungroup()

# First, identify the individuals that have problematic m -> *
bad_m <- filter(df_stages, stage == "m", !is.na(next_stage)) %>%
  select(POPNUM, IND_NUM)
df_2 <- anti_join(df_stages, bad_m)

# Now clean out the excess non-entries with NA -> NA
df_3 <- filter(df_2, !(is.na(stage) & is.na(next_stage)))

# Next, identify the missing 'm' entries, and discard the "recruit" observations (NA -> *),
# and the mortality row observations (m->NA). Don't forget that IND_NUM not unique.
df_4 <- df_3 %>%
  filter(!(stage == "m" & is.na(next_stage))) %>%
  filter(!(is.na(stage) & next_stage %in% c("a", "j", "p", "m"))) %>%
  group_by(POPNUM, IND_NUM) %>%
  mutate(
    first_year = first(year),
    last_year = last(year),
    recruited = ifelse(year == first_year, TRUE, FALSE),
    died = ifelse(year == last_year, year, NA_real_),
    lifespan = n()
  ) %>%
  mutate(next_stage = ifelse(is.na(died), next_stage, "m")) %>%
  filter(lifespan == (last_year - first_year + 1))

# first calculate fertility as next year's seedlings divided by this year's adults.
# this should probably be done by population and year
L_elto <- df_4 %>%
  group_by(POPNUM, year) %>%
  summarize(
    seedlings = sum(if_else(stage == "p" & recruited, 1, 0)),
    adults = sum(if_else(stage == "a", 1, 0))
  ) %>%
  ungroup() %>%
  group_by(POPNUM) %>%
  mutate(fertility = lead(seedlings) / adults) %>%
  right_join(df_4, by = c("POPNUM", "year")) %>%
  mutate(fertility = if_else(stage == "a", fertility, 0)) %>% # this fixes the few divide by zero issues
  ungroup()

# push it out to data/
usethis::use_data(L_elto, overwrite = TRUE)

# if changed deliberately, push to tests/testthat/one.rds as well
saveRDS(L_elto, "tests/testthat/one.rds")
