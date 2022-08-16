#Convert old format to new


# Packages
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)

# Import data
data <- read_csv("~/mathematicaproject/data/Set4_titation_data_old_format.csv")

data %>% select(contains("X")) -> x_vals
data %>% select(contains("Y")) -> y_vals

x_vals <- as_tibble(t(x_vals))
y_vals <- as_tibble(t(y_vals))

n_time_points <- dim(x_vals)[1]

names(x_vals) <- str_replace(names(x_vals), "V", "X")
names(y_vals) <- str_replace(names(y_vals), "V", "Y")

data_all <- bind_cols(x_vals, y_vals)
write_xlsx(data_all, "~/mathematicaproject/data/Set4_titration_data_new_format.xlsx",  col_names = TRUE)
