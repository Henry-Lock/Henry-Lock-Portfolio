# Henry Lock Project 1
library(ggplot2)
library(dplyr)
library(readr)

# Get all the data
setwd("/Users/henrylock/Library/Mobile Documents/com~apple~CloudDocs/R Stuff")
data <- read_table("Land_and_Ocean_summary.txt", 
                   comment = "%",
                   col_names = c("Year", "Anomaly_Air", "Unc_Air", 
                                 "FiveYr_Air", "FiveYrUnc_Air", 
                                 "Anomaly_Water", "Unc_Water", 
                                 "FiveYr_Water", "FiveYrUnc_Water"),
                   na = "NaN")

# Deal with NA values
data[data == "NaN"] <- NA
data <- type.convert(data, as.is = TRUE)

# Land (Air temp) plot
land_plot <- ggplot(data, aes(x = Year)) +
  geom_line(aes(y = Anomaly_Air), color = "red", alpha = 0.6) +
  geom_ribbon(aes(ymin = Anomaly_Air - Unc_Air, ymax = Anomaly_Air + Unc_Air), fill = "red", alpha = 0.2) +
  geom_line(aes(y = FiveYr_Air), color = "darkred", linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1850, 2050, by = 20)) +
  labs(title = "Land Average Temperature",
       y = "Temperature Anomaly (°C)", x = NULL,
        x = "Year") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# Ocean (Water temp) plot
ocean_plot <- ggplot(data, aes(x = Year)) +
  geom_line(aes(y = Anomaly_Water), color = "blue", alpha = 0.6) +
  geom_ribbon(aes(ymin = Anomaly_Water - Unc_Water, ymax = Anomaly_Water + Unc_Water), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = FiveYr_Water), color = "darkblue", linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1850, 2050, by = 20)) +
  labs(title = "Ocean Average Temperature",
       y = "Temperature Anomaly (°C)", x = "Year",
       x = "Year") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

land_plot
ocean_plot

##########################################

# Country by country average temps


# Column names
col_names <- c("Year", "Month", "Monthly_Anomaly", "Monthly_Uncertainty", 
                       "Annual_Anomaly", "Annual_Uncertainty", 
                       "FiveYr", "FiveYr_Unc", 
                       "TenYr", "TenYr_Unc", 
                       "TwentyYr", "TwentyYr_Unc")

# Denmark
denmark <- read_table("denmark-(europe)-TAVG-Trend.txt", 
                      skip = 626, 
                      na = "NaN",
                      col_names = col_names)
denmark_avg <- denmark |>
  filter(!is.na(FiveYr)) |>
  group_by(Year) |>
  summarise(
    FiveYr = mean(FiveYr, na.rm = TRUE)
  ) |>
  ungroup()

# UK
uk <- read_table("united-kingdom-TAVG-Trend.txt", 
                 skip = 627, 
                 na = "NaN",
                 col_names = col_names)

uk_avg <- uk |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  filter(!is.na(FiveYr)) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  ungroup()

# Chile
chile <- read_table("chile-TAVG-Trend.txt", skip = 626, na = "NaN", col_names = col_names)
chile_avg <- chile |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  ungroup()

# Canada
canada <- read_table("canada-TAVG-Trend.txt", skip = 70, na = "NaN", col_names = col_names)
canada_avg <- canada |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  ungroup()

# Morocco
morocco <- read_table("morocco-TAVG-Trend.txt", skip = 626, na = "NaN", col_names = col_names)
morocco_avg <- morocco |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  ungroup()

# Plot 5-year average for each country
best_countries_plot <- ggplot() +
  geom_line(data = denmark_avg |> mutate(Country = "Denmark"), aes(x = Year, y = FiveYr, color = Country), linewidth = 1.2) +
  geom_line(data = uk_avg |> mutate(Country = "United Kingdom"), aes(x = Year, y = FiveYr, color = Country), linewidth = 1.2) +
  geom_line(data = chile_avg |> mutate(Country = "Chile"), aes(x = Year, y = FiveYr, color = Country), linewidth = 1.2) +
  geom_line(data = canada_avg |> mutate(Country = "Canada"), aes(x = Year, y = FiveYr, color = Country), linewidth = 1.2) +
  geom_line(data = morocco_avg |> mutate(Country = "Morocco"), aes(x = Year, y = FiveYr, color = Country), linewidth = 1.2) +
  scale_color_manual(
    name = "Country",
    values = c(
      "Denmark" = "#E41A1C",
      "United Kingdom" = "#377EB8",
      "Chile" = "#984EA3",
      "Canada" = "#FF7F00",
      "Morocco" = "#4DAF4A"
    )
  ) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in 5 Countries Leading Climate Policy",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right")

best_countries_plot

#############################

# Iran
iran <- read_table("iran-TAVG-Trend.txt", skip = 626, na = "NaN", col_names = col_names)
iran_avg <- iran |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  mutate(Country = "Iran")

# Libya
libya <- read_table("libya-TAVG-Trend.txt", skip = 626, na = "NaN", col_names = col_names)
libya_avg <- libya |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  mutate(Country = "Libya")

# Saudi Arabia
saudi <- read_table("saudi-arabia-TAVG-Trend.txt", skip = 626, na = "NaN", col_names = col_names)
saudi_avg <- saudi |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  mutate(Country = "Saudi Arabia")

# United States
us <- read_table("united-states-TAVG-Trend.txt", skip = 70, na = "NaN", col_names = col_names)
us_avg <- us |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  mutate(Country = "United States")

# Yemen
yemen <- read_table("yemen-TAVG-Trend.txt", skip = 626, na = "NaN", col_names = col_names)
yemen_avg <- yemen |>
  filter(!is.na(FiveYr), Year >= 1850) |>
  group_by(Year) |>
  summarise(FiveYr = mean(FiveYr, na.rm = TRUE)) |>
  mutate(Country = "Yemen")

# Easier in one data
combined_data <- bind_rows(iran_avg, libya_avg, saudi_avg, us_avg, yemen_avg)

# Country colors
country_colors <- c(
  "Iran" = "#99000D",
  "Libya" = "#08519C",
  "Saudi Arabia" = "#54278F",
  "United States" = "#B35806",
  "Yemen" = "#006D2C"
)

# Plot worst countries
worst_countries_plot <- ggplot(combined_data, aes(x = Year, y = FiveYr, color = Country)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = country_colors, name = "Country") +
  scale_x_continuous(
    breaks = c(seq(1850, 2029, by = 20), 2030),
    limits = c(1850, 2030)
  ) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in 5 Countries with the worst Climate Policy",
    y = "Temperature Anomaly (°C)",
    x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

worst_countries_plot

###############################

# Ploting each individual country

denmark_plot <- ggplot(
  data = denmark_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#E41A1C", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in Denmark",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
denmark_plot

UK_plot <- ggplot(
  data = uk_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#377EB8", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in the UK",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
UK_plot

chile_plot <- ggplot(
  data = chile_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#984EA3", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in Chile",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
chile_plot

canada_plot <- ggplot(
  data = canada_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#FF7F00", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in Canada",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
canada_plot

morocco_plot <- ggplot(
  data = morocco_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#4DAF4A", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in Morocco",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
morocco_plot

###

iran_plot <- ggplot(
  data = iran_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#99000D", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in Iran",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
iran_plot

libya_plot <- ggplot(
  data = libya_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#08519C", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in Libya",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
libya_plot

saudi_plot <- ggplot(
  data = saudi_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#54278F", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in Saudi Arabia",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
saudi_plot

us_plot <- ggplot(
  data = us_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#B35806", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in the US",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
us_plot

yemen_plot <- ggplot(
  data = yemen_avg,
  aes(x = Year, y = FiveYr)
) +
  geom_line(color = "#006D2C", linewidth = 1.2) +
  scale_x_continuous(breaks = c(seq(1850, 2029, by = 20), 2030),
                     limits = c(1850, 2030)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(
    title = "Temperature Averages in Yemen",
    y = "Temperature Anomaly (°C)", x = "Year"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )
yemen_plot

#####

#Average warming table to plot bar graph

warming_df <- tibble(
  Country = c("Denmark", "United Kingdom", "Chile", "Canada", "Morocco",
              "Iran", "Libya", "United States", "Saudi Arabia", "Yemen"),
  Warming = c(2.86, 2.30, 0.96, 3.64, 2.71, 3.52, 2.54, 2.36, 2.50, 2.58),
  SE = c(0.23, 0.26, 0.37, 0.24, 0.23, 0.37, 0.63, 0.13, 0.55, 0.50),
  Group = c(rep("Climate Policy Leaders", 5), rep("Worst Climate Policy", 5))
) 
View(warming_df)

warming_plot <- ggplot(warming_df, aes(x = Country, y = Warming, fill = Group)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = Warming - SE, ymax = Warming + SE),
                width = 0.2) +
  facet_wrap(~Group, scales = "free_x", nrow = 1) +
  labs(
    title = "Historical Warming by Country",
    y = "Warming Since 1960 (°C / Century)",
    x = "Country"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
warming_plot
