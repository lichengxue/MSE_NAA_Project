#### Loading library ####
setwd("C:/Users/chengxue.li/Desktop/MSE_manuscript")
library(wham)
library(whamMSE)
library(ggplot2)
library(dplyr)
library(tidyr)

name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

#### Functions ####
create_plot <- function(data, use_n_years = 30, boxplot = FALSE,
                        ggtitle_text = "Model Performance", 
                        file_name = "Model_performance",
                        plot_width = 12, plot_height = 10, 
                        point_size = 3, line_thickness = 0.8, 
                        x_text_size = 12, y_text_size = 12,
                        axis_title_size = 15, strip_text_size = 12) {
  
  # Set up the ggplot object
  plot <- ggplot(data_summary, aes(x = factor(EM), y = Median, col = factor(EM), fill = factor(EM))) +
    # Facet with parsed labels
    facet_grid(Index ~ OM, scales = "free", labeller = label_parsed) +
    # Background shading for specific columns (1, 3, 5, 7, 9, 10)
    geom_rect(data = data.frame(xmin = c(0.5, 2.5, 4.5, 6.5, 8.5),
                                xmax = c(1.5, 3.5, 5.5, 7.5, 9.5)),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "lightgrey", alpha = 0.4, inherit.aes = FALSE) +
    scale_colour_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    scale_x_discrete(labels = labels) +
    labs(x = "Estimation Model", y = "") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = x_text_size),
      axis.text.y = element_text(size = y_text_size),
      axis.title = element_text(face = "bold", size = axis_title_size),
      strip.text = element_text(face = "bold", size = strip_text_size),
      legend.position = "none",
      title = element_text(face = "bold")
    ) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    ggtitle(ggtitle_text)
  
  # Add elements based on the boxplot parameter
  if (boxplot) {
    # Boxplot style: Box segment for IQR, median as horizontal line, whiskers as vertical lines
    plot <- plot +
      # Box for the IQR
      geom_rect(aes(xmin = as.numeric(EM) - 0.4, xmax = as.numeric(EM) + 0.4,
                    ymin = Q1, ymax = Q3, col = factor(EM)), fill = NA, size = line_thickness) +
      # Median line
      geom_segment(aes(x = as.numeric(EM) - 0.4, xend = as.numeric(EM) + 0.4,
                       y = Median, yend = Median, col = factor(EM)), size = line_thickness) +
      # Lower whisker (from min to Q1)
      geom_segment(aes(x = as.numeric(EM), xend = as.numeric(EM),
                       y = min_val, yend = Q1, col = factor(EM)), size = line_thickness) +
      # Upper whisker (from Q3 to max)
      geom_segment(aes(x = as.numeric(EM), xend = as.numeric(EM),
                       y = Q3, yend = max_val, col = factor(EM)), size = line_thickness)
  } else {
    # Point and IQR error bar style: Median as point and IQR as error bar
    plot <- plot +
      geom_point(aes(y = Median), size = point_size, shape = 19) +
      geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, size = line_thickness)
  }
  
  # Print the plot independently
  print(plot)
  
  # Save the plot
  ggsave(paste0(file_name, ".png"), 
         plot, width = plot_width, height = plot_height)
}

#### Diagnostics ####

# Define para.
n_oms = 12
n_ems = 1:10

# Load data
data <- read.csv("Parameter_Est.csv")

# Calculate the median and IQR
data_summary <- data %>%
  group_by(OM, Group, EM, Var) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75)
  ) %>%
  mutate(EM = factor(EM))

data_summary[data_summary$Median == 1,'Median'] = NA
data_summary[data_summary$Q1 == 1,'Q1'] = NA
data_summary[data_summary$Q3 == 1,'Q3'] = NA

# Create a named vector of expressions for the x-axis labels with subscripts
labels <- c(
  `1` = bquote(bold(PAN[NAA])),
  `2` = bquote(bold(PAN[noNAA])),
  `3` = bquote(bold(FAA[NAA])),
  `4` = bquote(bold(FAA[noNAA])),
  `5` = bquote(bold(SEP[NAA])),
  `6` = bquote(bold(SEP[noNAA])),
  `7` = bquote(bold(SpD[NAA])),
  `8` = bquote(bold(SpD[noNAA])),
  `9` = bquote(bold(SpE[NAA * "," * Est])),
  `10` = bquote(bold(SpE[NAA * "," * Fix]))
)


# Define labels for each subset of OMs
OM_labels_list <- list(
  c("OM 1" = bquote(bold(OM[high]~(F[A]))), "OM 2" = bquote(bold(OM[high]~(F[U]))), "OM 3" = bquote(bold(OM[high]~(F[C])))),
  c("OM 4" = bquote(bold(OM[low]~(F[A]))), "OM 5" = bquote(bold(OM[low]~(F[U]))), "OM 6" = bquote(bold(OM[low]~(F[C])))),
  c("OM 7" = bquote(bold(OM[high]~(F[A]))), "OM 8" = bquote(bold(OM[high]~(F[U]))), "OM 9" = bquote(bold(OM[high]~(F[C])))),
  c("OM 10" = bquote(bold(OM[low]~(F[A]))), "OM 11" = bquote(bold(OM[low]~(F[U]))), "OM 12" = bquote(bold(OM[low]~(F[C])))
  )
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- data_summary %>% filter(OM %in% names(OM_labels_list[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels_list[[i]]), labels = OM_labels_list[[i]])
  
  # Create the plot for the current subset
  p <- ggplot(om_data, aes(x = EM, y = Median, fill = Var, group = Var)) + # Add group = Var to dodge by each Var level
    # Background shading for specific columns (1, 3, 5, 7, 9, 10)
    geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 1
    geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 3
    geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 5
    geom_rect(aes(xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 7
    geom_rect(aes(xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 9
    geom_rect(aes(xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 10
    # Plot points and error bars with larger dodge width
    geom_point(position = position_dodge(width = 0.9), size = 3, shape = 19) + # Increase dodge width
    geom_errorbar(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.9), width = 0) + # Increase dodge width
    facet_wrap(OM ~ Group, scales = "free_y", ncol = 3, labeller = custom_labeller) +
    labs(
      title = paste("Model Parameter Estimates"),
      x = "Estimation Model", y = ""
    ) +
    scale_x_discrete(labels = labels) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      title = element_text(face = "bold")
    )
  
  # Print the plot independently
  print(p)
  
  # Save the plot independently with a unique filename
  ggsave(filename = paste0("Parameter_Estimates_OM_", i, ".png"), plot = p, width = 10, height = 7, units = "in")
}

#### 30 year performance ####
use.n.years <- 30
data <- read.csv(file = paste0("Variables_of_Interest_", use.n.years, "_years_relative.csv"))

OM_labels <- list(
  c("1" = bquote(bold(OM[high]~(F[A]))), "2" = bquote(bold(OM[high]~(F[U]))), "3" = bquote(bold(OM[high]~(F[C])))),
  c("4" = bquote(bold(OM[low]~(F[A]))), "5" = bquote(bold(OM[low]~(F[U]))), "6" = bquote(bold(OM[low]~(F[C])))),
  c("7" = bquote(bold(OM[high]~(F[A]))), "8" = bquote(bold(OM[high]~(F[U]))), "9" = bquote(bold(OM[high]~(F[C])))),
  c("10" = bquote(bold(OM[low]~(F[A]))), "11" = bquote(bold(OM[low]~(F[U]))), "12" = bquote(bold(OM[low]~(F[C])))
  )
)

labels <- c(
  `1` = bquote(bold(PAN[NAA])),
  `2` = bquote(bold(PAN[noNAA])),
  `3` = bquote(bold(FAA[NAA])),
  `4` = bquote(bold(FAA[noNAA])),
  `5` = bquote(bold(SEP[NAA])),
  `6` = bquote(bold(SEP[noNAA])),
  `7` = bquote(bold(SpD[NAA])),
  `8` = bquote(bold(SpD[noNAA])),
  `9` = bquote(bold(SpE[NAA * "," * Est])),
  `10` = bquote(bold(SpE[NAA * "," * Fix]))
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  
  filename = paste0("Performance_last_", use.n.years, "_years_relative_OM_", i)
  
  # Filter data for the current subset of OMs
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Filter to remove EM 10 and drop unused levels
  tmp <- tmp %>%
    filter(EM != 10) %>%
    droplevels()
  
  # Define custom colors
  # my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ADD8E6",
  #                "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
                 "#E74C3C", "#E74C3C", "#8E44AD")
  
  # Summarize the data
  data_summary <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  # Rename Index directly in data_summary with parsed expressions using bquote-like format
  data_summary <- data_summary %>%
    mutate(
      Index = case_when(
        Index == "Catch_R.1" ~ "bold(C[R1])",
        Index == "Catch_R.2" ~ "bold(C[R2])",
        Index == "Catch_Total" ~ "bold(C[G])",
        Index == "F_R.1" ~ "bold(F[R1])",
        Index == "F_R.2" ~ "bold(F[R2])",
        Index == "F_Total" ~ "bold(F[G])",
        Index == "SSB_S.1" ~ "bold(SSB[R1])",
        Index == "SSB_S.2" ~ "bold(SSB[R2])",
        Index == "SSB_Total" ~ "bold(SSB[G])"
      )
    )
  
  # Reorder Index to the desired level order
  data_summary <- data_summary %>%
    mutate(
      Index = factor(
        Index,
        levels = c("bold(C[R1])", "bold(C[R2])", "bold(C[G])",
                   "bold(F[R1])", "bold(F[R2])", "bold(F[G])",
                   "bold(SSB[R1])", "bold(SSB[R2])", "bold(SSB[G])")
      )
    )
  
  create_plot(data_summary, boxplot = F, file_name = filename)
}

#### Last 5 performance ####
use.n.years <- 5
data <- read.csv(file = paste0("Variables_of_Interest_", use.n.years, "_years_relative.csv"))
OM_labels <- list(
  c("1" = bquote(bold(OM[high]~(F[A]))), "2" = bquote(bold(OM[high]~(F[U]))), "3" = bquote(bold(OM[high]~(F[C])))),
  c("4" = bquote(bold(OM[low]~(F[A]))), "5" = bquote(bold(OM[low]~(F[U]))), "6" = bquote(bold(OM[low]~(F[C])))),
  c("7" = bquote(bold(OM[high]~(F[A]))), "8" = bquote(bold(OM[high]~(F[U]))), "9" = bquote(bold(OM[high]~(F[C])))),
  c("10" = bquote(bold(OM[low]~(F[A]))), "11" = bquote(bold(OM[low]~(F[U]))), "12" = bquote(bold(OM[low]~(F[C])))
  )
)

labels <- c(
  `1` = bquote(bold(PAN[NAA])),
  `2` = bquote(bold(PAN[noNAA])),
  `3` = bquote(bold(FAA[NAA])),
  `4` = bquote(bold(FAA[noNAA])),
  `5` = bquote(bold(SEP[NAA])),
  `6` = bquote(bold(SEP[noNAA])),
  `7` = bquote(bold(SpD[NAA])),
  `8` = bquote(bold(SpD[noNAA])),
  `9` = bquote(bold(SpE[NAA * "," * Est])),
  `10` = bquote(bold(SpE[NAA * "," * Fix]))
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  
  filename = paste0("Performance_last_", use.n.years, "_years_relative_OM_", i)
  
  # Filter data for the current subset of OMs
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Filter to remove EM 10 and drop unused levels
  tmp <- tmp %>%
    filter(EM != 10) %>%
    droplevels()
  
  # Define custom colors
  my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ADD8E6",
                 "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  
  my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
                 "#E74C3C", "#E74C3C", "#8E44AD")
  # Summarize the data
  data_summary <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  # Rename Index directly in data_summary with parsed expressions using bquote-like format
  data_summary <- data_summary %>%
    mutate(
      Index = case_when(
        Index == "Catch_R.1" ~ "bold(C[R1])",
        Index == "Catch_R.2" ~ "bold(C[R2])",
        Index == "Catch_Total" ~ "bold(C[G])",
        Index == "F_R.1" ~ "bold(F[R1])",
        Index == "F_R.2" ~ "bold(F[R2])",
        Index == "F_Total" ~ "bold(F[G])",
        Index == "SSB_S.1" ~ "bold(SSB[R1])",
        Index == "SSB_S.2" ~ "bold(SSB[R2])",
        Index == "SSB_Total" ~ "bold(SSB[G])"
      )
    )
  
  # Reorder Index to the desired level order
  data_summary <- data_summary %>%
    mutate(
      Index = factor(
        Index,
        levels = c("bold(C[R1])", "bold(C[R2])", "bold(C[G])",
                   "bold(F[R1])", "bold(F[R2])", "bold(F[G])",
                   "bold(SSB[R1])", "bold(SSB[R2])", "bold(SSB[G])")
      )
    )
  
  create_plot(data_summary, boxplot = F, file_name = filename)
}

#### First 5 performance ####
use.n.years <- 5
data <- read.csv(file = paste0("Variables_of_Interest_first_", use.n.years, "_years_relative.csv"))

OM_labels <- list(
  c("1" = bquote(bold(OM[high]~(F[A]))), "2" = bquote(bold(OM[high]~(F[U]))), "3" = bquote(bold(OM[high]~(F[C])))),
  c("4" = bquote(bold(OM[low]~(F[A]))), "5" = bquote(bold(OM[low]~(F[U]))), "6" = bquote(bold(OM[low]~(F[C])))),
  c("7" = bquote(bold(OM[high]~(F[A]))), "8" = bquote(bold(OM[high]~(F[U]))), "9" = bquote(bold(OM[high]~(F[C])))),
  c("10" = bquote(bold(OM[low]~(F[A]))), "11" = bquote(bold(OM[low]~(F[U]))), "12" = bquote(bold(OM[low]~(F[C])))
  )
)

labels <- c(
  `1` = bquote(bold(PAN[NAA])),
  `2` = bquote(bold(PAN[noNAA])),
  `3` = bquote(bold(FAA[NAA])),
  `4` = bquote(bold(FAA[noNAA])),
  `5` = bquote(bold(SEP[NAA])),
  `6` = bquote(bold(SEP[noNAA])),
  `7` = bquote(bold(SpD[NAA])),
  `8` = bquote(bold(SpD[noNAA])),
  `9` = bquote(bold(SpE[NAA * "," * Est])),
  `10` = bquote(bold(SpE[NAA * "," * Fix]))
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  
  filename = paste0("Performance_first_", use.n.years, "_years_relative_OM_", i)
  
  # Filter data for the current subset of OMs
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Filter to remove EM 10 and drop unused levels
  tmp <- tmp %>%
    filter(EM != 10) %>%
    droplevels()
  
  # Define custom colors
  my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ADD8E6",
                 "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  
  my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
                 "#E74C3C", "#E74C3C", "#8E44AD")
  
  # Summarize the data
  data_summary <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  # Rename Index directly in data_summary with parsed expressions using bquote-like format
  data_summary <- data_summary %>%
    mutate(
      Index = case_when(
        Index == "Catch_R.1" ~ "bold(C[R1])",
        Index == "Catch_R.2" ~ "bold(C[R2])",
        Index == "Catch_Total" ~ "bold(C[G])",
        Index == "F_R.1" ~ "bold(F[R1])",
        Index == "F_R.2" ~ "bold(F[R2])",
        Index == "F_Total" ~ "bold(F[G])",
        Index == "SSB_S.1" ~ "bold(SSB[R1])",
        Index == "SSB_S.2" ~ "bold(SSB[R2])",
        Index == "SSB_Total" ~ "bold(SSB[G])"
      )
    )
  
  # Reorder Index to the desired level order
  data_summary <- data_summary %>%
    mutate(
      Index = factor(
        Index,
        levels = c("bold(C[R1])", "bold(C[R2])", "bold(C[G])",
                   "bold(F[R1])", "bold(F[R2])", "bold(F[G])",
                   "bold(SSB[R1])", "bold(SSB[R2])", "bold(SSB[G])")
      )
    )
  
  create_plot(data_summary, boxplot = F, file_name = filename)
}

#### Pop Status ####
my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3","#ADD8E6",
               "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3","#000000")

my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
               "#E74C3C", "#E74C3C", "#8E44AD", "#000000")
# Calculate median and IQR over nsims for each EM and OM combination
calculate_medians_iqr <- function(temp) {
  temp %>%
    group_by(OM, EM, Index) %>%
    summarise(
      Median_Overfished = median(Overfished, na.rm = TRUE),
      Median_Overfishing = median(Overfishing, na.rm = TRUE),
      IQR_Overfished_Lower = quantile(Overfished, 0.25, na.rm = TRUE),
      IQR_Overfished_Upper = quantile(Overfished, 0.75, na.rm = TRUE),
      IQR_Overfishing_Lower = quantile(Overfishing, 0.25, na.rm = TRUE),
      IQR_Overfishing_Upper = quantile(Overfishing, 0.75, na.rm = TRUE)
    ) %>%
    ungroup()
}

temp <- read.csv("Population_status.csv")

# Corrected facet labels
facet_labels_index <- c(
  "Index_4" = "bold(P[R1])",
  "Index_5" = "bold(P[R2])",
  "Index_6" = "bold(P[G])"
)

custom_labeller <- labeller(
  Index = as_labeller(facet_labels_index, label_parsed),
  OM = label_parsed
)

# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(temp)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

ylim = c(0.6, 1.8)
xlim = c(0.4, 1.9)

xmin_rect <- max(0.5, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "Global",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "Global","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB[paste(.(40),"%")]))) +
    ylab(bquote(bold(F/F[paste(.(40),"%")]))) +
    ggtitle("Stock Status") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    xlim(xlim) + 
    ylim(ylim) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_OM_", i, ".png"), p, width = 15, height = 15, dpi = 150)
}

#### Regional Stock Status ####
temp <- read.csv("Population_status.csv")
# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(temp)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

ylim = c(0.2, 1.1)
xlim = c(0.3, 1.2)

xmin_rect <- max(0.5, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "S1",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "S1","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB[paste(.(40),"%")]))) +
    ylab(bquote(bold(F/F[paste(.(40),"%")]))) +
    ggtitle("Stock Status") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    xlim(xlim) + 
    ylim(ylim) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_S1_OM_", i, ".png"), p, width = 15, height = 15, dpi = 150)
}

#### Regional Stock Status ####
temp <- read.csv("Population_status.csv")
# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(temp)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

ylim = c(0.2, 1.1)
xlim = c(0.2, 1.1)

xmin_rect <- max(0.5, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "S2",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "S2","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB[paste(.(40),"%")]))) +
    ylab(bquote(bold(F/F[paste(.(40),"%")]))) +
    ggtitle("Stock Status") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    xlim(xlim) + 
    ylim(ylim) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_S2_OM_", i, ".png"), p, width = 15, height = 15, dpi = 150)
}

#### Catch variation relative ####
data <- read.csv("Catch_AAV_relative.csv")

OM_labels <- list(
  c("1" = bquote(bold(OM[high]~(F[A]))), "2" = bquote(bold(OM[high]~(F[U]))), "3" = bquote(bold(OM[high]~(F[C])))),
  c("4" = bquote(bold(OM[low]~(F[A]))), "5" = bquote(bold(OM[low]~(F[U]))), "6" = bquote(bold(OM[low]~(F[C])))),
  c("7" = bquote(bold(OM[high]~(F[A]))), "8" = bquote(bold(OM[high]~(F[U]))), "9" = bquote(bold(OM[high]~(F[C])))),
  c("10" = bquote(bold(OM[low]~(F[A]))), "11" = bquote(bold(OM[low]~(F[U]))), "12" = bquote(bold(OM[low]~(F[C])))
  ))

labels <- c(
  `1` = bquote(bold(PAN[NAA])),
  `2` = bquote(bold(PAN[noNAA])),
  `3` = bquote(bold(FAA[NAA])),
  `4` = bquote(bold(FAA[noNAA])),
  `5` = bquote(bold(SEP[NAA])),
  `6` = bquote(bold(SEP[noNAA])),
  `7` = bquote(bold(SpD[NAA])),
  `8` = bquote(bold(SpD[noNAA])),
  `9` = bquote(bold(SpE[NAA * "," * Est])),
  `10` = bquote(bold(SpE[NAA * "," * Fix]))
)

for (i in 1:4) {
  
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Filter to remove EM 10 and drop unused levels
  tmp <- tmp %>%
    filter(EM != 10) %>%
    droplevels()
  
  tmp <- tmp %>%
    mutate(
      Index = case_when(
        Index == "Catch.1" ~ "bold(C[R1])",
        Index == "Catch.2" ~ "bold(C[R2])",
        Index == "Catch.G" ~ "bold(C[G])"
      )
    )
  
  # Reorder Index to the desired level order
  tmp <- tmp %>%
    mutate(
      Index = factor(
        Index,
        levels = c("bold(C[R1])", "bold(C[R2])", "bold(C[G])")
      )
    )
  
  # Summarize the data
  tmp <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  plot <- ggplot(tmp, aes(x = factor(EM), y = Median, col = factor(EM), fill = factor(EM))) + 
    geom_point(size = 3, shape = 19) + 
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, size = 0.8) + 
    facet_grid(Index ~ OM, scales = "free", labeller = label_parsed) + 
    ggtitle(paste0("Average Annual Catch Variation (Relative Difference)")) + 
    geom_rect(data = data.frame(xmin = c(0.5, 2.5, 4.5, 6.5, 8.5),
                                xmax = c(1.5, 3.5, 5.5, 7.5, 9.5)),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "lightgrey", alpha = 0.4, inherit.aes = FALSE) + 
    scale_colour_manual(values = my_colors) + 
    scale_fill_manual(values = my_colors) + 
    scale_x_discrete(labels = labels) + 
    labs(x = "Estimation Model", y = "") + 
    theme_bw() +  # Set the background to white
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 15),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      title = element_text(face = "bold")
    ) + 
    geom_hline(yintercept = 0, col = "red", linetype = "dashed")
  
  print(plot)
  ggsave(paste0("Catch_variation_OM",i,"_relative.png"), plot, width = 10, height = 7, dpi = 150)
}

#### Probability overfishing #### 
data <- read.csv("Probability_Overfished_Overfishing_relative.csv")

data <- data %>% filter(Index == "Index_4" | Index == "Index_5" | Index == "Index_6")

# Corrected facet labels
facet_labels_index <- c(
  "Index_4" = "bold(P[R1])",
  "Index_5" = "bold(P[R2])",
  "Index_6" = "bold(P[G])"
)

for (i in 1:4) {
  
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Custom labeller function
  custom_labeller <- labeller(
    Index = as_labeller(facet_labels_index, label_parsed),
    OM = label_parsed
  )
  
  tmp <- tmp[tmp$EM != 10,]
  
  tmp <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  plot = ggplot(tmp, aes(x = factor(EM), y = Median, col = factor(EM), fill = factor(EM))) + 
    geom_point(size = 3, shape = 19) + 
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, size = 0.8) + 
    facet_grid(Index ~ OM, scales = "free", labeller = custom_labeller) + 
    ggtitle(paste0("Probability of Overfishing")) + 
    geom_rect(data = data.frame(xmin = c(0.5, 2.5, 4.5, 6.5, 8.5),
                                xmax = c(1.5, 3.5, 5.5, 7.5, 9.5)),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "lightgrey", alpha = 0.4, inherit.aes = FALSE) + 
    scale_colour_manual(values = my_colors) + 
    scale_fill_manual(values = my_colors) + 
    scale_x_discrete(labels = labels) + 
    labs(x = "Estimation Model", y = "") + 
    theme_bw() +  # Set the background to white
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 15),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      title = element_text(face = "bold")
    ) + 
    geom_hline(yintercept = 0, col = "red", linetype = "dashed")
  print(plot)
  ggsave(paste0("Probability_Overfishing_OM",i,".png"), plot, width = 10, height = 7, dpi = 150)
}

#### Trajectory Catch and SSB ####

data <- read.csv("quantities_trajectory.csv")
library(ggplot2)
library(dplyr)

sum_data <- data %>%
  group_by(OM, EM, Index, Years) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.4),
    Q3 = quantile(Value, 0.6)
  )

my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3","#ADD8E6",
               "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3","#000000")
my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
               "#E74C3C", "#E74C3C", "#8E44AD", "#000000")
index.name = unique(sum_data$Index)

sum_data <- sum_data[sum_data$EM != 10,]
sum_data <- sum_data[sum_data$Years >= 2023,]
sum_data$Years <- sum_data$Years-2022

sum_data$EM <- factor(sum_data$EM, levels = 1:10, labels = name_map)

sum_data[sum_data$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est","EM5_Fix"),"Group"] = "NAA"
sum_data[is.na(sum_data$Group),"Group"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# sum_data <- sum_data %>%
#   mutate(
#     Index = case_when(
#       Index == "Catch_1" ~ "bold(C[R1])",
#       Index == "Catch_2" ~ "bold(C[R2])",
#       Index == "Catch_G" ~ "bold(C[G])",
#       Index == "SSB_1" ~ "bold(SSB[S1])",
#       Index == "SSB_2" ~ "bold(SSB[S2])",
#       Index == "SSB_G" ~ "bold(SSB[G])"
#     )
#   )

for (name in unique(sum_data$Index)) {
  
  for (i in 1:4) {
    
    summary_data <- filter(sum_data, Index == name)
    
    summary_data <- summary_data %>% filter(OM %in% names(OM_labels[[i]]))
    
    # Set OM factor levels and labels for the current subset
    summary_data$OM <- factor(summary_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
    
    plot <- ggplot(summary_data, aes(x = Years, y = Median, color = as.factor(EM))) + 
      geom_line(size = 1) + 
      scale_color_manual(values = my_colors, labels = legend_labels) +  # Use legend_labels here
      facet_grid(OM ~ Group, scales = "free", labeller = custom_labeller) + 
      geom_ribbon(aes(ymin = Q1, ymax = Q3, fill = as.factor(EM)), color = NA, alpha = 0.1) + 
      labs(x = "Years", y = "", color = "EM", fill = "EM") + 
      ggtitle(name) + 
      theme_bw() + 
      theme(
        axis.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),
        legend.key.width = unit(0.5, "cm"),  # Increase width of legend keys
        legend.key.height = unit(1, "cm")
      ) + 
      guides(fill = "none") +  # Remove the fill legend
      geom_hline(yintercept = 0, col = "red", linetype = "dashed")
    print(plot)
    
    ggsave(paste0(name,"_OM",i,".png"), plot, width = 10, height = 10, dpi = 150) 
  }
}

#### Simulated Movement rate ####
data <- read.csv("move_demo.csv")

# Calculate mean and standard deviation for each scenario at each time step
summary_data <- data %>%
  group_by(OM, Movement, Years) %>%
  summarize(
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE)
  )

custom_labeller <- labeller(
  OM = label_parsed  # Only apply label_parsed to OM
)

summary_data$Years = summary_data$Years - 1992


for (i in 1:4) {
  
  tmp <- summary_data %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Plot the data with ggplot2
  plot <- ggplot(tmp, aes(x = Years, y = mean_value, color = as.factor(Movement))) +
    geom_line(size = 0.5) +  # Plot the mean line
    geom_point()+
    facet_wrap(~OM,ncol = 1, labeller = custom_labeller) + 
    #geom_ribbon(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value, fill = as.factor(Movement)), alpha = 0.2) +  # Add shaded area for SD
    labs(x = "Years", y = "", color = "Movement", fill = "Movement") +
    scale_colour_manual(values = my_colors) +
    ggtitle("Simulated Movement") +
    theme_bw() +
    theme(axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"))
  print(plot)
  ggsave(paste0("Movement",i,".png"), plot, width = 7.5, height = 7.5, dpi = 150)
}

#### Holistic Performance ####

# Catch and SSB last 10 years
data1 <- read.csv("Variables_of_Interest_5_years_relative.csv")

# Specify the indices you want to extract
indices_to_extract <- c("Catch_Total", "SSB_Total")

# Extract the specific indices
data1 <- data1 %>%
  filter(Index %in% indices_to_extract)
data1$Index <- ifelse(data1$Index == "Catch_Total","Long_term_Catch","Long_term_SSB")

data2 <- read.csv("Variables_of_Interest_first_5_years_relative.csv")
data2 <- data2 %>%
  filter(Index %in% indices_to_extract)
data2$Index <- ifelse(data2$Index == "Catch_Total","Short_term_Catch","Short_term_SSB")

data3 <- read.csv("Terminal_pop_status_overfished_relative.csv")
indices_to_extract <- c("Overfished_Global")
data3 <- data3 %>%
  filter(Index %in% indices_to_extract)
data3$Index <- ifelse(data3$Index == "Overfishing_Global","F/F40%","SSB/SSB40%")

data4 <- read.csv("Terminal_pop_status_overfishing_relative.csv")
indices_to_extract <- c("Overfishing_Global")
data4 <- data4 %>%
  filter(Index %in% indices_to_extract)
data4$Index <- ifelse(data4$Index == "Overfishing_Global","F/F40%","SSB/SSB40%")
data4$Value <- -data4$Value

data5 <- read.csv("Probability_overfished_overfishing_relative.csv")
indices_to_extract <- c("Index_3","Index_6")
data5 <- data5 %>%
  filter(Index %in% indices_to_extract)
data5$Index <- ifelse(data5$Index == "Index_3","P_noOverfished","P_noOverfishing")

data6 <- read.csv("Catch_AAV_relative.csv")
indices_to_extract <- c("Catch.G")
data6 <- data6 %>%
  filter(Index %in% indices_to_extract)
data6$Index <- ifelse(data6$Index == "Catch.G","AAV_Catch","AAV_Catch")
data6$Value <- -data6$Value

data <- rbind(data1, data2, data3, data4, data5, data6)
data <- data[data$EM != 10,]
data$EM <- factor(data$EM, levels = 1:9, labels = name_map[1:9])

scale_values <- function(x) {
  if (min(x) == max(x)) {
    return(rep(1, length(x)))
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}

data$nsim <- as.factor(data$nsim)
data$Value = data$Value * 100

sum_data <- data %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

sum_data <- sum_data %>%
  group_by(OM, EM, Index) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75)
  )

sum_data$EM = factor(sum_data$EM, levels = c("EM5_Est",
                                             "EM4_NAA",
                                             "EM3_NAA",
                                             "EM2_NAA",
                                             "EM1_NAA",
                                             "EM4_noNAA",
                                             "EM3_noNAA",
                                             "EM2_noNAA",
                                             "EM1_noNAA"))

position_dodge_val = position_dodge(width = 0.9)
custom_labeller <- labeller(
  OM = label_parsed
)

sum_data <- sum_data[sum_data$Index != 'P_noOverfished',]

sum_data$Index = factor(sum_data$Index, levels = c("Short_term_Catch",
                                                   "Short_term_SSB",
                                                   "Long_term_Catch",
                                                   "Long_term_SSB",
                                                   "P_noOverfishing",
                                                   "F/F40%",
                                                   "SSB/SSB40%",
                                                   "AAV_Catch"))

legend_labels <- c(
  "Short_term_Catch" = expression(Catch[short]),
  "Short_term_SSB" = expression(SSB[short]),
  "Long_term_Catch" = expression(Catch[long]),
  "Long_term_SSB" = expression(SSB[long]),
  "P_noOverfishing" = expression(P[overfishing]),
  "F/F40%" = expression(F / F[40 * "% (T.)"]),
  "SSB/SSB40%" = expression(SSB / SSB[40 * "% (T.)"]),
  "AAV_Catch" = expression(paste("Catch ", "Variation"))
)

OM_labels <- list(
  c("1" = bquote(bold(OM[high]~(F[A]))), "2" = bquote(bold(OM[high]~(F[U]))), "3" = bquote(bold(OM[high]~(F[C])))),
  c("4" = bquote(bold(OM[low]~(F[A]))), "5" = bquote(bold(OM[low]~(F[U]))), "6" = bquote(bold(OM[low]~(F[C])))),
  c("7" = bquote(bold(OM[high]~(F[A]))), "8" = bquote(bold(OM[high]~(F[U]))), "9" = bquote(bold(OM[high]~(F[C])))),
  c("10" = bquote(bold(OM[low]~(F[A]))), "11" = bquote(bold(OM[low]~(F[U]))), "12" = bquote(bold(OM[low]~(F[C])))
  ))

for (i in 1:4) {
  summary_data <- sum_data %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  summary_data$OM <- factor(summary_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ADD8E6",
                 "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  
  my_colors <- c("#2874A6", "#F39C12", "#27AE60", "#E74C3C", 
                 "#8E44AD", "#17A589", "#F4D03F", "#34495E")
  
  # Create bar graph with IQR
  plot <- ggplot(summary_data, aes(x = Median, y = EM, fill = Index, color = Index)) +
    geom_bar(stat = "identity", position = position_dodge_val, width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3),
                  position = position_dodge_val, width = 0.2, size = 0.5, alpha = 0.5) +  # IQR error bars
    facet_wrap(~OM, labeller = custom_labeller) +
    scale_fill_manual(values = my_colors, labels = legend_labels) +
    scale_color_manual(values = my_colors, labels = legend_labels) +
    scale_y_discrete(labels = c(
      "EM1_NAA" = expression(PAN[NAA]),
      "EM1_noNAA" = expression(PAN[noNAA]),
      "EM2_NAA" = expression(FAA[NAA]),
      "EM2_noNAA" = expression(FAA[noNAA]),
      "EM3_NAA" = expression(SEP[NAA]),
      "EM3_noNAA" = expression(SEP[noNAA]),
      "EM4_NAA" = expression(SpD[NAA]),
      "EM4_noNAA" = expression(SpD[noNAA]),
      "EM5_Est" = expression(SpE[NAA * "," * Est])
    )) +
    labs(title = "Holistic Model Performance", x = "Score", y = "Estimation Model", fill = "Index", color = "Index") +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      title = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"),
      legend.key.height = unit(2, "cm")
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 1.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 2.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 3.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 4.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 5.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 6.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 7.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 8.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 9.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 10.5, linetype = "dashed", size = 0.8) + 
    geom_rect(data = data.frame(ymin = c(0.5),
                                ymax = c(6.5)),
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
              fill = "lightgrey", alpha = 0.3, inherit.aes = FALSE) 
  
  print(plot)
  ggsave(paste0("Holistic_performance_OM_", i, ".png"), plot, width = 15, height = 15, dpi = 150)
}

for (i in 1:4) {
  summary_data <- sum_data %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  summary_data$OM <- factor(summary_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  my_colors <- c("#2874A6", "#F39C12", "#27AE60", "#E74C3C", 
                 "#8E44AD", "#17A589", "#F4D03F", "#34495E")
  
  # Create bar graph with IQR
  plot <- ggplot(summary_data, aes(x = Median, y = EM, fill = Index, color = Index)) +
    geom_bar(stat = "identity", position = position_dodge_val, width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3),
                  position = position_dodge_val, width = 0.2, size = 0.5, alpha = 0.5) +  # IQR error bars
    facet_wrap(~OM, labeller = custom_labeller) +
    scale_fill_manual(values = my_colors, labels = legend_labels) +
    scale_color_manual(values = my_colors, labels = legend_labels) +
    scale_y_discrete(labels = c(
      "EM1_NAA" = expression(PAN[NAA]),
      "EM1_noNAA" = expression(PAN[noNAA]),
      "EM2_NAA" = expression(FAA[NAA]),
      "EM2_noNAA" = expression(FAA[noNAA]),
      "EM3_NAA" = expression(SEP[NAA]),
      "EM3_noNAA" = expression(SEP[noNAA]),
      "EM4_NAA" = expression(SpD[NAA]),
      "EM4_noNAA" = expression(SpD[noNAA]),
      "EM5_Est" = expression(SpE[NAA * "," * Est])
    )) +
    labs(title = "Holistic Model Performance", x = "Score", y = "Estimation Model") +
    guides(fill = guide_legend(title = "Metric"),
           color = guide_legend(title = "Metric")) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      title = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"),
      legend.key.height = unit(2, "cm")
    ) +
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 1.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 2.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 3.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 4.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 5.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 6.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 7.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 8.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 9.5, linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 10.5, linetype = "dashed", size = 0.8) + 
    geom_rect(data = data.frame(ymin = c(0.5),
                                ymax = c(6.5)),
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
              fill = "lightgrey", alpha = 0.3, inherit.aes = FALSE) 
  
  print(plot)
  ggsave(paste0("Holistic_performance_OM_", i, ".png"), plot, width = 15, height = 15, dpi = 150)
}

n_oms = 12
# KOBE Plot
res1 <- NULL
res2 <- NULL

#### use a new method for the F reference points ####
for (om in 1:n_oms) {
  mods <- readRDS(paste0("OM_", om, ".RDS"))
  print(om)
  for (i in 1:length(mods)){
    for (j in 1:length(mods[[1]])){
      
      if(j <= 4){
        true <- mods[[i]][[j]]$om$rep$SSB[11:30,]
        true <- cbind(true,rowSums(true))
        tmp <- mods[[i]][[j]]$em_list[[1]]$SSB
        ssb.re = tmp/true[,3]-1
        ssb.re = cbind(mtrx = matrix(NA,20,2),ssb.re)
        
        true <- cbind(mods[[i]][[j]]$om$rep$NAA[1,1,11:30,1],mods[[i]][[j]]$om$rep$NAA[2,2,11:30,1])
        true <- cbind(true,rowSums(true))
        tmp <- matrix(mods[[i]][[j]]$em_list[[1]]$NAA[1,1,,1])
        rec.re = tmp/true[,3]-1
        rec.re = cbind(mtrx = matrix(NA,20,2),rec.re)
      } 
      
      if(j > 4 && j <= 6) {
        true <- mods[[i]][[j]]$om$rep$SSB[11:30,]
        true <- cbind(true,rowSums(true))
        tmp1 <- mods[[i]][[j]]$em_list[[1]][[1]]$SSB
        tmp2 <- mods[[i]][[j]]$em_list[[1]][[2]]$SSB
        tmp <- cbind(tmp1,tmp2,tmp1+tmp2)
        ssb.re = tmp/true-1
        
        true <- cbind(mods[[i]][[j]]$om$rep$NAA[1,1,11:30,1],mods[[i]][[j]]$om$rep$NAA[2,2,11:30,1])
        true <- cbind(true,rowSums(true))
        tmp1 <- matrix(mods[[i]][[j]]$em_list[[1]][[1]]$NAA[1,1,,1])
        tmp2 <- matrix(mods[[i]][[j]]$em_list[[1]][[2]]$NAA[1,1,,1])
        tmp <- cbind(tmp1,tmp2,tmp1+tmp2)
        rec.re = tmp/true-1
      }
      
      if(j > 6) {
        true <- mods[[i]][[j]]$om$rep$SSB[11:30,]
        true <- cbind(true,rowSums(true))
        tmp <- mods[[i]][[j]]$em_list[[1]]$SSB
        tmp <- cbind(tmp, rowSums(tmp))
        ssb.re = tmp/true-1
        
        true <- cbind(mods[[i]][[j]]$om$rep$NAA[1,1,11:30,1],mods[[i]][[j]]$om$rep$NAA[2,2,11:30,1])
        true <- cbind(true,rowSums(true))
        tmp1 <- matrix(mods[[i]][[j]]$em_list[[1]]$NAA[1,1,,1])
        tmp2 <- matrix(mods[[i]][[j]]$em_list[[1]]$NAA[2,2,,1])
        tmp <- cbind(tmp1,tmp2,tmp1+tmp2)
        rec.re = tmp/true-1
      }
      
      ssb.re <- data.frame(ssb.re)
      rec.re <- data.frame(rec.re)
      
      ssb.re$nsim <- i
      ssb.re$OM <- paste(om)
      ssb.re$EM <- paste(j)
      res1 <- rbind(res1, ssb.re)
      
      rec.re$nsim <- i
      rec.re$OM <- paste(om)
      rec.re$EM <- paste(j)
      res2 <- rbind(res2, rec.re)
    }
  }
}

res1[res1$EM == 1,3] = res1[res1$EM == 1,1]
res1[res1$EM == 1,1] = NA

res1[res1$EM == 2,3] = res1[res1$EM == 2,1]
res1[res1$EM == 2,1] = NA

res1[res1$EM == 3,3] = res1[res1$EM == 3,1]
res1[res1$EM == 3,1] = NA

res1[res1$EM == 4,3] = res1[res1$EM == 4,1]
res1[res1$EM == 4,1] = NA

res2[res2$EM == 1,3] = res2[res2$EM == 1,1]
res2[res2$EM == 1,1] = NA

res2[res2$EM == 2,3] = res2[res2$EM == 2,1]
res2[res2$EM == 2,1] = NA

res2[res2$EM == 3,3] = res2[res2$EM == 3,1]
res2[res2$EM == 3,1] = NA

res2[res2$EM == 4,3] = res2[res2$EM == 4,1]
res2[res2$EM == 4,1] = NA

write.csv(res1,"Bias_SSB.csv", row.names = F)
write.csv(res2,"Bias_Rec.csv", row.names = F)


library(ggplot2)
library(tidyr)
library(dplyr)
# Transform data into long format for ggplot
data <- res1 %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  mutate(
    Index = recode(Index,
                   "X1" = "SSB_1",
                   "X2" = "SSB_2",
                   "X3" = "SSB_G")
  )

data <- data %>%
  group_by(OM, EM, Index) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25, na.rm = T),
    Q3 = quantile(Value, 0.75, na.rm = T)
  )

# Define the new names corresponding to each integer value
name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

# Convert the numeric vector to a factor and then map the levels to the new names
data$EM <- factor(data$EM, levels = 1:10, labels = name_map)

# Create a named vector of expressions for the x-axis labels with subscripts
labels <- c(
  EM1_NAA = bquote(bold(EM1[NAA])), 
  EM1_noNAA = bquote(bold(EM1[noNAA])), 
  EM2_NAA = bquote(bold(EM2[NAA])), 
  EM2_noNAA = bquote(bold(EM2[noNAA])), 
  EM3_NAA = bquote(bold(EM3[NAA])), 
  EM3_noNAA = bquote(bold(EM3[noNAA])), 
  EM4_NAA = bquote(bold(EM4[NAA])), 
  EM4_noNAA = bquote(bold(EM4[noNAA])), 
  EM5_Est = bquote(bold(EM5_Est)),
  EM5_Fix = bquote(bold(EM5_Fix))
)

OM_labels_list <- list(
  c("OM 1" = bquote(bold(OM[high]~(F[A]))), "OM 2" = bquote(bold(OM[high]~(F[U]))), "OM 3" = bquote(bold(OM[high]~(F[C])))),
  c("OM 4" = bquote(bold(OM[low]~(F[A]))), "OM 5" = bquote(bold(OM[low]~(F[U]))), "OM 6" = bquote(bold(OM[low]~(F[C])))),
  c("OM 7" = bquote(bold(OM[high]~(F[A]))), "OM 8" = bquote(bold(OM[high]~(F[U]))), "OM 9" = bquote(bold(OM[high]~(F[C])))),
  c("OM 10" = bquote(bold(OM[low]~(F[A]))), "OM 11" = bquote(bold(OM[low]~(F[U]))), "OM 12" = bquote(bold(OM[low]~(F[C])))
  )
)

data$OM <- paste("OM",data$OM)

custom_labeller <- labeller(
  OM = label_parsed
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- data %>% filter(OM %in% names(OM_labels_list[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels_list[[i]]), labels = OM_labels_list[[i]])
  
  # Create the plot for the current subset
  p <- ggplot(om_data, aes(x = EM, y = Median, fill = Index, group = Index)) + # Add group = Var to dodge by each Var level
    # Background shading for specific columns (1, 3, 5, 7, 9, 10)
    geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 1
    geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 3
    geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 5
    geom_rect(aes(xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 7
    geom_rect(aes(xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 9
    geom_rect(aes(xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 10
    # Plot points and error bars with larger dodge width
    geom_point(position = position_dodge(width = 0.9), size = 3, shape = 19) + # Increase dodge width
    geom_errorbar(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.9), width = 0) + # Increase dodge width
    facet_wrap(OM ~ Index, scales = "fixed", ncol = 3, labeller = custom_labeller) +
    labs(
      title = paste("Relative Bias in SSB"),
      x = "Estimation Model", y = ""
    ) +
    scale_x_discrete(labels = labels) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      title = element_text(face = "bold")
    ) +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed")
  
  # Print the plot independently
  print(p)
  
  # Save the plot independently with a unique filename
  ggsave(filename = paste0("Bias_SSB_OM_", i, ".png"), plot = p, width = 10, height = 7, units = "in")
}


library(ggplot2)
library(tidyr)
library(dplyr)
# Transform data into long format for ggplot
data <- res2 %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  mutate(
    Index = recode(Index,
                   "X1" = "SSB_1",
                   "X2" = "SSB_2",
                   "X3" = "SSB_G")
  )

data <- data %>%
  group_by(OM, EM, Index) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25, na.rm = T),
    Q3 = quantile(Value, 0.75, na.rm = T)
  )

# Define the new names corresponding to each integer value
name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

# Convert the numeric vector to a factor and then map the levels to the new names
data$EM <- factor(data$EM, levels = 1:10, labels = name_map)

# Create a named vector of expressions for the x-axis labels with subscripts
labels <- c(
  EM1_NAA = bquote(bold(EM1[NAA])), 
  EM1_noNAA = bquote(bold(EM1[noNAA])), 
  EM2_NAA = bquote(bold(EM2[NAA])), 
  EM2_noNAA = bquote(bold(EM2[noNAA])), 
  EM3_NAA = bquote(bold(EM3[NAA])), 
  EM3_noNAA = bquote(bold(EM3[noNAA])), 
  EM4_NAA = bquote(bold(EM4[NAA])), 
  EM4_noNAA = bquote(bold(EM4[noNAA])), 
  EM5_Est = bquote(bold(EM5_Est)),
  EM5_Fix = bquote(bold(EM5_Fix))
)

OM_labels_list <- list(
  c("OM 1" = bquote(bold(OM[high]~(F[A]))), "OM 2" = bquote(bold(OM[high]~(F[U]))), "OM 3" = bquote(bold(OM[high]~(F[C])))),
  c("OM 4" = bquote(bold(OM[low]~(F[A]))), "OM 5" = bquote(bold(OM[low]~(F[U]))), "OM 6" = bquote(bold(OM[low]~(F[C])))),
  c("OM 7" = bquote(bold(OM[high]~(F[A]))), "OM 8" = bquote(bold(OM[high]~(F[U]))), "OM 9" = bquote(bold(OM[high]~(F[C])))),
  c("OM 10" = bquote(bold(OM[low]~(F[A]))), "OM 11" = bquote(bold(OM[low]~(F[U]))), "OM 12" = bquote(bold(OM[low]~(F[C])))
  )
)

data$OM <- paste("OM",data$OM)

custom_labeller <- labeller(
  OM = label_parsed
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- data %>% filter(OM %in% names(OM_labels_list[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels_list[[i]]), labels = OM_labels_list[[i]])
  
  # Create the plot for the current subset
  p <- ggplot(om_data, aes(x = EM, y = Median, fill = Index, group = Index)) + # Add group = Var to dodge by each Var level
    # Background shading for specific columns (1, 3, 5, 7, 9, 10)
    geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 1
    geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 3
    geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 5
    geom_rect(aes(xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 7
    geom_rect(aes(xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 9
    geom_rect(aes(xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 10
    # Plot points and error bars with larger dodge width
    geom_point(position = position_dodge(width = 0.9), size = 3, shape = 19) + # Increase dodge width
    geom_errorbar(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.9), width = 0) + # Increase dodge width
    facet_wrap(OM ~ Index, scales = "fixed", ncol = 3, labeller = custom_labeller) +
    labs(
      title = paste("Relative Bias in Recruitment"),
      x = "Estimation Model", y = ""
    ) +
    scale_x_discrete(labels = labels) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      title = element_text(face = "bold")
    ) +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed")
  
  # Print the plot independently
  print(p)
  
  # Save the plot independently with a unique filename
  ggsave(filename = paste0("Bias_Rec_OM_", i, ".png"), plot = p, width = 10, height = 7, units = "in")
}

# Probability of overfishing and overfished
res1 <- NULL
res2 <- NULL
use.n.years = 30

F.MSY <- rbind(c(0.2011254,0.2011254,0.4022508),c(0.2011254,0.2011254,0.4022508),c(0.2011254,0.2011254,0.4022508),
               c(0.2019615,0.2019615,0.4039231),c(0.2019615,0.2019615,0.4039231),c(0.2019615,0.2019615,0.4039231),
               c(0.2011254,0.2011254,0.4022508),c(0.2011254,0.2011254,0.4022508),c(0.2011254,0.2011254,0.4022508),
               c(0.2019615,0.2019615,0.4039231),c(0.2019615,0.2019615,0.4039231),c(0.2019615,0.2019615,0.4039231))

B.MSY <- rbind(c(31577.67, 20869.58, 52447.25), c(31577.67, 20869.58, 52447.25), c(31577.67, 20869.58, 52447.25),
               c(31889.43,20620.58,52510.01),c(31889.43,20620.58,52510.01),c(31889.43,20620.58,52510.01),
               c(31577.67, 20869.58, 52447.25),c(31577.67, 20869.58, 52447.25),c(31577.67, 20869.58, 52447.25),
               c(31889.43,20620.58,52510.01),c(31889.43,20620.58,52510.01),c(31889.43,20620.58,52510.01))

for (om in 1:n_oms) {
  mods <- readRDS(paste0("OM_", om, ".RDS"))
  print(om)
  for (i in 1:length(mods)){
    for (j in 1:length(mods[[1]])){
      
      tmp <- mods[[i]][[j]]$om$rep$SSB
      tmp <- cbind(tmp, rowSums(tmp))
      SSB_FXSPR <- B.MSY[om,]
      tmp <- data.frame(sweep(tmp,2,SSB_FXSPR, FUN = "/"))
      tmp <- tail(tmp, use.n.years)
      tmp <- apply(tmp, 2, function(x) sum(x < 1) / use.n.years)
      
      tmp <- data.frame(t(tmp))
      names(tmp) <- paste0("Index_",1:length(tmp))
      tmp$nsim <- i
      tmp$OM <- paste(om)
      tmp$EM <- paste(j)
      res1 <- rbind(res1, tmp)
      
      Fbar <- mods[[i]][[j]]$om$rep$Fbar
      Fbar <- cbind(Fbar, rowSums(Fbar))
      Fbar_XSPR = F.MSY[om,]
      tmp <- data.frame(sweep(Fbar,2,Fbar_XSPR, FUN = "/"))
      
      names(tmp) <- paste0("Index_",(length(tmp)+1:length(tmp)))
      tmp <- tail(tmp, use.n.years)
      tmp <- apply(tmp, 2, function(x) sum(x > 1) / use.n.years)
      tmp <- data.frame(t(tmp))
      tmp$nsim <- i
      tmp$OM <- paste(om)
      tmp$EM <- paste(j)
      res2 <- rbind(res2, tmp)
    }
  }
}

write.csv(res1, "Probability_overfished.csv",row.names = F)
write.csv(res2, "Probability_overfishing.csv",row.names = F)

res1 <- read.csv("Probability_overfished.csv")
res2 <- read.csv("Probability_overfishing.csv")
# Pivot the data and assign numeric labels 1-6 to the Index column
res <- rbind(
  pivot_longer(res1, cols = -c(nsim, OM, EM), names_to = "Index", values_to = "Value"),
  pivot_longer(res2, cols = -c(nsim, OM, EM), names_to = "Index", values_to = "Value")
)

write.csv(res, "Probability_overfished_overfishing.csv",row.names = F)

# ----------------------------------
# Population Status terminal year
use.n.years = 1
# ----------------------------------

# KOBE Plot
res1 <- NULL
res2 <- NULL

#### use a new method for the F reference points ####
for (om in 1:n_oms) {
  mods <- readRDS(paste0("OM_", om, ".RDS"))
  print(om)
  for (i in 1:length(mods)){
    for (j in 1:length(mods[[1]])){
      tmp <- mods[[i]][[j]]$om$rep$SSB
      tmp <- cbind(tmp,rowSums(tmp))
      SSB_FXSPR <- B.MSY[om,]
      tmp <- data.frame(sweep(tmp,2,SSB_FXSPR, FUN = "/"))
      tmp <- tail(tmp, 1)
      names(tmp) <- paste0("S", 1:ncol(tmp))
      tmp$nsim <- i
      tmp$OM <- paste0(om)
      tmp$EM <- paste0(j)
      res1 <- rbind(res1, tmp)
      
      Fbar <- mods[[i]][[j]]$om$rep$Fbar
      Fbar <- cbind(Fbar,rowSums(Fbar))
      Fbar_XSPR = F.MSY[om,]
      tmp <- data.frame(sweep(Fbar,2,Fbar_XSPR, FUN = "/"))
      names(tmp) <- paste0("S",1:ncol(tmp))
      tmp$nsim <- i
      tmp$OM <- paste(om)
      tmp$EM <- paste(j)
      tmp <- tail(tmp, 1)
      res2 <- rbind(res2, tmp)
    }
  }
}

temp1 <- pivot_longer(res1, cols = starts_with("S"), names_to = "Index", values_to = "Overfished")
temp2 <- pivot_longer(res2, cols = starts_with("S"), names_to = "Index", values_to = "Overfishing")

# Function to replace the highest S index with S_total
replace_largest_S <- function(temp) {
  # Extract the numeric part of the S indexes
  tmp <- as.numeric(gsub("S", "", temp$Index))
  
  # Find the maximum numeric value
  max_value <- max(tmp, na.rm = TRUE)
  
  # Replace the largest index with S_total
  temp$Index[temp$Index == paste0("S", max_value)] <- "Global" 
  
  return(temp)
}

# Apply the function to replace the largest S index in temp1$Index
temp1 <- replace_largest_S(temp1)
temp2 <- replace_largest_S(temp2)

temp <- cbind(temp1, temp2)
temp <- temp %>% select(unique(names(.))) 

temp1$Index <- paste0("Overfished_",temp1$Index)
temp2$Index <- paste0("Overfishing_",temp2$Index)
temp1$Value <- temp1$Overfished 
temp2$Value <- temp2$Overfishing 
temp1 <- temp1 %>% select(nsim,OM,EM,Index,Value)
temp2 <- temp2 %>% select(nsim,OM,EM,Index,Value)

write.csv(temp1,"Terminal_pop_status_overfished.csv", row.names = F)
write.csv(temp2,"Terminal_pop_status_overfishing.csv", row.names = F)


temp <- temp %>% select(-nsim)

n.col = length(unique(temp$EM))
my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3","#ADD8E6",
               "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3","#000000")


library(dplyr)
library(ggplot2)
library(tidyr)

# Calculate median and IQR over nsims for each EM and OM combination
calculate_medians_iqr <- function(temp) {
  temp %>%
    group_by(OM, EM, Index) %>%
    summarise(
      Median_Overfished = median(Overfished, na.rm = TRUE),
      Median_Overfishing = median(Overfishing, na.rm = TRUE),
      IQR_Overfished_Lower = quantile(Overfished, 0.25, na.rm = TRUE),
      IQR_Overfished_Upper = quantile(Overfished, 0.75, na.rm = TRUE),
      IQR_Overfishing_Lower = quantile(Overfishing, 0.25, na.rm = TRUE),
      IQR_Overfishing_Upper = quantile(Overfishing, 0.75, na.rm = TRUE)
    ) %>%
    ungroup()
}

write.csv(temp,"Population_status.csv",row.names = F)

#### Pop Status (Global) ####
my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3","#ADD8E6",
               "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3","#000000")

my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
               "#E74C3C", "#E74C3C", "#8E44AD", "#000000")
# Calculate median and IQR over nsims for each EM and OM combination
calculate_medians_iqr <- function(temp) {
  temp %>%
    group_by(OM, EM, Index) %>%
    summarise(
      Median_Overfished = median(Overfished, na.rm = TRUE),
      Median_Overfishing = median(Overfishing, na.rm = TRUE),
      IQR_Overfished_Lower = quantile(Overfished, 0.25, na.rm = TRUE),
      IQR_Overfished_Upper = quantile(Overfished, 0.75, na.rm = TRUE),
      IQR_Overfishing_Lower = quantile(Overfishing, 0.25, na.rm = TRUE),
      IQR_Overfishing_Upper = quantile(Overfishing, 0.75, na.rm = TRUE)
    ) %>%
    ungroup()
}

temp <- read.csv("Population_status.csv")

name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(temp)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

ylim = c(0.4, 1.8)
xlim = c(0.4, 2)

xmin_rect <- max(0.5, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "Global",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "Global","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  custom_labeller <- labeller(
    OM = label_parsed
  )
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB["MSY"]))) +
    ylab(bquote(bold(F/F["MSY"]))) +
    ggtitle("Terminal-year Status (Global)") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    xlim(xlim) + 
    ylim(ylim) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_G_OM_", i, ".png"), p, width = 15, height = 15, dpi = 150)
}


#### Pop Status S1 ####
temp <- read.csv("Population_status.csv")

# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(temp)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

ylim = c(0, 2.5)
xlim = c(0.4, 2.2)

xmin_rect <- max(0.4, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "S1",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "S1","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  custom_labeller <- labeller(
    OM = label_parsed
  )
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB["MSY"]))) +
    ylab(bquote(bold(F/F["MSY"]))) +
    ggtitle("Terminal-year Status (Region 1)") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    xlim(xlim) + 
    ylim(ylim) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_S1_OM_", i, ".png"), p, width = 15, height = 15, dpi = 150)
}


#### Pop Status S2 ####
temp <- read.csv("Population_status.csv")

# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(temp)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

ylim = c(0, 2.5)
xlim = c(0.4, 2.2)

xmin_rect <- max(0.5, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "S2",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "S2","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  custom_labeller <- labeller(
    OM = label_parsed
  )
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB["MSY"]))) +
    ylab(bquote(bold(F/F["MSY"]))) +
    ggtitle("Terminal-year Status (Region 2)") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    xlim(xlim) + 
    ylim(ylim) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_S2_OM_", i, ".png"), p, width = 15, height = 15, dpi = 150)
}

#### Pop Status (Global) ####
my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3","#ADD8E6",
               "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3","#000000")

my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
               "#E74C3C", "#E74C3C", "#8E44AD", "#000000")
# Calculate median and IQR over nsims for each EM and OM combination
calculate_medians_iqr <- function(temp) {
  temp %>%
    group_by(OM, EM, Index) %>%
    summarise(
      Median_Overfished = median(Overfished, na.rm = TRUE),
      Median_Overfishing = median(Overfishing, na.rm = TRUE),
      IQR_Overfished_Lower = quantile(Overfished, 0.25, na.rm = TRUE),
      IQR_Overfished_Upper = quantile(Overfished, 0.75, na.rm = TRUE),
      IQR_Overfishing_Lower = quantile(Overfishing, 0.25, na.rm = TRUE),
      IQR_Overfishing_Upper = quantile(Overfishing, 0.75, na.rm = TRUE)
    ) %>%
    ungroup()
}

#### Terminal year status (S1) Relative ####
data <- read.csv("Population_status.csv")

name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

# Calculate relative difference
for (OM in 1:12){
  data[data$OM == OM, 'nsim'] = rep(1:100, each = 30)
}

data <- data %>%
  group_by(OM, nsim, Index) %>%
  mutate(
    OFD = Overfished[EM == 10], 
    OFG = Overfishing[EM == 10],
    Rel_Diff_OFD = ifelse(EM != 10, (Overfished - OFD) / OFD, NA), # Calculate relative difference
    Rel_Diff_OFG = ifelse(EM != 10, (Overfishing - OFG) / OFG, NA) # Calculate relative difference
  ) %>%
  filter(EM != 10) %>% # Exclude EM 10 from results
  ungroup()

data <- select(data, -c(Overfished,Overfishing,OFD,OFG))
names(data)[5:6] <- c("Overfished","Overfishing")

# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(data)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

xmin_rect <- max(0.5, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "Global",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "Global","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  custom_labeller <- labeller(
    OM = label_parsed
  )
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB["MSY"]))) +
    ylab(bquote(bold(F/F["MSY"]))) +
    ggtitle("Terminal-Year Status (Global)") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_G_OM_", i, "_relative.png"), p, width = 15, height = 15, dpi = 150)
}

data <- read.csv("Population_status.csv")

name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

# Calculate relative difference
for (OM in 1:12){
  data[data$OM == OM, 'nsim'] = rep(1:100, each = 30)
}

data <- data %>%
  group_by(OM, nsim, Index) %>%
  mutate(
    OFD = Overfished[EM == 10], 
    OFG = Overfishing[EM == 10],
    Rel_Diff_OFD = ifelse(EM != 10, (Overfished - OFD) / OFD, NA), # Calculate relative difference
    Rel_Diff_OFG = ifelse(EM != 10, (Overfishing - OFG) / OFG, NA) # Calculate relative difference
  ) %>%
  filter(EM != 10) %>% # Exclude EM 10 from results
  ungroup()

data <- select(data, -c(Overfished,Overfishing,OFD,OFG))
names(data)[5:6] <- c("Overfished","Overfishing")

# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(data)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

xmin_rect <- max(0.5, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "S1",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "S1","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  custom_labeller <- labeller(
    OM = label_parsed
  )
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB["MSY"]))) +
    ylab(bquote(bold(F/F["MSY"]))) +
    ggtitle("Terminal-year Status (Region 1)") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_S1_OM_", i, "_relative.png"), p, width = 15, height = 15, dpi = 150)
}

#### Terminal year status (S2) Relative ####
data <- read.csv("Population_status.csv")

name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

# Calculate relative difference
for (OM in 1:12){
  data[data$OM == OM, 'nsim'] = rep(1:100, each = 30)
}

data <- data %>%
  group_by(OM, nsim, Index) %>%
  mutate(
    OFD = Overfished[EM == 10], 
    OFG = Overfishing[EM == 10],
    Rel_Diff_OFD = ifelse(EM != 10, (Overfished - OFD) / OFD, NA), # Calculate relative difference
    Rel_Diff_OFG = ifelse(EM != 10, (Overfishing - OFG) / OFG, NA) # Calculate relative difference
  ) %>%
  filter(EM != 10) %>% # Exclude EM 10 from results
  ungroup()

data <- select(data, -c(Overfished,Overfishing,OFD,OFG))
names(data)[5:6] <- c("Overfished","Overfishing")

# Apply the function to calculate the medians and IQRs
temp_summary <- calculate_medians_iqr(data)
temp_summary$EM <- factor(temp_summary$EM, levels = 1:10, labels = name_map)

xmin_rect <- max(0.5, xlim[1])
xmax_rect <- min(100, xlim[2])
ymin_rect <- max(-100, ylim[1])
ymax_rect <- min(1, ylim[2])

temp_summary <- temp_summary[temp_summary$Index == "S2",]
temp_summary[temp_summary$EM %in% c("EM1_NAA","EM2_NAA","EM3_NAA","EM4_NAA","EM5_Est", "EM5_Fix"), 'Index'] = "NAA"
temp_summary[temp_summary$Index == "S2","Index"] = "noNAA"

# Custom labels for legend, using parse expressions
legend_labels <- c(
  "EM1_NAA" = expression(PAN[NAA]),
  "EM1_noNAA" = expression(PAN[noNAA]),
  "EM2_NAA" = expression(FAA[NAA]),
  "EM2_noNAA" = expression(FAA[noNAA]),
  "EM3_NAA" = expression(SEP[NAA]),
  "EM3_noNAA" = expression(SEP[noNAA]),
  "EM4_NAA" = expression(SpD[NAA]),
  "EM4_noNAA" = expression(SpD[noNAA]),
  "EM5_Est" = expression(SpE[NAA * "," * Est]),
  "EM5_Fix" = expression(SpE[NAA * "," * Fix])
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- temp_summary %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  custom_labeller <- labeller(
    OM = label_parsed
  )
  
  # Create the plot
  p <- ggplot(data = om_data, aes(x = Median_Overfished, y = Median_Overfishing)) +
    facet_wrap(OM ~ Index, ncol = 2, labeller = custom_labeller) +
    geom_point(aes(color = EM), size = 6, shape = 19, stroke = 2) +
    geom_errorbarh(aes(xmin = IQR_Overfished_Lower, xmax = IQR_Overfished_Upper, color = EM), height = 0.05, alpha = 0.5) +
    geom_errorbar(aes(ymin = IQR_Overfishing_Lower, ymax = IQR_Overfishing_Upper, color = EM), width = 0.05, alpha = 0.5) +
    scale_color_manual(values = my_colors, labels = legend_labels) + # Use legend_labels here
    theme_bw() +
    xlab(bquote(bold(SSB/SSB["MSY"]))) +
    ylab(bquote(bold(F/F["MSY"]))) +
    ggtitle("Terminal-year Status (Region 2)") +
    theme(
      axis.text = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(face = "bold", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      legend.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      legend.key.width = unit(1, "cm"),  # Increase width of legend keys
      legend.key.height = unit(2, "cm"),  # Increase height of legend keys
      legend.spacing.x = unit(1, 'cm'),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8)
  
  # Print the plot
  print(p)
  
  # Save the plot sized for an A4 page
  ggsave(paste0("Population_Status_S2_OM_", i, "_relative.png"), p, width = 15, height = 15, dpi = 150)
}

#### Probability overfishing #### 
data <- read.csv("Probability_Overfished_Overfishing.csv")

# Corrected facet labels
facet_labels_index <- c(
  "Index_1" = "bold(P[SSB < SSB[MSY](R1)])",
  "Index_2" = "bold(P[SSB < SSB[MSY](R2)])",
  "Index_3" = "bold(P[SSB < SSB[MSY](G)])",
  "Index_4" = "bold(P[F > F[MSY](R1)])",
  "Index_5" = "bold(P[F > F[MSY](R2)])",
  "Index_6" = "bold(P[F > F[MSY](G)])"
)

facet_labels_index <- c(
  "Index_1" = "bold(P[SSB(R1)])",
  "Index_2" = "bold(P[SSB(R2)])",
  "Index_3" = "bold(P[SSB(G)])",
  "Index_4" = "bold(P[F(R1)])",
  "Index_5" = "bold(P[F(R2)])",
  "Index_6" = "bold(P[F(G)])"
)

# facet_labels_index <- c(
#   "Index_1" = "bold(atop(P[SSB < SSB[MSY]], S1))",
#   "Index_2" = "bold(atop(P[SSB < SSB[MSY]], S2))",
#   "Index_3" = "bold(atop(P[SSB < SSB[MSY]], G))",
#   "Index_4" = "bold(atop(P[F > F[MSY]], R1))",
#   "Index_5" = "bold(atop(P[F > F[MSY]], R2))",
#   "Index_6" = "bold(atop(P[F > F[MSY]], G))"
# )

for (i in 1:4) {
  
  labels <- c(
    `1` = bquote(bold(PAN[NAA])),
    `2` = bquote(bold(PAN[noNAA])),
    `3` = bquote(bold(FAA[NAA])),
    `4` = bquote(bold(FAA[noNAA])),
    `5` = bquote(bold(SEP[NAA])),
    `6` = bquote(bold(SEP[noNAA])),
    `7` = bquote(bold(SpD[NAA])),
    `8` = bquote(bold(SpD[noNAA])),
    `9` = bquote(bold(SpE[NAA * "," * Est])),
    `10` = bquote(bold(SpE[NAA * "," * Fix]))
  )
  
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Custom labeller function
  custom_labeller <- labeller(
    Index = as_labeller(facet_labels_index, label_parsed),
    OM = label_parsed
  )
  
  # tmp <- tmp[tmp$EM != 10,]
  
  tmp <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
                 "#E74C3C", "#E74C3C", "#8E44AD", "#000000")
  
  plot = ggplot(tmp, aes(x = factor(EM), y = Median, col = factor(EM), fill = factor(EM))) + 
    geom_point(size = 3, shape = 19) + 
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, size = 0.8) + 
    facet_grid(Index ~ OM, scales = "free", labeller = custom_labeller) + 
    ggtitle(expression(bold(paste("Probability of ", italic(SSB), " < ", italic(SSB)[MSY], " and ", italic(F), " > ", italic(F)[MSY])))) + 
    geom_rect(data = data.frame(xmin = c(0.5, 2.5, 4.5, 6.5, 8.5, 9.5),
                                xmax = c(1.5, 3.5, 5.5, 7.5, 9.5, 10.5)),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "lightgrey", alpha = 0.4, inherit.aes = FALSE) + 
    scale_colour_manual(values = my_colors) + 
    scale_fill_manual(values = my_colors) + 
    scale_x_discrete(labels = labels) + 
    labs(x = "Estimation Model", y = "") + 
    theme_bw() +  # Set the background to white
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 15),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      title = element_text(face = "bold")
    ) 
  # geom_hline(yintercept = 1, col = "red", linetype = "dashed")
  print(plot)
  ggsave(paste0("Probability_Overfished_Overfishing_OM",i,".png"), plot, width = 10, height = 10, dpi = 150)
}

#### Probability overfishing (Relative Difference) #### 
data <- read.csv("Probability_Overfished_Overfishing.csv")

data <- data %>% filter(Index == "Index_4" | Index == "Index_5" | Index == "Index_6")

# Load necessary library
library(dplyr)

# Calculate relative difference
data <- data %>%
  group_by(nsim, OM, Index) %>%
  mutate(
    Value_EM10 = Value[EM == 10], # Extract Value for EM 10
    Rel_Diff = ifelse(EM != 10, (Value - Value_EM10) / Value_EM10, NA) # Calculate relative difference
  ) %>%
  filter(EM != 10) %>% # Exclude EM 10 from results
  ungroup()

data <- select(data, -c(Value, Value_EM10))
names(data)[5] = "Value"

# Corrected facet labels
facet_labels_index <- c(
  "Index_4" = "bold(P[OFG][(R1)])",
  "Index_5" = "bold(P[OFG][(R2)])",
  "Index_6" = "bold(P[OFG][(G)])"
)

for (i in 1:4) {
  
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Custom labeller function
  custom_labeller <- labeller(
    Index = as_labeller(facet_labels_index, label_parsed),
    OM = label_parsed
  )
  
  tmp <- tmp[tmp$EM != 10,]
  
  tmp <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25, na.rm = T),
      Q3 = quantile(Value, 0.75, na.rm = T),
      min_val = boxplot.stats(Value[!is.na(Value)])$stats[1],
      max_val = boxplot.stats(Value[!is.na(Value)])$stats[5],
      .groups = 'drop'
    )
  
  plot = ggplot(tmp, aes(x = factor(EM), y = Median, col = factor(EM), fill = factor(EM))) + 
    geom_point(size = 3, shape = 19) + 
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, size = 0.8) + 
    facet_grid(Index ~ OM, scales = "free", labeller = custom_labeller) + 
    ggtitle(paste0("Probability of Overfishing (Relative Difference)")) + 
    geom_rect(data = data.frame(xmin = c(0.5, 2.5, 4.5, 6.5, 8.5),
                                xmax = c(1.5, 3.5, 5.5, 7.5, 9.5)),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "lightgrey", alpha = 0.4, inherit.aes = FALSE) + 
    scale_colour_manual(values = my_colors) + 
    scale_fill_manual(values = my_colors) + 
    scale_x_discrete(labels = labels) + 
    labs(x = "Estimation Model", y = "") + 
    theme_bw() +  # Set the background to white
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 15),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      title = element_text(face = "bold")
    ) + 
    geom_hline(yintercept = 0, col = "red", linetype = "dashed")
  print(plot)
  ggsave(paste0("Probability_Overfishing_OM",i,"_relative.png"), plot, width = 10, height = 7, dpi = 150)
}

#### New Holistic Performance ####

scale_values <- function(x) {
  if (min(x) == max(x)) {
    return(rep(1, length(x)))
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}

scale_values2 <- function(x) {
  if (min(x) == max(x)) {
    return(rep(1, length(x)))
  } else {
    return(1 - (x - min(x)) / (max(x) - min(x)))
  }
}

# Catch and SSB last 10 years
data1 <- read.csv("Variables_of_Interest_5_years.csv")

# Specify the indices you want to extract
indices_to_extract <- c("Catch_Total", "SSB_Total")

# Extract the specific indices
data1 <- data1 %>%
  filter(Index %in% indices_to_extract)
data1$Index <- ifelse(data1$Index == "Catch_Total","Long_term_Catch","Long_term_SSB")
data1 <- data1 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data2 <- read.csv("Variables_of_Interest_first_5_years.csv")
data2 <- data2 %>%
  filter(Index %in% indices_to_extract)
data2$Index <- ifelse(data2$Index == "Catch_Total","Short_term_Catch","Short_term_SSB")
data2 <- data2 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data3 <- read.csv("Terminal_pop_status_overfished.csv")
indices_to_extract <- c("Overfished_Global")
data3 <- data3 %>%
  filter(Index %in% indices_to_extract)
data3$Index <- ifelse(data3$Index == "Overfishing_Global","F/FMSY","SSB/SSBMSY")
data3 <- data3 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data4 <- read.csv("Terminal_pop_status_overfishing.csv")
indices_to_extract <- c("Overfishing_Global")
data4 <- data4 %>%
  filter(Index %in% indices_to_extract)
data4$Index <- ifelse(data4$Index == "Overfishing_Global","F/FMSY","SSB/SSBMSY")
data4 <- data4 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data5 <- read.csv("Probability_overfished_overfishing.csv")
indices_to_extract <- c("Index_3")
data5 <- data5 %>%
  filter(Index %in% indices_to_extract)
data5$Index <- ifelse(data5$Index == "Index_3","P_OFD","P_OFG")
data5 <- data5 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data6 <- read.csv("Probability_overfished_overfishing.csv")
indices_to_extract <- c("Index_6")
data6 <- data6 %>%
  filter(Index %in% indices_to_extract)
data6$Index <- ifelse(data6$Index == "Index_3","P_OFD","P_OFG")
data6 <- data6 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data7 <- read.csv("Catch_AAV.csv")
indices_to_extract <- c("Catch.G")
data7 <- data7 %>%
  filter(Index %in% indices_to_extract)
data7$Index <- ifelse(data7$Index == "Catch.G","AAV_Catch","AAV_Catch")
data7 <- data7 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data <- rbind(data1, data2, data3, data4, data5, data6, data7)
# data <- data[data$EM != 10,]
data$EM <- factor(data$EM, levels = 1:10, labels = name_map[1:10])

sum_data <- data %>%
  group_by(OM, EM, Index) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75)
  )

sum_data$EM = factor(sum_data$EM, levels = c("EM5_Fix",
                                             "EM5_Est",
                                             "EM4_NAA",
                                             "EM3_NAA",
                                             "EM2_NAA",
                                             "EM1_NAA",
                                             "EM4_noNAA",
                                             "EM3_noNAA",
                                             "EM2_noNAA",
                                             "EM1_noNAA"))

position_dodge_val = position_dodge(width = 0.9)
custom_labeller <- labeller(
  OM = label_parsed
)

sum_data$Index = factor(sum_data$Index, levels = rev(c("Long_term_Catch",
                                                       "Short_term_Catch",
                                                       "Long_term_SSB",
                                                       "Short_term_SSB",
                                                       "P_OFG",
                                                       "P_OFD",
                                                       "F/FMSY",
                                                       "SSB/SSBMSY",
                                                       "AAV_Catch")))

legend_labels <- c(
  "Short_term_Catch" = expression(Catch[short]),
  "Short_term_SSB" = expression(SSB[short]),
  "Long_term_Catch" = expression(Catch[long]),
  "Long_term_SSB" = expression(SSB[long]),
  "P_OFG" = expression(P[F<F[MSY]]),
  "P_OFD" = expression(P[SSB>SSB[MSY]]),
  "F/FMSY" = expression(OFG[terminal]),
  "SSB/SSBMSY" = expression(OFD[terminal]),
  "AAV_Catch" = expression(AACV)
)

# Calculate the sum of scores for each EM and OM
sum_data2 <- sum_data %>%
  group_by(OM, EM) %>%
  summarize(Total_Score = sum(Median, na.rm = TRUE),
            .groups = 'drop')
sum_data2$Total_Score <- round(sum_data2$Total_Score, 1)

for (i in 1:4) {
  summary_data <- sum_data %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  summary_data$OM <- factor(summary_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  summary_data2 <- sum_data2 %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  summary_data2$OM <- factor(summary_data2$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  my_colors <- c("#2874A6", "#F39C12", "#27AE60", "#E74C3C", 
                 "#8E44AD", "#17A589", "#F4D03F", "#34495E", "#b3b3b3")
  
  # Create bar graph with IQR
  plot <- ggplot(summary_data, aes(x = Median, y = EM, fill = Index, color = Index)) +
    geom_bar(stat = "identity", position = position_dodge_val, width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3),
                  position = position_dodge_val, width = 0.2, size = 0.5, alpha = 0.5) +  # IQR error bars
    facet_wrap(~OM, labeller = custom_labeller) +
    scale_fill_manual(values = my_colors, labels = legend_labels) +
    scale_color_manual(values = my_colors, labels = legend_labels) +
    scale_y_discrete(labels = c(
      "EM1_NAA"   = expression(PAN[NAA]),
      "EM1_noNAA" = expression(PAN[noNAA]),
      "EM2_NAA"   = expression(FAA[NAA]),
      "EM2_noNAA" = expression(FAA[noNAA]),
      "EM3_NAA"   = expression(SEP[NAA]),
      "EM3_noNAA" = expression(SEP[noNAA]),
      "EM4_NAA"   = expression(SpD[NAA]),
      "EM4_noNAA" = expression(SpD[noNAA]),
      "EM5_Est"   = expression(SpE[NAA * "," * Est]),
      "EM5_Fix"   = expression(SpE[NAA * "," * Fix])
    )) +
    labs(
      title = "Holistic Model Performance", 
      x = "Score", 
      y = "Estimation Model"
    ) +
    guides(
      fill  = guide_legend(title = "Metric"),
      color = guide_legend(title = "Metric")
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      title = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"),
      legend.key.height = unit(2, "cm")
    ) +
    geom_hline(yintercept = 0.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 1.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 2.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 3.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 4.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 5.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 6.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 7.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 8.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 9.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 10.5, linetype = "dashed", size = 0.8) + 
    geom_rect(data = data.frame(ymin = c(0.5),
                                ymax = c(6.5)),
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
              fill = "lightgrey", alpha = 0.3, inherit.aes = FALSE)  +
    geom_text(data = summary_data2,
              aes(x = 0.95, y = EM, label = paste(Total_Score)),
              inherit.aes = FALSE,
              size = 6,
              hjust = 0,
              color = "black")
  
  print(plot)
  ggsave(paste0("Holistic_OM_", i, ".png"), plot, width = 15, height = 15, dpi = 150)
}



#### Region 1 Holistic ####

# Catch and SSB last 10 years
data1 <- read.csv("Variables_of_Interest_5_years.csv")

# Specify the indices you want to extract
indices_to_extract <- c("Catch_R.1", "SSB_S.1")

# Extract the specific indices
data1 <- data1 %>%
  filter(Index %in% indices_to_extract)
data1$Index <- ifelse(data1$Index == "Catch_R.1","Long_term_Catch","Long_term_SSB")
data1 <- data1 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data2 <- read.csv("Variables_of_Interest_first_5_years.csv")
data2 <- data2 %>%
  filter(Index %in% indices_to_extract)
data2$Index <- ifelse(data2$Index == "Catch_R.1","Short_term_Catch","Short_term_SSB")
data2 <- data2 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data3 <- read.csv("Terminal_pop_status_overfished.csv")
indices_to_extract <- c("Overfished_S1")
data3 <- data3 %>%
  filter(Index %in% indices_to_extract)
data3$Index <- ifelse(data3$Index == "Overfishing_S1","F/FMSY","SSB/SSBMSY")
data3 <- data3 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data4 <- read.csv("Terminal_pop_status_overfishing.csv")
indices_to_extract <- c("Overfishing_S1")
data4 <- data4 %>%
  filter(Index %in% indices_to_extract)
data4$Index <- ifelse(data4$Index == "Overfishing_S1","F/FMSY","SSB/SSBMSY")
data4 <- data4 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data5 <- read.csv("Probability_overfished_overfishing.csv")
indices_to_extract <- c("Index_1")
data5 <- data5 %>%
  filter(Index %in% indices_to_extract)
data5$Index <- ifelse(data5$Index == "Index_1","P_OFD","P_OFG")
data5 <- data5 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data6 <- read.csv("Probability_overfished_overfishing.csv")
indices_to_extract <- c("Index_4")
data6 <- data6 %>%
  filter(Index %in% indices_to_extract)
data6$Index <- ifelse(data6$Index == "Index_1","P_OFD","P_OFG")
data6 <- data6 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data7 <- read.csv("Catch_AAV.csv")
indices_to_extract <- c("Catch.G")
data7 <- data7 %>%
  filter(Index %in% indices_to_extract)
data7$Index <- ifelse(data7$Index == "Catch.1","AAV_Catch","AAV_Catch")
data7 <- data7 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data <- rbind(data1, data2, data3, data4, data5, data6, data7)
# data <- data[data$EM != 10,]
data$EM <- factor(data$EM, levels = 1:10, labels = name_map[1:10])

sum_data <- data %>%
  group_by(OM, EM, Index) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75)
  )

sum_data$EM = factor(sum_data$EM, levels = c("EM5_Fix",
                                             "EM5_Est",
                                             "EM4_NAA",
                                             "EM3_NAA",
                                             "EM2_NAA",
                                             "EM1_NAA",
                                             "EM4_noNAA",
                                             "EM3_noNAA",
                                             "EM2_noNAA",
                                             "EM1_noNAA"))

position_dodge_val = position_dodge(width = 0.9)
custom_labeller <- labeller(
  OM = label_parsed
)

sum_data$Index = factor(sum_data$Index, levels = rev(c("Long_term_Catch",
                                                       "Short_term_Catch",
                                                       "Long_term_SSB",
                                                       "Short_term_SSB",
                                                       "P_OFG",
                                                       "P_OFD",
                                                       "F/FMSY",
                                                       "SSB/SSBMSY",
                                                       "AAV_Catch")))

legend_labels <- c(
  "Short_term_Catch" = expression(Catch[short]),
  "Short_term_SSB" = expression(SSB[short]),
  "Long_term_Catch" = expression(Catch[long]),
  "Long_term_SSB" = expression(SSB[long]),
  "P_OFG" = expression(P[F<F[MSY]]),
  "P_OFD" = expression(P[SSB>SSB[MSY]]),
  "F/FMSY" = expression(OFG[terminal]),
  "SSB/SSBMSY" = expression(OFD[terminal]),
  "AAV_Catch" = expression(AACV)
)

# Calculate the sum of scores for each EM and OM
sum_data2 <- sum_data %>%
  group_by(OM, EM) %>%
  summarize(Total_Score = sum(Median, na.rm = TRUE),
            .groups = 'drop')
sum_data2$Total_Score <- round(sum_data2$Total_Score, 1)

for (i in 1:4) {
  summary_data <- sum_data %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  summary_data$OM <- factor(summary_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  summary_data2 <- sum_data2 %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  summary_data2$OM <- factor(summary_data2$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  my_colors <- c("#2874A6", "#F39C12", "#27AE60", "#E74C3C", 
                 "#8E44AD", "#17A589", "#F4D03F", "#34495E", "#b3b3b3")
  
  # Create bar graph with IQR
  plot <- ggplot(summary_data, aes(x = Median, y = EM, fill = Index, color = Index)) +
    geom_bar(stat = "identity", position = position_dodge_val, width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3),
                  position = position_dodge_val, width = 0.2, size = 0.5, alpha = 0.5) +  # IQR error bars
    facet_wrap(~OM, labeller = custom_labeller) +
    scale_fill_manual(values = my_colors, labels = legend_labels) +
    scale_color_manual(values = my_colors, labels = legend_labels, guide = "none") + # remove duplicate legend
    scale_y_discrete(labels = c(
      "EM1_NAA"   = expression(PAN[NAA]),
      "EM1_noNAA" = expression(PAN[noNAA]),
      "EM2_NAA"   = expression(FAA[NAA]),
      "EM2_noNAA" = expression(FAA[noNAA]),
      "EM3_NAA"   = expression(SEP[NAA]),
      "EM3_noNAA" = expression(SEP[noNAA]),
      "EM4_NAA"   = expression(SpD[NAA]),
      "EM4_noNAA" = expression(SpD[noNAA]),
      "EM5_Est"   = expression(SpE[NAA * "," * Est]),
      "EM5_Fix"   = expression(SpE[NAA * "," * Fix])
    )) +
    labs(
      title = "Holistic Model Performance", 
      x = "Score", 
      y = "Estimation Model", 
      fill = "Metric"   # <- legend title updated here
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      title = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"),
      legend.key.height = unit(2, "cm")
    ) +
    geom_hline(yintercept = 0.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 1.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 2.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 3.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 4.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 5.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 6.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 7.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 8.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 9.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 10.5, linetype = "dashed", size = 0.8) + 
    geom_rect(data = data.frame(ymin = c(0.5),
                                ymax = c(6.5)),
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
              fill = "lightgrey", alpha = 0.3, inherit.aes = FALSE) +
    geom_text(data = summary_data2,
              aes(x = 0.95, y = EM, label = paste(Total_Score)),
              inherit.aes = FALSE,
              size = 6,
              hjust = 0,
              color = "black")
  
  print(plot)
  ggsave(paste0("Holistic_OM_R1_", i, ".png"), plot, width = 15, height = 15, dpi = 150)
}


#### Region 2 Holistic ####

# Catch and SSB last 10 years
data1 <- read.csv("Variables_of_Interest_5_years.csv")

# Specify the indices you want to extract
indices_to_extract <- c("Catch_R.2", "SSB_S.2")

# Extract the specific indices
data1 <- data1 %>%
  filter(Index %in% indices_to_extract)
data1$Index <- ifelse(data1$Index == "Catch_R.2","Long_term_Catch","Long_term_SSB")
data1 <- data1 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data2 <- read.csv("Variables_of_Interest_first_5_years.csv")
data2 <- data2 %>%
  filter(Index %in% indices_to_extract)
data2$Index <- ifelse(data2$Index == "Catch_R.2","Short_term_Catch","Short_term_SSB")
data2 <- data2 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data3 <- read.csv("Terminal_pop_status_overfished.csv")
indices_to_extract <- c("Overfished_S2")
data3 <- data3 %>%
  filter(Index %in% indices_to_extract)
data3$Index <- ifelse(data3$Index == "Overfishing_S2","F/FMSY","SSB/SSBMSY")
data3 <- data3 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values(Value))

data4 <- read.csv("Terminal_pop_status_overfishing.csv")
indices_to_extract <- c("Overfishing_S2")
data4 <- data4 %>%
  filter(Index %in% indices_to_extract)
data4$Index <- ifelse(data4$Index == "Overfishing_S2","F/FMSY","SSB/SSBMSY")
data4 <- data4 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data5 <- read.csv("Probability_overfished_overfishing.csv")
indices_to_extract <- c("Index_2")
data5 <- data5 %>%
  filter(Index %in% indices_to_extract)
data5$Index <- ifelse(data5$Index == "Index_2","P_OFD","P_OFG")
data5 <- data5 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data6 <- read.csv("Probability_overfished_overfishing.csv")
indices_to_extract <- c("Index_5")
data6 <- data6 %>%
  filter(Index %in% indices_to_extract)
data6$Index <- ifelse(data6$Index == "Index_2","P_OFD","P_OFG")
data6 <- data6 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data7 <- read.csv("Catch_AAV.csv")
indices_to_extract <- c("Catch.G")
data7 <- data7 %>%
  filter(Index %in% indices_to_extract)
data7$Index <- ifelse(data7$Index == "Catch.2","AAV_Catch","AAV_Catch")
data7 <- data7 %>%
  group_by(OM, Index, nsim) %>%
  mutate(Value = scale_values2(Value))

data <- rbind(data1, data2, data3, data4, data5, data6, data7)
# data <- data[data$EM != 10,]
data$EM <- factor(data$EM, levels = 1:10, labels = name_map[1:10])

sum_data <- data %>%
  group_by(OM, EM, Index) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25),
    Q3 = quantile(Value, 0.75)
  )

sum_data$EM = factor(sum_data$EM, levels = c("EM5_Fix",
                                             "EM5_Est",
                                             "EM4_NAA",
                                             "EM3_NAA",
                                             "EM2_NAA",
                                             "EM1_NAA",
                                             "EM4_noNAA",
                                             "EM3_noNAA",
                                             "EM2_noNAA",
                                             "EM1_noNAA"))

position_dodge_val = position_dodge(width = 0.9)
custom_labeller <- labeller(
  OM = label_parsed
)

sum_data$Index = factor(sum_data$Index, levels = rev(c("Long_term_Catch",
                                                       "Short_term_Catch",
                                                       "Long_term_SSB",
                                                       "Short_term_SSB",
                                                       "P_OFG",
                                                       "P_OFD",
                                                       "F/FMSY",
                                                       "SSB/SSBMSY",
                                                       "AAV_Catch")))

legend_labels <- c(
  "Short_term_Catch" = expression(Catch[short]),
  "Short_term_SSB" = expression(SSB[short]),
  "Long_term_Catch" = expression(Catch[long]),
  "Long_term_SSB" = expression(SSB[long]),
  "P_OFG" = expression(P[F<F[MSY]]),
  "P_OFD" = expression(P[SSB>SSB[MSY]]),
  "F/FMSY" = expression(OFG[terminal]),
  "SSB/SSBMSY" = expression(OFD[terminal]),
  "AAV_Catch" = expression(AACV)
)

# Calculate the sum of scores for each EM and OM
sum_data2 <- sum_data %>%
  group_by(OM, EM) %>%
  summarize(Total_Score = sum(Median, na.rm = TRUE),
            .groups = 'drop')
sum_data2$Total_Score <- round(sum_data2$Total_Score, 1)


for (i in 1:4) {
  summary_data <- sum_data %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  summary_data$OM <- factor(summary_data$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  summary_data2 <- sum_data2 %>% filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  summary_data2$OM <- factor(summary_data2$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  my_colors <- c("#2874A6", "#F39C12", "#27AE60", "#E74C3C", 
                 "#8E44AD", "#17A589", "#F4D03F", "#34495E", "#b3b3b3")
  
  # Create bar graph with IQR
  plot <- ggplot(summary_data, aes(x = Median, y = EM, fill = Index, color = Index)) +
    geom_bar(stat = "identity", position = position_dodge_val, width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3),
                  position = position_dodge_val, width = 0.2, size = 0.5, alpha = 0.5) +  # IQR error bars
    facet_wrap(~OM, labeller = custom_labeller) +
    scale_fill_manual(values = my_colors, labels = legend_labels) +
    scale_color_manual(values = my_colors, labels = legend_labels, guide = "none") + # remove duplicate legend
    scale_y_discrete(labels = c(
      "EM1_NAA"   = expression(PAN[NAA]),
      "EM1_noNAA" = expression(PAN[noNAA]),
      "EM2_NAA"   = expression(FAA[NAA]),
      "EM2_noNAA" = expression(FAA[noNAA]),
      "EM3_NAA"   = expression(SEP[NAA]),
      "EM3_noNAA" = expression(SEP[noNAA]),
      "EM4_NAA"   = expression(SpD[NAA]),
      "EM4_noNAA" = expression(SpD[noNAA]),
      "EM5_Est"   = expression(SpE[NAA * "," * Est]),
      "EM5_Fix"   = expression(SpE[NAA * "," * Fix])
    )) +
    labs(
      title = "Holistic Model Performance", 
      x = "Score", 
      y = "Estimation Model", 
      fill = "Metric"   # <- legend title updated here
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 20, face = "bold"),
      title = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20, face = "bold"),
      legend.key.height = unit(2, "cm")
    ) +
    geom_hline(yintercept = 0.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 1.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 2.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 3.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 4.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 5.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 6.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 7.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 8.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 9.5,  linetype = "dashed", size = 0.8) +
    geom_hline(yintercept = 10.5, linetype = "dashed", size = 0.8) + 
    geom_rect(data = data.frame(ymin = c(0.5),
                                ymax = c(6.5)),
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
              fill = "lightgrey", alpha = 0.3, inherit.aes = FALSE) +
    geom_text(data = summary_data2,
              aes(x = 0.95, y = EM, label = paste(Total_Score)),
              inherit.aes = FALSE,
              size = 6,
              hjust = 0,
              color = "black")
  
  print(plot)
  ggsave(paste0("Holistic_OM_R2_", i, ".png"), plot, width = 15, height = 15, dpi = 150)
}


#### Catch variation ####

data <- read.csv("Catch_AAV.csv")

OM_labels <- list(
  c("1" = bquote(bold(OM[high]~(F[A]))), "2" = bquote(bold(OM[high]~(F[U]))), "3" = bquote(bold(OM[high]~(F[C])))),
  c("4" = bquote(bold(OM[low]~(F[A]))), "5" = bquote(bold(OM[low]~(F[U]))), "6" = bquote(bold(OM[low]~(F[C])))),
  c("7" = bquote(bold(OM[high]~(F[A]))), "8" = bquote(bold(OM[high]~(F[U]))), "9" = bquote(bold(OM[high]~(F[C])))),
  c("10" = bquote(bold(OM[low]~(F[A]))), "11" = bquote(bold(OM[low]~(F[U]))), "12" = bquote(bold(OM[low]~(F[C])))
  ))

labels <- c(
  `1` = bquote(bold(PAN[NAA])),
  `2` = bquote(bold(PAN[noNAA])),
  `3` = bquote(bold(FAA[NAA])),
  `4` = bquote(bold(FAA[noNAA])),
  `5` = bquote(bold(SEP[NAA])),
  `6` = bquote(bold(SEP[noNAA])),
  `7` = bquote(bold(SpD[NAA])),
  `8` = bquote(bold(SpD[noNAA])),
  `9` = bquote(bold(SpE[NAA * "," * Est])),
  `10` = bquote(bold(SpE[NAA * "," * Fix]))
)

for (i in 1:4) {
  
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  tmp <- tmp %>%
    mutate(
      Index = case_when(
        Index == "Catch.1" ~ "bold(C[R1])",
        Index == "Catch.2" ~ "bold(C[R2])",
        Index == "Catch.G" ~ "bold(C[G])"
      )
    )
  
  # Reorder Index to the desired level order
  tmp <- tmp %>%
    mutate(
      Index = factor(
        Index,
        levels = c("bold(C[R1])", "bold(C[R2])", "bold(C[G])")
      )
    )
  
  # Summarize the data
  tmp <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  plot <- ggplot(tmp, aes(x = factor(EM), y = Median, col = factor(EM), fill = factor(EM))) + 
    geom_point(size = 3, shape = 19) + 
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, size = 0.8) + 
    facet_grid(Index ~ OM, scales = "free", labeller = label_parsed) + 
    ggtitle(paste0("Average Annual Catch Variation")) + 
    geom_rect(data = data.frame(xmin = c(0.5, 2.5, 4.5, 6.5, 8.5, 9.5),
                                xmax = c(1.5, 3.5, 5.5, 7.5, 9.5, 10.5)),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "lightgrey", alpha = 0.4, inherit.aes = FALSE) + 
    scale_colour_manual(values = my_colors) + 
    scale_fill_manual(values = my_colors) + 
    scale_x_discrete(labels = labels) + 
    labs(x = "Estimation Model", y = "") + 
    theme_bw() +  # Set the background to white
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 15),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      title = element_text(face = "bold")
    ) 
  
  print(plot)
  ggsave(paste0("Catch_variation_OM",i,".png"), plot, width = 10, height = 7, dpi = 150)
}

#### Functions ####
create_plot <- function(data, use_n_years = 30, boxplot = FALSE,
                        ggtitle_text = "Model Performance", 
                        file_name = "Model_performance",
                        plot_width = 12, plot_height = 10, 
                        point_size = 3, line_thickness = 0.8, 
                        x_text_size = 12, y_text_size = 12,
                        axis_title_size = 15, strip_text_size = 12) {
  
  # Set up the ggplot object
  plot <- ggplot(data_summary, aes(x = factor(EM), y = Median, col = factor(EM), fill = factor(EM))) +
    # Facet with parsed labels
    facet_grid(Index ~ OM, scales = "free", labeller = label_parsed) +
    # Background shading for specific columns (1, 3, 5, 7, 9, 10)
    geom_rect(data = data.frame(xmin = c(0.5, 2.5, 4.5, 6.5, 8.5,9.5),
                                xmax = c(1.5, 3.5, 5.5, 7.5, 9.5,10.5)),
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "lightgrey", alpha = 0.4, inherit.aes = FALSE) +
    scale_colour_manual(values = my_colors) +
    scale_fill_manual(values = my_colors) +
    scale_x_discrete(labels = labels) +
    labs(x = "Estimation Model", y = "") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = x_text_size),
      axis.text.y = element_text(size = y_text_size),
      axis.title = element_text(face = "bold", size = axis_title_size),
      strip.text = element_text(face = "bold", size = strip_text_size),
      legend.position = "none",
      title = element_text(face = "bold")
    ) +
    ggtitle(ggtitle_text)
  
  # Add elements based on the boxplot parameter
  if (boxplot) {
    # Boxplot style: Box segment for IQR, median as horizontal line, whiskers as vertical lines
    plot <- plot +
      # Box for the IQR
      geom_rect(aes(xmin = as.numeric(EM) - 0.4, xmax = as.numeric(EM) + 0.4,
                    ymin = Q1, ymax = Q3, col = factor(EM)), fill = NA, size = line_thickness) +
      # Median line
      geom_segment(aes(x = as.numeric(EM) - 0.4, xend = as.numeric(EM) + 0.4,
                       y = Median, yend = Median, col = factor(EM)), size = line_thickness) +
      # Lower whisker (from min to Q1)
      geom_segment(aes(x = as.numeric(EM), xend = as.numeric(EM),
                       y = min_val, yend = Q1, col = factor(EM)), size = line_thickness) +
      # Upper whisker (from Q3 to max)
      geom_segment(aes(x = as.numeric(EM), xend = as.numeric(EM),
                       y = Q3, yend = max_val, col = factor(EM)), size = line_thickness)
  } else {
    # Point and IQR error bar style: Median as point and IQR as error bar
    plot <- plot +
      geom_point(aes(y = Median), size = point_size, shape = 19) +
      geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, size = line_thickness)
  }
  
  # Print the plot independently
  print(plot)
  
  # Save the plot
  ggsave(paste0(file_name, ".png"), 
         plot, width = plot_width, height = plot_height)
}

#### Last 5 performance ####
use.n.years <- 5
data <- read.csv(file = paste0("Variables_of_Interest_", use.n.years, "_years.csv"))

OM_labels <- list(
  c("1" = bquote(bold(OM[high]~(F[A]))), "2" = bquote(bold(OM[high]~(F[U]))), "3" = bquote(bold(OM[high]~(F[C])))),
  c("4" = bquote(bold(OM[low]~(F[A]))), "5" = bquote(bold(OM[low]~(F[U]))), "6" = bquote(bold(OM[low]~(F[C])))),
  c("7" = bquote(bold(OM[high]~(F[A]))), "8" = bquote(bold(OM[high]~(F[U]))), "9" = bquote(bold(OM[high]~(F[C])))),
  c("10" = bquote(bold(OM[low]~(F[A]))), "11" = bquote(bold(OM[low]~(F[U]))), "12" = bquote(bold(OM[low]~(F[C])))
  )
)

labels <- c(
  `1` = bquote(bold(PAN[NAA])),
  `2` = bquote(bold(PAN[noNAA])),
  `3` = bquote(bold(FAA[NAA])),
  `4` = bquote(bold(FAA[noNAA])),
  `5` = bquote(bold(SEP[NAA])),
  `6` = bquote(bold(SEP[noNAA])),
  `7` = bquote(bold(SpD[NAA])),
  `8` = bquote(bold(SpD[noNAA])),
  `9` = bquote(bold(SpE[NAA * "," * Est])),
  `10` = bquote(bold(SpE[NAA * "," * Fix]))
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  
  filename = paste0("Performance_last_", use.n.years, "_years_OM_", i)
  
  # Filter data for the current subset of OMs
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Define custom colors
  my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ADD8E6",
                 "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  
  my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
                 "#E74C3C", "#E74C3C", "#8E44AD","#000000")
  # Summarize the data
  data_summary <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  # Rename Index directly in data_summary with parsed expressions using bquote-like format
  data_summary <- data_summary %>%
    mutate(
      Index = case_when(
        Index == "Catch_R.1" ~ "bold(C[R1])",
        Index == "Catch_R.2" ~ "bold(C[R2])",
        Index == "Catch_Total" ~ "bold(C[G])",
        Index == "F_R.1" ~ "bold(F[R1])",
        Index == "F_R.2" ~ "bold(F[R2])",
        Index == "F_Total" ~ "bold(F[G])",
        Index == "SSB_S.1" ~ "bold(SSB[R1])",
        Index == "SSB_S.2" ~ "bold(SSB[R2])",
        Index == "SSB_Total" ~ "bold(SSB[G])"
      )
    )
  
  # Reorder Index to the desired level order
  data_summary <- data_summary %>%
    mutate(
      Index = factor(
        Index,
        levels = c("bold(C[R1])", "bold(C[R2])", "bold(C[G])",
                   "bold(F[R1])", "bold(F[R2])", "bold(F[G])",
                   "bold(SSB[R1])", "bold(SSB[R2])", "bold(SSB[G])")
      )
    )
  
  create_plot(data_summary, boxplot = F, file_name = filename)
}

#### First 5 performance ####
use.n.years <- 5
data <- read.csv(file = paste0("Variables_of_Interest_first_", use.n.years, "_years.csv"))

OM_labels <- list(
  c("1" = bquote(bold(OM[high]~(F[A]))), "2" = bquote(bold(OM[high]~(F[U]))), "3" = bquote(bold(OM[high]~(F[C])))),
  c("4" = bquote(bold(OM[low]~(F[A]))), "5" = bquote(bold(OM[low]~(F[U]))), "6" = bquote(bold(OM[low]~(F[C])))),
  c("7" = bquote(bold(OM[high]~(F[A]))), "8" = bquote(bold(OM[high]~(F[U]))), "9" = bquote(bold(OM[high]~(F[C])))),
  c("10" = bquote(bold(OM[low]~(F[A]))), "11" = bquote(bold(OM[low]~(F[U]))), "12" = bquote(bold(OM[low]~(F[C])))
  )
)

labels <- c(
  `1` = bquote(bold(PAN[NAA])),
  `2` = bquote(bold(PAN[noNAA])),
  `3` = bquote(bold(FAA[NAA])),
  `4` = bquote(bold(FAA[noNAA])),
  `5` = bquote(bold(SEP[NAA])),
  `6` = bquote(bold(SEP[noNAA])),
  `7` = bquote(bold(SpD[NAA])),
  `8` = bquote(bold(SpD[noNAA])),
  `9` = bquote(bold(SpE[NAA * "," * Est])),
  `10` = bquote(bold(SpE[NAA * "," * Fix]))
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  
  filename = paste0("Performance_first_", use.n.years, "_years_OM_", i)
  
  # Filter data for the current subset of OMs
  tmp <- data %>%
    mutate(OM = factor(OM)) %>%
    filter(OM %in% names(OM_labels[[i]]))
  
  # Set OM factor levels and labels for the current subset
  tmp$OM <- factor(tmp$OM, levels = names(OM_labels[[i]]), labels = OM_labels[[i]])
  
  # Define custom colors
  my_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#ADD8E6",
                 "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  
  my_colors <- c("#2874A6", "#2874A6", "#F39C12", "#F39C12", "#27AE60", "#27AE60", 
                 "#E74C3C", "#E74C3C", "#8E44AD","#000000")
  
  # Summarize the data
  data_summary <- tmp %>%
    group_by(OM, EM, Index) %>%
    summarize(
      Median = median(Value),
      Q1 = quantile(Value, 0.25),
      Q3 = quantile(Value, 0.75),
      min_val = boxplot.stats(Value)$stats[1],
      max_val = boxplot.stats(Value)$stats[5],
      .groups = 'drop'
    )
  
  # Rename Index directly in data_summary with parsed expressions using bquote-like format
  data_summary <- data_summary %>%
    mutate(
      Index = case_when(
        Index == "Catch_R.1" ~ "bold(C[R1])",
        Index == "Catch_R.2" ~ "bold(C[R2])",
        Index == "Catch_Total" ~ "bold(C[G])",
        Index == "F_R.1" ~ "bold(F[R1])",
        Index == "F_R.2" ~ "bold(F[R2])",
        Index == "F_Total" ~ "bold(F[G])",
        Index == "SSB_S.1" ~ "bold(SSB[R1])",
        Index == "SSB_S.2" ~ "bold(SSB[R2])",
        Index == "SSB_Total" ~ "bold(SSB[G])"
      )
    )
  
  # Reorder Index to the desired level order
  data_summary <- data_summary %>%
    mutate(
      Index = factor(
        Index,
        levels = c("bold(C[R1])", "bold(C[R2])", "bold(C[G])",
                   "bold(F[R1])", "bold(F[R2])", "bold(F[G])",
                   "bold(SSB[S1])", "bold(SSB[S2])", "bold(SSB[G])")
      )
    )
  
  create_plot(data_summary, boxplot = F, file_name = filename)
}

#### Bias SSB and Rec Plot ####
n_oms = 12
# KOBE Plot
res1 <- NULL
res2 <- NULL

#### use a new method for the F reference points ####
for (om in 1:n_oms) {
  mods <- readRDS(paste0("OM_", om, ".RDS"))
  print(om)
  for (i in 1:length(mods)){
    for (j in 1:length(mods[[1]])){
      
      if(j <= 4){
        true <- mods[[i]][[j]]$om$rep$SSB[11:30,]
        true <- cbind(true,rowSums(true))
        tmp <- mods[[i]][[j]]$em_list[[1]]$SSB
        ssb.re = tmp/true[,3]-1
        ssb.re = cbind(mtrx = matrix(NA,20,2),ssb.re)
        
        true <- cbind(mods[[i]][[j]]$om$rep$NAA[1,1,11:30,1],mods[[i]][[j]]$om$rep$NAA[2,2,11:30,1])
        true <- cbind(true,rowSums(true))
        tmp <- matrix(mods[[i]][[j]]$em_list[[1]]$NAA[1,1,,1])
        rec.re = tmp/true[,3]-1
        rec.re = cbind(mtrx = matrix(NA,20,2),rec.re)
      } 
      
      if(j > 4 && j <= 6) {
        true <- mods[[i]][[j]]$om$rep$SSB[11:30,]
        true <- cbind(true,rowSums(true))
        tmp1 <- mods[[i]][[j]]$em_list[[1]][[1]]$SSB
        tmp2 <- mods[[i]][[j]]$em_list[[1]][[2]]$SSB
        tmp <- cbind(tmp1,tmp2,tmp1+tmp2)
        ssb.re = tmp/true-1
        
        true <- cbind(mods[[i]][[j]]$om$rep$NAA[1,1,11:30,1],mods[[i]][[j]]$om$rep$NAA[2,2,11:30,1])
        true <- cbind(true,rowSums(true))
        tmp1 <- matrix(mods[[i]][[j]]$em_list[[1]][[1]]$NAA[1,1,,1])
        tmp2 <- matrix(mods[[i]][[j]]$em_list[[1]][[2]]$NAA[1,1,,1])
        tmp <- cbind(tmp1,tmp2,tmp1+tmp2)
        rec.re = tmp/true-1
      }
      
      if(j > 6) {
        true <- mods[[i]][[j]]$om$rep$SSB[11:30,]
        true <- cbind(true,rowSums(true))
        tmp <- mods[[i]][[j]]$em_list[[1]]$SSB
        tmp <- cbind(tmp, rowSums(tmp))
        ssb.re = tmp/true-1
        
        true <- cbind(mods[[i]][[j]]$om$rep$NAA[1,1,11:30,1],mods[[i]][[j]]$om$rep$NAA[2,2,11:30,1])
        true <- cbind(true,rowSums(true))
        tmp1 <- matrix(mods[[i]][[j]]$em_list[[1]]$NAA[1,1,,1])
        tmp2 <- matrix(mods[[i]][[j]]$em_list[[1]]$NAA[2,2,,1])
        tmp <- cbind(tmp1,tmp2,tmp1+tmp2)
        rec.re = tmp/true-1
      }
      
      ssb.re <- data.frame(ssb.re)
      rec.re <- data.frame(rec.re)
      
      ssb.re$nsim <- i
      ssb.re$OM <- paste(om)
      ssb.re$EM <- paste(j)
      res1 <- rbind(res1, ssb.re)
      
      rec.re$nsim <- i
      rec.re$OM <- paste(om)
      rec.re$EM <- paste(j)
      res2 <- rbind(res2, rec.re)
    }
  }
}

res1[res1$EM == 1,3] = res1[res1$EM == 1,1]
res1[res1$EM == 1,1] = NA

res1[res1$EM == 2,3] = res1[res1$EM == 2,1]
res1[res1$EM == 2,1] = NA

res1[res1$EM == 3,3] = res1[res1$EM == 3,1]
res1[res1$EM == 3,1] = NA

res1[res1$EM == 4,3] = res1[res1$EM == 4,1]
res1[res1$EM == 4,1] = NA

res2[res2$EM == 1,3] = res2[res2$EM == 1,1]
res2[res2$EM == 1,1] = NA

res2[res2$EM == 2,3] = res2[res2$EM == 2,1]
res2[res2$EM == 2,1] = NA

res2[res2$EM == 3,3] = res2[res2$EM == 3,1]
res2[res2$EM == 3,1] = NA

res2[res2$EM == 4,3] = res2[res2$EM == 4,1]
res2[res2$EM == 4,1] = NA

write.csv(res1,"Bias_SSB.csv", row.names = F)
write.csv(res2,"Bias_Rec.csv", row.names = F)


library(ggplot2)
library(tidyr)
library(dplyr)

res1 = read.csv("Bias_SSB.csv")
res2 = read.csv("Bias_Rec.csv")


# Transform data into long format for ggplot
data <- res1 %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  mutate(
    Index = recode(Index,
                   "X1" = "SSB_1",
                   "X2" = "SSB_2",
                   "X3" = "SSB_G")
  )

data <- data %>%
  group_by(OM, EM, Index) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25, na.rm = T),
    Q3 = quantile(Value, 0.75, na.rm = T)
  )

# Define the new names corresponding to each integer value
name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

# Convert the numeric vector to a factor and then map the levels to the new names
data$EM <- factor(data$EM, levels = 1:10, labels = name_map)

# Create a named vector of expressions for the x-axis labels with subscripts
# Named vector of expressions for axis labels (keys = EM names, values = display labels)
labels <- c(
  EM1_NAA   = bquote(bold(PAN[NAA])),
  EM1_noNAA = bquote(bold(PAN[noNAA])),
  EM2_NAA   = bquote(bold(FAA[NAA])),
  EM2_noNAA = bquote(bold(FAA[noNAA])),
  EM3_NAA   = bquote(bold(SEP[NAA])),
  EM3_noNAA = bquote(bold(SEP[noNAA])),
  EM4_NAA   = bquote(bold(SpD[NAA])),
  EM4_noNAA = bquote(bold(SpD[noNAA])),
  EM5_Est   = bquote(bold(SpE[NAA * "," * Est])),
  EM5_Fix   = bquote(bold(SpE[NAA * "," * Fix]))
)


OM_labels_list <- list(
  c("OM 1" = bquote(bold(OM[high]~(F[A]))), "OM 2" = bquote(bold(OM[high]~(F[U]))), "OM 3" = bquote(bold(OM[high]~(F[C])))),
  c("OM 4" = bquote(bold(OM[low]~(F[A]))), "OM 5" = bquote(bold(OM[low]~(F[U]))), "OM 6" = bquote(bold(OM[low]~(F[C])))),
  c("OM 7" = bquote(bold(OM[high]~(F[A]))), "OM 8" = bquote(bold(OM[high]~(F[U]))), "OM 9" = bquote(bold(OM[high]~(F[C])))),
  c("OM 10" = bquote(bold(OM[low]~(F[A]))), "OM 11" = bquote(bold(OM[low]~(F[U]))), "OM 12" = bquote(bold(OM[low]~(F[C])))
  )
)

data$OM <- paste("OM",data$OM)

custom_labeller <- labeller(
  OM = label_parsed
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- data %>% filter(OM %in% names(OM_labels_list[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels_list[[i]]), labels = OM_labels_list[[i]])
  
  # Create the plot for the current subset
  p <- ggplot(om_data, aes(x = EM, y = Median, fill = Index, group = Index)) + # Add group = Var to dodge by each Var level
    # Background shading for specific columns (1, 3, 5, 7, 9, 10)
    geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 1
    geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 3
    geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 5
    geom_rect(aes(xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 7
    geom_rect(aes(xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 9
    geom_rect(aes(xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 10
    # Plot points and error bars with larger dodge width
    geom_point(position = position_dodge(width = 0.9), size = 3, shape = 19) + # Increase dodge width
    geom_errorbar(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.9), width = 0) + # Increase dodge width
    facet_wrap(OM ~ Index, scales = "fixed", ncol = 3, labeller = custom_labeller) +
    labs(
      title = paste("Relative Bias in SSB"),
      x = "Estimation Model", y = ""
    ) +
    scale_x_discrete(labels = labels) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      title = element_text(face = "bold")
    ) +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed")
  
  # Print the plot independently
  print(p)
  
  # Save the plot independently with a unique filename
  ggsave(filename = paste0("Bias_SSB_OM_", i, ".png"), plot = p, width = 10, height = 7, units = "in")
}


library(ggplot2)
library(tidyr)
library(dplyr)
# Transform data into long format for ggplot
data <- res2 %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  mutate(
    Index = recode(Index,
                   "X1" = "Rec_1",
                   "X2" = "Rec_2",
                   "X3" = "Rec_G")
  )

data <- data %>%
  group_by(OM, EM, Index) %>%
  summarize(
    Median = median(Value),
    Q1 = quantile(Value, 0.25, na.rm = T),
    Q3 = quantile(Value, 0.75, na.rm = T)
  )

# Define the new names corresponding to each integer value
name_map <- c("EM1_NAA", "EM1_noNAA", "EM2_NAA", "EM2_noNAA", 
              "EM3_NAA", "EM3_noNAA", "EM4_NAA", "EM4_noNAA", 
              "EM5_Est","EM5_Fix")

# Convert the numeric vector to a factor and then map the levels to the new names
data$EM <- factor(data$EM, levels = 1:10, labels = name_map)

# Create a named vector of expressions for the x-axis labels with subscripts
# Named vector of expressions for axis labels (keys = EM names, values = display labels)
labels <- c(
  EM1_NAA   = bquote(bold(PAN[NAA])),
  EM1_noNAA = bquote(bold(PAN[noNAA])),
  EM2_NAA   = bquote(bold(FAA[NAA])),
  EM2_noNAA = bquote(bold(FAA[noNAA])),
  EM3_NAA   = bquote(bold(SEP[NAA])),
  EM3_noNAA = bquote(bold(SEP[noNAA])),
  EM4_NAA   = bquote(bold(SpD[NAA])),
  EM4_noNAA = bquote(bold(SpD[noNAA])),
  EM5_Est   = bquote(bold(SpE[NAA * "," * Est])),
  EM5_Fix   = bquote(bold(SpE[NAA * "," * Fix]))
)


OM_labels_list <- list(
  c("OM 1" = bquote(bold(OM[high]~(F[A]))), "OM 2" = bquote(bold(OM[high]~(F[U]))), "OM 3" = bquote(bold(OM[high]~(F[C])))),
  c("OM 4" = bquote(bold(OM[low]~(F[A]))), "OM 5" = bquote(bold(OM[low]~(F[U]))), "OM 6" = bquote(bold(OM[low]~(F[C])))),
  c("OM 7" = bquote(bold(OM[high]~(F[A]))), "OM 8" = bquote(bold(OM[high]~(F[U]))), "OM 9" = bquote(bold(OM[high]~(F[C])))),
  c("OM 10" = bquote(bold(OM[low]~(F[A]))), "OM 11" = bquote(bold(OM[low]~(F[U]))), "OM 12" = bquote(bold(OM[low]~(F[C])))
  )
)

data$OM <- paste("OM",data$OM)

custom_labeller <- labeller(
  OM = label_parsed
)

# Loop through each subset of OMs and create a plot
for (i in 1:4) {
  # Filter data for the current subset of OMs
  om_data <- data %>% filter(OM %in% names(OM_labels_list[[i]]))
  
  # Set OM factor levels and labels for the current subset
  om_data$OM <- factor(om_data$OM, levels = names(OM_labels_list[[i]]), labels = OM_labels_list[[i]])
  
  # Create the plot for the current subset
  p <- ggplot(om_data, aes(x = EM, y = Median, fill = Index, group = Index)) + # Add group = Var to dodge by each Var level
    # Background shading for specific columns (1, 3, 5, 7, 9, 10)
    geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 1
    geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 3
    geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 5
    geom_rect(aes(xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 7
    geom_rect(aes(xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 9
    geom_rect(aes(xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.1) + # Column 10
    # Plot points and error bars with larger dodge width
    geom_point(position = position_dodge(width = 0.9), size = 3, shape = 19) + # Increase dodge width
    geom_errorbar(aes(ymin = Q1, ymax = Q3), position = position_dodge(width = 0.9), width = 0) + # Increase dodge width
    facet_wrap(OM ~ Index, scales = "fixed", ncol = 3, labeller = custom_labeller) +
    labs(
      title = paste("Relative Bias in Recruitment"),
      x = "Estimation Model", y = ""
    ) +
    scale_x_discrete(labels = labels) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      title = element_text(face = "bold")
    ) +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed")
  
  # Print the plot independently
  print(p)
  
  # Save the plot independently with a unique filename
  ggsave(filename = paste0("Bias_Rec_OM_", i, ".png"), plot = p, width = 10, height = 7, units = "in")
}
