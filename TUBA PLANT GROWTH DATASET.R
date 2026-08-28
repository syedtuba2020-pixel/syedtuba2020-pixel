# ============================================================
# PlantGrowth Dataset Analysis
# ============================================================


# 1. Load the PlantGrowth dataset

data("PlantGrowth")

cat("\n--- PlantGrowth Dataset ---\n")
print(PlantGrowth)

cat("\n--- First 6 Rows ---\n")
print(head(PlantGrowth))

cat("\n--- Structure of Dataset ---\n")
str(PlantGrowth)


# ============================================================
# 2. Mean plant weight for each treatment group
# ============================================================

mean_weight <- aggregate(weight ~ group,
                         data = PlantGrowth,
                         FUN = mean)

cat("\n--- Mean Plant Weight for Each Treatment ---\n")
print(mean_weight)


# ============================================================
# 3. Box Plot
# ============================================================

boxplot(weight ~ group,
        data = PlantGrowth,
        main = "Plant Weight for Different Treatments",
        xlab = "Treatment Group",
        ylab = "Plant Weight")


# ============================================================
# 4. Violin Plot
# ============================================================

# Install ggplot2 automatically if it is not already installed

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

violin_plot <- ggplot(PlantGrowth,
                      aes(x = group, y = weight)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15) +
  labs(title = "Plant Weight Distribution by Treatment",
       x = "Treatment Group",
       y = "Plant Weight") +
  theme_minimal()

print(violin_plot)


# ============================================================
# 5. Histogram
# ============================================================

hist(PlantGrowth$weight,
     main = "Histogram of Plant Weight",
     xlab = "Plant Weight",
     ylab = "Frequency",
     breaks = 8)


# Histogram for each treatment group

histogram_groups <- ggplot(PlantGrowth,
                           aes(x = weight)) +
  geom_histogram(bins = 6) +
  facet_wrap(~group) +
  labs(title = "Plant Weight Distribution in Each Treatment",
       x = "Plant Weight",
       y = "Frequency") +
  theme_minimal()

print(histogram_groups)


# ============================================================
# 6. One-Way ANOVA
# ============================================================

anova_model <- aov(weight ~ group,
                   data = PlantGrowth)

cat("\n--- One-Way ANOVA Results ---\n")
print(summary(anova_model))


# ============================================================
# 7. Tukey HSD Test
# ============================================================

cat("\n--- Tukey HSD Test ---\n")
print(TukeyHSD(anova_model))


# ============================================================
# 8. Extract ANOVA values
# ============================================================

anova_result <- summary(anova_model)[[1]]

f_value <- anova_result["group", "F value"]
p_value <- anova_result["group", "Pr(>F)"]

cat("\nF-value =", round(f_value, 3), "\n")
cat("P-value =", round(p_value, 4), "\n")


# ============================================================
# 9. Interpretation
# ============================================================

cat("\n--- Interpretation ---\n")

cat("The mean plant weight for the control group is",
    round(mean(PlantGrowth$weight[PlantGrowth$group == "ctrl"]), 3),
    ".\n")

cat("The mean plant weight for treatment 1 is",
    round(mean(PlantGrowth$weight[PlantGrowth$group == "trt1"]), 3),
    ".\n")

cat("The mean plant weight for treatment 2 is",
    round(mean(PlantGrowth$weight[PlantGrowth$group == "trt2"]), 3),
    ".\n\n")

cat("Treatment 2 has the highest average plant weight, while",
    "treatment 1 has the lowest average plant weight.\n")

cat("The box plot and violin plot show the differences in plant",
    "weight distribution among the three treatment groups.\n")

cat("The histogram shows the overall distribution of plant weights",
    "and how the observations are spread across the groups.\n\n")


# Interpret the ANOVA result automatically

if (p_value < 0.05) {
  
  cat("The ANOVA p-value is", round(p_value, 4),
      ", which is less than 0.05.\n")
  
  cat("Therefore, there is a statistically significant difference",
      "in mean plant weight among the treatment groups.\n")
  
} else {
  
  cat("The ANOVA p-value is", round(p_value, 4),
      ", which is greater than 0.05.\n")
  
  cat("Therefore, there is no statistically significant difference",
      "in mean plant weight among the treatment groups.\n")
}


# ============================================================
# 10. Conclusion
# ============================================================

cat("\n--- Conclusion ---\n")

if (p_value < 0.05) {
  
  cat("The analysis shows that treatment has a significant effect",
      "on plant weight.\n")
  
  cat("Treatment 2 produced the highest average plant weight,",
      "whereas treatment 1 produced the lowest average plant weight.\n")
  
  cat("The one-way ANOVA showed a statistically significant",
      "difference among the three groups (p < 0.05).\n")
  
  cat("Therefore, the type of treatment appears to influence",
      "plant growth.\n")
  
} else {
  
  cat("The analysis did not find a statistically significant",
      "difference in plant weight among the treatment groups.\n")
  
  cat("Therefore, there is not enough evidence to conclude that",
      "the treatments affected plant growth.\n")
}


# ============================================================
# End of Analysis
# ============================================================

