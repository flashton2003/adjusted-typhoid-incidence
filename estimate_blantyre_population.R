## Blantyre

# Data preparation
year <- c(2008, 2018)
population <- c(661256, 800264)

# Create a linear regression model
model <- lm(population ~ year)

# Generate a sequence of years from 2017 to 2024
future_years <- data.frame(year = 2017:2024)

# Predict the population for each year in the sequence
predicted_population <- predict(model, future_years)

# Combine years and predictions into a data frame
predictions <- data.frame(Year = future_years$year, Predicted_Population = predicted_population)

# Print the predictions
print(predictions)

# Optional: Save the predictions to a CSV file
write.csv(predictions, "population_predictions_2017_2024.csv", row.names = FALSE)

## Ndirande

year <- c(2016, 2018, 2022)

population <- c(97392, 102242, 102400)

# Create a linear regression model
model <- lm(population ~ year)

# Generate a sequence of years from 2017 to 2024
future_years <- data.frame(year = 2017:2024)

# Predict the population for each year in the sequence
predicted_population <- predict(model, future_years)

# Combine years and predictions into a data frame
predictions <- data.frame(Year = future_years$year, Predicted_Population = predicted_population)

# Print the predictions
print(predictions)

# Optional: Save the predictions to a CSV file
write.csv(predictions, "population_predictions_2017_2024.csv", row.names = FALSE)


