import pandas as pd
from scipy.interpolate import interp1d
from sklearn.metrics import mean_squared_error
from math import sqrt
import matplotlib.pyplot as plt

# File paths after re-upload
file_to_process = '~/Downloads/CD-smoothed-116-3.csv'
file_to_plot = '~/Downloads/preds.csv'

# Reload the data
df_to_process = pd.read_csv(file_to_process, delimiter='\t')
df_to_plot = pd.read_csv(file_to_plot)

# Sort the 'preds.csv' for interpolation
df_to_plot_sorted = df_to_plot.sort_values(by='WL')

# Interpolate preds.csv to match the wavelength points of CD-smoothed-116-3.csv
interp_function = interp1d(df_to_plot_sorted['WL'], df_to_plot_sorted.iloc[:, 1], kind='linear', fill_value='extrapolate')
preds_interpolated = interp_function(df_to_process['Wavelength'])

# Scale the 'Smoothed_CD' values by dividing by 10
df_to_process['Smoothed_CD'] = df_to_process['Smoothed_CD']

# Calculate MSE and RMSE with the scaled values
mse_scaled = mean_squared_error(df_to_process['Smoothed_CD'], preds_interpolated)
rmse_scaled = sqrt(mse_scaled)

# Plot the two datasets together
plt.figure(figsize=(10, 6))
plt.plot(df_to_process['Wavelength'], df_to_process['Smoothed_CD'], label='CD-smoothed-116-3 Scaled')
plt.plot(df_to_plot_sorted['WL'], df_to_plot_sorted.iloc[:, 1], label='Preds', linestyle='--')
plt.xlabel('Wavelength')
plt.ylabel('CD')
plt.title(f'CD-smoothed-116-3 vs. Preds (RMSE: {rmse_scaled:.2f})')
plt.legend()
plt.show()
