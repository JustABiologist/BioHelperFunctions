from scipy.interpolate import UnivariateSpline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def import_df(csv_path):
    # Load the CSV into a DataFrame
    df = pd.read_csv(csv_path, names=["Wavelength", "CD"])
    return df

def interpol_and_smooth(df):
    # Extract data
    x = np.array(df["Wavelength"])
    y = np.array(df["CD"])
    
    # Create a spline with a smoothing factor
    spline = UnivariateSpline(x, y)
    
    # Generate new x values
    x_new = np.linspace(191.0, 250)
    
    # Get smoothed y values
    y_smooth = spline(x_new)/10
    
    return x_new, y_smooth

df = import_df("/Users/floriangrun/Downloads/Default Dataset.csv")
x_new, y_smooth = interpol_and_smooth(df)

# Plotting
plt.plot(x_new, y_smooth)
plt.xlabel('Wavelength')
plt.ylabel('CD')
plt.title('Smoothed Interpolation')
plt.show()

# Save to CSV
output_df = pd.DataFrame({'Wavelength': x_new, 'Smoothed_CD': y_smooth})
output_df.to_csv('/Users/floriangrun/Downloads/CD-smoothed-116-inRightScale.csv', index=False, sep="\t")