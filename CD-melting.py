import argparse
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Define the function to fit
def uv_hairpin(t, mds, bds, mss, bss, delH, Tm):
    R = 8.314 / 4184
    exp_factor = np.exp(((1 / Tm) - (1 / (t + 273.16))) * delH / R)
    f = exp_factor / (1 + exp_factor)
    return (((mds * t) + bds) * f) + (((mss * t) + bss) * (1 - f))

# Define command line arguments
parser = argparse.ArgumentParser(description='Plot UV hairpin curve fit.')
parser.add_argument('data_files', nargs='+', help='Data files to plot')
parser.add_argument('--mds', type=float, default=30, help='Value for mds')
parser.add_argument('--bds', type=float, default=0.85, help='Value for bds')
parser.add_argument('--mss', type=float, default=0.0006, help='Value for mss')
parser.add_argument('--bss', type=float, default=0.95, help='Value for bss')
parser.add_argument('--delH', type=float, default=-40, help='Value for delH')
parser.add_argument('--Tm', type=float, default=340, help='Value for Tm')
args = parser.parse_args()

# Load the data and plot
fig, axs = plt.subplots(nrows=2, sharex=True)
for data_file in args.data_files:
    data = np.loadtxt(data_file)
    t, absorbance = data[:, 0], data[:, 1]
    #absorbance_norm = absorbance / np.max(absorbance)
    absorbance_norm = (absorbance - np.min(absorbance)) / (np.max(absorbance) - np.min(absorbance))
    popt, pcov = curve_fit(uv_hairpin, t, absorbance, p0=[args.mds, args.bds, args.mss, args.bss, args.delH, args.Tm])
    print("delH:", popt[4])
    print("Tm:", popt[5])
    delS = popt[4]/popt[5] 
    delG_25 = popt[4] - 298.15 * delS
    delG_37 = popt[4] - 310.15 * delS
    print("delS:", delS)
    print("delG_25:", delG_25)
    print("delG_37:", delG_37)
    popt_norm, pcov_norm = curve_fit(uv_hairpin, t, absorbance_norm, p0=[args.mds, args.bds, args.mss, args.bss, args.delH, args.Tm])
    delS = pcov[4]/pcov[5]
    axs[0].plot(t, absorbance, 'o', label=data_file)
    axs[0].plot(t, uv_hairpin(t, *popt), label=f'{data_file} fit')
    axs[1].plot(t, absorbance_norm, 'o', label=data_file)
    axs[1].plot(t, uv_hairpin(t, *popt_norm), label=f'{data_file} fit')
axs[0].set_ylabel('CD mdeg')
axs[1].set_ylabel('Fraction folded')
axs[1].set_xlabel('Temperature ($^\circ$C)')
axs[0].legend()
axs[1].legend()
axs[0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0].yaxis.set_minor_locator(AutoMinorLocator())
axs[1].xaxis.set_minor_locator(AutoMinorLocator())
axs[1].yaxis.set_minor_locator(AutoMinorLocator())




plt.show()

