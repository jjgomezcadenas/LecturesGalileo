import numpy as np
import pandas as pd
from scipy.integrate import quad
from  . system_of_units import *
import matplotlib.pyplot as plt


def gaussian(x, mean, sigma):
    """
    A gaussian function 
    """
    return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def integrate_gaussian_from_x0(mean, sigma, x0):
    """
    Integrate Gaussian from x0 to infinity

    """
    result, _ = quad(lambda x: gaussian(x, mean, sigma), x0, np.inf)
    return result


# graphics

def plot_gaussians(mean1, mean2, sigma1, sigma2, title):
    # Define the range of x values (energy in keV)
    x = np.linspace(2400, 2500, 1000)
    
    
    # Calculate the Gaussian curves
    y1 = gaussian(x, mean1, sigma1)
    y2 = gaussian(x, mean2, sigma2)
    
    # Plot the Gaussians
    plt.figure(figsize=(8, 4))
    plt.plot(x, y1, label=f'Bi-214 at {mean1} keV', color='blue')
    plt.plot(x, y2, label=f'bb2nu at {mean2} keV', color='red')
    
    # Add labels and legend
    plt.xlabel('Energy (keV)')
    plt.ylabel('Normalized Intensity')
    plt.title(title)
    plt.legend()
    
    # Show the plot
    plt.grid()
    plt.show()


def energy_window_eff(eg, fbb, fbi):
    plt.figure(figsize=(8, 4))
    plt.plot(eg, fbb, label=f'bb2u', color='blue')
    plt.plot(eg, fbi, label=f'bi-214', color='red')
       
    # Add labels and legend
    plt.xlabel('Energy (keV)')
    plt.ylabel('Normalized Intensity')
    plt.legend()
    
    # Show the plot
    plt.grid()
    plt.show()


def energy_window_fom(eg, fbb, fbi):
    plt.figure(figsize=(8, 4))
    plt.plot(eg, fbb/np.sqrt(fbi), label=f'bb2u', color='blue')
    #plt.plot(eg, fbi, label=f'bi-214', color='red')
       
    # Add labels and legend
    plt.xlabel('Energy (keV)')
    plt.ylabel('Normalized Intensity')
    plt.legend()
    
    # Show the plot
    plt.grid()
    plt.show()