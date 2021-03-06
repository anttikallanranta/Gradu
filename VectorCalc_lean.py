# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 11:11:13 2017
@author: Joni & Antti
"""

import pandas as pd
import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import collections
from pylab import figure, show, legend, ylabel, xlabel    



# Read data from .xlsx to Pandas dataFrame
fp = r"C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\Rakoilu-detail_karsittu.xlsx"
data = pd.read_excel(fp, sheet_name='Sheet1', na_values='', )  # index_col='ID')

######################################
# Define all functions
######################################

# Convert spherical coordinates to cartesian
def sphere2cartesian(array):
    nDecimal = 4
    r = array[0]
    theta = array[1]
    phi = array[2]
    x = r * np.round(np.sin(np.deg2rad(phi)), decimals=nDecimal) * np.round(np.cos(np.deg2rad(theta)),
                                                                            decimals=nDecimal)
    y = r * np.round(np.sin(np.deg2rad(phi)), decimals=nDecimal) * np.round(np.sin(np.deg2rad(theta)),
                                                                            decimals=nDecimal)
    z = r * np.round(np.cos(np.deg2rad(phi)), decimals=nDecimal)
    matrix = np.array([x, y, z])
    return matrix


# Calculate apparent dip
def apparentDip(R, No, Ot):
    AD = R - (np.dot(R, No) / (np.linalg.norm(No) ** 2)) * No  # Apparent dip vector
    ad = np.rad2deg(np.arccos(np.dot(AD, Ot) / (np.linalg.norm(AD) * np.linalg.norm(Ot))))  # Apparent dip angle
    return ad


def aDVector(R, No):
    AD = R - (np.dot(R, No) / (np.linalg.norm(No) ** 2)) * No
    return AD


def vectorAngle(array1, array2):
    alpha = np.rad2deg(np.arccos(np.dot(array1, array2) / (np.linalg.norm(array1) * np.linalg.norm(array2))))
    return alpha


# Returns the reciprocal of a value
def reciprocal(n):
    reciprocal = 1.0 / n
    return reciprocal


def _cdf(angles, x):
    return np.sum(angles <= x) / np.sum(~np.isnan(angles))


def make_rad_cdf_2(angles, n):
    cdf_out = np.zeros(n)
    cdf_range = np.linspace(0, np.pi, n)
    da = np.pi / (n-1)
    for i in range(n):
        cdf_out[i] = _cdf(np.pi*angles/180, (i+1)*da)
    return (cdf_out, cdf_range)


def make_deg_cdf_2(angles, n):
    cdf_out = np.zeros(n)
    cdf_range = np.linspace(0, 180, n)
    da = 180 / (n-1)
    for i in range(n):
        cdf_out[i] = _cdf(angles, (i+1)*da)
    return (cdf_out, cdf_range)


def make_rad_pdf_2(cdf, cdfr):
    n = len(cdf) - 1
    pdf_out = np.zeros(n)
    pdf_range = np.linspace(0, np.pi, n)
    for i in range(n-1):
        pdf_out[i] = (cdf[i+1] - cdf[i]) / (cdfr[i+1]-cdfr[i])
    return (pdf_out, pdf_range)


def make_deg_pdf_2(cdf, cdfr):
    n = len(cdf) - 1
    pdf_out = np.zeros(n)
    pdf_range = np.linspace(0, 180, n)
    for i in range(n-1):
        pdf_out[i] = (cdf[i+1] - cdf[i]) / (cdfr[i+1]-cdfr[i])
    return (pdf_out, pdf_range)


###################################


# Create the lists needed in the following loops // For fracture population direction and dip, and Fisher mean pole direction and dip
dirList = ['R1Dir', 'R2Dir', 'R3Dir', 'R4Dir']
dipList = ['R1Dip', 'R2Dip', 'R3Dip', 'R4Dip']
fishDirList = ['R1fDir', 'R2fDir', 'R3fDir', 'R4fDir']
fishDipList = ['R1fDip', 'R2fDip', 'R3fDip', 'R4fDip']

# Create DataFrame for Fisher mean pole direction and dip and kappa
fisherDF = pd.DataFrame(
    {'Domain': ['Area 1', 'Area 2', 'Area 3'], 'R1fDir': [128.08, 167.47, None], 'R2fDir': [19.35, 268.89, None],
     'R3fDir': [114.09, 107.31, None], 'R4fDir': [None, None, None],
     'R1fDip': [2.8, 5.47, None], 'R2fDip': [3.51, 1.15, None], 'R3fDip': [82.25, 75.07, None],
     'R4fDip': [None, None, None],
     'R1kappa': [1.22, 1.28, None], 'R2kappa': [1.60, 1.27, None], 'R3kappa': [30.36, 14.27, None],
     'R4kappa': [None, None, None]})

# Merge Fisher values into data
data = data.merge(fisherDF, how='outer', on='Domain')
data = data.sort_values(['Domain', 'Id'])
data = data.reset_index(drop=True)

##############################
# Begin loop
# Calculate alpha and rho
# Update the corresponding column
##############################

for index, row in data.iterrows():
    r = 1

    ##############################
    # Horizontal
    # Iterate over the data by row by row if outcrop normal = NaN -> Horizontal outcrop
    # Calculate scanline vectors
    # Calculate the angle alpha between scanline vector and fracture population normal
    # Calculate rho between scanline and Fisher mean pole
    # R = Fracture population dip vector
    ##############################

    if math.isnan(row['Outcrop normal']):
        for i in dirList:
            R = sphere2cartesian(np.array([r, row[i], 90]))
            alpha = vectorAngle(R, sphere2cartesian(np.array([r, row[i], row[dipList[dirList.index(i)]] + 180])))
            if math.isnan(alpha):
                pass
            else:
                data.loc[index, 'alpha{}'.format(i[:2])] = alpha
                if row[fishDirList[dirList.index(i)]] != None:
                    rho = vectorAngle(R, sphere2cartesian(
                        np.array([r, row[fishDirList[dirList.index(i)]], row[fishDipList[dirList.index(i)]] + 90])))
                    data.loc[index, 'rho{}'.format(i[:2])] = rho
                else:
                    pass

    ########################
    # Vertical
    # Get apparentDip for all fracture populations
    # Get apparentDipVector for all fracture populations
    # Calculate scanline vectors
    # Calculate the angle alpha between scanline and fracture population normal
    # Calculate rho between scanline and Fisher mean pole
    # No = Outcrop normal
    # Ot = Outcrop trend
    # R = Fracture population dip vector
    # S = Scanline vector
    ########################

    else:
        No = sphere2cartesian(np.array([r, row['Outcrop normal'], 90]))
        Ot = sphere2cartesian(np.array([r, row['Outcrop trend'], 90]))

        for i in dirList:
            R = sphere2cartesian(np.array([r, row[i], row[dipList[dirList.index(i)]] + 90]))
            ad = apparentDip(R, No, Ot)
            if math.isnan(ad):
                pass
            else:
                S = np.cross(aDVector(R, No), No)
                normal = sphere2cartesian((np.array([r, row[i], row[dipList[dirList.index(i)]] + 180])))
                alpha = vectorAngle(S, normal)
                data.loc[index, 'alpha{}'.format(i[:2])] = alpha
                if row[fishDirList[dirList.index(i)]] != None:
                    rho = vectorAngle(S, sphere2cartesian(
                        np.array([r, row[fishDirList[dirList.index(i)]], row[fishDipList[dirList.index(i)]] + 90])))
                    data.loc[index, 'rho{}'.format(i[:2])] = rho
                else:
                    pass

#########################
# C13
#########################


# DataFrame dictionary for the use of TC13, f_a_alpha and C13
uniqueAreas = data.Domain.unique()
DataFrameDict = {elem: pd.DataFrame for elem in uniqueAreas}
for key in DataFrameDict.keys():
    DataFrameDict[key] = data[:][data.Domain == key]
DataFrameDict = collections.OrderedDict(
    sorted(DataFrameDict.items()))  # Sort the Dictionary so everything works in order

# Create DataFrame for the C13 values
C13DF = pd.DataFrame(
    {'Population': ['R1', 'R2', 'R3', 'R4'], 'Area 1': [None, None, None, None], 'Area 2': [None, None, None, None],
     'Area 3': [None, None, None, None]})
C13DF = C13DF.set_index('Population')
alphaCols = ['alphaR1', 'alphaR2', 'alphaR3']

# Calculate C13 for each area // C13 as a function of area AND population! FIX
# Save conversion factors to a DataFrame

for key in DataFrameDict:
    tempFrame = DataFrameDict[key]
    tempFrame = tempFrame[alphaCols]
    
    for name in alphaCols:
        tempName = tempFrame[name].values
        tempName = tempName[~np.isnan(tempName)]
        tempName = np.sort(tempName)
        uncerts = np.zeros(len(tempName))+2
        alpha = 0.6
        pdf_range = np.linspace(0,180,10)
        dx = pdf_range[1]-pdf_range[0]
        print(tempName)
        cdf_angle, cdf_range = make_rad_cdf_2(tempName,10)
        pdf_angle, pdf_range = make_rad_pdf_2(cdf_angle, cdf_range)
        C13 = integrate.trapz(pdf_angle*np.abs(np.cos(pdf_range)), pdf_range)
        C13DF.loc[name[5:], key] = C13
        print(C13DF)
        outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
        fname = outfp + '\\' + key + '_' + name + '_combi' + '.png'

######################
        # Plotting
######################
        
        # create the general figure
        fig1 = figure()       
        # and the first axes using subplot populated with data 
        ax1 = fig1.add_subplot(111)
        line1 = ax1.plot(pdf_range, pdf_angle, 'o-')
        ylabel("PDF", color = 'b') 
        xlabel("Radians")
        # now, the second axes that shares the x-axis with the ax1
        ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
        line2 = ax2.plot(cdf_range, cdf_angle, 'xr-')
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        ylabel("CDF", color = 'r')        
        # for the legend, remember that we used two different axes so, we need 
        # to build the legend manually
#        legend((line1, line2), ("1", "2"))
        plt.title(key + ' ' + name, loc='center')
        plt.savefig(fname, dpi=300)
        show()

'''            
        plt.plot(pdf_range, pdf_angle)
        plt.plot(cdf_range, cdf_angle)
        plt.xlabel('Radians')
        plt.ylabel('CDF')
        plt.title(key + ' ' + name, loc='center')
        plt.savefig(fname, dpi=300)
        plt.show()
'''

'''
Normalisoi CDF ja PDF jakamalla arvot suurimmalla arvolla. Suurin = 1
Tee kaikille populaatioille PDF, syötä integraaliin
'''        
