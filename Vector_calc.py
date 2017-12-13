# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 11:11:13 2017

@author: Joni & Antti
"""

import pandas as pd
import numpy as np
import math
import scipy
import scipy.integrate as integrate
import scipy.stats.vonmises as vonmises
import scipy.stats  as st
import matplotlib.pyplot as plt
import collections
#from decimal import *
#import mplstereonet
#import sympy as sp

"""
Tähän skripti, joka poimii datan eism. excel-taulukosta dataframeksi ja/tai 
np.array -matriisiksi
"""

# Read data from .xlsx to Pandas dataFrame
fp = r"C:\Users\Antti\Desktop\Uni\Maisterivaihe\GeoGradu\Aineisto\Rakoilu-detail_karsittu.xlsx"
data = pd.read_excel(fp, sheet_name='Sheet1', na_values='', )  # index_col='ID')

'''
# Vector = np.array([r, theta, phi])  / r in meters, theta and phi in degrees
test = np.array([1, data['RR1aDir'].iloc[0], data['RR1aDip'].iloc[0]])
F = np.array([1, 45, 150])                                                                  # Spherical coordinates for fracture dip 
No = np.array([1, 90, 90])                                                                  # Spherical coordinates for outcrop normal 
Ot = np.array([1, 0, 90])                                                                   # Spherical coordinates for outcrop trend

F = sphere2cartesian(F)                                                                     # Fracture dip vector
No = sphere2cartesian(No)                                                                   # Outcrop normal vector
Ot = sphere2cartesian(Ot)                                                                   # Outcrop trend vector
AD = F - (np.dot(F, No) / (np.linalg.norm(No)**2)) * No                                     # Apparent dip vector
ad = np.rad2deg(np.arccos(np.dot(AD, Ot)/(np.linalg.norm(AD)* np.linalg.norm(Ot))))         # Apparent dip angle

# Columns for modifying BaseData-file 
data['alphaR1a'] = ''
data['alphaR1b'] = ''
data['alphaR2a'] = ''
data['alphaR2b'] = ''
data['alphaR3a'] = ''
data['alphaR3b'] = ''
data['alphaR4a'] = ''
data['alphaR4b'] = ''
data['alphaR4c'] = ''

dirList = ['RR1aDir','RR1bDir',
           'RR2aDir','RR2bDir',
           'RR3aDir','RR3bDir',
           'RR4aDir','RR4bDir','RR4cDir']
dipList = ['RR1aDip','RR1bDip',
           'RR2aDip','RR2bDip',
           'RR3aDip','RR3bDip',
           'RR4aDip','RR4bDip','RR4cDip']
'''


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


'''   
def wangFAalpha(alpha, rho, x, kappa):
    if alpha < rho:
        Rd = [rho-alpha, rho+alpha]
    else:
        Rd = [0,alpha+rho]
    f_a = (1/np.pi)*integrate.quad(lambda x: ((np.sin(alpha))/(np.sqrt(np.sin**2(x))*np.sin**2(rho)*rho-(np.cos(alpha)-np.cos(x)*np.cos(rho))**2))*
           ((kappa*np.exp**kappa*np.cos(x))/(np.exp**kappa-np.exp**(-kappa))), Rd[0], Rd[1]) 
    return f_a

def wangFC13(aArray,rhoArray,kappaArray):
    C13List = []
    for a in aArray:   
        c13 = (integrate.quad(lambda x: abs(np.cos(x))*wangFAalpha(a,rhoArray[aArray.index(a)],x,kappaArray[aArray.index(a)])), 0, np.pi)
        c13 = round(c13[0],2)
        C13List.append(c13)
    c13 = C13List.mean()
    return c13

def PDFVonmises(alpha, kappa):

def PDF():
    # Takes in array-like, returns ndarray
    st.rv_continuous(momtype=0,a=0,b=alpha)
'''
# Returns the reciprocal of a value
def reciprocal(n):
    reciprocal = 1.0 / n
    return reciprocal


# Returns theta for f_a_alpha
def TC13(aArray):
    # alphaArray = all alphas for an area as a dataFrame
    amin = min(aArray)  # amin = lowest value of alpha
    # Theta for C13
    tc13 = 1 + len(aArray) / np.sum(np.log(aArray / amin))
    return tc13


# Returns the value of the pdf of alpha at the point x
def f_a_alpha(aArray, x, TC13):
    # According to power-law
    a_min = min(aArray)  # a_min = lowest value of alpha
    fa = (TC13 - 1) / a_min * (x / a_min) ** TC13
    return fa


# Returns the conversion factor C13
def fC13(aArray):
    c13 = (integrate.quad(lambda x: abs(np.cos(x)) * f_a_alpha(aArray, x, TC13(aArray)), 0, np.pi))
    c13 = c13[0]
    c13 = reciprocal(c13)
    return c13


# Alternative, clearer way to express the integral (still gives tuples, integral and error margin)
def f(x):
    c13 = abs(np.cos(x)) * f_a_alpha(aArray, x, TC13(aArray))
    c13 = reciprocal(c13.item(0))
    return c13

def make_angle_pdf(angles,uncerts,alpha,pdf_range):
    pdf = np.zeros(len(pdf_range))
    spdf = np.zeros(len(pdf_range))
    for agei in range(len(angles)):
        for pdfagei in range(len(pdf_range)):
            pdf[pdfagei] = (1.0/(alpha*uncerts[agei]*np.sqrt(2.0*np.pi)))*np.exp(-0.5*((pdf_range[pdfagei]-angles[agei])/(alpha*uncerts[agei]))**2.0)
            spdf[pdfagei] += pdf[pdfagei]
    spdf /= float(len(angles))
    return spdf

def make_angle_cdf(spdf,dx):
    cdf = np.zeros(len(spdf))
    for agei in range(len(spdf)):
        if agei == 0:
            cdf[agei] = spdf[agei]
        else:
            cdf[agei] = cdf[agei-1] + (spdf[agei]+spdf[agei-1])/2.0 * dx
    return cdf


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
alphaCols = ['alphaR1', 'alphaR2', 'alphaR3', 'alphaR4']

# Calculate C13 for each area // C13 as a function of area AND population! FIX
# Save conversion factors to a DataFrame

for key in DataFrameDict:
    tempFrame = DataFrameDict[key]
    tempFrame = tempFrame[alphaCols]
    tempFrame[:].replace('', np.nan, inplace=True)  # Replace empty cells with NaN values
    columnNames = tempFrame.columns.values.tolist()

    for name in columnNames:
        tempName = tempFrame[name]  # Get one Series (fracture population) at a time
        tempName = tempName.dropna()  # Drop NaN values from the Series
        aArray = tempName.values  # Convert the Series into a numpy array
        aArray = aArray.astype(float)  # Convert values into float

        if len(aArray) > 1:  # Test if the array has more than 1 values (not possible to integrate)
            C13 = fC13(aArray)  # Calculate C13 for the array
#            C13 = integrate.quad(f,0,np.pi)
#            C13 = wangFC13(aArray)
            C13DF.loc[name[5:], key] = C13  # Save the C13 value into the DataFrame
            print(C13DF)

        else:
            pass


######################
# Plotting
######################

N = len(DataFrameDict['Area 1'])
area1 = DataFrameDict['Area 1']
p,x = np.histogram(area1['alphaR1'],bins=N)
x = x[:-1] + (x[1] - x[0])/2   # convert bin edges to centers
f = scipy.interpolate.UnivariateSpline(x, p, s=0)
#plt.plot(x, f(x))
#plt.show()
"""
hist = np.histogram(area1['alphaR1'], bins=N)
hist_dist = scipy.stats.rv_histogram(hist)
X = np.linspace(-5.0, 5.0, 100)
plt.title("PDF from Template")
plt.hist(data, normed=True, bins=100)
plt.plot(X, hist_dist.pdf(X), label='PDF')
plt.plot(X, hist_dist.cdf(X), label='CDF')
plt.show()
"""
area1 = area1['alphaR1'].values
uncerts = np.zeros(len(area1))+2
alpha = 0.6
pdf_range = np.linspace(0,180,1001)
dx = pdf_range[1]-pdf_range[0]
print(area1)
type(area1)
pdf_angle = make_angle_pdf(area1,uncerts,alpha,pdf_range)

cdf_angle = make_angle_cdf(pdf_angle,dx)
#plt.plot(pdf_range,pdf_angle)
plt.plot(pdf_range, cdf_angle)
plt.show()