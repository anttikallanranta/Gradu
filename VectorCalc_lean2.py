# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 11:11:13 2017
@author: Joni & Antti
"""

# %%
import pandas as pd
import numpy as np
import math
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import collections
from pylab import figure, show, legend, ylabel, xlabel    
import mplstereonet
import mplFunctions as mpl
import matplotlib.lines as mlines


# Read data from .xlsx to Pandas dataFrame
#fp = r"C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\Rakoilu-detail_karsittu.xlsx"
fp2 =  r"C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\Rakoilu-detail_karsittu_3.xlsx"
#data = pd.read_excel(fp, sheet_name='Sheet1', na_values='', )  # index_col='ID')
data2 = pd.read_excel(fp2, sheet_name='OnlyArea1', na_values='', )


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


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def vectorAngle(array1, array2):
    #alpha = np.rad2deg(np.arccos(np.dot(array1, array2) / (np.linalg.norm(array1) * np.linalg.norm(array2))))
    alpha = np.rad2deg(np.arccos(np.clip(np.dot(unit_vector(array1), unit_vector(array2)), -1.0, 1.0)))
    return alpha


# Returns the reciprocal of a value
def reciprocal(x):
    reciprocal = 1.0 / x
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
    for i in range(n):
        pdf_out[i] = (cdf[i+1] - cdf[i]) / (cdfr[i+1]-cdfr[i])
    return (pdf_out, pdf_range)


def plot_stereonet(strikes, dips, domain):
    bin_edges = np.arange(-5, 366, 10)
    number_of_strikes, bin_edges = np.histogram(strikes, bin_edges)
    number_of_strikes[0] += number_of_strikes[-1]
    half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
    two_halves = np.concatenate([half, half])

    fig = plt.figure(figsize=(16,8))
    
    ax = fig.add_subplot(121, projection='stereonet')

    ax.pole(strikes, dips, c='k', label='Poles of the Planes')
    ax.density_contourf(strikes, dips, measurement='poles', cmap='Reds', sigma = 10)
    ax.set_title('Density contour of the Poles', y=1.10, fontsize=15)
    ax.grid()
    
    ax = fig.add_subplot(122, projection='polar')
    
    ax.bar(np.deg2rad(np.arange(0, 360, 10)), two_halves,
           width=np.deg2rad(10), bottom=0.0, color='.8', edgecolor='k')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10))
    ax.set_rgrids(np.arange(1, two_halves.max() + 1, 10), angle=0, weight= 'black')
    ax.set_title('Rose Diagram of the Bearings', y=1.10, fontsize=15)
    fig.suptitle(domain)
#    fig.tight_layout()
#    show()


def weightByN(dataframe):
    weightcols =list(dataframe.columns.values)
    frameList = dataframe.values.tolist()
    list2 = []
    for idx, row in dataframe.iterrows():
        if row['n_fractures'] == 0:
            pass
        else:
            for i in range(row['n_fractures']):
                list2.append(frameList[idx])
    frame = pd.DataFrame(list2, columns=weightcols)
    return frame


def splitByAngle(dframe):
    testPoles = pd.DataFrame({'strike':[210,10,190,295,115], 'dip':[10,90,90,90,90]})

    for idx, row in dframe.iterrows():

        lon, lat = mpl.pole((row['strike']), row['dip'])
        for i, r in testPoles.iterrows():
            angle = np.degrees(mpl.angular_distance((lon,lat), (mpl.pole(r['strike'],r['dip'])), bidirectional=False))
            if (angle <= np.array([40])).all():
                dframe.loc[idx, 'population'] = int(i)
                break
            else:
                pass
    return dframe


def correctPop(dframe):
    for idx, row in dframe.iterrows():
        if row['population'] == 0:
            dframe.loc[idx, 'population'] = 1
            continue
        elif row['population'] == 1 or row['population'] == 2:
            dframe.loc[idx, 'population'] = 2
            continue
        elif row['population'] == 3 or row['population'] == 4:
            dframe.loc[idx, 'population'] = 3
            continue
        else:
            pass
    return dframe


def getSmallAlpha(S,normal):
    alpha1 = vectorAngle(S, normal)
    alpha2 = vectorAngle(S*-1, normal)
    alpha3 = vectorAngle(S, normal*-1)
    alpha4 = vectorAngle(S*-1, normal*-1)
    if alpha1 < alpha2 and alpha1 < alpha3 and alpha1 < alpha4:
        return alpha1
    elif alpha1 > alpha2 and alpha2 < alpha3 and alpha2 < alpha4:
        return alpha2
    elif alpha3 < alpha1 and alpha3 < alpha2 and alpha3 < alpha4:
        return alpha3
    else:
        return alpha4


def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]


#######################################################################3

##############################
# Loop over the data, row by row
# Calculate alpha - the solid angle between the scanline and fracture population normal
# Save alpha to the dataframe
##############################

# Begin loop
for index, row in data2.iterrows():
    # r = magnitude of the vector
    r = 1

    ##############################
    # Horizontal
    # Iterate over the data by row by row if outcrop normal = NaN -> Horizontal outcrop
    # Calculate scanline vectors
    # Calculate the angle alpha between scanline vector and fracture population normal
    # Calculate rho between scanline and Fisher mean pole
    # S = scanline vector
    # normal = fracture population normal
    ##############################

    if math.isnan(row['outcrop_normal']):
        # Stereonet coordinates for scanline vector
        Slon, Slat = mpl.line(row['dir'], 90)
        S = unit_vector(mpl.sph2cart(Slon, Slat))

        # Stereonet coordinates for fracture normal vector
        nlon, nlat = mpl.line(row['dir'], 90 - row['dip'])
        normal = unit_vector(mpl.sph2cart(nlon, nlat))

        # Calculate and save alpha
        alpha = np.degrees(mpl.angular_distance((Slon, Slat), (nlon, nlat), bidirectional=True))
        data2.loc[index, 'alpha'] = alpha

########################
# Vertical
# Get apparentDip for all fracture populations
# Get apparentDipVector for all fracture populations
# Calculate scanline vectors
# Calculate the angle alpha between scanline and fracture population normal
# Calculate rho between scanline and Fisher mean pole
# No = Outcrop normal
# Ot = Outcrop strike
# R = Fracture population dip vector
# ad = Apparent dip
# S = Scanline vector
# normal = Fracture plane normal
########################

    else:
        # Stereonet coordinates for outcrop normal vector
        Nolon, Nolat = mpl.pole(row['outcrop_normal'], 90)
        No = unit_vector(mpl.sph2cart(Nolon, Nolat))
        No = No.flatten()

        # Stereonet coordinates for outcrop strike vector
        Otlon, Otlat = mpl.pole(row['outcrop_strike'], 90)
        Ot = unit_vector(mpl.sph2cart(Otlon, Otlat))
        Ot = Ot.flatten()

        # Stereonet coordinates for fracture population dip vector
        Rlon, Rlat = mpl.pole(row['strike'], row['dip'])
        R = unit_vector(mpl.sph2cart(Rlon, Rlat))
        R = R.flatten()

        # Plunge and bearing for line of intersection of the outcrop and fracture planes - apparent dip
        ADplunge, ADbearing = mpl.plane_intersection(row['strike'], row['dip'], row['outcrop_strike'], 90)
        # Stereonet coordinates for apparent dip vector
        ADlon, ADlat = mpl.line(ADplunge, ADbearing)
        #ADstrike, ADdip = mpl.geographic2pole(ADlon, ADlat)

        # Stereonet coordinates for scanline vector
        Slon, Slat = mpl.line(row['outcrop_strike'], ADplunge-90)
        S = unit_vector(mpl.sph2cart(Slon, Slat))
        S = S.flatten()

        # Stereonet coordinates for fracture population normal vector
        nlon, nlat = mpl.line(row['dir'], 90-row['dip'])
        normal = unit_vector(mpl.sph2cart(nlon, nlat))
        normal = normal.flatten()

        # Append fracture normals to arrays
        data2.loc[index, 'lon'] = nlon
        data2.loc[index, 'lat'] = nlat

        # Calculate and save alpha
        alpha = np.degrees(mpl.angular_distance((Slon, Slat), (nlon, nlat), bidirectional = True))
        data2.loc[index, 'alpha'] = alpha



#########################
# Stereoplots and division
# Weighting by density
# Convert population column into int
# Sort by outcrop id
#########################

# Weighting the data by density. Approximate all outcrops as 10 m radius circular disks
data2 = weightByN(data2)

# Classify the data into populations by comparing the poles with a set of test poles
data2 = splitByAngle(data2)

# Replace all NaN values with 999
data2[['population']] = data2[['population']].fillna(999).astype(int)
# Convert population column into int
data2[['population']] = data2[['population']].astype(int)
# Sort data by outcrop index
data2 = data2.sort_values('Id', ascending = True)
# Reclassify the populations, taking repoles into consideration
data2 = correctPop(data2)


# Create stereoplot and roseplot
#grouped = data2.groupby('Domain')
#for idx, group in grouped:
strikes = data2.as_matrix(columns = ['strike'])
dips = data2.as_matrix(columns = ['dip'])
domain = 'Area 1'
plot_stereonet(strikes,dips,domain)
outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
fname = outfp + '\\' + domain + '_' + 'stereo_rose' + '.png'
plt.savefig(fname, dpi=300)


#########################
# C13
#########################

# Empty DataFrame for C13
C13DF2 = pd.DataFrame(columns = ['Population', 'Area 1']) #, 'Area 2': [None], 'Area 3': [None]})
C13DF2 = C13DF2.set_index('Population')

# Group data by population
data2 = data2.sort_values('population', ascending = True)
grouped = data2.groupby('population')

for idx,group in grouped:
    if idx == 999:
        continue

    # Select the column containing alpha values
    tempFrame = group['alpha']
    # Convert pandas series into numpy array
    tempName = tempFrame.values
    # Delete NaN values from the array
    tempName = tempName[~np.isnan(tempName)]
    # Sort the array into ascending order
    tempName = np.sort(tempName)

    #uncerts = np.zeros(len(tempName)) + 2???
    #alpha = 0.6???
    #dx = pdf_range[1] - pdf_range[0]???

    # Calculate Cumulative distribution function values and its range for values of alpha
    cdf_angle, cdf_range = make_rad_cdf_2(tempName, len(tempName) - 1)
    # Calculate Probability Distribution Function values and its range for values of CDF
    pdf_angle, pdf_range = make_rad_pdf_2(cdf_angle, cdf_range)

    # Calculate the conversion factor C13 by integrating over the PDF values
    C13 = reciprocal(integrate.trapz(pdf_angle * np.abs(np.cos(pdf_range)), pdf_range))
    # Save C13 in a dataframe
    C13DF2.loc[idx, 'Area 1'] = C13

    # Print the dataframe for viewing
    print(C13DF2)

    ######################
    # Fisher Stats
    ######################

    # Select stereoplot longitude column as numpy array
    lons = group['lon'].as_matrix()
    # Select stereoplot latitude column as numpy array
    lats = group['lat'].as_matrix()

    # Remove NaN values
    lons = lons[~np.isnan(lons)]
    lats = lats[~np.isnan(lats)]

    # Calculate Fisher statistics
    mean_vec, (r_value, angle, kappa) = mpl.fisher_stats(lons, lats, conf=95)
    meanStrike, meanDip = mpl.geographic2pole(mean_vec[0], mean_vec[1])
    print('Mean_vec strike ', meanStrike, ' dip ', meanDip )
    print('r_value', r_value)
    print('angle', angle)
    print('kappa', kappa)


######################
        # Plotting
######################
    # Create output folder path and filename
    outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
    fname = outfp + '\\' + 'Area 1' + '_' + str(idx) + '_test' + '.png'

    # Create the general figure
    fig1 = figure()
    # Create the first axes using subplot populated with data
    ax1 = fig1.add_subplot(111)

    # Convert radians into degrees for visualization
    pdf_range = np.rad2deg(pdf_range)
    cdf_range = np.rad2deg(cdf_range)

    # Prepare for masking arrays - 'conventional' arrays won't do it
    pdf_angle = np.ma.array(pdf_angle)

    # Mask values where pdf_angle < 0
    pdf_angle_masked = np.ma.masked_where(pdf_angle <= 0, pdf_angle)
    pdf_range_masked = np.ma.masked_where(pdf_angle <= 0, pdf_range)

    # Choose values where pdf_range > 92 degrees, for visualization
    pdf_angle_masked = pdf_angle_masked[pdf_range_masked < 92]
    pdf_range_masked = pdf_range_masked[pdf_range_masked < 92]

    # State the width of the bars
    width = 0.5

    # Plot the first part, PDF as a bar chart
    pdfBar = ax1.bar(pdf_range_masked, pdf_angle_masked, width, color ='b')
    # Set first y axis label and x axis label
    ylabel("Probability Density Function", color = 'b')
    xlabel("Alpha")

    # Add a second axes that shares the x axis with the first one
    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)

    # Make PDF and CDF arrays the same size
    cdf_range = cdf_range[0:-1]
    cdf_angle = cdf_angle[0:-1]

    # Mask values where pdf_range > 92 degrees, for visualization
    cdf_range_masked = np.ma.masked_where(pdf_range > 92, cdf_range)
    cdf_angle_masked = np.ma.masked_where(pdf_range > 92, cdf_angle)

    # Plot the second part, CDF as a line plot
    cdfLine = ax2.plot(cdf_range_masked, cdf_angle_masked, 'r-',)
    # Set second y axis to the right
    ax2.yaxis.tick_right()
    # Set second y axis label position and label
    ax2.yaxis.set_label_position("right")
    ylabel("Cumulative Density Function", color = 'r')

    # Build the legend
    cdfLineLegend = mlines.Line2D([], [], color='r', label='CDF')
    legend((pdfBar, cdfLineLegend), ("PDF", "CDF"))
    plt.title('Area 1' + '_R' + str(idx), loc='center')
    plt.savefig(fname, dpi=300)






#Normalisoi CDF ja PDF jakamalla arvot suurimmalla arvolla. Suurin = 1
#Tee kaikille populaatioille PDF, syötä integraaliin

#fishDirList = ['R1fDir', 'R2fDir', 'R3fDir', 'R4fDir']
#fishDipList = ['R1fDip', 'R2fDip', 'R3fDip', 'R4fDip']

# Create DataFrame for Fisher mean pole direction and dip and kappa
# FisherDF = pd.DataFrame(
#    {'Domain': ['Area 1', 'Area 2', 'Area 3'], 'R1fDir': [128.08, 167.47, None], 'R2fDir': [19.35, 268.89, None],
#     'R3fDir': [114.09, 107.31, None], 'R4fDir': [None, None, None],
#     'R1fDip': [2.8, 5.47, None], 'R2fDip': [3.51, 1.15, None], 'R3fDip': [82.25, 75.07, None],
#     'R4fDip': [None, None, None],
#     'R1kappa': [1.22, 1.28, None], 'R2kappa': [1.60, 1.27, None], 'R3kappa': [30.36, 14.27, None],
#     'R4kappa': [None, None, None]})

# Merge Fisher values into data
#data = data.merge(fisherDF, how='outer', on='Domain')
#data = data.sort_values(['Domain', 'Id'])
#data = data.reset_index(drop=True)


# DataFrame dictionary for the use of TC13, f_a_alpha and C13
#uniqueAreas = data2.Domain.unique()
#DataFrameDict = {elem: pd.DataFrame for elem in uniqueAreas}
#for key in DataFrameDict.keys():
#    DataFrameDict[key] = data2[:][data2.Domain == key]
#DataFrameDict = collections.OrderedDict(
#    sorted(DataFrameDict.items()))  # Sort the Dictionary so everything works in order

