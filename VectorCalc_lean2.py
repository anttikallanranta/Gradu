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
import matplotlib
import matplotlib.pyplot as plt
import collections
from pylab import figure, show, legend, ylabel, xlabel    
import mplstereonet
import mplFunctions as mpl
import analysis
import stereonet_math
import matplotlib.lines as mlines
import openpyxl
import more_itertools as itertools

# Read data from .xlsx to Pandas dataFrame
fp2 =  r"C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\Rakoilu-detail_karsittu_4.xlsx"
data2 = pd.read_excel(fp2, sheet_name='OnlyArea1', na_values='', )
dataBeni = pd.read_excel(fp2, sheet_name='BeniOnlyArea1', na_values='', )


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
    alpha = np.rad2deg(np.arccos(np.clip(np.dot(unit_vector(array1), unit_vector(array2)), -1.0, 1.0)))
    return alpha


# Returns the reciprocal of a value
def reciprocal(x):
    reciprocal = 1.0 / x
    return reciprocal


def _cdf(angles, x):
    return np.sum(angles <= x) / np.sum(~np.isnan(angles))


def make_rad_cdf(angles, n):
    cdf_out = np.zeros(n)
    cdf_range = np.linspace(0, np.pi, n)
    da = np.pi / (n-1)
    for i in range(n):
        cdf_out[i] = _cdf(np.pi*angles/180, (i+1)*da)
    return (cdf_out, cdf_range)


def make_rad_pdf(cdf, cdfr):
    n = len(cdf) - 1
    pdf_out = np.zeros(n)
    pdf_range = np.linspace(0, np.pi, n)
    for i in range(n-1):
        pdf_out[i] = (cdf[i+1] - cdf[i]) / (cdfr[i+1]-cdfr[i])
    return (pdf_out, pdf_range)


def make_deg_cdf(angles, n):
    cdf_out = np.zeros(n)
    cdf_range = np.linspace(0, 180, n)
    da = 180 / (n-1)
    for i in range(n):
        cdf_out[i] = _cdf(angles, (i+1)*da)
    return (cdf_out, cdf_range)


def make_deg_pdf(cdf, cdfr):
    n = len(cdf) - 1
    pdf_out = np.zeros(n)
    pdf_range = np.linspace(0, 180, n)
    for i in range(n):
        pdf_out[i] = (cdf[i+1] - cdf[i]) / (cdfr[i+1]-cdfr[i])
    return (pdf_out, pdf_range)


def plot_stereonet(strikes, dips, domain, sigma, *args, **kwargs):

    bin_edges = np.arange(-5, 366, 10)
    number_of_strikes, bin_edges = np.histogram(strikes, bin_edges)
    number_of_strikes[0] += number_of_strikes[-1]
    half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
    two_halves = np.concatenate([half, half])

    fig = plt.figure(figsize=(16,8))
    
    ax = fig.add_subplot(121, projection='stereonet')

    ax.pole(strikes, dips, c='k', label='Poles of the Planes')
    ax.density_contourf(strikes, dips, measurement='poles', cmap='Reds', sigma = sigma)
    ax.set_title('Density contour of the Poles', y=1.10, fontsize=15)
    ax.grid()
    
    ax = fig.add_subplot(122, projection='polar')
    
    ax.bar(np.deg2rad(np.arange(0, 360, 10)), two_halves,
           width=np.deg2rad(10), bottom=0.0, color='0.8', edgecolor='k')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10))
    ax.set_rgrids(np.arange(1, two_halves.max() + 1, 10), angle=0, weight= 'black')
    ax.set_title('Rose Diagram of the Bearings', y=1.10, fontsize=15)
    fig.suptitle(domain)


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


def splitByAngle(dframe, cPoles):
    testPoles = cPoles

    # Calculate the minimum angle between the cluster centers
    minangle = np.array([180])
    for idx, row in testPoles.iterrows():
        lon, lat = mpl.pole((row['strike']), row['dip'])
        for i, r in testPoles.iterrows():
            angle = np.degrees(mpl.angular_distance((lon,lat), (mpl.pole(r['strike'],r['dip'])), bidirectional=False))
            if angle > 1 and angle < minangle:
                minangle = angle

    # Divide minimum angle by 2 to avoid overlap
    minangle = minangle/2
    print('Cluster solid angle: ' + str(minangle))

    # Iterate over both DataFrames and classify poles to clusters
    for idx, row in dframe.iterrows():
        lon, lat = mpl.pole((row['strike']), row['dip'])

        for i, r in testPoles.iterrows():
            angle = np.degrees(mpl.angular_distance((lon,lat), (mpl.pole(r['strike'],r['dip'])), bidirectional=False))
            if (angle <= minangle).all():
                dframe.loc[idx, 'population'] = int(i+1)
                break
            else:
                pass

    return dframe


def checkPoles(dataframe):
    # Check if some of the created population poles are part of the same population. I.e. on the other side of the stereonet
    copyFrame = dataframe
    copyFrame['opposite'] = False
    for idx, row in dataframe.iterrows():
        rowlon, rowlat = mpl.pole(row['strike'], row['dip'])
        for i, r in dataframe.iterrows():
            ilon, ilat = mpl.pole(r['strike'], r['dip'])
            angle = np.degrees(mpl.angular_distance((rowlon, rowlat), (ilon, ilat), bidirectional = False))
            print(idx,angle)
            if angle >= 130:
                copyFrame.loc[i,'opposite'] = True
            else:
                continue
    return copyFrame


def correctPop(dframe):
    for idx, row in dframe.iterrows():
        if row['population'] == 0:
            dframe.loc[idx, 'population'] = 1
            continue
        elif row['population'] == 2 or row['population'] == 3:
            dframe.loc[idx, 'population'] = 2
            continue
        elif row['population'] == 4 or row['population'] == 5:
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


def calcAlpha(dataframe):
##############################
# Loop over the data, row by row
# Calculate alpha - the solid angle between the scanline and fracture population normal
# Save alpha to the dataframe
##############################

    # Begin loop
    for index, row in dataframe.iterrows():
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

            # Append fracture normals to arrays
            dataframe.loc[index, 'lon'] = nlon
            dataframe.loc[index, 'lat'] = nlat

            # Calculate and save alpha
            alpha = np.degrees(mpl.angular_distance((Slon, Slat), (nlon, nlat), bidirectional=True))
            dataframe.loc[index, 'alpha'] = alpha

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
            dataframe.loc[index, 'lon'] = nlon
            dataframe.loc[index, 'lat'] = nlat

            # Calculate and save alpha
            alpha = np.degrees(mpl.angular_distance((Slon, Slat), (nlon, nlat), bidirectional = True))
            dataframe.loc[index, 'alpha'] = alpha

    return dataframe


def calcPlotC13(dataframe, type):
    # Empty DataFrame for C13
    C13DF = pd.DataFrame(columns=['population', 'Area 1_' + type])
    C13DF = C13DF.set_index('population')

    # Empty DataFrame for the Fisher values
    clist = ['meanStrike', 'meanDip', 'r_value', 'angle', 'kappa']
    fisherDF = pd.DataFrame(columns = clist)

    # Empty DataFrames for the CDFs
    cdfAngleDF = pd.DataFrame()
    cdfRangeDF = pd.DataFrame()

    # Group data by population
    dataframe = dataframe.sort_values('population', ascending=True)
    grouped = dataframe.groupby('population')

    for idx, group in grouped:
        if idx == 999:
            continue
        print('IDX ', idx)
        # Select the column containing alpha values
        tempFrame = group['alpha']
        # Convert pandas series into numpy array
        tempName = tempFrame.values
        # Delete NaN values from the array
        tempName = tempName[~np.isnan(tempName)]
        # Sort the array into ascending order
        tempName = np.sort(tempName)

        # Calculate Cumulative distribution function values and its range for values of alpha
        cdf_angle, cdf_range = make_rad_cdf(tempName, len(tempName) - 1)
        # Calculate Probability Distribution Function values and its range for values of CDF
        pdf_angle, pdf_range = make_rad_pdf(cdf_angle, cdf_range)

        print(cdf_angle)

        # Calculate the conversion factor C13 by integrating over the PDF values
        C13 = reciprocal(integrate.trapz(pdf_angle * np.abs(np.cos(pdf_range)), pdf_range))
        # Save C13 in a dataframe
        C13DF.loc[idx, 'Area 1_' + type] = C13

        # Print the dataframe for viewing
        print(C13DF)

        ##########################
        # Fisher Stats
        ##########################

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
        print('Mean density', group['meanDens'].mean())
        print('Mean_vec strike ', meanStrike, ' dip ', meanDip)
        print('r_value', r_value)
        print('angle', angle)
        print('kappa', kappa)
        C13DF.loc[idx, 'kappa'] = kappa

        # Create a DF for the Fisher values
        fisherTemp = pd.DataFrame([[meanStrike, meanDip, r_value, angle, kappa]], columns = clist)
        fisherDF = fisherDF.append(fisherTemp)


        ######################
        # Plotting
        ######################
        # Create output folder path and filename
        outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
        fname = outfp + '\\' + 'Area 1_' + str(idx) + '_' + type + '.png'

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
        pdfBar = ax1.bar(pdf_range_masked, pdf_angle_masked, width, color='b')
        # Set first y axis label and x axis label
        ylabel("Probability Density Function", color='b')
        xlabel("Alpha")

        # Add a second axes that shares the x axis with the first one
        ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)

        # Make PDF and CDF arrays the same size
        cdf_range = cdf_range[:-1]
        cdf_angle = cdf_angle[:-1]

        # Mask values where pdf_range > 92 degrees, for visualization
        cdf_range_masked = np.ma.masked_where(pdf_range > 92, cdf_range)
        cdf_angle_masked = np.ma.masked_where(pdf_range > 92, cdf_angle)

        # Plot the second part, CDF as a line plot
        cdfLine = ax2.plot(cdf_range_masked, cdf_angle_masked, 'r-', )
        # Set second y axis to the right
        ax2.yaxis.tick_right()
        # Set second y axis label position and label
        ax2.yaxis.set_label_position("right")
        ylabel("Cumulative Density Function", color='r')

        # Build the legend
        cdfLineLegend = mlines.Line2D([], [], color='r', label='CDF')
        unid = group.Id.unique()
        legend((pdfBar, cdfLineLegend), ("PDF", "CDF, n=" + str(len(unid))))
        plt.title('Area 1' + '_R' + str(idx), loc='center')
        plt.savefig(fname, dpi=300)

        #KORJAA
        cdfAngleDF = pd.concat([cdfAngleDF, pd.Series(cdf_angle.tolist())], ignore_index=True)
        cdfRangeDF = pd.concat([cdfRangeDF, pd.Series(cdf_range.tolist())], ignore_index=True)


    fisherDF = fisherDF.reset_index(drop=True)

    return  C13DF, fisherDF, cdfAngleDF, cdfRangeDF

def groupByDensity(dframe):
    meanSeries = dframe['dens'].groupby(dframe['population']).mean()
    meanFrame = pd.DataFrame({'meanDens': meanSeries})
    meanFrame = meanFrame.assign(population=meanFrame.index.tolist())
    dframe = pd.merge(dframe, meanFrame, how='outer', on='population')

    return dframe

def popByDensity2(dframe, centers):
    meanSeries = dframe['dens'].groupby(dframe['population']).mean()
    meanFrame = pd.DataFrame({'meanDens':meanSeries})
    centerSeries = pd.Series(centers)
    centerSeries.index +=1
    meanFrame = meanFrame.assign(centers=centerSeries)
    meanFrame = meanFrame.assign(population=meanFrame.index.tolist())
    meanFrame = meanFrame.loc[meanFrame['population'] != 999, ['population', 'meanDens', 'centers']]
    meanFrame = meanFrame.sort_values('meanDens', ascending = False)
    meanFrame = meanFrame.reset_index(drop = True)
    meanFrame = meanFrame.assign(corrPop=meanFrame.index.tolist())
    meanFrame['corrPop'] += 1

    centers = meanFrame['centers'].tolist()

    dframe = pd.merge(dframe,meanFrame, how='outer', on='population')
    dframe['population'] = dframe['corrPop']
    dframe['population'] = dframe['population'].fillna(999).astype(int)
    dframe['population'] = dframe['population'].astype(int)
    dframe = dframe.drop(['corrPop'], axis=1)
    dframe = dframe.sort_values(['Id', 'population'])

    return dframe, centers


def calcP32(dframe, c13frame):
    c13frame.columns = ['C13', 'kappa']
    c13frame = c13frame.assign(population=c13frame.index.tolist())
    newframe = pd.merge(dframe, c13frame, on='population', how='outer')

    newframe['P32'] = newframe['dens']*newframe['C13']

    return newframe


#########################################################################
# Main script
# Calculate alpha
# Stereoplots and division
# Weighting by density
# Calculating the centers of n clusters
# Extracting cluster strike/dip into a DataFrame
# Convert population column into int
# Sort by outcrop id
#########################################################################

# Calculate alpha for each fracture population on every outcrop
data2 = calcAlpha(data2)
dataBeni = calcAlpha(dataBeni)

# Weighting the data by density. Approximate all outcrops as 10 m radius circular disks
data2 = weightByN(data2)
data2 = data2.sort_values('Id', ascending = True)

dataBeni = weightByN(dataBeni)
dataBeni = dataBeni.sort_values('Id', ascending = True)

# Calculate the centers of n clusters
strikes = data2['strike'].as_matrix()
dips = data2['dip'].as_matrix()
centers1 = analysis.kmeans(strikes, dips, num=5, bidirectional=True, measurement='poles', tolerance=1e-15)

cstrikes = np.zeros(1)
cdips = np.zeros(1)
counter = 0

for i in centers1:
    cstrike, cdip = stereonet_math.geographic2pole(i[0],i[1])
    cstrikes = np.append(cstrikes,cstrike)
    cdips = np.append(cdips,cdip)
    counter += 1

cstrikes = cstrikes[1:]
cdips = cdips[1:]

# Classify the data into populations by comparing the poles with a set of test poles
cTest = pd.DataFrame({'strike':cstrikes.tolist(), 'dip':cdips.tolist()})
data2 = splitByAngle(data2, cTest)

# Replace all NaN values with 999
data2[['population']] = data2[['population']].fillna(999).astype(int)
# Convert population column into int
data2[['population']] = data2[['population']].astype(int)
# Sort data by outcrop index
data2 = data2.sort_values('Id', ascending = True)
# Reclassify the populations from densest to sparsest (R1-Rn)
# Reclassify the cluster centers to match the new order
data2, centers1 = popByDensity2(data2, centers1)


############################################
# Plot some stereonets
############################################

# Group data by population
grouped = data2.groupby('population')
counter = 0

# Plot all clustered populations on stereonets
for idx, group in grouped:
    s = group['strike'].as_matrix()
    d = group['dip'].as_matrix()
    plot_stereonet(s, d, 'Area 1 clustered population ' + str(idx), 3)
    outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
    fname = outfp + '\\' + 'Area 1_' + 'kmeans_population ' + str(idx) + '_stereo_rose' + '.png'
    plt.savefig(fname, dpi=300)
    counter += 1

# Plot all expert approach populations on stereonets
grouped = dataBeni.groupby('population')
counter = 0

for idx, group in grouped:
    s = group['strike'].as_matrix()
    d = group['dip'].as_matrix()
    plot_stereonet(s, d, 'Area 1 expert population ' + str(idx), 3)
    outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
    fname = outfp + '\\' + 'Area 1_' + 'beni_population ' + str(idx) + '_stereo_rose' + '.png'
    plt.savefig(fname, dpi=300)
    counter += 1

# Plot clustered centers on stereonets
counter = 0

for i in centers1:
    cstrike, cdip = stereonet_math.geographic2pole(i[0],i[1])
    plot_stereonet(cstrike, cdip, 'Area 1 population ' + str(i), 3)
    counter += 1
    print(counter)
    outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
    fname = outfp + '\\' + 'Area 1_' + 'kmeans_' + str(counter) + '_stereo_rose' + '.png'
    print(fname)
    plt.savefig(fname, dpi=300)

# Plot all centers on a stereonet
plot_stereonet(cstrikes,cdips,'Clustered centers',3)
outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
fname = outfp + '\\' + 'Area 1_' + 'kmeans_all'  + '_stereo_rose' + '.png'
plt.savefig(fname, dpi=300)

# Plot all fracture poles on a stereonet
plot_stereonet(strikes,dips,'All clustered fracture poles',10)
outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
fname = outfp + '\\' + 'Area 1_' + 'stereo_rose' + '.png'
plt.savefig(fname, dpi=300)

# Plot all Expert approach poles on a stereonet
beniStrikes = dataBeni['strike'].as_matrix()
beniDips = dataBeni['dip'].as_matrix()
plot_stereonet(beniStrikes,beniDips,'All Bence fracture poles',10)
outfp = r'C:\Users\antti\Documents\Uni\GeoGradu\Aineisto\PythonScripts'
fname = outfp + '\\' + 'Area 1_' + 'Beni' + '.png'
plt.savefig(fname, dpi=300)

#########################
# C13
# P32
# Fisher stats
# Save to Excel
#########################
myC13, myfisherDF, myCDFAngles, myCDFRanges= calcPlotC13(data2, 'kmeans')

dataBeni = groupByDensity(dataBeni)
beniC13, bfisherDF, aaa, bbbb= calcPlotC13(dataBeni, 'Beni')

data2 = calcP32(data2, myC13)
dataBeni = calcP32(dataBeni, beniC13)

data2 = data2.drop_duplicates()
dataBeni = dataBeni.drop_duplicates()

writer = pd.ExcelWriter(outfp + '\\' +'C13.xlsx')
myC13.to_excel(writer, 'Kmeans_C13')
beniC13.to_excel(writer, 'HardSectoring_C13')
data2.to_excel(writer, 'HardSectoring_origdata')
dataBeni.to_excel(writer, 'HardSectoring_groupByDensity')
myfisherDF.to_excel(writer, 'Kmeans_Fisherstats')
bfisherDF.to_excel(writer, 'HardSectoring_Fisherstats')
myCDFAngles.to_excel(writer,'Kmeans_CDF_angles')
myCDFRanges.to_excel(writer,'Kmeans_CDF_range')
aaa.to_excel(writer,'HardSectoring_CDF_angles')
bbbb.to_excel(writer,'HardSectoring_CDF_range')
writer.save()

grouped = data2.groupby('population')
for idx, group in grouped:
    fname = outfp + '\\' + 'Cluster' + str(idx) + '.csv'
    group.to_csv(fname, index=False, encoding='utf-8')

grouped =dataBeni.groupby('population')
for idx, group in grouped:
    fname = outfp + '\\' + 'Beni' + str(idx) + '.csv'
    group.to_csv(fname, index=False, encoding='utf-8')