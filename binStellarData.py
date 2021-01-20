# R is resolving power/Resolution
# How finely are you sampling an optical element

# Sampling is unit specific: .001 microns, or .001 mm spacing
# Resolving power is constant in a fraction of a unit: lambda / delta lambda
# Really think of it as delta lambda / lambda, inverted for easier numbers (they get bigger)

# We want to keep R = lambda / delta lambda = 5000

# Start with the smallest wavelength value
# Start at 2 microns
# bin (average) the values from 1.9996 to 2.0004

# 

import configparser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from statistics import mean

def findEndIndices():
    pass

def readmodel(
    filename,
    cut_range=False,
    rangemin=None,  # in micron
    rangemax=None,  # in micron
    bin_data=False,
    ngridpoints=None,
    resolvingPower=5000,
    topValues=None,
    cwValues=None,
    start_time=None
):
    if (cut_range and rangemin is None) | (cut_range and rangemax is None):
        raise TypeError("If cut_range is True, cutmin and cutmax must not be None")
    if (
        (bin_data and rangemin is None)
        | (rangemin and rangemin is None)
        | (rangemin and ngridpoints is None)
    ):
        raise TypeError("If grid_data is True, cutmin and cutmax must not be None")

    model = pd.read_csv(
        filename, skiprows=7, delim_whitespace=True, names=["wavelength", "flux"]
    )

    # lets convert to micron
    model.wavelength /= 1.0e4
    
    # cut the wavelength values tot he specified cut range (.15-20.5 microns)
    if cut_range:
        model = model.loc[
            (model.wavelength >= rangemin) & (model.wavelength <= rangemax)
        ]

    ultimaLam=0.0
    lastLam=0.0
    lam = 0.0
    ultimaVal=0.0
    lastVal=0.0
    val = 0.0
    meanVal=0.0
    numValues = 0
    l=0
    bottom = 0.0
    outputY = []
    outputYstring = []
    cwValuesString = []

    lambdaIndex = 0
    while lambdaIndex < len(model.wavelength):
        lam = model.wavelength.values[lambdaIndex]
        val = model.flux.values[lambdaIndex]

        if lastLam == 0:
            lastLam = lam - 1e-6
            ultimaLam = lastLam - 1e-6
            lastVal = val
            ultimaVal = lastVal
            meanVal=0.0
            numberOfValues=0
            bottom = cwValues[0] - (topValues[0] - cwValues[0])
        if lam < bottom:
            numberOfValues = 0
            meanVal = 0.0

        # print("Lambda Index = ", lambdaIndex)
        # print("current l = ", l)
        # print("Current Top = ", topValues[l])
        # print("Bottom = ", bottom)
        # print("Current Lam = ", lam)
        # print("Last Lam = ", lastLam)
        # print("Ultimate Lam = ", ultimaLam)
        # print("Current Val = ", val)
        # print("Last Val = ", lastVal)
        # print("Ultimate Val = ", ultimaVal)
        # print("Mean Val = ", meanVal)
        # print("Number of Vals = ", numberOfValues)
        while l < len(cwValues) and topValues[l] < lam:
            if numberOfValues > 1:
                temp = meanVal/numberOfValues
                outputY.append(temp)
                outputYstring.append("{:.7e}".format(temp))
            elif topValues[l] > lastLam:
                temp = (val  - lastVal)*(cwValues[l] - lastLam)/( lam - lastLam) + lastVal
                outputY.append(temp)
                outputYstring.append("{:.7e}".format(temp))
            elif topValues[l] > ultimaLam:
                temp = (lastVal - ultimaVal) * (cwValues[l] - ultimaLam) / (lastLam - ultimaLam) + ultimaVal
                outputY.append(temp)
                outputYstring.append("{:.7e}".format(temp))
            else:
                outputY.append(ultimaVal)
                outputYstring.append("{:.7e}".format(ultimaVal))
            bottom = topValues[l]
            cwValuesString.append("{:.9e}".format(cwValues[l]))
            l += 1
            numberOfValues = 0
            meanVal = 0.0

        ultimaLam = lastLam
        lastLam = lam
        ultimaVal = lastVal
        lastVal = val
        meanVal += lastVal
        numberOfValues += 1
        lambdaIndex += 1

    if start_time != None:
        print("Time after binning")
        print("--- %s seconds ---" % (time.time() - start_time))
    binnedModel = pd.DataFrame(list(zip(cwValues, outputY)), columns =['wavelength', 'flux'])
    binnedModelStrings = pd.DataFrame(list(zip(cwValuesString, outputYstring)), columns =['wavelength', 'flux'])

    # This code plots the difference in wavelength values of the raw data. Shows a resolution across the model
    wavelengthDifference = []
    wavelengthResolution = []
    frequencyList = []
    frequencyDiff = []
    c = 2.99e14
    count = 0
    for value in binnedModel.wavelength:
        if count == len(binnedModel.wavelength) - 1 or count == 0:
            wavelengthDifference.append(float('nan'))
            wavelengthResolution.append(float('nan'))
            frequencyDiff.append(float('nan'))
            frequencyList.append(float('nan'))
            # print("End of file")
        else:
            frequency = c/value
            diffWave = binnedModel.wavelength[count + 1] - value
            diffFreq = (frequency) - (c/binnedModel.wavelength[count + 1])
            wavelengthDifference.append(diffWave)
            wavelengthResolution.append(value/diffWave)
            frequencyList.append(frequency)
            frequencyDiff.append(diffFreq)

            # print(diffWave)
            # print(value/diffWave)
            # print(diffFreq)
            # if diffWave > 0.05:
            #     print(model.wavelength[count])
            #     print("Greater than 0.05")
            #     print(diffWave)
        count += 1
    
    fig = plt.figure()
    plt.yscale('log')
    plt.plot(binnedModel.wavelength, wavelengthDifference)
    plt.title("Wavelength Spacing of NextGen Data")
    plt.ylabel("Delta Wavelength")
    plt.xlabel("Wavelength (um)")
    plt.savefig('./BinnedNextGenModels/VariabilityGraphs/wavelengthSpacingOfRawData.png', bbox_inches='tight')
    # plt.show()
    plt.close("all")
    print("WavelengthDiffPlot done")

    fig = plt.figure()
    # plt.yscale('log')
    plt.ticklabel_format(useOffset=False, style='plain')
    plt.ylim(4999, 5001)
    plt.plot(binnedModel.wavelength, wavelengthResolution)
    plt.title("Wavelength Resolution of NextGen Data")
    plt.ylabel("Wavelength/Delta Wavelength (Wavelength Resolving Power)")
    plt.xlabel("Wavelength (um)")
    plt.savefig('./BinnedNextGenModels/VariabilityGraphs/wavelengthResolvingPower.png', bbox_inches='tight')
    # plt.show()
    plt.close("all")
    print("WavelengthResolvingPowerPlot done")

    fig = plt.figure()
    plt.yscale('log')
    plt.plot(frequencyList, frequencyDiff)
    plt.title("Delta Frequency of NextGen Data")
    plt.ylabel("Delta Frequency")
    plt.xlabel("Frequency (um/s)")
    plt.savefig('./BinnedNextGenModels/VariabilityGraphs/frequencyDiff.png', bbox_inches='tight')
    # plt.show()
    plt.close("all")
    print("FrequencyDiffPlot done")

    return binnedModelStrings


if __name__ == "__main__":
    configParser = configparser.RawConfigParser()
    while True:
        try:
            fileName = input("Config File Path ./Config/")
            configParser.read_file(open("./Config/%s" % fileName))
            break
        except FileNotFoundError:
            print("There is no file by that name, please try again.")

    # The file names for the wavelength/flux values of the star, spots, and faculae
    phot_model_file = configParser.get('Star', 'phot_model_file')
    phot_model_file = phot_model_file.strip('"') # configParser adds extra "" that I remove
    spot_model_file = configParser.get('Star', 'spot_model_file')
    spot_model_file = spot_model_file.strip('"')
    fac_model_file = configParser.get('Star', 'fac_model_file')
    fac_model_file = fac_model_file.strip('"')

    resolvingPower = int(configParser.get('Star', 'resolvingPower'))
    binnedWavelengthMin = float(configParser.get('Star', 'binnedWavelengthMin'))
    binnedWavelengthMax = float(configParser.get('Star', 'binnedWavelengthMax'))

    topValues = []
    cwValues = []
    CW = binnedWavelengthMin
    while CW < binnedWavelengthMax:
        deltaLambda = CW / resolvingPower
        topValue = CW + (deltaLambda / 2)
        topValues.append(topValue)
        cwValues.append(CW)
        CW = CW + deltaLambda


    start_time = time.time()
    # Returns a dataframe of the binned wavelength/flux pairs as formatted strings
    binnedStarspectrum = readmodel(
        phot_model_file,
        cut_range=True,
        rangemin=.1999,  # in micron, based on MIRECLE constraints
        rangemax=20.5,  # in micron, based on  MIRECLE constraints
        bin_data=True,
        ngridpoints=3000,
        resolvingPower = resolvingPower,
        topValues = topValues,
        cwValues=cwValues,
        start_time=start_time
        )
    
    print("Time after binning and plotting")
    print("--- %s seconds ---" % (time.time() - start_time))

    binnedStarspectrumCSV = binnedStarspectrum.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
    print("Type of CSV = ", type(binnedStarspectrumCSV))
    file = open(r'./BinnedNextGenModels/binned3000StellarModel.txt','w')
    file.write(binnedStarspectrumCSV)
    file.close()

    print("\n")
    print("Time after binning and plotting and txt output")
    print("--- %s seconds ---" % (time.time() - start_time))

    binnedSpotspectrum = readmodel(
        spot_model_file,
        cut_range=True,
        rangemin=.1999,  # in micron, based on MIRECLE constraints
        rangemax=20.5,  # in micron, based on MIRECLE constraints
        bin_data=True,
        ngridpoints=3000,
        resolvingPower=5000,
        topValues = topValues,
        cwValues=cwValues,
        )

    binnedSpotspectrumCSV = binnedSpotspectrum.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
    print("Type of CSV = ", type(binnedStarspectrumCSV))
    file = open(r'./BinnedNextGenModels/binned2600StellarModel.txt','w')
    file.write(binnedSpotspectrumCSV)
    file.close()

    binnedFaculaespectrum = readmodel(
        fac_model_file,
        cut_range=True,
        rangemin=.1999, # in micron, based on MIRECLE constraints
        rangemax=20.5, # in micron, based on MIRECLE constraints
        bin_data=True,
        ngridpoints=3000,
        resolvingPower=5000,
        topValues = topValues,
        cwValues=cwValues,
        )

    binnedFaculaespectrumCSV = binnedFaculaespectrum.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
    print("Type of CSV = ", type(binnedStarspectrumCSV))
    file = open(r'./BinnedNextGenModels/binned3100StellarModel.txt','w')
    file.write(binnedFaculaespectrumCSV)
    file.close()