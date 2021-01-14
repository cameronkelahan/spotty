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
from statistics import mean

def findEndIndices():
    pass

def readmodel(
    filename,
    cut_range=False,
    rangemin=None,  # in micron
    rangemax=None,  # in micron
    grid_data=False,
    ngridpoints=None,
    resolvingpower=5000,
):
    if (cut_range and rangemin is None) | (cut_range and rangemax is None):
        raise TypeError("If cut_range is True, cutmin and cutmax must not be None")
    if (
        (grid_data and rangemin is None)
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

    # Start with the smallest wavelength value
    # To achieve a resolving power of 5000, divide the current lambda by 5000 to get the delta lambda value
    # Once the flux value of that lambda is calculated, add the delta lambda to the lambda to get new center wavelength
    
    # Start at 2 microns
    # 2 is the center, with .0004 as the delta lambda
    # 1) bin (average) the values from 1.9996 to 2.0004

    # Next central wavelength is current center (2) + current delta lambda (.0004)
    # 2) 2.0004 is the next center, with .00040008 as the next delta lambda
        # - lower = 1.99999992, uppper = 2.00080008

    # position = 0
    # for wavelength in model.wavelength:
    #     if wavelength >= 1.999:
    #         print("Wavelength = ", wavelength, "; position = ", position)
    #     position += 1

    # Create an empty dictionary to store the binned wavelength/flux pairs
    binnedValuesDict = {}
    # Create an empty dictionary to store the binned wavelength/flux pairs as formatted strings
    binnedStringValuesDict={}

    centerWavelength = .2
    wavelengthCounter = 0
    prevFirst = 0
    while wavelengthCounter < len(model.wavelength):
        wavelength = model.wavelength.values[wavelengthCounter]

        # finds the first center wavelength
        if wavelength >= centerWavelength:
            # Calculate the deltalambda of the current center wavelength
            deltaLambda = centerWavelength / resolvingpower

            # find the high and low value of the bin for the center wavelength
            lowValue = centerWavelength - deltaLambda
            highValue = centerWavelength + deltaLambda
            
            fluxValuesInRange = []
            
            # run through each wavelength looking for the values in this bin's range
            # These variables keep track of the smallest difference between the high or low limit of the bin
            # and the wavelength closest to them in the dataset
            # Theya re used when there are no wavelengths from the dataset that fall in the bin range
            smallestDiff = 9999999999999999999
            closestWavelengthIndex = -1
            waveIndex = prevFirst
            first = True

            while True:
                value = model.wavelength.values[waveIndex]
                currentDiff = min(abs(lowValue - value), abs(highValue - value))
                if currentDiff < smallestDiff:
                    smallestDiff = currentDiff
                    closestWavelengthIndex = waveIndex
                if value >= lowValue and value <= highValue:
                    # Add the appropriate flux values to the list
                    fluxValuesInRange.append(model.flux.values[waveIndex])
                    if first:
                        prevFirst = waveIndex
                        first = False
                elif value > highValue:
                    # Check to see if the list of flux values in range is empty
                    if not fluxValuesInRange:
                        fluxValuesInRange.append(model.flux.values[closestWavelengthIndex])
                    break
                waveIndex += 1
            
            # calcualte the bin's average
            binAverageFlux = mean(fluxValuesInRange)
            binAverageFluxString = "{:.7e}".format(binAverageFlux)
            
            # add the wavelength/flux pair to the dictionary
            binnedValuesDict[centerWavelength] = binAverageFlux
            centerWavelengthString = "{:.9e}".format(centerWavelength)
            binnedStringValuesDict[centerWavelengthString] = binAverageFluxString

            # calculate the new center wavelength by adding it to the current delta lambda
            centerWavelength = centerWavelength + deltaLambda

        # print(wavelength)
        if centerWavelength > 20:
            break

        wavelengthCounter += 1

    # Create a list of wavelength, flux pairs of the binnned data
    # binnedData = []
    # for key in binnedValuesDict:
    #     listPair = [key, binnedValuesDict[key]]
    #     binnedData.append(listPair)

    print(type(binnedValuesDict.items()))
    binnedModel = pd.DataFrame(list(binnedValuesDict.items()), columns=['wavelength', 'flux'])
    binnedModelStrings = pd.DataFrame(list(binnedStringValuesDict.items()), columns=['wavelength', 'flux'])
    print(binnedModel)
    print(binnedModel.columns)
    print("done")

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
            print("End of file")
        else:
            frequency = c/value
            diffWave = binnedModel.wavelength[count + 1] - value
            diffFreq = (frequency) - (c/binnedModel.wavelength[count + 1])
            wavelengthDifference.append(diffWave)
            wavelengthResolution.append(value/diffWave)
            frequencyList.append(frequency)
            frequencyDiff.append(diffFreq)

            print(diffWave)
            print(value/diffWave)
            print(diffFreq)
            if diffWave > 0.05:
                print(model.wavelength[count])
                print("Greater than 0.05")
                print(diffWave)
        count += 1
    
    fig = plt.figure()
    plt.yscale('log')
    plt.plot(binnedModel.wavelength, wavelengthDifference)
    plt.title("Wavelength Spacing of NextGen Data")
    plt.ylabel("Delta Wavelength")
    plt.xlabel("Wavelength (um)")
    plt.savefig('./spotty/BinnedNextGenModels/VariabilityGraphs/wavelengthSpacingOfRawData.png', bbox_inches='tight')
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
    plt.savefig('./spotty/BinnedNextGenModels/VariabilityGraphs/wavelengthResolvingPower.png', bbox_inches='tight')
    # plt.show()
    plt.close("all")
    print("WavelengthResolvingPowerPlot done")

    fig = plt.figure()
    plt.yscale('log')
    plt.plot(frequencyList, frequencyDiff)
    plt.title("Delta Frequency of NextGen Data")
    plt.ylabel("Delta Frequency")
    plt.xlabel("Frequency (um/s)")
    plt.savefig('./spotty/BinnedNextGenModels/VariabilityGraphs/frequencyDiff.png', bbox_inches='tight')
    # plt.show()
    plt.close("all")
    print("FrequencyDiffPlot done")

    return binnedModelStrings


if __name__ == "__main__":
    configParser = configparser.RawConfigParser()
    while True:
        try:
            fileName = input("Config File Path ./Config/")
            configParser.read_file(open("./spotty/Config/%s" % fileName))
            break
        except FileNotFoundError:
            print("There is no file by that name, please try again.")

    # The file names for the wavelength/flux values of the star, spots, and faculae
    phot_model_file = configParser.get('ProxCen', 'phot_model_file')
    phot_model_file = phot_model_file.strip('"') # configParser adds extra "" that I remove
    spot_model_file = configParser.get('ProxCen', 'spot_model_file')
    spot_model_file = spot_model_file.strip('"')
    fac_model_file = configParser.get('ProxCen', 'fac_model_file')
    fac_model_file = fac_model_file.strip('"')

    # Returns a dataframe of the binned wavelength/flux pairs as formatted strings
    binnedStarspectrum = readmodel(
        phot_model_file,
        cut_range=True,
        rangemin=.1999,  # in micron, based on MIRECLE constraints
        rangemax=20.5,  # in micron, based on  MIRECLE constraints
        bin_data=True,
        ngridpoints=3000,
        resolvingpower=5000,
        )
    
    binnedStarspectrumCSV = binnedStarspectrum.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
    print("Type of CSV = ", type(binnedStarspectrumCSV))
    file = open(r'./spotty/BinnedNextGenModels/binned3000StellarModel.txt','w')
    file.write(binnedStarspectrumCSV)
    file.close()

    binnedSpotspectrum = readmodel(
        spot_model_file,
        cut_range=True,
        rangemin=.1999,  # in micron, based on MIRECLE constraints
        rangemax=20.5,  # in micron, based on MIRECLE constraints
        bin_data=True,
        ngridpoints=3000,
        resolvingpower=5000,
        )

    binnedSpotspectrumCSV = binnedSpotspectrum.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
    print("Type of CSV = ", type(binnedStarspectrumCSV))
    file = open(r'./spotty/BinnedNextGenModels/binned2600StellarModel.txt','w')
    file.write(binnedSpotspectrumCSV)
    file.close()

    binnedFaculaespectrum = readmodel(
        fac_model_file,
        cut_range=True,
        rangemin=.1999, # in micron, based on MIRECLE constraints
        rangemax=20.5, # in micron, based on MIRECLE constraints
        bin_data=True,
        ngridpoints=3000,
        resolvingpower=5000,
        )

    binnedFaculaespectrumCSV = binnedFaculaespectrum.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
    print("Type of CSV = ", type(binnedStarspectrumCSV))
    file = open(r'./spotty/BinnedNextGenModels/binned3100StellarModel.txt','w')
    file.write(binnedFaculaespectrumCSV)
    file.close()