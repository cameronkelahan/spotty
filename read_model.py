import matplotlib.pyplot as plt
import pandas as pd
import time

def readmodel(
    filename,
    starName,
    cut_range=False,
    rangemin=None,  # in micron
    rangemax=None,  # in micron
    bin_data=False,
    ngridpoints=None,
    resolvingPower=5000,
    topValues=None,
    cwValues=None,
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
    
    # cut the wavelength values tot he specified cut range (.1999-20.5 microns)
    if cut_range:
        model = model.loc[
            (model.wavelength >= rangemin) & (model.wavelength <= rangemax)
        ]

    if bin_data:
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

        start_time = time.time()

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

        print("Time after binning")
        print("--- %s seconds ---" % (time.time() - start_time))
        model = pd.DataFrame(list(zip(cwValues, outputY)), columns =['wavelength', 'flux'])
        modelStrings = pd.DataFrame(list(zip(cwValuesString, outputYstring)), columns =['wavelength', 'flux'])

    # This code plots the difference in wavelength values of the raw data. Shows a resolution across the model
    wavelengthDifference = []
    wavelengthResolution = []
    frequencyList = []
    frequencyDiff = []
    c = 2.99e14
    count = 0
    for value in model.wavelength:
        if count == len(model.wavelength) - 1 or count == 0:
            wavelengthDifference.append(float('nan'))
            wavelengthResolution.append(float('nan'))
            frequencyDiff.append(float('nan'))
            frequencyList.append(float('nan'))
            # print("End of file")
        else:
            frequency = c/value
            diffWave = model.wavelength[count + 1] - value
            diffFreq = (frequency) - (c/model.wavelength[count + 1])
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
    plt.plot(model.wavelength, wavelengthDifference)
    plt.title("Wavelength Spacing of NextGen Data")
    plt.ylabel("Delta Wavelength")
    plt.xlabel("Wavelength (um)")
    if bin_data:
        plt.savefig('./BinnedNextGenModels/VariabilityGraphs/wavelengthSpacingOfRawData.png', bbox_inches='tight')
    else:
        plt.savefig('./%s/VariabilityGraphs/wavelengthSpacingOfRawData.png' % starName, bbox_inches='tight')
    # plt.show()
    plt.close("all")
    print("WavelengthDiffPlot done")

    fig = plt.figure()
    # plt.yscale('log')
    plt.ticklabel_format(useOffset=False, style='plain')
    plt.ylim(4999, 5001)
    plt.plot(model.wavelength, wavelengthResolution)
    plt.title("Wavelength Resolution of NextGen Data")
    plt.ylabel("Wavelength/Delta Wavelength (Wavelength Resolving Power)")
    plt.xlabel("Wavelength (um)")
    if bin_data:
        plt.savefig('./BinnedNextGenModels/VariabilityGraphs/wavelengthResolvingPower.png', bbox_inches='tight')
    else:
        plt.savefig('./%s/VariabilityGraphs/wavelengthResolvingPower.png' % starName, bbox_inches='tight')
    # plt.show()
    plt.close("all")
    print("WavelengthResolvingPowerPlot done")

    fig = plt.figure()
    plt.yscale('log')
    plt.plot(frequencyList, frequencyDiff)
    plt.title("Delta Frequency of NextGen Data")
    plt.ylabel("Delta Frequency")
    plt.xlabel("Frequency (um/s)")
    if bin_data:
        plt.savefig('./BinnedNextGenModels/VariabilityGraphs/frequencyDiff.png', bbox_inches='tight')
    else:
        plt.savefig('./%s/VariabilityGraphs/frequencyDiff.png' % starName, bbox_inches='tight')
    # plt.show()
    plt.close("all")
    print("FrequencyDiffPlot done")

    return model, modelStrings