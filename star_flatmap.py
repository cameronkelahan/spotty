import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
from numba import jit, njit
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from statistics import mean
import matplotlib
import configparser
import os

# Temporary:
import sys
np.set_printoptions(threshold=sys.maxsize)

# This program generates the spots on the star and plots the star's surface as a rectangle
# It takes into account the spherical nature and distorts the spots as appropriate on the 2d projection
# It saves the star's surface as an array to be loaded in later to create the round, hemisphere graphs

interactive_plot = True
if not interactive_plot:
    matplotlib.use("Agg")
    plt.interactive(False)
    matplotlib.rcParams["figure.dpi"] = 300

# constants
k = 0.0172020989  # Gaussian Grav Constant
grav = k ** 2  # work in units of AU, DAYS and solar masses
rsun_to_AU = 0.00465047
gm = 2.959122083e-4
rEarth = 6370.0  # km
rEarth_to_rSun = 0.0091577

def readmodel(
    filename,
    cut_range=False,
    rangemin=None,  # in micron
    rangemax=None,  # in micron
    bin_data=False,
    resolvingPower=None,
):
    if (cut_range and rangemin is None) | (cut_range and rangemax is None):
        raise TypeError("If cut_range is True, cutmin and cutmax must not be None")
    if (
        (bin_data and rangemin is None)
        | (rangemin and rangemin is None)
        | (rangemin and resolvingPower is None)
    ):
        raise TypeError("If bin_data is True, cutmin and cutmax must not be None")

    model = pd.read_csv(
        # # For the binned NextGenModels
        # filename, skiprows=1, delim_whitespace=True, names=["wavelength", "flux"]
        # For the original btNextGenStellarModels
        filename, skiprows=7, delim_whitespace=True, names=["wavelength", "flux"]
    )

    # lets convert to micron
    model.wavelength /= 1.0e4
    if cut_range:
        model = model.loc[
            (model.wavelength >= rangemin) & (model.wavelength <= rangemax)
        ]

    # Bin the nextGen data to the specified resolving power
    if bin_data:
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
                deltaLambda = centerWavelength / resolvingPower

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
    plt.ylim(resolvingPower - 1, resolvingPower + 1)
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

    return binnedModel, binnedModelStrings

class Spotmodel:
    def __init__(
        self,
        spotCoverage,
        spotNumber,
        facCoverage,
        facNumber,
        starName,
    ):
        self.spotCoverage = spotCoverage
        self.spotNumber = spotNumber
        self.facCoverage = facCoverage
        self.facNumber = facNumber
        self.starName = starName

    def generate_spots(self, randomSeed=None):
        if randomSeed is not None:
            np.random.seed(randomSeed)
        # Why is radius here chosen to be 1?
        surface_area = 4.0 * np.pi * 1.0 ** 2
        # Picks spotnumber of radius'. all between 0 and 1
        spot_radius = np.random.random_sample(self.spotNumber)
        fac_radius = np.random.random_sample(self.facNumber)

        total_spot_coverage = np.sum(np.pi * spot_radius ** 2)
        total_fac_coverage = np.sum(np.pi * fac_radius ** 2)

        spotNormalization = surface_area * self.spotCoverage / total_spot_coverage
        facNormalization = surface_area * self.facCoverage / total_fac_coverage

        true_spot_radius = spot_radius * spotNormalization ** 0.5
        true_fac_radius = fac_radius * facNormalization ** 0.5

        # Limits latitude to between 60 and -60? Based on Butterfly effect?
        spotLat = -60 + 120 * np.random.random_sample(self.spotNumber)
        facLat = -60 + 120 * np.random.random_sample(self.facNumber)
        # Limits Longitude to between -180 and 180
        spotLon = -180 + 360 * np.random.random_sample(self.spotNumber)
        facLon = -180 + 360 * np.random.random_sample(self.facNumber)

        surface_map = self.generate_flat_surface_map(true_spot_radius, spotLon, spotLat,
                                                     true_fac_radius, facLon, facLat,)
        # print("Type of surface_map = ", type(surface_map))

        flat_image = surface_map.flatten()
        length = len(flat_image)
        print("Length = ", length)
        spot_pixels = np.where(flat_image == 1)
        print(len(spot_pixels[0]))
        fac_pixels = np.where(flat_image == 2)
        print(len(fac_pixels[0]))
        # summ = sum(flat_image)
        surfaceMapSpotCoverage = len(spot_pixels[0]) / length
        surfaceMapFacCoverage = len(fac_pixels[0]) / length

        print("SpotCoverage = ", surfaceMapSpotCoverage)
        print("FacCoverage = ", surfaceMapFacCoverage)

        if os.path.isfile('./%s/surfaceMapInfo.txt' % self.starName):
            os.remove('./%s/surfaceMapInfo.txt' % self.starName)
        f = open('./%s/surfaceMapInfo.txt' % self.starName, 'a')
        f.write('Spot Coverage Percentage = {:.2f}%\n'.format(surfaceMapSpotCoverage * 100))
        f.write('Fac Coverage Percentage = {:.2f}%\n'.format(surfaceMapFacCoverage * 100))
        f.close()

        plt.close("all")
        return surface_map

    # @njit
    def generate_flat_surface_map(self, spot_radii, spotLon, spotLat, fac_radii, facLon, facLat):
        # we create an image using matplotlib (!!)
        fig = plt.figure(figsize=[5.00, 2.5], dpi=1200)
        proj = ccrs.PlateCarree()
        ax = plt.axes(projection=proj, fc="r")
        canvas = FigureCanvas(fig)
        plt.gca().set_position([0, 0, 1, 1])

        ax.set_global()
        ax.outline_patch.set_linewidth(0.0)
        ax.set_extent([-180, 180, -90, 90])

        # loop through each spot, adding it to the image
        # tissot assume the sphere is earth, so multiply by radius of earth
        for spot in range(self.spotNumber):
            add_spots = ax.tissot(
                rad_km=spot_radii[spot] * rEarth,
                lons=spotLon[spot],
                lats=spotLat[spot],
                n_samples=1000,
                fc="k",
                alpha=1,
            )
        
        canvas.draw()
        buf = canvas.buffer_rgba()
        surface_map_image = np.asarray(buf)
        # print(surface_map_image[1000][1000])
        # 0 = photosphere
        # 1 = spot
        # 2 = fac
        surface_map = np.where(surface_map_image[:, :, 0] == 255, 0, 1)

        # fac_map = np.where(surface_map_image[:, :, 2] == 255)
        # surface_map[fac_map] = 2
        # print(surface_map)

        plt.savefig('./%s/TempFlatMap.png' % self.starName)

        for fac in range(self.facNumber):
            add_facs = ax.tissot(
                rad_km=fac_radii[fac] * rEarth,
                lons=facLon[fac],
                lats=facLat[fac],
                n_samples=1000,
                fc="b",
                alpha=1,
            )
        
        canvas.draw()
        buf = canvas.buffer_rgba()
        surface_map_image = np.asarray(buf)
        # 0 = photosphere
        # 1 = spot
        # 2 = fac
        fac_map = np.where(surface_map_image[:, :, 2] == 255, 2, 0)
        for row in range(len(fac_map)):
            for col in range(len(fac_map[row])):
                if fac_map[row][col] == 2:
                    surface_map[row][col] = 2
        # print(surface_map)

        
        # Save and show the surface map values (1 or 0) that create the red/black map image (rectangular shape)
        # NOTE: the plt.show() function does not scale well with this plot, must view from file
        plt.savefig('./%s/FlatMap.png' % self.starName)
        # plt.show()
        return surface_map