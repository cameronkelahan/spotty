import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
from scipy import stats
from scipy.signal import savgol_filter
from numba import jit, njit
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib as mpl
import imageio
import statistics
import configparser
import ast
from star_spectra import Spectramodel
from star_flatmap import readmodel
import math
from matplotlib import ticker, cm
import copy
from numpy import ma

# This program

class Starmodel(Spectramodel):
    def __init__(
        self,
        teffStar,
        starspectrum,
        spotspectrum,
        Rotstar,
        surface_map,
        flux_min,
        flux_max,
        dictionary,
    ):
        # Inherit parent class's methods and attributes
        super().__init__(
        teffStar,
        starspectrum,
        spotspectrum,
        Rotstar,
        surface_map,
        flux_min,
        flux_max,
        )

        self.teffStar = teffStar
        self._set_model_spectra(starspectrum, spotspectrum)
        self.Rotstar = Rotstar
        self.surface_map = surface_map
        self.dictionary = dictionary

    def buildEmptyDict(self):
        for value in self.modelspectra.wavelength:
            self.dictionary[value] = []
        print("Dict Built")


    def build_variability_dictionary(self):

        # Plot values will be the sumflux values for this phase
        plot_vals = self.modelspectra.sumflux

        # This adds the current phases flux value to a list corresponding to each wavelength value
        # Creates a 3000 length list of wavelengths, and each wavelength has a list of 360 length
        # of the fluxes from each phase
        index = 0
        for value in self.modelspectra.wavelength:
            self.dictionary[value].append(plot_vals[index])
            index += 1

    # x is 360 in length, one for each phase image
    # y is 3000 in length, one for each wavelength
    # z is the variability value for each x, y, pair (360, 3000)
    # Plot the color contour map
    def plot_color_variability(self, x, y, z):

        posIndexList = []
        negIndexList = []
        PPMList = []
        # absMin = 0

        # Plot a graph of the following equation/array
        waveCount = 0
        for wavelengthValueList in z:
            waveAvg = statistics.mean(wavelengthValueList)
            z[waveCount] = (wavelengthValueList / waveAvg) - 1
            
            # Record the location of each positive and negative value in the current phase array
            negValuesIndex = np.where(z[waveCount] < 0)
            posValuesIndex = np.where(z[waveCount] > 0)
            zeroValuesIndex = np.where(z[waveCount] == 0)

            negIndexList.append(negValuesIndex)
            posIndexList.append(posValuesIndex)

            # Change all 0 values to nan to avoid -inf after log10
            for value in zeroValuesIndex:
                z[waveCount][value] = float('nan')

            # Change all negative fractional changes to positive so log10 can be taken
            for value in negValuesIndex:
                z[waveCount][value] = z[waveCount][value] * -1
            
            tinyValueIndex = np.where(z[waveCount] < 2e-14)

            # Change all 0 values to nan to avoid -inf after log10
            for value in tinyValueIndex:
                z[waveCount][value] = float('nan')

            # convert to ppm
            z[waveCount] = z[waveCount] * 1e6
            PPMList.append(z[waveCount])

            # Take the log10 of the Fractional Change values
            z[waveCount] = np.log10(z[waveCount])

            # Change originally negative values back to negative
            for index in negValuesIndex:
                z[waveCount][index] = z[waveCount][index] * -1
                PPMList[waveCount][index] = PPMList[waveCount][index] * -1

            # minZ = min(z[waveCount])
            # print("min(z[phaseCount]) = ", minZ)
            # if minZ < absMin:
            #     absMin = minZ

            # if not math.isnan(minZ):
            #     minZ = math.floor(minZ)
            #     z[phaseCount] = z[phaseCount] - minZ
            #     for value in negValuesIndex:
            #         z[phaseCount][value] *= -1
            
            waveCount += 1

        # absMin = math.floor(absMin)

        # phaseCount = 0
        # for phase in z:
        #     # z[phaseCount] = z[phaseCount] - absMin
        #     for value in negIndexList[phaseCount]:
        #         z[phaseCount][value] *= -1
        #     phaseCount += 1

        X, Y = np.meshgrid(x, y)
        cmap = mpl.cm.Spectral
        # PPM fractional change color graph
        plt.contourf(X, Y, PPMList, levels=25, cmap=cmap)
        cbar = plt.colorbar()
        cbar.set_label('\N{GREEK CAPITAL LETTER DELTA}F (ppm)')
        plt.title('Fractional Change in Intensity for Each Wavelength')
        plt.xlabel('Phase (0-360)')
        plt.ylabel('Wavelength (0-20 um)')
        plt.savefig('./ProxCen/VariabilityGraphs/contourFracChange.png', dpi=200)
        plt.show()
        plt.close("all")

        # logarithmic PPM fractional change color graph
        plt.contourf(X, Y, z, levels=25, cmap=cmap) #  locator=locator

        # cmap = mpl.cm.cool
        # norm = mpl.colors.Normalize(vmin=-6, vmax=6)
        # cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap))

        cbar = plt.colorbar()
        # cbar = plt.colorbar(ticks=[-5, -4, -3, -2, -1, 0, 1, 2, 3, 4], drawedges=True) # , spacing='proportional'

        cbar.set_ticklabels(['-10\N{SUPERSCRIPT FOUR}\u22C5\N{SUPERSCRIPT EIGHT}',
                             '-10\N{SUPERSCRIPT THREE}\u22C5\N{SUPERSCRIPT SIX}',
                             '-10\N{SUPERSCRIPT TWO}\u22C5\N{SUPERSCRIPT FOUR}',
                             '-10\N{SUPERSCRIPT ONE}\u22C5\N{SUPERSCRIPT TWO}',
                             '0',
                             '10\N{SUPERSCRIPT ONE}\u22C5\N{SUPERSCRIPT TWO}',
                             '10\N{SUPERSCRIPT TWO}\u22C5\N{SUPERSCRIPT FOUR}',
                             '10\N{SUPERSCRIPT THREE}\u22C5\N{SUPERSCRIPT SIX}',
                             '10\N{SUPERSCRIPT FOUR}\u22C5\N{SUPERSCRIPT EIGHT}'])  # vertically oriented colorbar
        cbar.set_label('\N{GREEK CAPITAL LETTER DELTA}F (ppm)')
        
        plt.xlabel('Phase (0-360)')
        plt.ylabel('Wavelength (0-20 um)')

        plt.title('Logarithmic Fractional Change in Intensity for Each Wavelength')
        plt.savefig('./ProxCen/VariabilityGraphs/contourLogFracChange.png', dpi=200)
        plt.show()
        print("Done")

if __name__ == "__main__":
    
    configParser = configparser.RawConfigParser()
    while True:
        try:
            fileName = input("Config File Path ./Config/")
            configParser.read_file(open("./Config/%s" % fileName))
            break
        except FileNotFoundError:
            print("There is no file by that name, please try again.")
    
    # # Temperature values for the star, spots, and faculae
    teffStar = int(configParser.get('ProxCen', 'teffStar'))
    # teffSpot = int(configParser.get('ProxCen', 'teffSpot'))
    # teffFac = int(configParser.get('ProxCen', 'teffFac'))

    # The file names for the wavelength/flux values of the star, spots, and faculae
    phot_model_file = configParser.get('ProxCen', 'phot_model_file')
    phot_model_file = phot_model_file.strip('"') # configParser adds extra "" that I remove
    spot_model_file = configParser.get('ProxCen', 'spot_model_file')
    spot_model_file = spot_model_file.strip('"')
    fac_model_file = configParser.get('ProxCen', 'fac_model_file')
    fac_model_file = fac_model_file.strip('"')

    Rotstar = float(configParser.get('ProxCen', 'Rotstar'))

    # Name of the .npy array which contains the flat surface map
    surface_map_file = configParser.get('HemiMap', 'surface_map_file')
    surface_map_file = surface_map_file.strip('"')
    surface_map = np.load(surface_map_file)

    # Total number of exposures to be taken
    num_exposures = int(configParser.get('HemiMap', 'num_exposures'))
    # Time *in days) between exposures
    time_between_exposures = float(configParser.get('HemiMap', 'time_between_exposures'))

    # Number of wavelength bins to be graphed
    num_bins = int(configParser.get('Spectra', 'num_bins'))
    # List of the Wavelength bins' lower and upper bounds
    list_bins = ast.literal_eval(configParser.get('Spectra', 'list_bins'))
    # Tells the program whether or not to normalize the binned flux values
    normalized_bins = configParser.get('Spectra', 'normalized_bins')
    if normalized_bins == "True" or normalized_bins == "true" or normalized_bins == "TRUE":
        normalized_bins = True
        print("Normalized = ", normalized_bins)
    else: normalized_bins = False

    # The chosen wavelength to have its variability examined
    focused_wavelength = float(configParser.get('Variability', 'focused_wavelength'))

    # hasPlanet = configParser.get('Planet', 'hasPlanet')
    # Rplanet = configParser.get('Planet', 'Rplanet')
    # Pplanet = configParser.get('Planet', 'Pplanet')
    # Impactplanet = configParser.get('Planet', 'Impactplanet')
    # TransitDuration = configParser.get('Planet', 'TransitDuration')
    
    # Calculates what percent the time between exposures is compared to the stellar rotation period
    # This is used to calculate the change in phase between images
    exposure_time_percent = time_between_exposures * (100 / Rotstar)
    deltaPhase = exposure_time_percent / 100

    
    starspectrum = readmodel(
        phot_model_file,
        cut_range=True,
        rangemin=0,  # in micron, based on MIRECLE constraints
        rangemax=20.000,  # in micron, based on  MIRECLE constraints
        grid_data=True,
        ngridpoints=3000,
    )

    spotspectrum = readmodel(
        spot_model_file,
        cut_range=True,
        rangemin=0,  # in micron, based on MIRECLE constraints
        rangemax=20.000,  # in micron, based on MIRECLE constraints
        grid_data=True,
        ngridpoints=3000,
    )

    faculaespectrum = readmodel(
        fac_model_file,
        cut_range=True,
        rangemin=0, # in micron, based on MIRECLE constraints
        rangemax=20.000, # in micron, based on MIRECLE constraints
        grid_data=True,
        ngridpoints=3000,
    )

    # Create a HemiMap class
    SM = Starmodel(
        teffStar,
        starspectrum,
        spotspectrum,
        Rotstar,
        surface_map,
        flux_min = min(spotspectrum.flux),
        flux_max = max(starspectrum.flux),
        dictionary = {},
    )

    SM.buildEmptyDict()

    phase = 0
    count = 0
    x = []
    y = []
    
    photflux_list = []
    smoothed_photflux_list = []
    sumflux_list = []
    smoothed_sumflux_list = []
    first_bool = True
    wavelength_position = 0

    color_plot_values = []

    # This for loop is for each hemisphere image
    for value in range(num_exposures):
        hemi_map = np.load('./ProxCen/HemiMapImages+Arrays/hemiArray_%d' % count, allow_pickle=True)

        # Calculates the star's combined spectrum for the given hemisphere and adds a column to the SM's
        # model info called sumflux
        SM.calculate_combined_spectrum(hemi_map)

        SM.build_variability_dictionary()

        phase += deltaPhase
        print("Phase = ", phase)
        if (phase > 1):
            phase = phase % 1
        count += 1

        del SM.modelspectra['sumflux']

    # Creates a list that is 3000 values in length, each entry containing a list of 360 length
    for wavelength in SM.modelspectra.wavelength:
        color_plot_values.append(SM.dictionary[wavelength])

    y = SM.modelspectra.wavelength
    x = list(range(360))

    SM.plot_color_variability(x, y, color_plot_values)
    color_plot_values = np.asarray(color_plot_values)
    # print(color_plot_values.size)

    index = 0
    for wavelength in SM.modelspectra.wavelength:
        pass
        

    np.save('./ProxCen/VariabilityGraphs/contourPlotValues.npy', color_plot_values)
    # temp = np.load('./ProxCen/VariabilityGraphs/contourPlotValues.npy')
    # print("Temp = ", temp)
    # print(type(temp))
    # print(temp.shape)
    print("Fin")
