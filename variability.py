import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
from scipy import stats
from scipy.signal import savgol_filter
from numba import jit, njit
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib
import imageio
import statistics
import configparser
import ast
from star_spectra import Spectramodel
from star_flatmap import readmodel
import math

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

    # Smooth all of the flux values in photflux and sumflux through the use of a moving average smoothing technique.
    # For every wavelength, look at the flux values of the 2 neighboring wavelengths in each direction.
    # Recomupute the current flux value to be the average of those 5 values.
    def reduceDataSmoothing(self):
        # Lists of smoothed flux values, both potoflux and sumflux
        smoothedPhotFlux = []
        smoothedSumFlux = []

        wavelengthIndex = 0
        for wavelength in self.modelspectra.wavelength:
            
            # Look to the flux values for the wavelenghts 2 in front and 2 behind the current wavelength
            startIndex = wavelengthIndex - 2
            stopIndex = wavelengthIndex + 2

            if startIndex < 0:
                startIndex = 0
            if stopIndex > len(self.modelspectra.wavelength) - 1:
                stopIndex = len(self.modelspectra.wavelength) - 1

            sumValuePhotFlux = 0
            sumValueSumFlux = 0
            startStopDiff = stopIndex - startIndex

            while startIndex <= stopIndex:
                sumValuePhotFlux += self.modelspectra.photflux[startIndex]
                sumValueSumFlux += self.modelspectra.sumflux[startIndex]
                startIndex += 1
            sumValuePhotFlux /= startStopDiff + 1
            sumValueSumFlux /= startStopDiff + 1

            smoothedPhotFlux.append(sumValuePhotFlux)
            smoothedSumFlux.append(sumValueSumFlux)

            wavelengthIndex += 1
        
        # Insert 2 new columns to the model's dataframe, one for smoothed photoplux and one for smoothed sumflux
        self.modelspectra.insert(4, 'smoothPhotFlux', smoothedPhotFlux)
        self.modelspectra.insert(5, 'smoothSumFlux', smoothedSumFlux)

    def plot_variability(self, count, color_plot_values, phase0, phase180):
        # fig = plt.figure()

        # index = 0
        # for value in self.modelspectra.photflux:
        #     print(self.modelspectra.photflux[index])
        #     print(self.modelspectra.sumflux[index])
        #     if self.modelspectra.photflux[index] == 0.0:
        #         self.modelspectra.photflux[index] = 1e-54
        #     if self.modelspectra.sumflux[index] == 0.0:
        #         self.modelspectra.sumflux[index] = 1e-54
        #     index += 1

        # Geronimo's first equation.
        B = (self.modelspectra.sumflux/self.modelspectra.photflux) - 1.0
        zeroIndex = np.where(B == 0)
        tinyValueIndex = np.where((B*-1) < 1e-10)
        for index in zeroIndex:
            B[index] = float('nan')
        for index in tinyValueIndex:
            B[index] = float('nan')
        x = self.modelspectra.wavelength
        x = np.asarray(x)
        x = list(x.flat)

        plt.plot(x, B)
        plt.title('Variability Between Reference and SumFlux')
        plt.ylabel('(SumFlux[L] / ReferenceFlux[L]) - 1')
        plt.xlabel('Wavelength (microns)')
        plt.xlim(-1, 21)

        plt.ylim(-.2, 0.01)
        plt.savefig("./ProxCen/SmoothedVariabilityGraphs/variableFlux_%d" % count, bbox_inches='tight', dpi=200)
        # plt.show()

        plt.close("all")

        #########################################
        # Plot 2: plotting the difference between the regular variability and the smoothed variability
        fig = plt.figure()

        # Adds a column of smoothed star flux (photflux) and a column of smoothed sumflux to modeldata DataFrame
        self.reduceDataSmoothing()

        # Geronimo's second equation
        # Variance between the smoothed versions of photflux and sumflux
        A = (self.modelspectra.smoothSumFlux/self.modelspectra.smoothPhotFlux) - 1.0

        # Difference between the regular variance between the photflux/sumflux (B) and the smoothed variance between
        # the smooth photflux/sumflux (A)
        # TLDR: Difference in variance of regular and smooth flux values
        plot_vals = B - A

        if(count==0):
            phase0 = plot_vals
        if(count==180):
            phase180=plot_vals

        # This adds a list of each wavelengths variability to the color plot list for the current
        # phase
        # Should have a shape of (360, 3000) after all phases run through
        color_plot_values.append(plot_vals)
        
        # plt.plot(self.modelspectra.wavelength, self.modelspectra.smoothPhotFlux)
        # plt.plot(self.modelspectra.wavelength, self.modelspectra.photflux, alpha=0.5)

        plt.plot(x, plot_vals)
        plt.title('B - A')
        plt.savefig('./ProxCen/SmoothedVariabilityGraphs/normVsSmoothVariability_%d.png' % count, bbox_inches='tight', dpi=200)
        # plt.show()
        plt.close("all")

        return phase0, phase180

        ####### OLD CODE THAT PLOTS BINNED VALUES (DATA REDUCTION) ########

        # fig = plt.figure()

        # # plt.title("High Resolution Variability")
        # plt.xlabel('Wavelength (microns)')
        # # plt.ylabel('')

        # # Convolve the x and y values of the plot down to 21 points, 0-20 microns inclusive with a dL of 1 micron
        # xlow = np.linspace(0, 20, 21)
        # print("xlow = ", xlow)
        # ylowSpec, ylowRef = self.reduceDataBinning(xlow)
        # # ylowRef = savgol_filter(self.modelspectra.sumflux, 21, 4)
        # print("ylowRef = ", ylowRef)
        # print("ylowSpec = ", ylowSpec)
        # # ylowSpec

        # plt.plot(xlow, ylowRef)
        # plt.savefig("./ProxCen/VariabilityGraphs/lowResRef.png")
        # plt.close("all")

        # fig = plt.figure()
        # plt.plot(xlow, ylowSpec)
        # plt.savefig("./ProxCen/VariabilityGraphs/lowResSpec.png")
        # plt.close("all")

        # fig = plt.figure()

        # A = []
        # count = 0
        # for value in ylowSpec:
        #     temp = (ylowSpec[count]/ylowRef[count]) - 1
        #     A.append(temp)
        # plt.plot(xlow, A)
        # plt.savefig("./ProxCen/VariabilityGraphs/highRes_%d" % count)
        # plt.show()

        # plt.close("all")

    # Plot the smoothed data over the original data to show difference
    def plot_smoothing_data_reduction(self, count):
        fig = plt.figure()

        plt.plot(self.modelspectra.wavelength, self.modelspectra.smoothPhotFlux, label='SmoothedPhotFlux')
        plt.plot(self.modelspectra.wavelength, self.modelspectra.photflux, alpha=0.5, label='PhotFlux')

        plt.title('Variability Between Normal Reference (photflux) and Smoothed Reference')
        plt.savefig('./ProxCen/SmoothedVariabilityGraphs/smoothedPhotFluxSpectrum_%d.png' % count)
        # plt.show()
        print("Done")

        plt.close("all")

    # Not yet implemented
    # Allows for focusing on one wavelength value, to see how it changes
    def find_focused_wavelength(self, photflux_list, smoothed_photflux_list, sumflux_list,
                                smoothed_sumflux_list, first_bool, focused_wavelength,
                                wavelength_position, count):
        if first_bool:
            wavelength_position = 0
            while first_bool:
                # if wavelength_position == 749 or wavelength_position == 750 or wavelength_position == 751 or wavelength_position == 752:
                #     print(value)
                if self.modelspectra.wavelength[wavelength_position + 1] >= focused_wavelength:
                    less_than_diff = focused_wavelength - self.modelspectra.wavelength[wavelength_position]
                    more_than_diff = focused_wavelength - self.modelspectra.wavelength[wavelength_position + 1]
                    if abs(less_than_diff) < abs(more_than_diff):
                        focused_wavelength = self.modelspectra.wavelength[wavelength_position]
                    else:
                        focused_wavelength = self.modelspectra.wavelength[wavelength_position + 1]

                    first_bool = False
                wavelength_position += 1
        
        photflux_list.append(self.modelspectra.photflux[wavelength_position])
        smoothed_photflux_list.append(self.modelspectra.smoothPhotFlux[wavelength_position])
        sumflux_list.append(self.modelspectra.sumflux[wavelength_position])
        smoothed_sumflux_list.append(self.modelspectra.smoothSumFlux[wavelength_position])

        return (photflux_list, smoothed_photflux_list, sumflux_list, smoothed_sumflux_list, 
               first_bool, focused_wavelength, wavelength_position)

    # x is 3000 in length, one for each wavelength
    # y is 360 in length, one for each ohase image
    # z is the variability value for each x, y, pair (360, 3000)
    # Plot the color contour map
    def plot_color_variability(self, x, y, z):

        X, Y = np.meshgrid(x, y)
        plt.contourf(X, Y, z)
        cbar = plt.colorbar()

        cbar.set_label('B - A')
        plt.ylabel('Phase (0-360)')
        plt.xlabel('Wavelength (0-20 um)')

        plt.title('Variability Between Smoothed Reference and Smoothed SumFlux')
        plt.savefig('./ProxCen/SmoothedVariabilityGraphs/contourVariability.png')
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
    )

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
    phase0 = []
    phase180 = []

    # This for loop is for each hemisphere image
    for value in range(num_exposures):
        hemi_map = np.load('./ProxCen/HemiMapImages+Arrays/hemiArray_%d' % count, allow_pickle=True)

        # Calculates the star's combined spectrum for the given hemisphere and adds a column to the SM's
        # model info called sumflux
        SM.calculate_combined_spectrum(hemi_map)

        phase0, phase180 = SM.plot_variability(count, color_plot_values, phase0, phase180)
        SM.plot_smoothing_data_reduction(count)

        # (photflux_list, smoothed_photflux_list, sumflux_list, smoothed_sumflux_list,
        # first_bool, focused_wavelength, wavelength_position) = SM.find_focused_wavelength(photflux_list,
        #                                                             smoothed_photflux_list,
        #                                                             sumflux_list,
        #                                                             smoothed_sumflux_list,
        #                                                             first_bool, focused_wavelength,
        #                                                             wavelength_position, count)

        phase += deltaPhase
        print("Phase = ", phase)
        if (phase > 1):
            phase = phase % 1
        count += 1

        del SM.modelspectra['sumflux']
        del SM.modelspectra['smoothPhotFlux']
        del SM.modelspectra['smoothSumFlux']

    x = SM.modelspectra.wavelength
    y = phase0 - phase180
    plt.plot(x, y)
    plt.savefig("./ProxCen/VariabilityGraphs/phase0Phase180Comparison.png", dpi=200, bbox_inches='tight')
    plt.show()
    plt.close("all")

    x = SM.modelspectra.wavelength
    y = list(range(360))
    SM.plot_color_variability(x, y, color_plot_values)
    # color_plot_values = np.asarray(color_plot_values)
    # print(color_plot_values.size)

    # np.save('./ProxCen/SmoothedVariabilityGraphs/contourPlotValues.npy', color_plot_values)
    # temp = np.load('./ProxCen/VariabilityGraphs/contourPlotValues.npy')
    # print("Temp = ", temp)
    # print(type(temp))
    # print(temp.shape)
    print("Fin")
