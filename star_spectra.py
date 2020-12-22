import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statistics
import configparser
import ast
from star_hemimap import Hemimodel
from star_flatmap import readmodel

# This program creates plots of the star's spectrum and light curve as it rotates

class Spectramodel(Hemimodel):
    def __init__(
        self,
        teffStar,
        starspectrum,
        spotspectrum,
        Rotstar,
        surface_map,
        inclincation,
        flux_min,
        flux_max,
    ):
        # Inherit parent class's methods and attributes
        super().__init__(
            teffStar,
            Rotstar,
            surface_map,
            inclination,
        )

        self.teffStar = teffStar
        self._set_model_spectra(starspectrum, spotspectrum)
        self.Rotstar = Rotstar
        self.surface_map = surface_map
        self.inclination = inclination
        self.flux_min = flux_min
        self.flux_max = flux_max

    # Create a a Dataframe called modelspectra which houses the wavelength, photoflux, and spotflux values
    def _set_model_spectra(self, starspectrum, spotspectrum):
        if not np.all(starspectrum.wavelength == spotspectrum.wavelength):
            raise ValueError("The star and spot spectra should be on the same wavelength scale")
        data = {'wavelength': starspectrum.wavelength, 'photflux': starspectrum.flux, 'spotflux': spotspectrum.flux}
        self.modelspectra = pd.DataFrame(data)


    # Calculate the percentage of the current phase's surface taken up by the spots
    def calculate_coverage(self, hemi_map, ignore_planet=False):
        flat_image = hemi_map[~hemi_map.mask].flatten()
        total_size = flat_image.shape[0]
        photo = np.where(flat_image == 0)[0].shape[0]
        spot = np.where(flat_image == 1)[0].shape[0]
        planet  = np.where(flat_image == 2)[0].shape[0]

        if ignore_planet:
           total_size_mod = total_size - planet
        else:
            total_size_mod = total_size

        photo_frac =  photo / total_size_mod
        spot_frac = spot / total_size_mod
        return photo_frac, spot_frac

    # Calculate the total output flux value of the current phase based on the percentage of surface area taken up by spots
    def calculate_combined_spectrum(self, hemi_map):
        photo_frac, spot_frac = self.calculate_coverage(hemi_map,)
        self.modelspectra.insert(3, 'sumflux', (self.modelspectra.photflux * photo_frac) + (self.modelspectra.spotflux * spot_frac))

    # Plot the spectra of the linearly combined spot and photo flux values based on the current phase
    def plot_combined_spectrum(self, count):
        fig = plt.figure()

        plt.xlim(-1, 21)

        ymin = self.flux_min * 0.05
        ymin = self.flux_min - ymin
        ymax = self.flux_max * 0.05
        ymax = self.flux_max + ymax
        plt.ylim(ymin, ymax)

        x = self.modelspectra.wavelength
        y = self.modelspectra.sumflux
        
        plt.plot(x, y)
        plt.xlabel("Wavelength (Microns)")
        plt.ylabel("Flux (ERG/CM2/S/A)")
        plt.title("Linearly Combined Star and Spot Flux Values")
        
        plt.savefig("./ProxCen/SpectrumGraphs/combinedHemiSpectrum_%d" % count)
        # plt.show()
        plt.close("all")

    def normalize(self, y):
        ymin = min(y)
        ymax = max(y)

        count = 0
        for value in y:
            y[count] = (y[count] - ymin)/(ymax - ymin)
            count += 1

        return y

    # Added functionality to plot specified wavelength bins (e.g. 5-7.5 microns)
    def plot_binned_spectrum(self, list_bins, count, normalized_bins):
        fig = plt.figure()

        plt.xlim(-1, 21)

        # These lists contain the lists of x and y values for each bin's plots
        list_of_x_value_lists = []
        list_of_y_value_lists = []

        # Keeps track of which bin's data is being created
        # User can specify any number of bins, greater than 0
        current_bin = 0
        for bin in list_bins:
            bins_x = []
            bins_y = []

            index = 0
            for wavelength in self.modelspectra.wavelength.values:
                # print("Current Wavelength = ", wavelength)
                
                # If the current wavelength is in the bin, append the wavelength to the x-values list
                # and append the sumflux value tot he y_calues list
                if wavelength >= list_bins[current_bin][0] and wavelength <= list_bins[current_bin][1]:
                    bins_x.append(wavelength)
                    bins_y.append(self.modelspectra.sumflux[index])
                
                index += 1
            
            if normalized_bins:
                # Normalize the values within a bin to be between 0 and 1
                bins_y = self.normalize(bins_y)

            # After all the appropraite wavelengths and sumflux's have been added, add the lists of x and y values
            # for this bin to the master list of x and y values for all of the bins.
            list_of_x_value_lists.append(bins_x)
            list_of_y_value_lists.append(bins_y)

            current_bin += 1
        
        if normalized_bins:
            plt.ylim(0, 1)
        else:
            ymin = self.flux_min * 0.05
            ymin = self.flux_min - ymin

            ymax = self.flux_max * 0.05
            ymax = self.flux_max + ymax

            plt.ylim(ymin, ymax)

        plot_count = 0
        for bin_list in range(len(list_of_x_value_lists)):
            label = 'Bin Values: %.3f - %.3f' % (list_bins[plot_count][0], list_bins[plot_count][1])
            plt.plot(list_of_x_value_lists[plot_count], list_of_y_value_lists[plot_count], alpha=0.5, label=label)
            plt.legend(loc='upper right')
            plot_count += 1
        if normalized_bins:
            plt.savefig("./ProxCen/BinnedSpectraNorm/binnedNormHemiSpectrum_%d.png" % count)
        else:
            plt.savefig("./ProxCen/BinnedSpectra/binnedHemiSpectrum_%d.png" % count)
        # plt.show()
        plt.close("all")

    # Plot a time-series light curve as the star rotates
    def plot_light_curve(self, count, x, y):
        # Calculate the mean value for the current bin, for the current star hemisphere
        avg_sumflux = statistics.mean(self.modelspectra.sumflux)

        # The X-Axis represents the number of the image taken. Presents itself as a time series, from firt image
        # to last
        x.append(count)

        # The y value is the flux value of the current hemisphere image
        y.append(avg_sumflux)

        plt.plot(x, y)
        plt.yscale("log")
        plt.ylabel("Flux (ERG/CM2/S/A)")
        plt.xlabel("Phase (0-360)")
        plt.savefig('./ProxCen/LightCurves/LightCurve_%d' % count, bbox_inches='tight')
        # plt.show()
        plt.close("all")

        return x, y

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

    # Rotation of star in days
    Rotstar = float(configParser.get('ProxCen', 'Rotstar'))

    # Rotation of star in days
    inclination = float(configParser.get('ProxCen', 'Inclination'))

    # Name of the .npy array which contains the flat surface map
    surface_map_file = configParser.get('HemiMap', 'surface_map_file')
    surface_map_file = surface_map_file.strip('"')
    surface_map = np.load(surface_map_file)

    # Total number of exposures to be takn
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
    else: normalized_bins = False
    
    # Calculates what percent the time between exposures is compared to the stellar rotation period
    # Used to calculate the change in phase between images
    exposure_time_percent = time_between_exposures * (100 / Rotstar)
    # print("Exposure Time percent = ", exposure_time_percent)
    deltaPhase = exposure_time_percent / 100
    # print("Delta Phase = ", deltaPhase)
    
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
    SM = Spectramodel(
        teffStar,
        starspectrum,
        spotspectrum,
        Rotstar,
        surface_map,
        inclination,
        flux_min = min(spotspectrum.flux),
        flux_max = max(starspectrum.flux),
    )

    phase = 0
    count = 0
    x = []
    y = []
    # This for loop is for each hemisphere image/each phase
    for value in range(num_exposures):
        # Loads in the hemiMap information for the current phase
        hemi_map = np.load('./ProxCen/HemiMapImages+Arrays/hemiArray_%d' % count, allow_pickle=True)

        # Calculates the star's combined spectrum for the given hemisphere and adds a column to the SM's
        # model info called sumflux
        SM.calculate_combined_spectrum(hemi_map)
        SM.plot_combined_spectrum(count)

        # This keeps track of which bin is being looked at currently
        list_index = 0
        # Plot the flux values of the chose n binned wavelengths
        # SM.plot_binned_spectrum(list_bins, count, normalized_bins)

        x, y = SM.plot_light_curve(count, x, y)

        # Advance the ohase; a.k.a. rotate the star
        phase += deltaPhase
        print("Phase = ", phase)
        if (phase > 1):
            phase = phase % 1
        count += 1

        # remove this phase's sumflux column of the modelspectra data frame of the SM object to make way for the next phase
        del SM.modelspectra['sumflux']