import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statistics
import configparser
import ast
from star_hemimap import Hemimodel
from star_flatmap import readmodel

# This program creates plots of the star's spectrum and light curve as it rotates

class Spectramodel():
    def __init__(
        self,
        teffStar,
        starName,
        starspectrum,
        spotspectrum,
        faculaespectrum,
        Rotstar,
        surface_map,
        inclination,
        flux_min,
        flux_max,
    ):

        self.teffStar = teffStar
        self.starName = starName
        self._set_model_spectra(starspectrum, spotspectrum, faculaespectrum)
        self.Rotstar = Rotstar
        self.surface_map = surface_map
        self.inclination = inclination
        self.flux_min = flux_min
        self.flux_max = flux_max

    # Create a a Dataframe called modelspectra which houses the wavelength, photoflux, and spotflux values
    def _set_model_spectra(self, starspectrum, spotspectrum, faculaespectrum):
        if not np.all(starspectrum.wavelength == spotspectrum.wavelength) or not np.all(starspectrum.wavelength == faculaespectrum.wavelength):
            raise ValueError("The star, spot, and faculae spectra should be on the same wavelength scale")
        data = {'wavelength': starspectrum.wavelength, 'photflux': starspectrum.flux, 'spotflux': spotspectrum.flux, 'facflux': faculaespectrum.flux}
        self.modelspectra = pd.DataFrame(data)


    # Calculate the percentage of the current phase's surface taken up by the spots
    def calculate_coverage(self, hemi_map, ignore_planet=False):
        # add where the 3rd dimension index is certain color
        flat_image = hemi_map[~hemi_map.mask].flatten()
        total_size = flat_image.shape[0]
        photo = np.where(flat_image == 0)[0].shape[0]
        spot = np.where(flat_image == 1)[0].shape[0]
        fac = np.where(flat_image == 2)[0].shape[0]

        # planet  = np.where(flat_image == 2)[0].shape[0]

        if ignore_planet:
        #    total_size_mod = total_size - planet
            pass
        else:
            total_size_mod = total_size

        photo_frac =  photo / total_size_mod
        spot_frac = spot / total_size_mod
        fac_frac = fac / total_size_mod
        return photo_frac, spot_frac, fac_frac

    # Calculate the total output flux value of the current phase based on the percentage of surface area taken up by spots
    def calculate_combined_spectrum(self, hemi_map):
        photo_frac, spot_frac, fac_frac = self.calculate_coverage(hemi_map,)
        self.modelspectra.insert(4, 'sumflux', (self.modelspectra.photflux * photo_frac) + (self.modelspectra.spotflux * spot_frac) + (self.modelspectra.facflux * fac_frac))

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
        
        plt.savefig("./%s/SpectrumGraphs/combinedHemiSpectrum_%d" % (self.starName, count))
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
            plt.savefig("./%s/BinnedSpectraNorm/binnedNormHemiSpectrum_%d.png" % (self.starName, count))
        else:
            plt.savefig("./%s/BinnedSpectra/binnedHemiSpectrum_%d.png" % (self.starName, count))
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
        plt.savefig('./%s/LightCurves/LightCurve_%d' % (self.starName, count), bbox_inches='tight')
        # plt.show()
        plt.close("all")

        return x, y