import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
from numba import jit, njit
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib
import configparser

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
    grid_data=False,
    ngridpoints=None,
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
    if cut_range:
        model = model.loc[
            (model.wavelength >= rangemin) & (model.wavelength <= rangemax)
        ]
    
    # # This code plots the difference in wavelength values of the raw data. Shows aaresolution across the model
    # wavelengthDifference = []
    # count = 0
    # for value in model.wavelength:
    #     if count == len(model.wavelength) - 1:
    #         wavelengthDifference.append(0.001)
    #         print("End of file")
    #     else:
    #         diff = model.wavelength[count + 1] - value
    #         wavelengthDifference.append(diff)
    #         print(diff)
    #         if diff > 0.05:
    #             print(model.wavelength[count])
    #             print("Greater than 0.05")
    #             print(diff)
    #     count += 1
    
    # fig = plt.figure()
    # plt.plot(model.wavelength, wavelengthDifference)
    # plt.title("Wavelength Spacing of NextGen Data")
    # plt.ylabel("Delta Wavelength")
    # plt.xlabel("Wavelength (um)")
    # plt.savefig('./ProxCen/VariabilityGraphs/wavelengthSpacingOfRawData.png')
    # # plt.show()
    # print("Plot done")

    if grid_data:
        wavenew = np.linspace(rangemin, rangemax, ngridpoints)
        fluxnew = np.interp(wavenew, model.wavelength, model.flux)
        data = {"wavelength": wavenew, "flux": fluxnew}
        model = pd.DataFrame(data=data)
    return model

class Spotmodel:
    def __init__(
        self,
        spotcoverage,
        spotnumber,
    ):
        self.spotcoverage = spotcoverage
        self.spotnumber = spotnumber

    def generate_spots(self, randomSeed=None):
        if randomSeed is not None:
            np.random.seed(randomSeed)
        # Why is radius here chosen to be 1?
        surface_area = 4.0 * np.pi * 1.0 ** 2
        # Picks spotnumber of radius'. all between 0 and 1
        spot_radius = np.random.random_sample(self.spotnumber)
        print("Spot Radius= ", spot_radius)
        total_coverage = np.sum(np.pi * spot_radius ** 2)
        print("Total Coverage = ", total_coverage)
        normalization = surface_area * self.spotcoverage / total_coverage
        print("Normalization = ", normalization)
        true_radius = spot_radius * normalization ** 0.5
        print("True_Radius = ", true_radius)

        # Limits latitude to between 60 and -60? Based on Butterfly effect?
        lat = -60 + 120 * np.random.random_sample(self.spotnumber)
        # Limits Longitude to between -180 and 180
        lon = -180 + 360 * np.random.random_sample(self.spotnumber)

        surface_map = self.generate_flat_surface_map(true_radius, lon, lat,)
        print("Type of surface_map = ", type(surface_map))

        plt.close("all")
        return surface_map

    # @njit
    def generate_flat_surface_map(self, spot_radii, lon, lat):
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
        for spot in range(self.spotnumber):
            add_spots = ax.tissot(
                rad_km=spot_radii[spot] * rEarth,
                lons=lon[spot],
                lats=lat[spot],
                n_samples=1000,
                fc="k",
                alpha=1,
            )
        canvas.draw()
        buf = canvas.buffer_rgba()
        surface_map_image = np.asarray(buf)
        # 0 = photosphere
        # 1 = spot
        # 2 = planet
        surface_map = np.where(surface_map_image[:, :, 0] == 255, 0, 1)
        
        # Save and show the surface map values (1 or 0) that create the red/black map image (rectangular shape)
        # NOTE: the plt.show() function does not scale well with this plot, must view from file
        plt.savefig('./ProxCen/FlatMap.png')
        # plt.show()
        return surface_map

if __name__ == "__main__":
    
    configParser = configparser.RawConfigParser()
    while True:
        try:
            fileName = input("Config File Path ./Config/")
            configParser.read_file(open("./Config/%s" % fileName))
            break
        except FileNotFoundError:
            print("There is no file by that name, please try again.")
    
    spotcoverage = int(configParser.get('ProxCen', 'spotcoverage'))
    spotnumber = int(configParser.get('ProxCen', 'spotnumber'))

    # Turn spot coverage into a percentage
    spotcoverage /= 100

    SM = Spotmodel(
        spotcoverage,
        spotnumber,
    )

    surface_map = SM.generate_spots()

    # Saves the numpy array version of the flat surface map to be loaded into other programs
    np.save('./ProxCen/flatMap.npy', surface_map)
