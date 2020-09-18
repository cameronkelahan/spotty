import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
from scipy import stats
from numba import jit, njit
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib

interactive_plot = True
if not interactive_plot:
    matplotlib.use("Agg")
    plt.interactive(False)
    matplotlib.rcParams["figure.dpi"] = 300

##  todo: use astropy units
        # non-zero impact parameters
        # ingress/egress/grazing transits

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
        filename, header=7, delim_whitespace=True, names=["wavelength", "flux"]
    )
    # lets convert to micron
    model.wavelength /= 1.0e4
    if cut_range:
        model = model.loc[
            (model.wavelength >= rangemin) & (model.wavelength <= rangemax)
        ]
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
        starspectrum,
        spotspectrum,
        Mstar,
        Rstar,
        Rotstar,
        hasPlanet=False,
        Rplanet=None,
        Pplanet=None,
        Impactplanet=None,
        TransitDuration=None,
    ):
        self._set_model_spectra(starspectrum, spotspectrum)
        self.spotcoverage = spotcoverage
        self.spotnumber = spotnumber
        self.Mstar = Mstar
        self.Rstar = Rstar
        self.Rotstar = Rotstar
        self.hasPlanet = hasPlanet
        if self.hasPlanet:
            self.Rplanet = Rplanet
            self.Pplanet = Pplanet
            self.Impactplanet = Impactplanet
            self.TransitDuration = TransitDuration
        else:
            self.Rplanet = None
            self.Pplanet = None
            self.Impactplanet = None
            self.TransitDuration = None
        # I'll worry about this part later
        if int(self.Impactplanet) != 0:
            raise NotImplementedError
        self._calculate_planet_values()

    def _calculate_planet_values(self,):
        if self.hasPlanet:
            self.rprs = self.Rplanet / self.Rstar * rEarth_to_rSun
        else:
            self.rprs = None
    
    def _set_model_spectra(self, starspectrum, spotspectrum):
        if not np.all(starspectrum.wavelength == spotspectrum.wavelength):
            raise ValueError("The star and spot spectra should be on the same wavelength scale")
        data = {'wavelength': starspectrum.wavelength, 'photflux': starspectrum.flux, 'spotflux': spotspectrum.flux}
        self.modelspectra = pd.DataFrame(data)

    def observations(self, numvisits=1, obsPerVisit=4, intsPerOrbit=10):

        for visit in range(numvisits):
            self._single_observation()

    @njit
    def _single_observation(self):
        pass

    def generate_spots(self, randomSeed=None):
        if randomSeed is not None:
            np.random.seed(randomSeed)
        surface_area = 4.0 * np.pi * 1.0 ** 2
        spot_radius = np.random.random_sample(self.spotnumber)
        total_coverage = np.sum(np.pi * spot_radius ** 2)
        normalization = surface_area * self.spotcoverage / total_coverage
        true_radius = spot_radius * normalization ** 0.5

        lat = -60 + 120 * np.random.random_sample(self.spotnumber)
        lon = -180 + 360 * np.random.random_sample(self.spotnumber)

        surface_map = self.generate_flat_surface_map(true_radius, lon, lat,)
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
        return surface_map

    def generate_hemisphere_map(self, surface_map, phase):
        # phase is between 0 and 1
        lon = phase * 360
        if np.abs(lon - 180) < 0.01:
            lon += 0.01  # needs a litle push at 180 degrees
        image_lon_min = -180 + lon
        image_lon_max = 180 + lon
        proj = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)
        fig = plt.figure(figsize=(5, 5), dpi=100, frameon=False)
        ax = plt.gca(projection=proj, fc="r")
        ax.outline_patch.set_linewidth(0.0)
        hemi_map = ax.imshow(
            surface_map,
            origin="upper",
            transform=ccrs.PlateCarree(),
            extent=[image_lon_min, image_lon_max, -90, 90],
            interpolation="none",
            regrid_shape=3000
        ).get_array()
        plt.close("all")
        return hemi_map

    def add_planet_to_image(self, hemi_map, transit_fraction):
        # this makes the image grid larger to account for ingress/egress
    

        # no grazing eclipses for the time being
        if (transit_fraction < -1) or (transit_fraction > -1):
            return hemi_map
        center = np.asarray(hemi_map.shape) // 2

        # center[1] is both the location of the center and the radius
        planet_location_y = np.round(
            (self.Impactplanet * center[1]) + center[1]
        ).astype(int)

        # for the time being transit_fraction is for a central transit
        planet_location_x = np.round(transit_fraction * hemi_map.shape[0]).astype(int)

        # we are going to use a square planet for numerical reasons
        Astar = np.pi * center[0] ** 2
        Aplanet = self.rprs ** 2 * Astar
        square_radius = np.ceil(Aplanet ** 0.5 * 0.5).astype(int)
        hemi_map[
            planet_location_y - square_radius : planet_location_y + square_radius,
            planet_location_x - square_radius : planet_location_x + square_radius,
        ] = 2
        return hemi_map
    
    def calculate_transit_fraction(self, ):



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
    
    def calculate_spectrum(self, hemi_map):
        photo_frac, spot_frac = self.calculate_coverage(hemi_map,)
        self.modelspectra.sumflux = (self.modelspectra.photflux * photo_frac) + (self.modelspectra.spotflux * spot_frac)



if __name__ == "__main__":
    phot_model_file = "T3500_g5.0_solar.txt"
    spot_model_file = "T3000_g5.0_solar.txt"
    starspectrum = readmodel(
        phot_model_file,
        cut_range=True,
        rangemin=0.350,  # in micron
        rangemax=2.000,  # in micron
        grid_data=True,
        ngridpoints=3000,
    )
    spotspectrum = readmodel(
        spot_model_file,
        cut_range=True,
        rangemin=0.350,  # in micron
        rangemax=2.000,  # in micron
        grid_data=True,
        ngridpoints=3000,
    )

    spotcoverage = 0.05
    spotnumber = 20
    Mstar = 0.4
    Rstar = 0.4
    Rotstar = 5
    hasPlanet = True
    Rplanet = 4
    Pplanet = 10
    Impactplanet = 0
    TransitDuration = 3.0
    SM = Spotmodel(
        spotcoverage,
        spotnumber,
        starspectrum,
        spotspectrum,
        Mstar,
        Rstar,
        Rotstar,
        hasPlanet,
        Rplanet,
        Pplanet,
        Impactplanet,
        TransitDuration,
    )

    surface_map = SM.generate_spots()
    phase = 0
    hemi_map = SM.generate_hemisphere_map(surface_map, phase)
    transit_fraction = 0.1
    hemi_map_with_planet = SM.add_planet_to_image(hemi_map, transit_fraction)
    SM.calculate_spectrum(hemi_map_with_planet)

