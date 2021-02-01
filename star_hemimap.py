import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import configparser

# This program creates and saves the images and arrays of the hemishpere maps produced
# as the star rotates.
# Loads in the flatmap model that has the spots on it.

class Hemimodel():
    def __init__(
        self,
        teffStar,
        Rotstar,
        surface_map,
        inclination,
    ):
        self.teffStar = teffStar
        self.Rotstar = Rotstar
        self.surface_map = surface_map
        self.inclination = inclination
    
    def generate_hemisphere_map(self, phase, count):
        # phase is between 0 and 1
        lon = phase * 360
        if np.abs(lon - 180) < 0.01:
            lon += 0.01  # needs a litle push at 180 degrees
        image_lon_min = -180 + lon
        image_lon_max = 180 + lon
        proj = ccrs.Orthographic(
            central_longitude=0.0, central_latitude=self.inclination)
        fig = plt.figure(figsize=(5, 5), dpi=100, frameon=False)

        ax = plt.gca(projection=proj, fc="r")
        ax.outline_patch.set_linewidth(0.0)
        hemi_map = ax.imshow(
            self.surface_map,
            origin="upper",
            transform=ccrs.PlateCarree(),
            extent=[image_lon_min, image_lon_max, -90, 90],
            interpolation="none",
            regrid_shape=100
            # # Optional regrid_shape of 100 runs significantly faster
            # regrid_shape=100
        ).get_array()

        plt.title("Hemishpere Map: Prox Cen\nT=%dK; Log g=5; Met=0" % self.teffStar)

        # Saves each hemishpere map image to file
        # Number of saved images will be equal to num_exposures value in config file
        plt.savefig("./ProxCen/HemiMapImages+Arrays/hemiMap_%d.png" % count)
        # plt.show()
        plt.close("all")

        # Save the pickle numpy array version of the hemisphere map to load in other programs later
        hemi_map.dump('./ProxCen/HemiMapImages+Arrays/hemiArray_%d' % count)

        return hemi_map