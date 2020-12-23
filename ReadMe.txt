The user interacts with the input parameters through a config file. There is an example config file in the Config folder called "proxCen.cfg" with information for the star Proxima Centauri.
- This is also where the user inputs how many images are to be taken, across what time period. The program then calulates how often an image is to be taken.
- Stellar/Spot temperature definition
- size, rotation, radius info of the star
- number of sunspots and what percent of the surface they cover
- can specify an inclination of the star

The use of this code has an order.

1)  If 'FlatMap.npy' already exists in the current star's folder (ex: ProxCen), skip this step. If not, the user must run star_flatmap.py in order to generate the 'flat surface map' of the star. This is necessary for the next step.

2) If the "HemiMapImages+Arrays" folder inside the current star's folder is already populated with .png and .txt files, skip this step. Otherwise, run this program to create the Hemishpere Images of the star and save the arrays of that information. The number
of hemisphere map images and arrays created will equal the 'num_exposures' value in the chosen config file.
- **Note: Current setup creates hemispheres with gridshape of 3000x3000. This take a long time to run. You can lower the 'regrid_shape' variable in 'generate_hemisphere_map' to create hemisphere maps quicker, but at a lower resolution.

3) Once the hemisphere maps/arrays have been created, the rest of the programs can be run in any order, as they all
only depend on the hemisphere arrays.

- star_spectra.py:
    This program will plot the flux spectrum of the current phase's hemisphere. It will create a new plot every time
    the phase increments (a.k.a. the star rotates and a new 'image' is taken). It will also create a plot depicting
    the star's light curve as it rotates.

- variability.py
    This program creates 2 plots based on two equations
        1) # Geronimo's first equation.
        B = (self.modelspectra.sumflux/self.modelspectra.photflux) - 1.0

            - This is the star's photo flux divided by the linearly combined flux, yielding a value close to 1, then it
              is subtracted by 1 to obtain a very small number close to 0.

        2) # Geronimo's second equation
           # Variance between the smoothed versions of photflux and sumflux
        A = (self.modelspectra.smoothSumFlux/self.modelspectra.smoothPhotFlux) - 1.0

            - This equation is the same as the first equation, but it uses the smoothed combined flux and the smoothed
              photo flux.
    
    A graph depicting equation 1 is created, then a graph depicting B - A is created.
    
    **Update: Also creates a graph showing the difference in flux of each wavelength at phases 0 and 180.
    
- variabilityWoutSmoothing
    - This program creates two colormaps which show each wavelength's variability across all phases/images taken.
    - The colors on the map represents the fractional change in parts per million or log(parts per million) respectively

- make_gif.py
    This program will create .gif's of the images created from the other programs. The user can specify which .gif's to
    make.
