The use of this code has an order.

1)  If 'FlatMap.npy' already exists in the ProxCen folder, skip this step. The user must run star_flatmap.py
in order to generate the 'flat surface map' of the star. This is necessary for the next step.

2) If the "HemiMapImages+Arrays" folder is already populated with .png and .txt files, skip this step. Otherwise,
run this program to create the Hemishpere Images of the star and save the arrays of that information. The number
of hemisphere map images and arrays created will equal the 'num_exposures' value in the chosen config file.

3) Once the hemisphere arrays have been created, the rest of the programs can be run in any prder, as they all
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

- make_gif.py
    This program will create .gif's of the images created from the other programs. The user can specify which .gif's to
    make.