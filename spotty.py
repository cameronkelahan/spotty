import configparser
import numpy as np
import os
import star_flatmap
import star_hemimap
import star_spectra
import read_model

##########################################################################################################
# 1) Create a configparser object to read in the info from the config file
#   - User inputs name of config file; found in Config folder
if __name__ == "__main__":
    
    configParser = configparser.RawConfigParser()
    while True:
        try:
            fileName = input("Config File Path ./Config/")
            configParser.read_file(open("./Config/%s" % fileName))
            break
        except FileNotFoundError:
            print("There is no file by that name, please try again.")
    
##########################################################################################################
    # 2) If there is no folder with the star's name as given in the config file, make one
    #   - All files produced will be saved into this directory, based on this star's config file
    starName = configParser.get('Star', 'starName')

    if not os.path.isdir('./%s/' % starName):
        os.mkdir('./%s/' % starName)

    # Read in the information of the star from the config file
    spotCoverage = int(configParser.get('Star', 'spotCoverage'))
    # Turn spot coverage into a percentage
    spotCoverage = spotCoverage/100
    spotNumber = int(configParser.get('Star', 'spotNumber'))
    
    facCoverage = int(configParser.get('Star', 'facCoverage'))
    # Turn facCoverage into a percentage
    facCoverage = facCoverage/100
    facNumber = int(configParser.get('Star', 'facNumber'))

    # If a flat 2D map of a star and its spots has not been created, make one
    #   - This means there is no star model yet
    if not os.path.isfile('./%s/FlatMap.png' % starName):
        # Create a 2D spotmmodel from the star_flatmap.py file
        SpotModel = star_flatmap.Spotmodel(
            spotCoverage,
            spotNumber,
            facCoverage,
            facNumber,
            starName,
        )
        
        # Generate the spots on the star
        surface_map = SpotModel.generate_spots()

        # Convert the 1's and 0's in the ndarray, which store the spot locations, to a smaller data type
        surface_map = surface_map.astype(np.int8)

        # Saves the numpy array version of the flat surface map to be loaded while creating the hemisphere views
        np.save('./%s/flatMap.npy' % starName, surface_map)

##########################################################################################################
    # 3) With a FlatMap.npy file created, the program then creates hemisphere views of the star
    #   - The number of hemisphere views created are based on the information in the conig file
    #       - Number of images taken, time between images taken, and the star's rotation

    # Temperature values for the star, spots, and faculae
    teffStar = int(configParser.get('Star', 'teffStar'))
    teffSpot = int(configParser.get('Star', 'teffSpot'))
    teffFac = int(configParser.get('Star', 'teffFac'))

    generateHemispheres = configParser.getboolean('Star', 'generateHemispheres')
    print(generateHemispheres)
    print(type(generateHemispheres))

    # Rotation of star in days
    Rotstar = float(configParser.get('Star', 'Rotstar'))

    # Inclination of the star
    inclination = int(configParser.get('Star', 'Inclination'))

    # Name of the .npy array which contains the flat surface map
    surface_map_file = configParser.get('HemiMap', 'surface_map_file')
    surface_map_file = surface_map_file.strip('"')
    surface_map = np.load(surface_map_file)

    # Total number of exposures to be takn
    num_exposures = int(configParser.get('HemiMap', 'num_exposures'))
    # Time (in days) between exposures
    time_between_exposures = float(configParser.get('HemiMap', 'time_between_exposures'))

    # Calculates what percent the time between exposures is compared to the stellar rotation period
    # Used to calculate the change in phase between images
    exposure_time_percent = time_between_exposures * (100 / Rotstar)
    print("Exposure Time percent = ", exposure_time_percent)
    deltaPhase = exposure_time_percent / 100
    print("Delta Phase = ", deltaPhase)

    # If the hemisphere maps have not been ccreated yet, create them
    if generateHemispheres:
        # Create a HemiMap object
        HM = star_hemimap.Hemimodel(
            teffStar,
            Rotstar,
            surface_map,
            inclination,
        )

        # Begin at phase = 0
        # Count keeps track of which hemisphere map image is being looked at currently
        phase = 0
        count = 0
        for value in range(num_exposures):
            hemi_map = HM.generate_hemisphere_map(phase, count)
            phase += deltaPhase
            print("Phase = ", phase)
            if (phase > 1):
                phase = phase % 1
            count += 1
        
##########################################################################################################
    # 4) Begin analyzing the flux output from the star

    # The file names for the wavelength/flux values of the star, spots, and faculae
    phot_model_file = configParser.get('Star', 'phot_model_file')
    phot_model_file = phot_model_file.strip('"') # configParser adds extra "" that I remove
    spot_model_file = configParser.get('Star', 'spot_model_file')
    spot_model_file = spot_model_file.strip('"')
    fac_model_file = configParser.get('Star', 'fac_model_file')
    fac_model_file = fac_model_file.strip('"')

    # Boolean which says whether or not the user wishes to bin the data to a certain resolving power
    binData = configParser.getboolean('Star', 'binData')
    # The resolving power to bin to
    resolvingPower = int(configParser.get('Star', 'resolvingPower'))
    binnedWavelengthMin = float(configParser.get('Star', 'binnedWavelengthMin'))
    binnedWavelengthMax = float(configParser.get('Star', 'binnedWavelengthMax'))

    topValues = []
    cwValues = []
    CW = binnedWavelengthMin
    while CW < binnedWavelengthMax:
        deltaLambda = CW / resolvingPower
        topValue = CW + (deltaLambda / 2)
        topValues.append(topValue)
        cwValues.append(CW)
        CW = CW + deltaLambda

    starspectrum, starspectrumString = read_model.readmodel(
        phot_model_file,
        starName,
        cut_range=True,
        rangemin=0.1999,  # in micron, based on MIRECLE constraints
        rangemax=20.5,  # in micron, based on  MIRECLE constraints
        bin_data=binData,
        ngridpoints=3000,
        resolvingPower=resolvingPower,
        topValues=topValues,
        cwValues=cwValues,
    )

    if binData:
        starspectrumStringCSV = starspectrumString.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
        file = open(r'./BinnedNextGenModels/binned%sStellarModel.txt' % teffStar,'w')
        file.write(starspectrumStringCSV)
        file.close()

    spotspectrum, spotspectrumString = read_model.readmodel(
        spot_model_file,
        starName,
        cut_range=True,
        rangemin=0.1999,  # in micron, based on MIRECLE constraints
        rangemax=20.5,  # in micron, based on  MIRECLE constraints
        bin_data=binData,
        ngridpoints=3000,
        resolvingPower=resolvingPower,
        topValues=topValues,
        cwValues=cwValues,
    )

    if binData:
        binnedSpotspectrumCSV = spotspectrumString.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
        file = open(r'./BinnedNextGenModels/binned%sStellarModel.txt' % teffSpot,'w')
        file.write(binnedSpotspectrumCSV)
        file.close()

    faculaespectrum, faculaespectrumString = read_model.readmodel(
        fac_model_file,
        starName,
        cut_range=True,
        rangemin=0.1999,  # in micron, based on MIRECLE constraints
        rangemax=20.5,  # in micron, based on  MIRECLE constraints
        bin_data=binData,
        ngridpoints=3000,
        resolvingPower=resolvingPower,
        topValues=topValues,
        cwValues=cwValues,
    )

    if binData:
        binnedFaculaespectrumCSV = faculaespectrumString.to_csv(index=False, header=['WAVELENGTH (MICRONS)','FLUX (ERG/CM2/S/A)'], sep=' ')
        file = open(r'./BinnedNextGenModels/binned%sStellarModel.txt' % teffFac,'w')
        file.write(binnedFaculaespectrumCSV)
        file.close()

        # Create a Spectramodel class
    SM = star_spectra.Spectramodel(
        teffStar,
        starName,
        starspectrum,
        spotspectrum,
        faculaespectrum,
        Rotstar,
        surface_map,
        inclination,
        flux_min = min(spotspectrum.flux),
        flux_max = max(faculaespectrum.flux),
    )

    phase = 0
    count = 0
    x = []
    y = []

    hemi_fluxes = [[SM.modelspectra.wavelength], []]

    # This for loop is for each hemisphere image/each phase
    for value in range(num_exposures):
        # Loads in the hemiMap information for the current phase
        hemi_map = np.load('./%s/HemiMapImages+Arrays/hemiArray_%d' % (starName, count), allow_pickle=True)

        # Calculates the star's combined spectrum for the given hemisphere and adds a column to the SM's
        # model info called sumflux
        SM.calculate_combined_spectrum(hemi_map)
        SM.plot_combined_spectrum(count)

        hemi_fluxes[1].append(SM.modelspectra.sumflux)

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

    hemi_fluxes = np.asarray(hemi_fluxes)
    np.save('./%s/SumfluxArraysByPhase/wavelengthsAndSumfluxValuesByPhase.npy' % starName, hemi_fluxes, allow_pickle=True)

    loaded = np.load('./%s/SumfluxArraysByPhase/wavelengthsAndSumfluxValuesByPhase.npy' % starName, allow_pickle=True)
    print("Done")