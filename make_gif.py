import imageio
import configparser
from PIL import Image

# This program takes plot images of the star as it rotates and turns them into a .gif
# Must comment/uncomment the files you wish to turn into .gif's

def make_gif(num_exposures):
    
    filenamesHemi = []
    filenamesSpectrum = []
    filenamesLightCurve = []
    filenamesBinnedSpectrum = []
    filenamesNormalizedBinnedSpectrum = []
    filenamesVariability = []
    filenamesSmoothedVariability = []
    
    for value in range(num_exposures):
        filenamesHemi.append("./ProxCen/HemiMapImages+Arrays/hemiMap_%d.png" % value)
        filenamesSpectrum.append("./ProxCen/SpectrumGraphs/combinedHemiSpectrum_%d.png" % value)
        filenamesLightCurve.append("./ProxCen/LightCurves/LightCurve_%d.png" % value)
        filenamesBinnedSpectrum.append("./ProxCen/BinnedSpectra/binnedHemiSpectrum_%d.png" % value)
        filenamesNormalizedBinnedSpectrum.append("./ProxCen/BinnedSpectraNorm/binnedNormHemiSpectrum_%d.png" % value)
        filenamesVariability.append("./ProxCen/VariabilityGraphs/variableFlux_%d.png" % value)
        filenamesSmoothedVariability.append("./ProxCen/SmoothedVariabilityGraphs/normVsSmoothVariability_%d.png" % value)

    # Creates a .gif of all the hemisphere phase images
    imagesHemi = []
    imagesLC = []
    combinedImages = []
    count = 0
    for filenameHemi in filenamesHemi:
        imagesHemi.append(imageio.imread(filenameHemi))
        imagesLC.append(imageio.imread(filenamesLightCurve[count]))
        
        #Read the two images
        hemi = Image.open(filenameHemi)
        # hemi.show()
        lc = Image.open(filenamesLightCurve[count])
        # lc.show()
        #resize, first image
        # image1 = image1.resize((426, 240))
        hemi_size = hemi.size
        lc_size = lc.size
        lc = lc.resize((500, 500))
        lc_size = lc.size
        new_image = Image.new('RGB',(2*hemi_size[0], hemi_size[1]), (250,250,250))
        new_image.paste(hemi,(0,0))
        new_image.paste(lc,(hemi_size[0],0))
        combinedImages.append(new_image)
        # new_image.save("images/merged_image.jpg","JPEG")
        # new_image.show()
        count += 1

    imageio.mimsave("./ProxCen/HemiMapImages+Arrays/hemiMapMovie.gif", imagesHemi)
    imageio.mimsave("./ProxCen/LightCurves/lightCurveMovie.gif", imagesLC)
    imageio.mimsave("./ProxCen/Hemi+LightCurve/HemiLightCurveMovie.gif", combinedImages)

    # Creates a .gif of all the combined spectra images
    images = []
    for filename in filenamesSpectrum:
        images.append(imageio.imread(filename))
    imageio.mimsave("./ProxCen/SpectrumGraphs/combinedHemiSpectrum.gif", images)

    # # Creates a .gif of the light curve of the star as it rotates
    # images = []
    # for filename in filenamesLightCurve:
    #     images.append(imageio.imread(filename))
    # imageio.mimsave("./ProxCen/LightCurves/lightCurveMovie.gif", images)

    # Creates a .gif of all the binned spectra plots
    # images = []
    # for filename in filenamesBinnedSpectrum:
    #     images.append(imageio.imread(filename))
    # imageio.mimsave("./ProxCen/BinnedSpectra/binnedSpectraMovie.gif", images)

    # Creates a .gif of the normalized (0-1) binned spectra values
    # images = []
    # for filename in filenamesNormalizedBinnedSpectrum:
    #     images.append(imageio.imread(filename))
    # imageio.mimsave("./ProxCen/BinnedSpectraNorm/binnedNormSpectraMovie.gif", images)

    # Creates a .gif of the variable flux plots
    # images = []
    # for filename in filenamesVariability:
    #     images.append(imageio.imread(filename))
    # imageio.mimsave("./ProxCen/VariabilityGraphs/variableFluxMovie.gif", images)

    # Creates a .gif of the smoothed variable flux plots
    images = []
    for filename in filenamesSmoothedVariability:
        images.append(imageio.imread(filename))
    imageio.mimsave("./ProxCen/SmoothedVariabilityGraphs/variableSmoothedFluxMovie.gif", images)

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
    
    num_exposures = int(configParser.get('HemiMap', 'num_exposures'))

    make_gif(num_exposures)