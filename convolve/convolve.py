from astropy.io import fits
from radio_beam import Beam
from astropy.wcs import WCS
from radio_beam import Beams
from astropy.convolution import convolve
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import Paths.Paths as paths

Path = paths.filepaths()


def convolveb3b6(b3data, b6data, outdir, label):
    fitsdata_b6 = fits.open(b6data)
    image_b6 = fitsdata_b6[0].data[0][0]
    fitsdata_b3 = fits.open(b3data)
    image_b3 = fitsdata_b3[0].data[0][0]
    
    hdrNB6 = fits.getheader(b6data)  
    hdrNB3 = fits.getheader(b3data)  
    wcsNB6 = WCS(hdrNB6,naxis=2)
    wcsNB3 = WCS(hdrNB3,naxis=2)

    my_beamNB3 = Beam.from_fits_header(hdrNB3)
    my_beamNB6 = Beam.from_fits_header(hdrNB6)

    scaleNB6 = wcsNB6.proj_plane_pixel_scales()[0]
    scaleNB3 = wcsNB3.proj_plane_pixel_scales()[0]

    beamsN =  Beams(beams=[my_beamNB3,my_beamNB6])

    common_beam = beamsN.common_beam()
    print(common_beam.major)
    
    area_rat_B3 = (common_beam.sr/my_beamNB3.sr).value

    if area_rat_B3!=1:
        print('convolution', b3data)
        kernelB3 = common_beam.deconvolve(my_beamNB3).as_kernel(scaleNB3)
        conv_B3 = convolve(image_b3, kernelB3,preserve_nan=True)
        #conv_B3 = conv_B3 * area_rat_B3
        if common_beam.major.unit=='arcsec':
            hdrNB3['BMAJ'] = common_beam.major.value/3600
            hdrNB3['BMIN'] = common_beam.minor.value/3600
            hdrNB3['BPA'] = common_beam.pa.value
        else:
            hdrNB3['BMAJ'] = common_beam.major.value
            hdrNB3['BMIN'] = common_beam.minor.value
            hdrNB3['BPA'] = common_beam.pa.value
        fits.writeto(outdir+'/%s_B3_conv.fits'%label, conv_B3, hdrNB3, overwrite = True)
    
    area_rat_B6 = (common_beam.sr/my_beamNB6.sr).value
    if area_rat_B6 != 1:
        print('convolution', b6data)
        kernelB6 = common_beam.deconvolve(my_beamNB6).as_kernel(scaleNB6)
        conv_B6 = convolve(image_b6, kernelB6,preserve_nan=True)  
        #conv_B6 = conv_B6 * area_rat_B6
        if common_beam.major.unit=='arcsec':
            hdrNB6['BMAJ'] = common_beam.major.value/3600
            hdrNB6['BMIN'] = common_beam.minor.value/3600
            hdrNB6['BPA'] = common_beam.pa.value
        else:
            hdrNB3['BMAJ'] = common_beam.major.value
            hdrNB3['BMIN'] = common_beam.minor.value
            hdrNB3['BPA'] = common_beam.pa.value
        fits.writeto(outdir+'/%s_B6_conv.fits'%label, conv_B6, hdrNB6, overwrite = True)

convolveb3b6(Path.w51e_b3_tt0,Path.w51e_b6_tt0,'/orange/adamginsburg/w51/TaehwaYoo/w51e_b6_imaging_2025/','w51e')
convolveb3b6(Path.w51n_b3_tt0,Path.w51n_b6_tt0,'/orange/adamginsburg/w51/TaehwaYoo/w51n_b6_imaging_2025/','w51n')