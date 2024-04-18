from astropy.convolution import convolve
from astropy.modeling import functional_models
from astropy.convolution import Gaussian1DKernel
from astropy import constants as consts
from scipy.interpolate import interp1d

"""
in the two following functions we find the disperion tables used for each object/spectrum
we can then download them by code or from the website:
https://hst-crds.stsci.edu/context_table/hst_1069.pmap
under the tab 'disptab'

in this website you can find many details about COS observations
note that the disptab files are not based on the LP of observation. for example the same
dispersion table could be use for objectA observed in LP1 and objectB in LP3

the dispersion coefficient is different depending on the CENWAVE and grating (settings), like 
the numbers that we download their corresponding files in HST website. However after the reduction
the delta-lambda in COS-1dsum spectra or my spectra is same through out the spectrum  
"""
cwd = Path(".")
datadir = cwd / "data"
outputdir = cwd / "output"
plotsdir = cwd / "output" / "plots"

# Make the directories if they don't exist
datadir.mkdir(exist_ok=True), outputdir.mkdir(exist_ok=True), plotsdir.mkdir(exist_ok=True)

fuvHeader0 = fits.getheader('../lcx220rxq_x1d.fits', ext=0)

param_dict = {}  # Make a dict to store what you find here

for hdrKeyword in [
    "DETECTOR",
    "OPT_ELEM",
    "LIFE_ADJ",
    "CENWAVE",
    "DISPTAB",
]:  # Print out the relevant values
    try:  # For DISPTAB
        value = fuvHeader0[hdrKeyword].split("$")[1]  # Save the key/value pairs to the dictionary
        param_dict[hdrKeyword] = value                # DISPTAB needs the split here
    except:  # For other params
        value = fuvHeader0[hdrKeyword]  # Save the key/value pairs to the dictionary
        param_dict[hdrKeyword] = value
    print(f"{hdrKeyword} = {value}") 


def fetch_files(det, grating, lpPos, cenwave, disptab):
    """
    Given all the inputs: (detector, grating, LP-POS, cenwave, dispersion table,) this will download both
    the LSF file and Disptab file you should use in the convolution and return their paths.
    Returns:
    LSF_file_name (str): filename of the new downloaded LSF file
    disptab_path (str): path to the new downloaded disptab file
    """
    COS_site_rootname = (
        "https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/cos/"
        "performance/spectral-resolution/_documents/"
    )  # Link to where all the files live - split into 2 lines
    if det == "NUV":  # Only one file for NUV
        LSF_file_name = "nuv_model_lsf.dat"
    elif det == "FUV":  # FUV files follow a naming pattern
        LSF_file_name = f"aa_LSFTable_{grating}_{cenwave}_LP{lpPos}_cn.dat"

    LSF_file_webpath = COS_site_rootname + LSF_file_name  # Where to find file online
    urllib.request.urlretrieve(
        LSF_file_webpath, str(datadir / LSF_file_name)
    )  # Where to save file to locally
    print(f"Downloaded LSF file to {str(datadir/ LSF_file_name)}")

    # And we'll need to get the DISPTAB file as well
    disptab_path = str(datadir / disptab)
    urllib.request.urlretrieve(
        f"https://hst-crds.stsci.edu/unchecked_get/references/hst/{disptab}",
        disptab_path,
    )
    print(f"Downloaded DISPTAB file to {disptab_path}")

    return LSF_file_name, disptab_path



def read_lsf(filename):
    # This is the table of all the LSFs: called "lsf"
    # The first column is a list of the wavelengths corresponding to the line profile, so we set our header accordingly
    if "nuv_" in filename:  # If its an NUV file, header starts 1 line later
        ftype = "nuv"

    else:  # assume its an FUV file
        ftype = "fuv"
    hs = 0
    lsf = Table.read(filename, format="ascii", header_start=hs)

    # This is the range of each LSF in pixels (for FUV from -160 to +160, inclusive)
    # middle pixel of the lsf is considered zero ; center is relative zero
    pix = np.arange(len(lsf)) - len(lsf) // 2  # integer division to yield whole pixels

    # the column names returned as integers.
    lsf_wvlns = np.array([int(float(k)) for k in lsf.keys()])

    return lsf, pix, lsf_wvlns

# lsf, pix, lsf_wvlns = read_lsf('aa_LSFTable_G130M_1291_LP3_cn.dat')


def get_disp_params(disptab, cenwave, segment, x=[]):
    """
    Helper function to redefine_lsf(). Reads through a DISPTAB file and gives relevant\
    dispersion relationship/wavelength solution over input pixels.
    Parameters:
    disptab (str): Path to your DISPTAB file.
    cenwave (str): Cenwave for calculation of dispersion relationship.
    segment (str): FUVA or FUVB?
    x (list): Range in pixels over which to calculate wvln with dispersion relationship (optional).
    Returns:
    disp_coeff (list): Coefficients of the relevant polynomial dispersion relationship
    wavelength (list; if applicable): Wavelengths corresponding to input x pixels 
    """
    with fits.open(disptab) as d:
        wh_disp = np.where(
            (d[1].data["cenwave"] == cenwave)
            & (d[1].data["segment"] == segment)
            & (d[1].data["aperture"] == "PSA")
        )[0]
        disp_coeff = d[1].data[wh_disp]["COEFF"][0] # 0 is needed as this returns nested list [[arr]]
        d_tv03 = d[1].data[wh_disp]["D_TV03"]  # Offset from WCA to PSA in Thermal Vac. 2003 data
        d_orbit = d[1].data[wh_disp]["D"]  # Current offset from WCA to PSA

    delta_d = d_tv03 - d_orbit

    if len(x):  # If given a pixel range, build up a polynomial wvln solution pix -> λ
        wavelength = np.polyval(p=disp_coeff[::-1], x=np.arange(16384))
        return disp_coeff, wavelength
    else:  # If x is empty:
        return disp_coeff


def redefine_lsf(lsf_file, cenwave, disptab, detector="FUV"):
    """
    Helper function to convolve_lsf(). Converts the LSF kernels in the LSF file from a fn(pixel) -> fn(λ)\
    which can then be used by convolve_lsf() and re-bins the kernels.
    Parameters:
    lsf_file (str): path to your LSF file
    cenwave (str): Cenwave for calculation of dispersion relationship
    disptab (str): path to your DISPTAB file
    detector (str): FUV or NUV?
    Returns:
    new_lsf (numpy.ndarray): Remapped LSF kernels.
    new_w (numpy.ndarray): New LSF kernel's LSF wavelengths.
    step (float): first order coefficient of the FUVA dispersion relationship; proxy for Δλ/Δpixel.
    """

    if detector == "FUV":
        xfull = np.arange(16384)

        # Read in the dispersion relationship here for the segments
        ### FUVA is simple
        disp_coeff_a, wavelength_a = get_disp_params(disptab, cenwave, "FUVA", x=xfull)
        ### FUVB isn't taken for cenwave 1105, nor 800:
        if (cenwave != 1105) & (cenwave != 800):
            disp_coeff_b, wavelength_b = get_disp_params(
                disptab, cenwave, "FUVB", x=xfull)
        elif cenwave == 1105:
            # 1105 doesn't have an FUVB so set it to something arbitrary and clearly not real:
            wavelength_b = [-99.0, 0.0]

        # Get the step size info from the FUVA 1st order dispersion coefficient
        step = disp_coeff_a[1]

        # Read in the lsf file
        lsf, pix, w = read_lsf(lsf_file)

        # take median spacing between original LSF kernels
        deltaw = np.median(np.diff(w))

        lsf_array = [np.array(lsf[key]) for key in lsf.keys()]
        if (deltaw < len(pix) * step * 2):  # resamples if the spacing of the original LSF wvlns is too narrow
            # this is all a set up of the bins we want to use
            # The wvln difference between kernels of the new LSF should be about twice their width
            new_deltaw = round(len(pix) * step * 2.0)  
            new_nw = (int(round((max(w) - min(w)) / new_deltaw)) + 1)  # nw = number of LSF wavelengths
            new_w = min(w) + np.arange(new_nw) * new_deltaw  # new version of lsf_wvlns

            # populating the lsf with the proper bins
            new_lsf = np.zeros((len(pix), new_nw))  # empty 2-D array to populate
            for i, current_w in enumerate(new_w):
                dist = abs(current_w - w)  # Find closest original LSF wavelength to new LSF wavelength
                lsf_index = np.argmin(dist)
                orig_lsf_wvln_key = lsf.keys()[lsf_index]  # column name corresponding to closest orig LSF wvln
                new_lsf[:, i] = np.array(lsf[orig_lsf_wvln_key])  # assign new LSF wvln the kernel of the closest original lsf wvln
        else:
            new_lsf = lsf
            new_w = w
        return new_lsf, new_w, step

    elif detector == "NUV":
        xfull = np.arange(1024)
        # Read in the dispersion relationship here for the segments
        disp_coeff_a, wavelength_a = get_disp_params(disptab, cenwave, "NUVA", x=xfull)
        disp_coeff_b, wavelength_b = get_disp_params(disptab, cenwave, "NUVB", x=xfull)
        disp_coeff_c, wavelength_c = get_disp_params(disptab, cenwave, "NUVC", x=xfull)

        # Get the step size info from the NUVB 1st order dispersion coefficient
        step = disp_coeff_b[1]

        # Read in the lsf file
        lsf, pix, w = read_lsf(lsf_file)

        # take median spacing between original LSF kernels
        deltaw = np.median(np.diff(w))

        lsf_array = [np.array(lsf[key]) for key in lsf.keys()]

        # this section is a set up of the new bins we want to use:
        new_deltaw = round(len(pix) * step * 2.0)  # The wvln difference between kernels of the new LSF should be about twice their width
        new_nw = (int(round((max(w) - min(w)) / new_deltaw)) + 1)  # nw = number of LSF wavelengths
        new_w = min(w) + np.arange(new_nw) * new_deltaw  # new version of lsf_wvlns

        # populating the lsf with the proper bins
        new_lsf = np.zeros((len(pix), new_nw))  # empty 2-D array to populate
        for i, current_w in enumerate(new_w):
            dist = abs(current_w - w)  # Find closest original LSF wavelength to new LSF wavelength
            lsf_index = np.argmin(dist)
            orig_lsf_wvln_key = lsf.keys()[lsf_index]  # column name corresponding to closest orig LSF wvln
            new_lsf[:, i] = np.array(lsf[orig_lsf_wvln_key])  # assign new LSF wvln the kernel of the closest original lsf wvln
        return new_lsf, new_w, step

