# mira-sparco-multi
More verstaile mira-sparco plugin

SPARCO (Semi-Parametric Approach for Reconstruction of Chromatic Objects) is
an approach to reconstruct chromatic images from optical interferometric data
by adding a geometrical model to the reconstructed image.
This plugin is designed to work with [MiRA2](https://github.com/emmt/MiRA) directly from the command line.

## Installation

Firs make sure [MiRA2](https://github.com/emmt/MiRA) is installed.

Either install it by copying it to MIRA_HOME or use the option
`-plugin=/the/path/to/the/plugin/mira2_plugin_sparcomulti.i`


## Usage

add the following options to mira2:

`-sparco_model=`      a list the model(s) you want to use (See Sect. Implemented Models), for example: `-sparco_model=star,star,UD,bg`

`-sparco_flux=`  a list of flux ratios for the model(s). Must be smaller than 1.

`-sparco_w0=`         central wavelengths where the flux ratios are effective

`-sparco_spectrum=`  Type of spectral behaviour for **the reconstructed image and model(s)**. Can be either: `pow` for a power law, `BB` for a blackbody or `spectrum` for a spectrum specified in an ascii file (TO BE IMPLEMENTED).

`-sparco_index=` if at least one `sparco_spectrum` is a `pow` then you have to specify the spectral index list here

`-sparco_temp=` if at least one `-sparco_spectrum=BB` then the specify the black body temperature here

`-sparco_params=` if at least one of your `sparco_type` is a UD, then specify here the UD diameter

`-sparco_xy=` Specify here the (x,y)-shifts (in mas) of the different models w.r.t. the reconstructed image (for ex. `-sparco_xy=0,0,10,-5` for x1=0, y1=0, x2=10, y2=-5).

`-sparco_file` if at leats one `-sparco_spectrum=spectrum`, specify here the path and na;e to the ASCII files defining the spectrum of your objects (TO BE IMPLEMENTED).

`-sparco_image` Specify here the fits file have the image of the models (TO BE IMPLEMENTED).

## Implemented Models

* **star:** It adds a point source.

* **UD:** It adds a Uniform Disk.

Needs `-sparco_params=UD` where UD is the disk diameter in mas.

* **bg**  It adds an over-resolved component (V=0).

* **image**  It adds the image computed from specified fits-file (TO BE IMPLEMENTED)

Needs `-sparco_image=file` where file is the fits-file of the image

## Example

* If you want to reconstruct a ring (at 1500K) around an evolved binary where the primary (flux ratio at 1.65mum is 60%) is marginally resolved (UD=0.5mas) and the secondary has a temperature of 2000K, has a 3% of flux ratio at 1.65mum and is shifted by d(ra)=-0.7mas and d(dec)=-0.6mas, then write:

`ymira -plugin=sparcomulti -your_mira_parameters -sparco_model=UD,star -sparco_w0=1.65 -sparco_flux=0.6,0.03 -sparco_spectrum=BB,pow,BB -sparco_index=-4 -sparco_temp=1500,4000 -sparco_params=0.5 -sparco_xy=0,0,-0.7,-0.6 -initial=Dirac your_data.fits outputfile.fits`

where `-your-mira-parameters` could be `-pixelsize=0.5mas -fov=80mas -regul=hyperbolic -mu=3e2 -tau=5e-5 -ftol=0 -gtol=0 -maxeval=1000 -overwrite -save_visibilities -save_initial -normalization=1 -min=0 -verb=10`

## Contact

If you have problems or you want to implement a specific model to be used with SPARCO,
please drop me an email at: jacques.kluska@kuleuven.be

## References

* Kluska, J. et al.: *"SPARCO : a semi-parametric approach for image reconstruction of chromatic objects. Application to young stellar objects"* in Astronomy & Astrophysics, Volume 564, id.A80, 11 pp. [DOI](https://ui.adsabs.harvard.edu/link_gateway/2014A&A...564A..80K/doi:10.1051/0004-6361/201322926)

## Examples of astrophysical results using SPARCO

* Kluska, J. et al.: *"A disk asymmetry in motion around the B[e] star MWC158"* in
    Astronomy & Astrophysics, Volume 591, id.A82, 15 pp. [DOI](https://ui.adsabs.harvard.edu/link_gateway/2016A&A...591A..82K/doi:10.1051/0004-6361/201527924)

* Hillen, M., Kluska, J. et al.: *"Imaging the dust sublimation front of a circumbinary disk"* in Astronomy & Astrophysics, Volume 588, id.L1, 6 pp. [DOI](https://ui.adsabs.harvard.edu/link_gateway/2016A&A...588L...1H/doi:10.1051/0004-6361/201628125)
