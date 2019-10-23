/*
 * mira2_plugin_sparco-multi.i
 *
 * Implement a more verstaile version of SPARCO plugin to mira2
 *
 */

//MIRA_PLUGDIR = "~/apps/src/mira2_plugin_sparco";


func mira_plugin_sparcomulti_init(nil) {
  /*DOCUMENT mira_plugin_sparco-multi_init(nil);

  Initialisation of the plugin from the parameters read in the
  command line.
  */

  inform, "Loading \"sparcomulti\" plugin...";

  SPARCO_options = _lst("\nSparco specific options",
    _lst("sparco_model", "star", "NAME", OPT_STRING_LIST,
         "Name of the SPARCO model used"),
    _lst("sparco_params", [], "VALUES", OPT_REAL_LIST,
         "Parameters used with SPARCO"),
    _lst("sparco_w0", [], "VALUE", OPT_REAL,
         "Central wavelength (in microns) for SPARCO"),
    _lst("sparco_image", [], "NAME", OPT_STRING,
         "Name of the fits file that are needed if image types specified"),
    _lst("sparco_spectrum", [], "NAME", OPT_STRING_LIST,
         "Type of spectrum that model has/have"),
    _lst("sparco_file", [], "NAME", OPT_STRING_LIST,
         "Name of the ascii file where the spectrum(a) is/are"),
    _lst("sparco_index", [], "VALUE", OPT_REAL_LIST,
         "Spectral index(es)"),
    _lst("sparco_temp", [], "VALUE", OPT_REAL_LIST,
         "Temperature(s) for blackbody"),
    _lst("sparco_xy", [], "VALUE", OPT_REAL_LIST,
         "x and y shifts (mas) of the models"),
    _lst("sparco_flux", [], "VALUE", OPT_REAL_LIST,
         "flux ratios of the models at w0")
         );

  plugin = mira_new_plugin(options = SPARCO_options,
                parse_options = parse_options,
                tweak_visibilities=tweak_visibilities,
                tweak_gradient=tweak_gradient,
                add_keywords=add_keywords,
                add_extensions=add_extensions,
                read_keywords=read_keywords
                 );

  inform, "Plugin \"sparcomulti\" ready.";

  return plugin;
}

func parse_options(plugin, opt)
/* DOCUMENT parse_options(plugin, opt)

  plugin is a hash_table with following arguments:

  plugin.model = name of the model in SPARCO. It can be:
                       "star", "binary", "UD".
  plugin.params = vector of necessary parameters for the
                         models
  plugin.w0 = central wavelenghts for computation of chromaticity

  plugin.image = a fits-file with the image that will be used as a
                        SPARCO model

*/

{
  local sparco, params, w0, image;

  inform, "parse_options() \"sparcomulti\" ...";

  h_set, opt, flags=opt.flags | MIRA_KEEP_WAVELENGTH;

  /* Read SPARCO settings */
  sparco = opt.sparco_model;
  params = opt.sparco_params;
  w0 = opt.sparco_w0;
  image = opt.sparco_image;
  type =   opt.sparco_type;
  file = opt.sparco_file;
  temp = opt.sparco_temp;
  index = opt.sparco_index;
  shift = opt.sparco_xy;
  flux = opt.sparco_flux;
  spectrum = opt.sparco_spectrum;

  if (!is_void(sparco)) {
    /* use upper case */
    sparco = strupper(sparco);

    if (!is_void(spectrum)) {
      spectrum = strupper(spectrum);
    }

    if (!is_void(w0)) {
      if (w0 > 0.1) { // assume it is micron
        w0 *= 1e-6;
      }
    }

    /* How many models? */
    nmods = numberof(sparco);
    nstar = numberof(where(sparco == "STAR"));
    nUD   = numberof(where(sparco == "UD"));
    nbg   = numberof(where(sparco == "BG"));
    nimage = numberof(where(sparco == "IMAGE"));

    if (nmods != nstar + nUD + nimage + nbg) {
      throw, "The list of names is not correct can only be one of: \"star\", \"UD\", \"bg\" or \"image\" ";
    }

    /* Test the right number of parameters given the models */
    if (numberof(image) != nimage) {
      throw, "The same number of image files need to be specified than the number of image models (n= %i )", nimage;
    }
    if (numberof(shift) == 0) {
      inform, "All the models are centered at (0,0).";
      shift = array(0.0, (nmods-nbg)*2);
    } else if (numberof(shift) != 2*(nmods-nbg) | numberof(shift) != 2*nmods) {
      throw, "Each model (except bg) should have a specified the relative shift in the x and y direction w.r.t. to the reconstructed image";
    } else if (numberof(shift) == 2*(nmods-nbg)) {
      shift2 = [];
      for (i=1; i<=nmods; ++i) {
        if (sparco(i) == "BG") {
          grow, shift2, [0,0];
        } else {
          idx = 2*i-1;
          grow, shift2, shift(idx:idx+1);
        }
      }
      shift = shift2;
    }

    if (numberof(flux) != nmods) {
      throw, "Each model should have its flux ratio. The flux ratio of the reconstructed image will be 1-f_models."
    }

    if (sum(flux) >= 1 | sum(flux) < 0) {
      throw, "The sum of the relative fluxes should be lower than one and positive (currently %f)", sum(flux);
    }
    
    if (numberof(spectrum) != nmods+1 ) {
      throw, "Each model and the reconstrcuted image should have a relative spectrum. It can be either \"pow\", \"BB\" or \"spectrum\" "
    }

    print("n params: ", numberof(params));
    print("nUD: ", nUD);

    if (numberof(params) != nUD ) {
      throw, "The number of parameters does not match the number of \"UD\" models";
    }

    p=0;
    for (i=1; i<=nmods; ++i) {
      if (flux(i) < 0) {
        throw, "fluxes should be positive (negative fluxes not implemented)";
      }
      if (sparco(i) == "UD") {
        p += 1;
        if (params(p) <= 0) {
          throw, "The UD diameter should be positive";
        }
        inform, "Model %i is a %s with a %f mas diameter, with a flux ratio of %f and a %s spectrum.", i, sparco(i), params(p), flux(i), spectrum(i+1) ;
      } else {
        inform, "Model %i is a %s, with a flux ratio of %f and a %s spectrum.", i, sparco(i), flux(i), spectrum(i+1) ;
      }
    }
  } else {
    throw, "no models specified...";
  }

  h_set, plugin, model=sparco,
                 params=params,
                 spectrum=spectrum,
                 flux=flux,
                 w0=w0,
                 image=image,
                 index=index,
                 temp=temp,
                 file=file,
                 shift=shift,
                 type=type,
                 nmods=nmods,
                 nimages=nimage,
                 nbg=nbg,
                 nstar=nstar,
                 nUD=nUD,
                 specinit=1n,
                 imageinit=1n,
                 visinit=1n;

  inform, "parse_options() \"sparcomulti\" done.";
}

func _get_spectrum(type, w, w0, index=, temp=, file=)
/* DOCUMENT _get_spectrum(type, w, w0, index=, temp=, file=);

  DO NOT USE OUTSIDE THE SPARCO PLUGIN

  Defines the spectrum of a component with respect of
  what is specified
  w    = the wavelengths at which the spectrum will be computed
         (meters)
  w0   = central wavelength (meters)
  type = "pow",  "BB",        "spectrum"
  index= powerlaw index
  temp = black body temperature (K)
  file = ascii file with the spectrum

*/
{
  local star;
  if (type == "POW") {
    star = (w/w0)^index;
  } else if (type == "BB") {
    star = _BB(w,temp)/_BB(w0,temp);
  } else {
    throw, "Spectrum not implemented yet";
  }

  return star
}

func readSpectrum(file, w)
/* DOCUMENT readSpectrum(master, file);

  Reads the ascii file and interpolates the spectrum to
  measured wavelengths

 */
{}

func readImage(file, master)
/* DOCUMENT readImage(master, file);

  Reads the fits file and reshape the image to fit the
  reconstructed one.

 */
{}

func tweak_visibilities (master, vis)
/* DOCUMENT tweak_complex_visibilities (master, vis);

  Takes the complex visibilities from the image and
  adds the selected sparco model to return complex
  visibilities of image + SPARCO model

*/
{
 plugin = mira_plugin(master);
 model = plugin.model;
 spectrum = plugin.spectrum;
 index = plugin.index;
 temp = plugin.temp;
 file = plugin.file;
 nmods = plugin.nmods;
 shift = plugin.shift;
 params = plugin.params;

 /* compute and store spectra */
 if (plugin.specinit) {
   w = mira_model_wave(master);
   nw = numberof(w);
   w0 = plugin.w0;
   spectra = array(double, nw, nmods+1);
   ii = it = is = 0;
   for (i=1; i<=plugin.nmods+1; ++i) {
     if (spectrum(i) == "POW") {
       ii +=1
       spectra(,i) = _get_spectrum(spectrum(i), w, w0, index=index(ii));
     } else if (spectrum(i) == "BB") {
       it +=1;
       spectra(,i) = _get_spectrum(spectrum(i), w, w0, temp=temp(it));
     } else {
      is +=1;
      spectra(,i) = readSpectrum(file(is), w);
    }
 }
 h_set, plugin, spectra = spectra, specinit = 0n;

 /* Read the images */
 if (plugin.nimages>0 & plugin.imageinit) {
   nimg = mira_image_size(master, 1);
   images = array(double, n, n, plugin.nimages);
   throw, "model made of image wait for Jacques to implement it..."
   h_set, plugin, images=images;
                  imageinit=0n;
   }
 }

 /* Compute the model complex visibilities */
 if (plugin.visinit) {
   w = mira_model_wave(master);
   visibilities = array(double, 2, numberof(w), nmods)
   ip = 0;
   for (i=1; i<=nmods; ++i) {
     is = (2*i)-1;
     if (model(i) == "STAR") {
       visibilities(,,i) = mira_sparco_star(master, shift(is:is+1));
     } else if (model(i) == "UD") {
       ip +=1;
       visibilities(,,i) = mira_sparco_UD(master, shift(is:is+1), params(ip));
     } else if (model(i) == "BG" ) {
       visibilities(,,i) = mira_sparco_bg(master);
     } else if (model(i) == "IMAGE") {
       visibilities(,,i) = mira_sparco_image(master, shift(is:is+1));
     } else {
       throw, " %s not implemented yet.", model(i)
     }
   }

   h_set, plugin, visibilities=visibilities,
                  visinit=0n;
 }

 vis = mira_sparco_vis(master, vis);

 return vis;
}

func tweak_gradient (master, grd)
/* DOCUMENT tweak_complex_gradient (master, grd);

  Takes the complex gradient on the image and
  normalise it to take into account SPARCO

  SEE ALSO: mira_sparco_star, mira_sparco_UD, mira_sparco_binary.
*/
{
  plugin = mira_plugin(master);
  /* Gradient modification for SPARCO */
  spectra = plugin.spectra;
  flux = plugin.flux;
  nmods = plugin.nmods;

  fimg0 = 1 - sum(flux);
  ftot = fimg = fimg0 * spectra(,1);
  for (i=1; i<=nmods; ++i) {
    ftot += flux(i) * spectra(,i+1);
  }

  grd_re = grd(1,..);
  grd_im = grd(2,..);

  grd_re *= ftot / fimg;
  grd_im *= ftot / fimg;

  grd = [grd_re, grd_im];
  grd = transpose(grd);

  return grd;
}


func mira_sparco_vis (master, vis)
/* DOCUMENT mira_sparco_vis (master, vis);

  Linearly adds all the requested models to the current
  reconstructed image.

*/
{
  plugin = mira_plugin(master);
  spectra = plugin.spectra;
  visibilities = plugin.visibilities;
  flux = plugin.flux;


  fimg0 = 1 - sum(flux);
  ftot = fimg = fimg0 * spectra(,1);

  vis_re = fimg * vis(1,..);
  vis_im = fimg * vis(2,..);

  for (i=1; i<=nmods; ++i) {
    vis_re += flux(i) * spectra(,i+1) * visibilities(1,,i);
    vis_im += flux(i) * spectra(,i+1) * visibilities(2,,i);
    ftot += flux(i) * spectra(,i+1);
  }

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  return vis;
}


func mira_sparco_star(master, shift)
  /* DOCUMENT mira_sparco_star(master, vis);

     Compute the total complex visibilities of a shifted point source

     SEE ALSO: mira_sparco_UD, mira_sparco_binary.
   */
{
  plugin = mira_plugin(master);
  w = mira_model_wave(master);
  u = mira_model_u(master) / w;
  v = mira_model_v(master) / w;
  xbin = shift(1) * MIRA_MAS;
  ybin = shift(2) * MIRA_MAS;

  vis_re = cos( -2*pi*(xbin*u + ybin*v) );
  vis_im = sin( -2*pi*(xbin*u + ybin*v) );

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  return vis;
}

func mira_sparco_bg(master)
  /* DOCUMENT mira_sparco_bg(master);

     Compute the total complex visibilities of a background

     SEE ALSO: mira_sparco_UD, mira_sparco_binary.
   */
{
  plugin = mira_plugin(master);
  w = mira_model_wave(master);

  vis_re = vis_im = array(0., dimsof(w));

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  return vis;
}


func mira_sparco_UD(master, shift, UD)
/* DOCUMENT mira_sparco_UD(master, vis);

   Compute the total complex visibilities of a Uniform Disk

   SEE ALSO: mira_sparco_star.
 */
{
  plugin = mira_plugin(master);
  w = mira_model_wave(master);
  u = mira_model_u(master)/w;
  v = mira_model_v(master)/w;
  xbin = shift(1) * MIRA_MAS;
  ybin = shift(2) * MIRA_MAS;
  B = abs(u,v);

  if (UD==0.) {
    V_UD = array(1., dimsof(u));
  } else {
    UD *= MIRA_MAS;
    V_UD = 2*bessj1(pi * B * UD) / ( pi * B * UD);
  }
  vis_re = V_UD * cos( -2*pi*(xbin*u + ybin*v) );
  vis_im = sin( -2*pi*(xbin*u + ybin*v) );

  vis = [vis_re, vis_im];
  vis = transpose(vis);

  return vis;
}


func mira_sparco_image(master, vis)
  /* DOCUMENT mira_sparco_imageBB(master, vis);

     Compute the total complex visibilities by adding a predefined image (im0)
     the reconstructed image using the im0-to-total flux ratio (fim0) and
     the Black Body temperatures of the two images (T0, Timg)
     fim0 = fim0 * BB(T0, lambda) / BB(T0, lam0)
     fd = (1-fim0) * BB(Timg, lambda) / BB(Timg, lam0)
     Vtot = fd*Vimg + fim0*Vimg0
     Vtot /= fd +fim0

     SEE ALSO: mira_sparco_star, mira_sparco_binary.
   */
{
  local vis, vis_re, vis_im, vis_amp, vis_phi, fs0, denv, B;

  vis_re = vis(1,..);
  vis_im = vis(2,..);

  plugin = mira_plugin(master);
  fim0 = plugin.params(1);
  T0 = plugin.params(2);
  Tim = plugin.params(3);
  w = mira_model_wave(master);
  w0 = plugin.w0;
  img0 = plugin.image;
  BB = []; //TODO Find a function for blackbody


  fim = fim0 * plugin.star_spectrum;
  fd = (1-fim0) * _BB(Tim, w) / _BB(Tim, w0);
  ftot = fim + fd;

  Vimg0 = master.xform(img0);
  Vimg0_re = Vimg0(1,..);
  Vimg0_im = Vimg0(2,..);

  vis_re = vis_re * fd + fim * Vim0_re;
  vis_im = vis_im * fd + fim * Vim0_im;

  vis_re /= ftot;
  vis_im /= ftot;

  vis = [vis_re, vis_im];
  vis = transpose(vis);

//  h_set, master, model_vis_re = vis_re, model_vis_im = vis_im, model_vis = vis;

  return vis;
}




func _BB(lambda, T)
    /* DOCUMENT _BB(lambda, T);

       DESCRIPTION
       Computatiof black body radiation function with a given temperature T
       at a specific wavelength lambda.
       (Inspired from yocoAstroBBodyLambda)

       PARAMETERS
       - lambda : wavelength (m)
       - T      : temperature (K)

       RETURN VALUES
       Return the energy radiated (W/m2/m)
    */
{
    local c, h, kb;
    c = 2.99792453e8;
    h = 6.626070040e-34;
    kb = 1.38064852e-23;
    mask = abs(h*c / (kb*T*lambda) ) > 700;
    if(numberof(where(mask))==0)
        flambda = 2*h*c^2 / lambda^5 / (exp(h*c / (kb*T*lambda)) - 1);
    else
    {
        flambda = lambda;
        flambda(where(mask)) = 0.0;
        flambda(where(!mask)) = 2*h*c^2 / lambda^5 / (exp(h*c / (kb*T*lambda)) - 1);
    }

    return flambda;
}

func add_keywords (master, fh)
/* DOCUMENT add_keywords (master, fh);

    This function adds the sparco plugin add_keywords
    to the finale saved fits file.

    SEE ALSO: add_extension.
*/
{

  plugin = mira_plugin(master);
  nmods = plugin.nmods;
  fits_set, fh, "SWAVE0",  plugin.w0,  "SPARCO: Central wavelength (m) for chromatism";
  fits_set, fh, "SNMODS",  nmods,  "SPARCO: Number of models used in SPARCO";
  ip = ii = it = 0;
  im0 = 1-sum(plugin.flux);
  fits_set, fh, "SFLU0",  im0,  "SPARCO: Flux ratio of the image";
  fits_set, fh, "SPEC0", plugin.spectrum(1),  "SPARCO: spectrum associated with the reconstructed image";
  if (plugin.spectrum(1)=="POW") {
    ii+=1;
    fits_set, fh, "SIDX0", plugin.index(ii),  "SPARCO: spectral index of the model";
  } else if (plugin.spectrum(1)=="BB") {
    it+=1;
    fits_set, fh, "STEM0", plugin.temp(it),  "SPARCO: black body temperature of the model";
  }
  for (i=1; i<=nmods; ++i) {
    modelid = swrite(format="SMOD%d", i);
    fluxid = swrite(format="SFLU%d", i);
    specid = swrite(format="SPEC%d", i);
    indexid = swrite(format="SIDX%d", i);
    tempid = swrite(format="STEM%d", i);
    xid = swrite(format="SDEX%d", i);
    yid = swrite(format="SDEY%d", i);
    is = 2*i-1;
    modelid = swrite(format="SMOD%d", i);
    fits_set, fh, modelid,  plugin.model(i),  "SPARCO: Model used";
    fits_set, fh, fluxid,  plugin.flux(i),  "SPARCO: Flux ratio of the model";
    if (plugin.model(i)=="UD") {
      ip += 1;
      paramid = swrite(format="SPAR%d", i);
      fits_set, fh, paramid,  plugin.params(ip),  "SPARCO: UD diameter (mas)";
    }
    fits_set, fh, specid, plugin.spectrum(i+1),  "SPARCO: spectrum associated with model";
    if (strupper(plugin.spectrum(i+1))=="POW") {
      ii += 1;
      fits_set, fh, indexid, plugin.index(ii),  "SPARCO: spectral index of the model";
    } else if (plugin.spectrum(i+1)=="BB") {
      it+=1
      fits_set, fh, tempid, plugin.temp(it),  "SPARCO: black body temperature of the model";
    }
    fits_set, fh, xid,  plugin.shift(is),  "SPARCO: RA shift of the model (mas)";
    fits_set, fh, yid,  plugin.shift(is+1),  "SPARCO: DEC shift of the model (mas)";
  }
//  fits_set, fh, "SWAVE0",  plugin.w0,  "Central wavelength (mum) for chromatism";
}

func add_extension (master, fh)
/* DOCUMENT add_extension (master, fh);

    This function adds the sparco plugin add_keywords
    to the finale saved fits file.

    SEE ALSO: add_keywords.
*/
{
  plugin = mira_plugin(master);


 //FIXME: add stellar spectrum if used by sparco
 // fits_new_hdu, fh, "IMAGE", "SPARCO adds an extension";
 // fits_write_array, fh, random(128,128);
 // fits_pad_hdu, fh;

}


func read_keywords (tab, fh)
/* DOCUMENT read_keywords (master, fh);

    This function reads sparco plugin keywords
    to the finale saved fits file.
    SEE ALSO: write_keyword.
*/
{
  inform, "read_keywords() \"sparcomulti\" ...";

  nmods     = mira_get_fits_integer(fh, "SNMODS");
  sparco_w0 = mira_get_fits_real(   fh, "SWAVE0");

  if (is_void(nmods)) {
    nmods = 0;
  }
  for (i=0; i<=nmods; ++i) {
    sparco_spectrum_i =  mira_get_fits_string(fh, swrite(format="SPEC%d", i));
    sparco_spectrum = _(sparco_spectrum,sparco_spectrum_i);

    if (sparco_spectrum_i =="POW") {
      sparco_index = _(sparco_index, mira_get_fits_real(fh, swrite(format="SIDX%d", i)));
    } else if (sparco_spectrum_i  =="BB") {
      sparco_temp = _(sparco_temp, mira_get_fits_real(fh, swrite(format="STEM%d", i)));
    }

    if(i>0){
      sparco_flux = _(sparco_flux, mira_get_fits_real(fh, swrite(format="SFLU%d", i)));
      sparco_model_i = mira_get_fits_string(fh, swrite(format="SMOD%d", i));
      if (sparco_model_i=="UD") {
        sparco_params = _(sparco_params, mira_get_fits_real(fh, swrite(format="SPAR%d", i)));
      }
      sparco_model= _(sparco_model,sparco_model_i);
      sparco_xy= _(sparco_xy, mira_get_fits_real(fh, swrite(format="SDEX%d", i)), mira_get_fits_real(fh, swrite(format="SDEY%d", i)));
    }
  }

  print("sparco_w0 :", sparco_w0);
  print("sparco_flux :", sparco_flux);
  print("sparco_model :", sparco_model);
  print("sparco_xy :", sparco_xy);
  print("sparco_spectrum :", sparco_spectrum);
  print("sparco_index :", sparco_index);
  print("sparco_temp :", sparco_temp);
  print("sparco_params :", sparco_params);

  /* Set SPARCO settings */
  h_set, tab, sparco_w0=sparco_w0,
    sparco_flux=sparco_flux,
    sparco_model=sparco_model,
    sparco_xy=sparco_xy,
    sparco_spectrum=sparco_spectrum,
    sparco_index=sparco_index,
    sparco_temp=sparco_temp,
    sparco_params=sparco_params ;

  inform, "read_keywords() \"sparcomulti\" done.";

  return tab;
}
