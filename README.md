# JWSTmorph

This package, working with python3, allows to calculate morphological parameters of galaxies, namely Gini, M20, concentration, asymmetry, shape asymmetry, smoothness, maximum radius (Rmax), and clumpiness. It is specifically adapted to be used on multiwavelength images from JWST-NIRCAM.
The code can be already used on a small sample of images simulated with MIRAGE (included in the package). For any problem contact me at antonello.calabro@inaf.it.


<strong>Python packages required :</strong>
```diff
matplotlib,skimage,pandas,astLib,astropy,photutils,statmorph,scipy
```

<strong>Instructions :</strong>
To run the code simply type:  <br/> 
```diff
python main_super_easy_from_catalog.py
```
<br/>

<strong>Input of the code </strong> (user defined input parameters, to be modified directly in the main code)  :
* <strong>which_subset :</strong> write your favourite label to distinguish your subset of galaxies (the name of the subset is included in the name of final table saved).
* <strong>segmap_detection_type :</strong> 'automatic' = to define the object, it takes automatically all segmentation regions in a central box whose size should be defined for each galaxy in the input catalog. If nothing is detected in this central box, an empty segmentation map is returned ; 'indices' = the user should write for each galaxy the corresponding index in the segmentation map (in the file /data/galaxy_segmap_indices.txt), and then perform a first run to check the right index and update that file. [N.B. currently working only with the automatic setup]
* <strong>use_statmorph :</strong> True = also calculates the morphological parameters with the *statmorph* package (https://statmorph.readthedocs.io/en/latest/) (Rodriguez-Gomez et al. 2019). ; False = the program skips *statmorph* (-9 is returned for all the output parameters starting with the prefix 'SM_' .
* <strong>segmentation_type :</strong> =  can be 'Pawlik' [DEFAULT] = the segmentation map is derived following Pawlik et al. (2016). Other possible choices are: 'petrosian_based', 'photutils', and 'square'. [currently fully tested with the default option only]

In the section 'INPUT CATALOG', the program looks for a catalog named *'assembled_catalog.cat'* in the data folder, which includes the list of galaxies to process, with their ID, redshift, and boxsize for the initial detection. Then, it looks for fits images in the data folder with the following structure name:  *'ID_'+str(IDgal)+'_'+band+'.fits'*, where *IDgal* is an integer, and *band* is a string specifying the photometric band (can be 'f090w', 'f115w', 'f150w', 'f200w', 'f277w', 'f356w', or 'f444w'). The number of bands, the pixel size of the images, and the FWHM resolution for each filter can be modified in the first part of the main code.
Image cutouts centered on each objects should be created with your favourite tool prior to running this program. In this package, example cutouts are given for all the 7 different bands listed above.

<br/><br/>

<strong>Output of the code </strong> (column name and brief description) :
* <strong>c1 :</strong> clumpiness parameter, derived as explained in Calabrò et al. (2019) https://doi.org/10.1051/0004-6361/201935778 . The pixel coordinates of the galactic nuclei to remove should be specified in the file file_nuclei.txt . No nucleus is removed if the table contains all zeros.
* <strong>Rmax :</strong> maximum radius of the galaxy R<sub>max</sub>, following the definition of Pawlik et al. (2016)
* <strong>gini :</strong> gini parameter, calculated following the definitions of Abraham et al. (2003) and Lotz et al. (2004)
* <strong>m20 :</strong> M<sub>20</sub> parameter, following the definition of Lotz et al. (2004)
* <strong>C :</strong> concentration parameter
* <strong>A1 :</strong> asymmetry parameter
* <strong>shapeasy :</strong> shape asymmetry, as defined in Pawlik et al. (2016)
* <strong>S :</strong> smoothness parameter, as in Conselice et al. (2003)
* <strong>SNpixel1 :</strong> average S/N per pixel of the galaxy
* <strong>SM_gini, SM_m20, SM_C, SM_A, SM_S, SM_SNpixel SM_Rpetro :</strong> Gini, M<sub>20</sub>, concentration, asymmetry, smoothness parameters, average S/N per pixel and petrosian radius calculated with *statmorph*.
* <strong>SM_flag :</strong> *statmorph* flag (0= good fit, 1= indicates a problem with the basic morphological measurements, such as a discontinuous Gini segmentation map).

<br/><br/><br/>

<strong>Example of output images for a set of simulated clump configurations as described in Treu, Calabrò, Castellano, et al. (2022) :</strong>
![ID_400_distclumps_0.4_mag_27.0.png](https://github.com/Anthony96/JWSTmorph/blob/2cf8bf4613d1b1c55bb182015f4fff3d77edee78/results/ID_400_distclumps_0.4_mag_27.0.png)
![400_clumpmap_distclumps_0.4_mag_27.0.png](https://github.com/Anthony96/JWSTmorph/blob/2cf8bf4613d1b1c55bb182015f4fff3d77edee78/results/400_clumpmap_distclumps_0.4_mag_27.0.png)

The latter figure is a composition of 6 simulated images for the same configuration in different filters (F090W, F200W, F277W, F444W), including also F090W and F277W bands, convolved respectively to the psf resolution of F200W and F444W. All the stamps are 2'' on a side. 
![image_cutouts_1_distclumps_0.2_mag_27.0](https://github.com/Anthony96/JWSTmorph/blob/87e2e9940e5fabee68811f2b4056932b56f04ff6/results/image_cutouts_1_distclumps_0.2_mag_27.0.png)
