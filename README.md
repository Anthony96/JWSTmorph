# JWSTmorph

This package, working with python3, allows to calculate morphological parameters of galaxies, namely Gini, M20, concentration, asymmetry, shape asymmetry, smoothness, maximum radius, and clumpiness. It is specifically designed to be applied on new JSWT NIRCAM images.
Even though it is working (it can be applied on a small sample of images simulated with MIRAGE and included in the package), it is under constant improvement. For any problem contact me at antonello.calabro@inaf.it.


<strong>Python packages required :</strong>
```diff
matplotlib,skimage,pandas,astLib,astropy,photutils,statmorph,scipy
```

<strong>Instructions :</strong>
To run the code simply type:  <br/> 
```diff
python main_super_easy.py
```
<br/>

<strong>Input of the code </strong> (user defined input parameters, to be modified directly in the main code)  :
* <strong>segmap_detection_type :</strong> 'automatic' = automatically chooses the central segmentation region of the segmentation map to define the object (if nothing is detected in the central 10 pixel size box, an empty map is returned) ; 'indices' = the user should write for each galaxy the corresponding index in the segmentation map (in the file /data/galaxy_segmap_indices.txt), and then perform a first run to check the right index and update that file.
* <strong>use_statmorph :</strong> True = also calculates the morphological parameters with the *statmorph* package (https://statmorph.readthedocs.io/en/latest/) (Rodriguez-Gomez et al. 2019). ; False = the program skips *statmorph* (-9 is returned for all the parameters starting with the prefix 'SM_' .
* <strong>which_subset :</strong> write your favourite label to distinguish your subset of galaxies (the name of the subset is included in the name of final table saved).
* <strong>segmentation_type :</strong> =  can be 'Pawlik' [DEFAULT] = the segmentation map is derived following Pawlik et al. (2016). Other possible choices are: 'petrosian_based', 'photutils', and 'square'.

In the section 'INPUT CATALOG', the program looks for a catalog named *'assembled_catalog.cat'* in the data folder, which includes the list of galaxy IDs to process. Then it searches for fits images in the data folder with the following structure name:  *'ID_'+str(IDgal)+'_'+band+'.fits'*, where *IDgal* is an integer, and *band* is a string specifying the photometric band (can be 'f090w', 'f115w', 'f150w', 'f200w', 'f277w', 'f356w', or 'f444w').
Image cutouts centered on each objects should be created with your favourite tool prior to running this program. In the example, image cutouts of 3'' box size are given for all the 7 different bands listed above.

<br/><br/>

<strong>Output of the code </strong> (column name and brief description) :
* <strong>c1 :</strong> clumpiness parameter, derived as explained in Calabrò et al. (2019) https://doi.org/10.1051/0004-6361/201935778 . No nucleus is removed in the calculation.
* <strong>Rmax :</strong> maximum radius of the galaxy R<sub>max</sub>, following the definition of Pawlik et al. (2016)
* <strong>gini :</strong> gini parameter, calculated following the definitions of Abraham et al. (2003) and Lotz et al. (2004)
* <strong>m20 :</strong> M<sub>20</sub> parameter, following the definition of Lotz et al. (2004)
* <strong>C :</strong> concentration parameter
* <strong>A1 :</strong> asymmetry parameter
* <strong>shapeasy :</strong> shape asymmetry, as defined in Pawlik et al. (2016)
* <strong>S :</strong> smoothness parameter, as in Conselice et al. (2003)
* <strong>SNpixel2 :</strong> S/N per pixel
* <strong>SM_gini, SM_m20, SM_C, SM_A, SM_S, SM_SNpixel SM_Rpetro :</strong> Gini, M<sub>20</sub>, concentration, asymmetry, smoothness parameters, S/N per pixel and petrosian radius calculated with *statmorph*.
* <strong>SM_flag :</strong> *statmorph* flag (0= good fit, 1= indicates a problem with the basic morphological measurements, such as a discontinuous Gini segmentation map).

<br/><br/><br/>

<strong>Example of output images (Galaxy ID 104, band F090W) :</strong>
![ID_104_band_f090w_ok](https://user-images.githubusercontent.com/12728781/178595197-63f5971b-51e8-43c5-960f-f278e7245674.png)
![104_clumpmap_band_f090w](https://user-images.githubusercontent.com/12728781/178595024-9d1030fe-c87a-4061-befe-66da95611c0a.png)

