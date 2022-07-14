# JWSTmorph

This package, working with python3, allows to calculate morphological parameters of galaxies, namely Gini, M20, concentration, asymmetry, shape asymmetry, smoothness, maximum radius, and clumpiness. It is specifically designed to be applied on new JSWT NIRCAM images.
Even though it is working (it can be applied on a small sample of images simulated with MIRAGE and included in the package), it is under constant improvement. For any problem contact me at antonello.calabro@inaf.it.

<strong>Python packages required :</strong>
matplotlib,skimage,pandas,astLib,astropy,photutils,statmorph,scipy

<strong>Instructions :</strong>
To run the code simply type: <i> <font size="10"> python main_super_easy.py </i> </font> <br/>
make_final_plots_results.py is for post-processing <br/><br/><br/>

<strong>Output of the code </strong> (column name and brief description)  :
* <strong>c1 :</strong> clumpiness parameter, derived as explained in Calabrò et al. (2019) https://doi.org/10.1051/0004-6361/201935778 . No nucleus is removed in the calculation.
* <strong>Rmax :</strong> maximum radius of the galaxy R<sub>max</sub>, following the definition of Pawlik et al. (2016)
* <strong>gini :</strong> gini parameter, calculated following the definitions of Abraham et al. (2003) and Lotz et al. (2004)
* <strong>m20 :</strong> M<sub>20</sub> parameter, following the definition of Lotz et al. (2004)
* <strong>C :</strong> concentration parameter
* <strong>A1 :</strong> asymmetry parameter
* <strong>shapeasy :</strong> shape asymmetry, as defined in Pawlik et al. (2016)
* <strong>S :</strong> smoothness parameter, as in Conselice et al. (2003)
* <strong>SNpixel2 :</strong> S/N per pixel
* <strong>SM_gini, SM_m20, SM_C, SM_A, SM_S, SM_SNpixel SM_Rpetro :</strong> Gini, M<sub>20</sub>, concentration, asymmetry, smoothness parameters, S/N per pixel and petrosian radius calculated with *statmorph* (https://statmorph.readthedocs.io/en/latest/) (Rodriguez-Gomez et al. 2019).
* <strong>SM_flag :</strong> *statmorph* flag (0= good fit, 1= indicates a problem with the basic morphological measurements, such as a discontinuous Gini segmentation map).

<br/><br/><br/>

<strong>Example of output images (Galaxy ID 104, band F090W) :</strong>
![ID_104_band_f090w_ok](https://user-images.githubusercontent.com/12728781/178595197-63f5971b-51e8-43c5-960f-f278e7245674.png)
![104_clumpmap_band_f090w](https://user-images.githubusercontent.com/12728781/178595024-9d1030fe-c87a-4061-befe-66da95611c0a.png)

