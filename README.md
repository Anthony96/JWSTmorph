# JWSTmorph

This package, working with python3, allows to calculate morphological parameters of galaxies, namely Gini, M20, concentration, asymmetry, shape asymmetry, smoothness, maximum radius, and clumpiness. It is specifically designed to be applied on new JSWT NIRCAM images.
Even though it is working (it can be applied on a small sample of images simulated with MIRAGE and included in the package), it is under constant improvement. For any problem contact me at antonello.calabro@inaf.it.

<strong>Python packages required :</strong>
matplotlib,skimage,pandas,astLib,astropy,photutils,statmorph,scipy

<strong>Instructions :</strong>
To run the code simply type: <i> <font size="10"> python main.py </i> </font>


make_final_plots_results.py is for post-processing

<strong>Example of output images (Galaxy ID 104, band F090W) :</strong>
![ID_104_band_f090w_ok](https://user-images.githubusercontent.com/12728781/178595197-63f5971b-51e8-43c5-960f-f278e7245674.png)
![104_clumpmap_band_f090w](https://user-images.githubusercontent.com/12728781/178595024-9d1030fe-c87a-4061-befe-66da95611c0a.png)

