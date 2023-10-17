# Activating-Amplifiers

This github contains sample codes for the following manuscript: "Design and mathematical analysis of activating transcriptional amplifiers that enable modular temporal control in synthetic juxtacrine circuits"
https://doi.org/10.1016/j.synbio.2023.09.008

This code runs on CompuCell3D v3.7.8 which can be found at the organization's website, along with instructions. 
Instructions can also be found at the page for the original model: https://github.com/lmorsut/Lam_Morsut_GJSM

The repository here contains example codes for the above manuscript. 

I have included the file: 
NP_TempTemp which is the code for the Ncad Pcad structure with temporarily amplified gray cells (Ncad) and temporarily amplified orange cells (Pcad).
Parameters can be changed to model the other tested structures in this study by following the changes made in the Supplementary Table provided with the manuscript.
For example, to obtain the E-cadherin based structures change the adhesion to the E-cadherin parameters and the circuit equations as desired.
I have also included a "Multipole Demo" file that has unecessary comments and unused parameters removed as a template.

I have also included the code for the synthetic immunotherapy experiments with the heterogenous tumor with pink, red, and orange cells.
The description of these cells can be found in the manuscript.
The file is SbTempAmpHML. All the other experiments can be derived from the code by changing the tested parameter.
A template flie is provided as well: "SBCART Demo"

Please feel free to contact me for questions or help: calvin.lam.k@gmail.com
