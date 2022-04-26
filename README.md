# bbgs_elmer
Elmer/Ice code for the Bering-Bagley Glacier System Model

Provided by Thomas Trantow as required by paper submission to Computer & Geosciecnes and the Journal of Geophysical Research Earth Surfaces. See the relevant papers titled "Sensitivity of Glacier Elevation Analysis and Numerical Modeling to CryoSat-2 SIRAL Retracking Techniques" (2020) and "Evolution of a Surge Cycle of the Bering-Bagley Glacier System from Observations and Numerical Modeling" (2022) respectively.

For 2020 C&G paper:
An example Solver Input File (SIF) to run the BBGS model in Elmer/Ice is provided (crev\_BBGS\_C2\_swath\_IAMG\_20181126.sif) along with code the BBGS-specific User Functions (USF\_Bering.f90). Code was written by Thomas Trantow and adapted from open source code provided by the Elmer/Ice community (see http://elmerice.elmerfem.org/courses-tutorials).

For 2022 JGR paper:
All relevant code for this project is found in the folder titled elmer_code_JGR. Privided is the SIF file for the 20-year quiescent phase experiment (full_BBGS_s11_veit_swath_bedv7_accumlars_20211005.sif), the SIF file for the surge-wave experiment (crev_BBGS_surgewave_paper_20220128.sif), the input surface DEM (s11_VHswath_wfront_dem4model.dat), the input bedrock DEM (bed_JPL_v7_fourptsall.dat) and the Bering Glacier specific user functions (USF_Bering.f90) which includes definition of the surge-wave friction law introduced in the paper that requiers the definition of the flowline (flowline_BBGS.dat).
