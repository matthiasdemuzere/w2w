---
title: "W2W: A Python package that injects WUDAPT's Local Climate Zone information in WRF"
tags:
  - Python
  - WUDAPT
  - WRF
  - Local Climate Zones
  - Urban canopy parameters
  - Climate modelling
authors:
    - name: Matthias Demuzere^[corresponding author]
      orcid: 0000-0003-3237-4077
      affiliation: 1 #
    - name: Daniel Argüeso
      orcid: 0000-0002-4792-162X
      affiliation: 2
    - name: Andrea Zonato
      orcid: 0000-0002-9174-1618
      affiliation: 3
affiliations:
    - name: Urban Climatology Group, Department of Geography, Ruhr-University Bochum, Bochum, Germany
      index: 1
    - name: Physics Department, University of the Balearic Islands, Palma, Spain
      index: 2
    - name: Atmospheric Physics Group, Department of Civil, Environmental and Mechanical Engineering, University of Trento, Trento, Italy
      index: 3
date: 10 August 2021
bibliography: paper.bib
---

# Summary
An important objective of WUDAPT, the World Urban Database and Acces Portals Tools community project, is 1) to acquire and make accessible coherent and consistent information on form and function of urban morphology relevant to climate weather, and environment studies, and 2) to provide tools that extract relevant urban parameters and properties for models and model applications at appropriate scales for various climate, weather, environment, and urban planning purposes [@Ching2018].

The Python-based WUDAPT-to-WRF (`W2W`) package is developed in this context, and translates Local Climate Zone (LCZ) maps into urban canopy parameters readable by WRF, the community "Weather Research and Forecasting" model [@Skamarock2021]. It is the successor of the Fortran-based `W2W` package developed by @Brousse2016 and @Martilli2016, and provides an improved, more simple, and more efficient procedure to use LCZ information in WRF. Some important changes include a direct manipulation of the geogrid files (without the creation of temporary files), and the use of average LCZ-based urban morphological parameters instead of assigning them to the modal LCZ class.

# Statement of need
Since the pioneering work of @Brousse2016 and @Martilli2016, the level-0 WUDAPT information, the Local Climate Zone maps, have been used increasingly in WRF.

We expect this trend to continue, because of two recent developments: 1) the creation of city-wide LCZ maps is now easier than ever with the launch of the LCZ Generator web application [@Demuzere2021], and 2) WRF versions > 4.3 [@Skamarock2021] are able to ingest 10 or 11 built classses (corresponding to WUDAPT's LCZs) by default, whereas previous WRF versions required manual code changes (see @Martilli2016, @Zonato2021a and @Zonato2021b for more information).

Because of these developments, an improved, Python-based, WUDAPT-to-WRF (`W2W`) routine is presented here, so as to make the translation of LCZ-based parameters better and more simple.

# Initial data requirements
In order to use the tool, two input files are required:

1. A **geo_em.d0X** (.nc) file for the inner WRF model domain in which one would like to use the LCZ-based information. This file can be produced by WRF's geogrid.exe tool as part of the WRF Preprocessing System (WPS), without additional modifications of the standard procedure.

2.  A **Local Climate Zone map** (.tif) file that is slightly bigger than the domain extent of the geo_em.d0X.nc file. There are a number of ways to obtain an LCZ map for your region of interest (ROI):

   * Extract your ROI from the continental-scale LCZ maps for Europe [@Demuzere2019] or the United States [@Demuzere2020] (see [here](https://www.wudapt.org/lcz-maps/) for more info).
   * Check if your ROI is already covered by the many LCZ maps available in the [submission table](https://lcz-generator.rub.de/submissions) of the LCZ Generator.
   * Use the [LCZ Generator](https://lcz-generator.rub.de/) to make an LCZ map for your ROI. See also [here](https://www.wudapt.org/create-lcz-classification/) for more information.


# Workflow
The goal of the Python-based `W2W` tool is to obtain a WRF domain file (*geo_em.d0X.nc*) that contains the built LCZ classes and their corresponding urban canopy parameters relevant for all urban parameterizations embedded in WRF: the single layer urban canopy model (Noah/SLUCM, @Kusaka2001), the Building Environment Parameterization (BEP, @Martilli2002), and BEP+BEM (Building Energy Model, @Salamanca2010).

To get to that point, a number of sequential steps are required:

* _Step 1: Remove the default urban land cover_

The default urban land cover from MODIS is replaced with the dominant surrounding vegetation category, as done in @Li2020. This procedure affects WRF's parameters LU_INDEX, LANDUSEF and GREENFRAC. LU_INDEX is selected as the dominant category from the $nlus$ (default = 45) nearest grid points (excluding ocean, urban and lakes). GREENFRAC is calculated as the mean over all grid points with that dominant vegetation category among the $nlus$ nearest points. For each grid point, if LANDUSEF had any percentage of urban, it is set to zero and the percentage is added to the dominant vegetation category assigned to that grid point.

Resulting output: **geo_em.d0X_NoUrban.nc**

* _Step 2: Define the LCZ-based urban extent_

LCZ-based impervious fraction values (FRC_URB2D, available from `LCZ_UCP_default.csv`) are assigned to the original 100 m resolution LCZ map, and are aggregated to the WRF resolution. Areas with FRC_URB2D < 0.2 ($frc$) are currently considered non-urban. This choice has been made to avoid the use of the urban schemes in areas where the majority of the landuse is vegetated, since the impact of the impervious surfaces is low. The FRC_URB2D field is also used to mask all other urban parameter fields, so that their extent is consistent.

Resulting output: **geo_em.d0X_LCZ_extent.nc**

* _Step 3: Introduce modal built LCZ classes_

For each WRF grid cell, the mode of the underlying built LCZ classes is added to LU_INDEX (numbered from 31-41). See [here](https://ral.ucar.edu/sites/default/files/public/product-tool/urban-canopy-model/WRF_urban_update_Readme_file_WRF4.3.pdf) for more info. Note that the `W2W` routine by default considers LCZ classes 1-10 as built classes ($bc$). In some cases, also LCZ E (or 15 - Bare rock or paved) can be considered as a built LCZ class, as it might reflect large asphalt surfaces such as big parking lots or airstrips. In that case, the user must make sure the $bc$ argument is set appropriately.

* _Step 4: Assign urban canopy parameters_

Two procedures are followed when assigning the various urban canopy parameters to the LCZ map and translating this information onto WRF's grid:

**Procedure 1**: **Morphological parameters** are assigned directly to the high-resolution LCZ map, and are afterwards aggregated to the lower-resolution WRF grid. As a result, the method produces a unique urban morphology parameter value for each WRF grid cell. This was found to be more efficient in reproducing urban boundary layer features, especially in the outskirts of the city [@Zonato2020], and is in line with the [WUDAPT-to-COSMO](https://github.com/matthiasdemuzere/WUDAPT-to-COSMO) routine [@Varentsov2020].

Morphological urban canopy parameter values are provided in `LCZ_UCP_default.csv`, and are generally based on values provided in @Stewart2012 and @Stewart2014. In addition:

* Building width (BW), available in `LCZ_UCP_default.csv`, is taken from `URBPARM_LCZ.TBL` (stored in WRF's run/ folder).
* While `URBPARM_LCZ.TBL` also has values on street width, `W2W` derives street width from the mean building height (MH_URB2D) and the Height-to-Width ratio (H2W), to have these fields consistent.
* Plan (LP_URB2D), frontal (LF_URB2D) and total (LB_URB2D) area indices are based on formulas in @Zonato2020.
* HI_URB2D is obtained by fitting a bounded normal distribution to the minimum (MH_URB2D_MIN), mean (MH_URB2D), and maximum (MH_URB2D_MAX) building height, as provided in `LCZ_UCP_default.csv`. The building height standard deviation is also required, and is approximated as (MH_URB2D_MAX - MH_URB2D_MIN) / 4.
* For computational efficiency, HI_URB2D values lower than 5% were set to 0 after resampling, the remaining HI_URB2D percentages are re-scaled to 100%.


**Procedure 2**: In line with the former Fortran-based `W2W` procedure, **radiative and thermal parameters** are assigned to the modal LCZ class that is assigned to each WRF grid cell (see _Step 3_). These parameter values are not stored in the netcdf output, but are read from `URBPARM_LCZ.TBL` and assigned automatically to the modal LCZ class when running the model.


* _Step 5: Adjust global attributes_

In a final step, some global attributes are adjusted in the resulting netcdf files:

* NBUI_MAX is added as a global attribute, reflecting the maximum amount of HI_URB2D classes that are not 0 across the model domain. This paramater can be used when compiling WRF, to optimize memory storage.
* NUM_LAND_CAT is set to 41, to reflect the addition of 10 (or 11) built LCZ classes. This is not only done for the highest resolution domain file (e.g. d04), but also for **all of its lower-resolution parent domain files (e.g. d01, d02, d03)**. As such, make sure these files are also available in the input data directory.

Resulting output: **geo_em.d0X_LCZ_params.nc**


# Integration in WRF's preprocessing
 The current tool is designed to work with the geo_em.d0X files produced by geogrid.exe, which is available in the WRF Preprocessing System (WPS). WPS needs to be at a version >3.8, in order to incorporate the urban geometrical parameters in the `URB_PARAM` matrix [@Glotfelty2013]. The user should run geogrid.exe using its default settings, which will provide the various geo_em.d0X.nc files containing the static data fields. No additional variables are required, neither in the namelist.wps nor within the GEOGRID.TBL table. The `W2W` tool (\autoref{fig:w2w_workflow}) reads the standard geo_em.d0X.nc files (for all the domains) and produces the aforementioned **geo_em.d0X_LCZ_params.nc** files. The user should then simply rename these files to the standard name for each of the domains (e.g. geo_em.d04_LCZ_params.nc to geo_em.d04.nc), which will serve as input to the metgrid.exe module (\autoref{fig:w2w_workflow}).

![Modified workflow to set-up and run a WRF simulations including urban parameters derived from LCZs using W2W.\label{fig:w2w_workflow}](w2w_workflow.jpg)



# Potential use cases
The files provided as output by `W2W` allow a wide range of applications, including - but not limited to - addressing the impact of:

* urbanization, by running WRF with the default geo_em.d0X.nc and the geo_em.d0X_NoUrban.nc files (see for example @Li2020 and @Hirsch2021).
* an improved urban land cover extent description, by running WRF with the default geo_em.d0X.nc and the geo_em.d0X_LCZ_extent.nc files (similar to for example @Bhati2018 and @Mallard2018).
* a more detailed (LCZ-based) urban description, by running WRF with the default geo_em.d0X.nc and the geo_em.d0X_LCZ_params.nc files (see for example @Brousse2016, @Hammerberg2018, @Molnar2019, @Wong2019, @Patel2020, @Zonato2020, @Ribeiro2021, @Hirsch2021 and @Patel2021).


# Important notes
* The LCZ-based urban canopy parameter values provided in `LCZ_UCP_default.csv` and `URBPARM_LCZ.TBL` are universal and generic, and might not be appropriate for your ROI. If available, please adjust the urban canopy parameters values according to the characteristics of your ROI.
* It is advised to use this tool with urban parameterization options BEP or BEP+BEM (`sf_urban_physics = 2 or 3`, respectively). In case you use this tool with the SLUCM model (`sf_urban_physics = 1`), make sure your lowest model level is above the highest building height. If not, real.exe will provide the following error message: `ZDC + Z0C + 2m is larger than the 1st WRF level - Stop in subroutine urban - change ZDC and Z0C`.
* It is advised to use WRF versions > 4.3, that are able to ingest 10 or 11 built classses (corresponding to WUDAPT's LCZs) by default [@Skamarock2021], and WPS versions > 3.8, in order to incorporate the urban geometrical parameters in the `URB_PARAM` matrix [@Glotfelty2013].  


# Acknowledgements
We acknowledge contributions and support from Alberto Martilli, Alejandro Rodriguez Sanchez and Oscar Brousse.

# References
