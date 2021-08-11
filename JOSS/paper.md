---
title: "W2W: A Python package that injects WUDAPT's Local Climate Zone information in WRF"
tags:
  - Python
  - Local Climate Zones
  - WRF
  - WUDAPT
  - Urban canopy parameters
  - Climate modelling 
authors:
  - name: Matthias Demuzere^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-3237-4077
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Daniel Argüeso # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-4792-162X
    affiliation: 2
  - name: Andrea Zonato
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
An important objective of WUDAPT, the World Urban Database and Acces Portals Tools community project, is to 1) to acquire and make accessible coherent and consistent information on form and function of urban morphology relevant to climate weather, and environment studies on a worldwide basis, and 2) to provide tools that extract relevant urban parameters and properties for models and for model applications at appropriate scales for various climate, weather, environment, and urban planning purposes [@Ching2018, @Ching2019]. 

The Python-based WUDAPT-to-WRF (`W2W`) package is developed in this context, and translates Local Climate Zone (LCZ) maps into urban canopy parameters readable by WRF, the community "Weather Research and Forecasting" model. It is the successor of the Fortran-based `W2W` package developed by @Brousse2016 and @Martilli2016, and provides a more simple, efficient and improved procedure to use LCZ information in WRF.   

# Statement of need
Since the pioneering work of @Brousse2016 and @Martilli2016, the level-0 WUDAPT information, the Local Climate Zone maps, have been used increasingly in WRF. We expect this trend to continue, because of two recent developments: 1) the creation of city-wide LCZ maps is now easier than ever with the online LCZ Generator [@Demuzere2021], and 2) as of spring 2021, the new version 4.3 of WRF [@Skamarock2021] is able to ingest 10 or 11 built classses (corresponding to WUDAPT's LCZs) by default, whereas previous versions required manual WRF code changes by the user (see @Martilli2016, @Zonato2021a and @Zonato2021b for more information). Because of these developments, we decided to simultaneously built an improved, Python-based, WUDAPT-to-WRF (`W2W`) routine, to make the translation of LCZ-based parameters better and more simple. 

# Initial data requirements
In order to use the tool, two input files are required: 

1. A **geo_em.d0X** (.nc) file for the inner WRF model domain in which you would like to use the LCZ-based information. This file can be produced by WRF's geogrid.exe tool as part of the WRF Preprocessing System (WPS). ** @ANDREA: does a user needs to use specific settings here to create this file?? Please extend this section if needed.** 

2.  A **Local Climate Zone map** (.tif) file that is slightly bigger than the domain of the geo_em.d0X.nc file. There are a number of ways to obtain an LCZ map for your region of interest: 

   * Extract your domain from the continental-scale LCZ maps for Europe [@Demuzere2019] or the United States [@Demuzere2020] (see [here](https://www.wudapt.org/lcz-maps/) for more info).
   * Check if your region of interest is already covered by the many LCZ maps available in the [submission table](https://lcz-generator.rub.de/submissions) of the LCZ Generator.
   * Use the [LCZ Generator](https://lcz-generator.rub.de/) to make an LCZ map for your region of interest. See also [here](https://www.wudapt.org/create-lcz-classification/) for more information.


# Workflow
The goal of the Python-based `W2W` tool is to obtain a WRF domain file (*geo_em.d0X.nc*) that contains the built LCZ classes and their corresponding urban canopy parameters (see TABLE XX) relevant for all urban parameterizations embedded in WRF: the single layer urban canopy model Noah/SLUCM (@Kusaka2001), the Building Environment Parameterization (BEP, @Martilli2002), and BEP+BEM (Building Energy Model, @Salamanca2010). 

MAKE A TABLE WITH ALL PARAMETERS, including abbrevation, long name, unit, type, source, etc ...

To get to that point, a number of sequential steps are followed:

* _Step 1: Remove the default urban land cover_

The default urban land cover from MODIS is replaced with the dominant surrounding vegetation category, as is done in @Li2020. This procedure affects WRF's variables LU_INDEX (land use index), LANDUSEF (land use fraction) and GREENFRAC (vegetation fraction). LU_INDEX is selected as the dominant category from the $nlus$ (default = 45) nearest grid points (excluding ocean, urban and lakes). LANDUSEF and GREENFRAC are calculated as the mean over all grid points with that category among the $nlus$ nearest points. @DANIEL: CORRECT??

Resulting output: **geo_em.d0X_NoUrban.nc**

* _Step 2: Define the LCZ-based urban extent_

LCZ-based impervious fraction (FRC_URB2D) values are assigned to the original 100 m resolution LCZ map, and are aggregated to the WRF resolution. Areas with FRC_URB2D < .2 ($frc$) are currently considered non-urban @ANDREA - ADD SMALL SENTENCE TO STATE WHY THAT IS. The FRC_URB2D field is also used to mask all other urban fields, so that they are consistent.

Resulting output: **geo_em.d0X_LCZ_extent.nc**

* _Step 3: Introduce modal built LCZ classes_

For each WRF grid cell, the mode of the underlying built LCZ classes is added to LU_INDEX, numbered from 31-41. See [here](https://ral.ucar.edu/sites/default/files/public/product-tool/urban-canopy-model/WRF_urban_update_Readme_file_WRF4.3.pdf) for more info. Note that the `W2W` routine by default considers LCZ classes 1-10 as built classes ($bc$). Sometimes, also LCZ E (or 15 - Bare rock or paved) can be considered as a built LCZ classes, as it might reflect large asphalt surfaces such as big parking lots or airstrips. In that case, make sure to set argument $bc$ appropriately.

* Step 4: Assign urban canopy parameters

Two pathways are followed when assigning the various urban canopy parameters to the Local Climate Zone Map:

  * Pathway 1: **Morphological** parameters are assigned directly to the high-resolution LCZ map, and are afterwards aggregated to the lower-resolution WRF grid. In this way, the method produces a unique value of the different urban morphology parameters for each WRF grid cell. This was found to be more efficient in reproducing urban boundary layer features, especially in the outskirts of the city [@Zonato2020], and is in line with the [WUDAPT-to-COSMO](https://github.com/matthiasdemuzere/WUDAPT-to-COSMO) routine [@Varentsov2020].

    Morphological urban canopy parameter values are provided in `LCZ_UCP_default.csv` (available in the github repository), and are generally based on @Stewart2012 and @Stewart2014. Building width (BW) is taken from URBPARM_LCZ.TBL (available in WRF's run/ folder). And while URBPARM_LCZ.TBL also has values on street width, `W2W` derives street width from the mean building height (MH_URB2D) and the Height-to-Width ratio (H2W), to have these fields consistent.

    In addition:
    * Plan (LP_URB2D), frontal (LF_URB2D) and total (LB_URB2D) area indices are based on formulas in @Zonato2020.
    * HI_URB2D is obtained by fitting a bounded normal distribution to the minimum, mean, maximum and standard deviation of the building heights (BH), as provided in `LCZ_UCP_default.csv`. In addition, for computational efficiency, values lower than 5% were set to 0 after resampling, the remaining HI_URB2D percentages are re-scaled to 100%.

  * Pathway 2: In line with the former Fortran-based `W2W` procedure, **radiative and thermal parameters** are assigned to the modal LCZ class that is assigned to each WRF grid cell. These parameter values are not stored in the netcdf output, but are read from URBPARM_LCZ.TBL and assigned automatically to the modal LCZ class when running the model. 

Resulting output: **geo_em.d0X_LCZ_params.nc**

# Run the tool
With respect to the WRF pre-processing chain?


# Potential use cases


# Things to keep in mind (come up with better section title!!)
* NBUI_MAX is added as a global attribute in the netcdf file, reflecting the maximum amount of HI_URB2D classes that are not 0 across the model domain. This paramater can be used during compilation to optimize memory storage.
* best to use with BEP or BEP+BEM, because of the building heights / lowest model layer
* replace generic LCZ-based UCP values with site-specific ones when available
* Important to have good quality LCZ map, if not: garbage in, garbage out.
* netcdf4/hdf5 compilation?

# Acknowledgements
We acknowledge contributions and support from Alberto Martilli, Alejandro Rodriguez Sanchez and Oscar Brousse.

# References