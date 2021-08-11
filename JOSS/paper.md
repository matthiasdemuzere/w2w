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

The Python-based WUDAPT-to-WRF (W2W) package is developed in this context, and translates Local Climate Zone (LCZ) maps into urban canopy parameters readable by WRF, the community "Weather Research and Forecasting" model. It is the successor of the Fortran-based W2W package developed by @Brousse2016 and @Martilli2016, and provides a more simple, efficient and improved procedure to use LCZ information in WRF.   

# Statement of need
Since the pioneering work of @Brousse2016 and @Martilli2016, the level-0 WUDAPT information, the Local Climate Zone maps, have been used increasingly in WRF. We expect this trend to continue, because of two recent developments: 1) the creation of city-wide LCZ maps is now easier than ever with the online LCZ Generator [@Demuzere2021], and 2) as of spring 2021, the new version 4.3 of WRF [@Skamarock2021] is able to ingest 11 urban classses (corresponding to WUDAPT's LCZs) by default, whereas previous versions required manual WRF code changes by the user (see @Martilli2016, @Zonato2021a and @Zonato2021b for more information). Because of these developments, we decided to simultaneously built an improved, Python-based, WUDAPT-to-WRF (W2W) routine, to make the translation of LCZ-based parameters better and more simple. 

# Initial data requirements
In order to use the tool, two input files are required: 

1. A **geo_em.d0X.nc file** for the inner WRF model domain in which you would like to use the LCZ-based information. This file can be produced by WRF's geogrid.exe tool as part of the WRF Preprocessing System (WPS). ** @ANDREA: does a user needs to use specific settings here to create this file?? Please extend this section if needed.** 

2.  A **Local Climate Zone map** that is slightly bigger than the domain of the geo_em.d0X.nc file. There are a number of ways to obtain an LCZ map for your region of interest: 

   * Extract your domain from the continental-scale LCZ maps for Europe [@Demuzere2019] or the United States [@Demuzere2020] (see [here](https://www.wudapt.org/lcz-maps/) for more info).
   * Check if your region of interest is already covered by the many LCZ maps available in the [submission table](https://lcz-generator.rub.de/submissions) of the LCZ Generator.
   * Use the [LCZ Generator](https://lcz-generator.rub.de/) to make an LCZ map for your region of interest. 


# General workflow
The goal of the Python-based W2W tool is to obtain a WRF domain file (*geo_em.d0X.nc*) that contains the urban LCZ classes and their corresponding urban canopy parameters relevant for all urban parameterizations embedded in WRF: the single layer urban canopy model Noah/SLUCM (@Kusaka2001), the Building Environment Parameterization (BEP, @Martilli2002), and BEP+BEM (Building Energy Model, @Salamanca2010). To get to that point, the following three general steps are followed, which are partly inspired by the work of @Li2020:

* Step 1: 
* Step 2: 
* Step 3:


# Urban canopy parameter assignment
MAKE A TABLE WITH ALL PARAMETERS, including abbrevation, long name, unit, type, source, etc ...


Two pathways are followed when assigning the various urban canopy parameters to the Local Climate Zone Map (\autoref{fig:workflow}):

* Pathway 1: **Morphological** parameters are assigned directly to the high-resolution LCZ map, and are only afterwards aggregated to the lower-resolution WRF grid. In this way, the method produces a unique value of the different urban morphology parameters for each WRF grid cell. This was found to be more efficient in reproducing urban boundary layer features, especially in the outskirts of the city [@Zonato2020], and is in line with the [WUDAPT-to-COSMO](https://github.com/matthiasdemuzere/WUDAPT-to-COSMO) routine [@Varentsov2020]. 
* Pathway 2: In line with the former Fortran-based W2W procedure, **radiative and thermal parameters** are assigned to the modal LCZ class that is assigned to each WRF grid cell. 

![General workflow. The example maps are derived from the sample data for Zaragoza (Spain), available in the github repository \label{fig:workflow}](workflow.png){ width=100% }

As before, the LCZ-based urban canopy parameters generally follow the values provided by [Stewart and Oke (2012)](http://doi.org/10.1175/BAMS-D-11-00019.1) and [Stewart et al. (2014)](http://doi.org/10.1002/joc.3746).


# Potential use case


# Things to keep in mind (come up with better section title!!)
* best to use with BEP or BEP+BEM, because of the building heights / lowest model layer
* replace generic LCZ-based UCP values with site-specific ones when available
* Important to have good quality LCZ map, if not: garbage in, garbage out.





# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](wudapt_logo.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](wudapt_logo.png){ width=20% }

# Acknowledgements
We acknowledge contributions and support from Alberto Martilli, Alejandro Rodriguez Sanchez and Oscar Brousse.

# References