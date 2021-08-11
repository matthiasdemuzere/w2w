---
title: "W2W: A Python package that makes WUDAPT's Local Climate Zone information available in WRF"
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
An important objective of WUDAPT, the World Urban Database and Acces Portals Tools community project, is to 1) to acquire and make accessible coherent and consistent information on form and function of urban morphology relevant to climate weather, and environment studies on a worldwide basis, and 2) to provide tools that extract relevant urban parameters and properties for models and for model applications at appropriate scales for various climate, weather, environment, and urban planning purposes [@Ching2018]. 

This python WUDAPT-to-WRF (W2W) package is developed in this context, and translates Local Climate Zone (LCZ) maps into urban canopy parameters readable by WRF, the community "Weather Research and Forecasting" model. It is the successor of the Fortran-based W2W package developed by @Brousse2016 and @Martilli2016, and provides a more simple, efficient and improved procedure to use LCZ information in WRF.   

# Statement of need
Since the pioneering work of @Brousse2016 and Martilli2016, the level-0 WUDAPT information, the Local Climate Zone maps, have been used increasingly in WRF. We expect this trend to continue, because of two recent developments: 1) the creation of city-wide LCZ maps is now easier than ever with the online LCZ Generator [@Demuzere2021], and 2) as of spring 2021, the new version 4.3 of WRF [@Skamarock2021] is able to ingest 11 urban classses (corresponding to WUDAPT's LCZs) by default, whereas previous versions required manual WRF code changes by the user (see @Martilli2016, @Zonato2021a and @Zonato2021b for more information). Because of these developments, we decided to simultaneously built an improved WUDAPT-to-WRF (W2W) routine, to make the translation of LCZ-based parameters better and more simple. 

# Workflow
As before, the LCZ-based urban canopy parameters generally follow the values provided by [Stewart and Oke (2012)](http://doi.org/10.1175/BAMS-D-11-00019.1) and [Stewart et al. (2014)](http://doi.org/10.1002/joc.3746).

The procedure in this new `w2w` tool is different from the former tool. Morphological parameters are assigned directly to the high-resolution LCZ map, and only afterwards aggregated to the WRF grid. In this way, the method produces a unique value of the different urban morphology parameters for each model cell. This was found to be more efficient in reproducing urban boundary layer features, especially in the outskirts of the city [(Zonato et al., 2020)](https://doi.org/10.1016/j.uclim.2020.100584), and is in line with the [WUDAPT-to-COSMO](https://github.com/matthiasdemuzere/WUDAPT-to-COSMO) routine [(Varentsov et al., 2020)](https://www.mdpi.com/2073-4433/11/12/1349). Other radiative and thermal parameters are for now still assigned to the modal LCZ class. More details on the procedure and its assumptions will soon be available [here](https://github.com/matthiasdemuzere/w2w#how-to-cite).  

Yet because of its complex workflow, and  complex procedure, and Note that this tool was first assigning the LCZ mode to each WRF grid cell, and only afterwards assigning corresponding morphological, radiative and thermal properties to this modal LCZ class. This is done differently in w2w, see below. 

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

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