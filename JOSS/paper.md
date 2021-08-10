---
title: 'W2W: A Python package that makes WUDAPT's Local Climate Zone information available in WRF'
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

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
# Test compilation here: https://whedon.theoj.org/

---

# Summary

An important objective of the WUDAPT (World Urban Database and Acces Portals Tools) community project, is to generate urban canopy information and provide the (open-source) tools to facilitate urban-focused modelling studies CHING2018. As a part of that, BROUSSE2016 pioneered the integration of Local Climate Zone maps and their corresponding urban canopy parameters into WRF, the community "Weather Research and Forecasting" model. 

As of spring 2021, WRF v4.3.x is able to ingest LCZ information by default (previous versions required manual WRF code changes by the user). See more details on "Updates of WRF-urban in WRF 4.3: Local Climate Zones, Mitigation Strategies, building materials permeability and new buildings drag coefficient" here. Because of this, we decided to simultaneously built an improved WUDAPT-to-WRF routine, to make the translation of LCZ-based parameters better and more simple.

# Statement of need
Their original guide and code on how to use WUDAPT information into WRF (originally designed for WRF v3.2) is available here. Note that this tool was first assigning the LCZ mode to each WRF grid cell, and only afterwards assigning corresponding morphological, radiative and thermal properties to this modal LCZ class. This is done differently in w2w, see below.


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
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements
We acknowledge contributions and support from Alberto Martilli, Alejandro Rodriguez Sanchez and Oscar Brousse.

# References