w2w.py
======
A WUDAPT-to-WRF python tool that translates [WUDAPT's](www.wudapt.org) Local Climate Zone information into [WRF](https://github.com/wrf-model/WRF).


Context
-------
An	important	objective	of	WUDAPT, the World Urban Database and Acces Portals Tools community project, is to generate urban canopy information and provide the (open-source) tools to facilitate urban-focused modelling studies ([Ching et al., 2018](http://journals.ametsoc.org/doi/10.1175/BAMS-D-16-0236.1)). 

Since the work of [Brousse et al. (2016)](10.1016/j.uclim.2016.04.001), the level-0 WUDAPT information, the Local Climate Zone maps, have been used increasingly in [WRF](https://github.com/wrf-model/WRF), the community “Weather Research and Forecasting” model. Their original guide and code on how to use WUDAPT information into WRF (originally designed for WRF v3.2) is available [here](https://www.wudapt.org/wudapt-to-wrf/).

As of spring 2021, [WRF v4.3.x](https://github.com/wrf-model/WRF/releases/tag/v4.3) is able to ingest LCZ information by default (previous versions required manual WRF code changes by the user). See more details on "*Updates of WRF-urban in WRF 4.3: Local Climate Zones, Mitigation Strategies, building materials permeability and new buildings drag coefficient*" [here](https://ral.ucar.edu/sites/default/files/public/product-tool/urban-canopy-model/WRF_urban_update_Readme_file_WRF4.3.pdf). Because of this, we decided to simultaneously built an improved WUDAPT-to-WRF routine, to make the translation of LCZ-based parameters better and more simple. As before, the LCZ-based urban canopy parameters generally follow the values provided by [Stewart and Oke (2012)](http://doi.org/10.1175/BAMS-D-11-00019.1) and [Stewart et al. (2014)](http://doi.org/10.1002/joc.3746).


Python environment
-------
Create a virtual python environment, for example with anaconda: 
1. If you do not have anaconda on your system yet, please install it first (information [here](https://docs.conda.io/en/latest/miniconda.html))
2. Clone the w2w repository in your folder of choice, and enter it:
```sh   
> git clone https://github.com/matthiasdemuzere/w2w.git
> cd w2w 
```

3. Then, create the virtual environment from within this repository:
```sh   
> conda env create -f w2w.yml
```
4. Activate this new environment:
```sh   
> conda activate w2w
```

Requirements
-------
1. A **geo_em.d0X.nc file** (produced by WRF's WPS geoegrid.exe), for the inner WRF model domain in which you would like to use the LCZ-based information. 


2. A **Local Climate Zone map** (lcz.tif) that is slightly bigger than the domain of the geo_em.d0**X**.nc file. There are a number of ways to obtain an LCZ map for your region of interest:
   
   * Extract your domain from the continental-scale LCZ maps for Europe ([Demuzere et al., 2019](https://doi.org/10.1371/journal.pone.0214474)) or the United States ([Demuzere et al., 2020](https://doi.org/10.1038/s41597-020-00605-z)). For more information, see [here](https://www.wudapt.org/lcz-maps/). Make sure you use the version ending with *_epsg4326.tif*.
   * Check if your region of interest is already covered by the many LCZ maps available in the [LCZ Generator submission table](https://lcz-generator.rub.de/submissions).
   * Use the [LCZ Generator](https://lcz-generator.rub.de/) to make an LCZ map for your region of interest. In case the geo_em.d0**X**.nc domain is larger than ~ 2.5 x 2.5°, the LCZ Generator will fail. In that case, please contact [Matthias Demuzere](mailto:matthias.demuzere@rub.de) for support. 

**Important notes**: 
* Also the geo_em.d0[**0 to X**].nc file(s) of the parent domain(s) should be available in the same INPUT_DIRECTORY. This is needed because the `w2w.py` routine will check whether `NUM_LAND_CAT` is set to 41 in all these parent domain files. If that is not the case, this will be fixed.
* Your LCZ .tif and geo_em*.d0X.nc files should live in the same directory.
* Also, this directory should be writeable by the user.
* In case you use an LCZ map produced by the LCZ Generator, make sure to use `-lcz_band 1`. As such, the best-quality gaussian filtered LCZ map will be used in the process (see [Demuzere et al. (2021)](https://doi.org/10.3389/fenvs.2021.637455) for more info).

Run the tool
-------

1. Check out its help:
```sh
python w2w.py -h 
```
2. Try with the provided sample:
```sh
python w2w.py ./sample_data lcz_zaragoza.tif geo_em.d04.nc       
```

3. Deploy using your own data:
```sh
python w2w.py INPUT_DIRECTORY YOUR_LCZ.TIF YOUR_GEO_EM.d0X.NC       
```


License
-------
The project is licensed under the MIT license.


How to cite
-------
Not available yet. 

Contributing
-------

Contributions to `w2w` are welcome! This is how:

- **Bugs:** If you find a bug, please report it by opening an issue. if possible, please attach the error, code, version, and other details. 

- **Fixing Issues:** If you want to contributte by fixing an issue, please check the issues: contributions are welcome for open issues with labels :code:`bug` and :code:`help wanted`.

- **Enhancement:** New features and modules are welcome! You can check the issues: contributions are welcome for open issues with labels :code:`enhancement` and :code:`help wanted`.

Artists
-------
* [Matthias Demuzere](https://github.com/matthiasdemuzere): Lead developer
* [Daniel Argueso](https://github.com/dargueso): WRF expert 
* [Andrea Zonato](https://github.com/andreazonato): WUDAPT-WRF expert and liaison

Credits
-------
* We appreciate the feedback and suggestions provided by [Alberto Martilli](https://github.com/albertomartilli) and [Oscar Brousse](https://github.com/oscarbrousse), lead developers of the original Wudapt-to-wrf fortran package. 
* Thanks to [Alberto Martilli](https://github.com/albertomartilli) and Alejandro Rodriguez Sanchez for allowing us to use of their Zaragoza case-study files as sample data. 

