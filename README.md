w2w.py
======
A python tool that translates and ingests [WUDAPT's](www.wudapt.org) Local Climate Zone information into [WRF](https://github.com/wrf-model/WRF).


## Context
To add.

## Python environment
Create a virtual python environment, for example with anaconda: 
1. If you not have anaconda on your system yet, please install it first (information [here](https://docs.conda.io/en/latest/miniconda.html))
2. Then, create the virtual environment from within this repository: `conda env create -f w2w.yml`
3. Activate this new environment: `conda activate w2w`

## Requirements
- A geo_em.d0**X**.nc file (produced by WRF's WPS geoegrid.exe), for your inner domain you would like to use with LCZ information. 
  
**Note**: also the geo_em.d0[**0 to X**].nc files of the parent domain(s) should be available in the same INPUT_DIRECTORY. As the `w2w.py` routine will make sure `NUM_LAND_CAT` is set to 41 in all these files.

- A Local Climate Zone map (lcz.tif) that is slightly bigger than the domain of the geo_em.d0**X**.nc file. There are a number of ways to obtain an LCZ map for your region of interest:
  * Extract your domain from the continental-scale LCZ maps for Europe ([Demuzere et al., 2019](https://doi.org/10.1371/journal.pone.0214474)) or the United States ([Demuzere et al., 2020](https://doi.org/10.1038/s41597-020-00605-z)). For more information, see [here](https://www.wudapt.org/lcz-maps/). Make sure you use the version _epsg4326.tif.
  * Check if your region of interest is already covered by the many LCZ maps available in the [LCZ Generator submission table](https://lcz-generator.rub.de/submissions).
  * Use the [LCZ Generator](https://lcz-generator.rub.de/) to make an LCZ map for your region of interest. 

**Important**: In case the geo_em.d0**X**.nc domain is larger than ~ 2.5 x 2.5°, the LCZ Generator will fail. In that case, please contact [Matthias Demuzere](mailto:matthias.demuzere@rub.de) for support. 

- the lcz .tif file and geo_em*.d0X.nc should live in the same directory. Also, this directory should be writeable by the user.

## Run the tool
```
python w2w.py INPUT_DIRECTORY lcz.tif geo_em.d0**X**.nc       
```

Try with the provided sample:
```
python w2w.py ./sample_data lcz.tif geo_em.d0**X**.nc       
```

More information:
```python
python w2w.py -h 
```

## License
The project is licensed under the MIT license.


## How to cite
Not available yet. 


