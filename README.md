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
- A Local Climate Zone map (lcz.tif) that is slightly bigger than the domain of the geo_em.d0**X**.nc file. The [LCZ Generator](https://lcz-generator.rub.de/) can be used for this purpose. 

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


