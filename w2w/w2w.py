import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import pandas as pd
import rioxarray
import rasterio
import xarray as xr
from rasterio.warp import reproject, Resampling
from scipy.stats import mode, truncnorm
import os, sys
import argparse
from argparse import RawTextHelpFormatter
import traceback
from typing import Dict
from scipy import stats

if sys.version_info < (3, 8):  # pragma: no cover (>=py38)
    import importlib_metadata
    import importlib_resources
else:  # pragma: no cover (<py38)
    import importlib.metadata as importlib_metadata
    import importlib.resources as importlib_resources

def main(argv=None):

    ''' Add WUDAPT info to WRF's '''

    parser = argparse.ArgumentParser(
        description="PURPOSE: Add LCZ-based info to WRF geo_em.d0X.nc\n \n"
                    "OUTPUT:\n"
                    "- *_NoUrban.nc: MODIS Urban replaced by surrounding natural LC\n"
                    "- *_LCZ_extent.nc: LCZ urban extent implemented, no LCZ UCPs yet\n"
                    "- *_LCZ_params.nc: LCZ urban extent + UPC parameter values\n"
                    "- *_d0X_41.nc: Parent domain files reflecting 41 Land categories",
        formatter_class=RawTextHelpFormatter
    )

    # Required arguments
    parser.add_argument(type=str, dest='io_dir',
                        help='Directory that contains geo_em.d0X.nc and LCZ.tif file',
                        )
    parser.add_argument(type=str, dest='lcz_file',
                        help='LCZ map file name',
                        )
    parser.add_argument(type=str, dest='wrf_file',
                        help='WRF geo_em* file name',
                        )

    # Additional arguments
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=f'%(prog)s {importlib_metadata.version("w2w")}',
    )
    parser.add_argument('-b', '--built-lcz',
                        nargs='+',
                        metavar='',
                        type=int,
                        dest='built_lcz',
                        help='LCZ classes considered as urban '
                             '(DEFAULT: 1 2 3 4 5 6 7 8 9 10)',
                        default=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    parser.add_argument('-l', '--lcz-band',
                        metavar='',
                        type=int,
                        dest='LCZ_BAND',
                        help='Band to use from LCZ file (DEFAULT: 0). '
                             'For maps produced with LCZ Generator, use 1',
                        default=0)
    parser.add_argument('-f', '--frc-threshold',
                        metavar='',
                        type=float,
                        dest='FRC_THRESHOLD',
                        help='FRC_URB2D treshold value to assign pixel as urban '
                             '(DEFAULT: 0.2)',
                        default=0.2)
    parser.add_argument('-n', '--npix-nlc',
                        metavar='',
                        type=int,
                        dest='NPIX_NLC',
                        help='Number of pixels to use for sampling neighbouring '
                             'natural land cover (DEFAULT: 45)',
                        default=45)
    args = parser.parse_args(argv)

    # Define output and tmp file(s), the latter is removed when done.
    dst_nu_file = os.path.join(
        args.io_dir,
        args.wrf_file.replace('.nc','_NoUrban.nc')
    )
    dst_gridinfo = os.path.join( # TMP file, will be removed
        args.io_dir,
        args.wrf_file.replace('.nc','_gridinfo.tif')
    )
    dst_lcz_extent_file = os.path.join(
        args.io_dir,
        args.wrf_file.replace('.nc', '_LCZ_extent.nc')
    )
    dst_lcz_params_file = os.path.join(
        args.io_dir,
        args.wrf_file.replace('.nc', '_LCZ_params.nc')
    )

    # Put all information in info dictionary
    info = {
        'io_dir': args.io_dir,
        'src_file': os.path.join(args.io_dir, args.lcz_file),
        'dst_file': os.path.join(args.io_dir, args.wrf_file),
        'dst_nu_file': dst_nu_file,
        'dst_gridinfo': dst_gridinfo,
        'dst_lcz_extent_file': dst_lcz_extent_file,
        'dst_lcz_params_file': dst_lcz_params_file,
        'BUILT_LCZ': args.built_lcz,
    }

    # Execute the functions
    print("Check if LCZ domain extends WRF domain in all directions?")
    check_lcz_wrf_extent(
        info=info,
    )
    print("")

    print("Replace WRF MODIS urban LC with surrounding natural LC")
    wrf_remove_urban(
        info=info,
        NPIX_NLC=args.NPIX_NLC,
    )
    print("")

    print("Create temporary WRF grid .tif file for resampling")
    create_wrf_gridinfo(
        info=info,
    )
    print("")

    print("+ FRC_URB2D, alter LU_INDEX, GREENFRAC and LANDUSEF")
    frc_mask = add_frc_lu_index_2_wrf(
        info=info,
        LCZ_BAND=args.LCZ_BAND,
        FRC_THRESHOLD=args.FRC_THRESHOLD,
        LCZ_NAT_MASK=True,
    )
    print("")

    print("+ LCZ-based UCP values into WRF's URB_PARAM")
    nbui_max = add_urb_params_to_wrf(
        info=info,
        LCZ_BAND=args.LCZ_BAND,
    )
    print("")

    print("Create LCZ-based urban extent file (excluding other LCZ-based info).")
    create_extent_file(
        info=info,
        frc_mask=frc_mask,
    )
    print("")

    print("Expanding land categories of parent domain(s) to 41:")
    expand_land_cat_parents(
        info=info,
    )
    print("")
    print("******************************")
    print(f"Set nbui_max to {nbui_max} during compilation, "
          "in order to optimize memory storage.")
    print("******************************")
    print("")

    print("Start sanity check and clean-up ...")
    checks_and_cleaning(
        info=info,
    )
    print("")
    print("********* All done ***********")


def check_lcz_wrf_extent(info: Dict[str, str]) -> None:

    # Read the data
    lcz = rasterio.open(info['src_file'])
    wrf = xr.open_dataset(info['dst_file'])

    # Get bounding box coordinates
    lcz_xmin, lcz_ymin, lcz_xmax, lcz_ymax = lcz.bounds
    wrf_xmin, wrf_ymin, wrf_xmax, wrf_ymax = \
        float(wrf.XLONG_M.min()), float(wrf.XLAT_M.min()), \
        float(wrf.XLONG_M.max()), float(wrf.XLAT_M.max())

    # Evaluate and throw error if wrf not within LCZ domain
    if not (wrf_xmin > lcz_xmin ) & (wrf_xmax < lcz_xmax ) & \
           (wrf_ymin > lcz_ymin) & (wrf_ymax < lcz_ymax):

        print("ERROR: LCZ domain should be larger than WRF domain "
              "in all directions.")
        print(f"LCZ bounds  (xmin, ymin, xmax, ymax): "
              f"{lcz_xmin, lcz_ymin, lcz_xmax, lcz_ymax}")
        print(f"WRF bounds  (xmin, ymin, xmax, ymax): "
              f"{wrf_xmin, wrf_ymin, wrf_xmax, wrf_ymax}")

        sys.exit()
    else:
        print("OK - LCZ domain is covering WRF domain")

def wrf_remove_urban(
        info,
        NPIX_NLC,
):

    '''Remove MODIS urban extent from geo_em*.nc file'''

    # Make a copy of original dst file
    dst_data = xr.open_dataset(info['dst_file'])

    # Read the relevant parameters
    luse = dst_data.LU_INDEX.squeeze()
    luf = dst_data.LANDUSEF.squeeze()
    greenf = dst_data.GREENFRAC.squeeze()
    lat = dst_data.XLAT_M.squeeze()
    lon = dst_data.XLONG_M.squeeze()
    newluse=luse.values.copy()
    newluf=luf.values.copy()
    newgreenf=greenf.values.copy()

    # Convert urban to surrounding natural characteristics
    for i in dst_data.south_north:
        for j in dst_data.west_east:
            if luse.isel(south_north=i,west_east=j) == 13:
                dis = (lat.where((luse!=13) &
                                 (luse!=17) &
                                 (luse!=21)
                                 )-lat.isel(south_north=i,west_east=j))**2 +\
                      (lon.where((luse!=13) &
                                 (luse!=17) &
                                 (luse!=21)
                                 )-lon.isel(south_north=i,west_east=j))**2

                disflat = dis.stack(gridpoints=('south_north','west_east'))\
                    .reset_index('gridpoints').drop_vars(['south_north','west_east'])
                aux = luse.where(dis<disflat.sortby(disflat)
                                 .isel(gridpoints=NPIX_NLC),drop=True)
                m = stats.mode(aux.values.flatten(), nan_policy="omit")[0]
                newluse[i, j] = int(m)

                auxg = greenf.where(dis<disflat.sortby(disflat)
                                    .isel(gridpoints=NPIX_NLC),drop=True)\
                    .where(aux==newluse[i,j]).mean(dim=['south_north','west_east'])
                newgreenf[:,i,j]=auxg

            if luf.isel(south_north=i,west_east=j,land_cat=12)>0.:
                dis = (lat.where((luf.isel(land_cat=12)==0.) &
                                 (luf.isel(land_cat=16)==0.) &
                                 (luf.isel(land_cat=20)==0.)
                                 )-lat.isel(south_north=i,west_east=j))**2 +\
                      (lon.where((luf.isel(land_cat=12)==0.) &
                                 (luf.isel(land_cat=16)==0.) &
                                 (luf.isel(land_cat=20)==0.)
                                 )-lon.isel(south_north=i,west_east=j))**2

                disflat = dis.stack(gridpoints=('south_north','west_east'))\
                    .reset_index('gridpoints').drop_vars(['south_north','west_east'])
                aux = luse.where(dis<disflat.sortby(disflat)
                                 .isel(gridpoints=NPIX_NLC),drop=True)
                m = stats.mode(aux.values.flatten(), nan_policy="omit")[0]
                newlu = int(m) - 1
                #newlu = int(mode(aux.values.flatten())[0])-1
                newluf[newlu,i,j]+=luf.isel(south_north=i,west_east=j,land_cat=12).values
                newluf[12,i,j]=0.


    dst_data.LU_INDEX.values[0,:]=newluse[:]
    dst_data.LANDUSEF.values[0,:]=newluf[:]
    dst_data.GREENFRAC.values[0,:]=newgreenf[:]

    # Save to final _lcz_params file
    if os.path.exists(info['dst_nu_file']):
        os.remove(info['dst_nu_file'])
    dst_data.to_netcdf(info['dst_nu_file'])


# Make WRF grid info available for Resampler (tmp file)
def create_wrf_gridinfo(
        info,
):

    # Read  gridded WRF data
    dst_data = xr.open_dataset(info['dst_nu_file'])

    # Create simpler WRF grid target.
    da_lu = xr.Dataset(
        {'LU_INDEX': (['y', 'x'],  dst_data['LU_INDEX'][0, :, :].values)},
        coords={'y': dst_data.XLAT_M.values[0, :, 0],
                'x': dst_data.XLONG_M.values[0, 0, :]}
    )

    # Add projection information as attributes, save and read back in.
    da_lu.rio.write_crs("epsg:4326", inplace=True)
    da_lu.rio.to_raster(info['dst_gridinfo'])

    return 0

def _ucp_resampler(
        info,
        ucp_key,
        RESAMPLE_TYPE,
        LCZ_BAND,
        **kwargs,
):

    '''Helper function to resample lcz ucp data to WRF grid'''

    # Read the look-up table
    ucp_table = pd.read_csv(
        importlib_resources.open_text('w2w.resources', 'LCZ_UCP_lookup.csv'),
        sep=',', index_col=0
    ).iloc[:17, :]

    # Read gridded data: LCZ and WRF grid
    src_data = xr.open_rasterio(info['src_file'])[LCZ_BAND, :, :]
    dst_grid = xr.open_rasterio(info['dst_gridinfo'])

    # Get Look-up for FRC_values
    if ucp_key in ['LB_URB2D', 'LF_URB2D', 'LP_URB2D']:

        # Following Zonato et al (2020)
        # and building width values from URB_PARAM.LCZ_TBL
        # Street width extracted from S02012 Building heighht and H2W.
        SW = ucp_table['MH_URB2D'] / ucp_table['H2W']
        LAMBDA_P = ucp_table['BW'] / (ucp_table['BW'] + SW)
        LAMBDA_F = 2 * ucp_table['MH_URB2D'] / (ucp_table['BW'] + SW)
        LAMBDA_B = LAMBDA_P + LAMBDA_F

        if ucp_key == 'LB_URB2D':
            lookup = LAMBDA_B.loc[info['BUILT_LCZ']]
        elif ucp_key == 'LP_URB2D':
            lookup = LAMBDA_P.loc[info['BUILT_LCZ']]
        elif ucp_key == 'LF_URB2D':
            lookup = LAMBDA_F.loc[info['BUILT_LCZ']]

    elif ucp_key == 'STDH_URB2D':
        lookup = ((ucp_table['MH_URB2D_MAX']-
                  ucp_table['MH_URB2D_MIN'])/4).loc[info['BUILT_LCZ']]
    else:
        lookup = ucp_table[ucp_key].loc[info['BUILT_LCZ']]

    # Get mask of selected built LCZs
    lcz_urb_mask = xr.DataArray(
        np.in1d(src_data, info['BUILT_LCZ']).reshape(src_data.shape),
        dims=src_data.dims, coords=src_data.coords
    )

    # Get LCZ class values only.
    lcz_arr = src_data.values

    # Set LCZ classes not in BUILT_LCZ to 0
    lcz_arr[~lcz_urb_mask] = 0

    # Make replacer object to map UCP values on LCZ class values
    replacer = np.zeros((max(info['BUILT_LCZ']) + 1,), object)
    replacer[lookup.index.values] = lookup
    lcz_data = np.array(replacer[lcz_arr], dtype='float')

    # Store into dataarray for resampling
    lcz_data_da = xr.Dataset(
        {'band': (['y', 'x'], lcz_data)},
        coords={'y': src_data.y.values, 'x': src_data.x.values},
        attrs={'transform': src_data.transform, 'crs': src_data.crs}
    ).to_array()

    # Info: https://rasterio.readthedocs.io/en/latest/api/rasterio.warp.html?highlight=reproject(#rasterio.warp.reproject
    ucp_2_wrf = reproject(
        lcz_data_da,
        dst_grid,
        src_transform=lcz_data_da.attrs['transform'],
        src_crs=lcz_data_da.attrs['crs'],
        dst_transform=dst_grid.attrs['transform'],
        dst_crs=dst_grid.attrs['crs'],
        resampling=Resampling[RESAMPLE_TYPE])[0]

    # In case of FRC_URB2D, filter for too low values
    if 'FRC_THRESHOLD' in kwargs.keys():
        ucp_2_wrf = ucp_2_wrf.where(
            ucp_2_wrf > kwargs['FRC_THRESHOLD'],
            0
        )

    ## In case nans occur, set to zero
    ucp_2_wrf.values[0, np.isnan(ucp_2_wrf[0, :, :])] = 0

    return ucp_2_wrf


def _hgt_resampler(
        info,
        RESAMPLE_TYPE,
        LCZ_BAND,
):

    '''Helper function to resample lcz ucp data to WRF grid'''

    # Read the look-up table
    ucp_table = pd.read_csv(
        importlib_resources.open_text('w2w.resources', 'LCZ_UCP_lookup.csv'),
        sep=',', index_col=0
    ).iloc[:17, :]

    # Read gridded data: LCZ and WRF grid
    src_data = xr.open_rasterio(info['src_file'])[LCZ_BAND, :, :]
    dst_grid = xr.open_rasterio(info['dst_gridinfo'])

    # Get Look-up for HGT values
    lookup_nom = ucp_table['BW'].loc[info['BUILT_LCZ']] ** 2 \
                 * ucp_table['MH_URB2D'].loc[info['BUILT_LCZ']]
    lookup_denom = ucp_table['BW'].loc[info['BUILT_LCZ']] ** 2

    # Get mask of selected built LCZs
    lcz_urb_mask = xr.DataArray(
        np.in1d(src_data, info['BUILT_LCZ']).reshape(src_data.shape),
        dims=src_data.dims, coords=src_data.coords
    )

    # Get LCZ class values only.
    lcz_arr = src_data.values

    # Set LCZ classes not in BUILT_LCZ to 0
    lcz_arr[~lcz_urb_mask] = 0

    # Make replacer object for nominator
    replacer_nom = np.zeros((max(info['BUILT_LCZ']) + 1,), object)
    replacer_nom[lookup_nom.index.values] = lookup_nom
    dataLcz_nom = np.array(replacer_nom[lcz_arr], dtype='float')

    # Make replacer object for denominator
    replacer_denom = np.zeros((max(info['BUILT_LCZ']) + 1,), object)
    replacer_denom[lookup_denom.index.values] = lookup_denom
    dataLcz_denom = np.array(replacer_denom[lcz_arr], dtype='float')

    # Store into dataarray for resampling
    lcz_data_da_nom = xr.Dataset(
        {'band': (['y', 'x'], dataLcz_nom)},
        coords={'y': src_data.y.values, 'x': src_data.x.values},
        attrs={'transform': src_data.transform, 'crs': src_data.crs}
    ).to_array()
    lcz_data_da_denom = xr.Dataset(
        {'band': (['y', 'x'], dataLcz_denom)},
        coords={'y': src_data.y.values, 'x': src_data.x.values},
        attrs={'transform': src_data.transform, 'crs': src_data.crs}
    ).to_array()

    # Get the aggregated values on WRF grid - nominator
    ucp_2_wrf_nom = reproject(
        lcz_data_da_nom,
        dst_grid,
        src_transform=lcz_data_da_nom.attrs['transform'],
        src_crs=lcz_data_da_nom.attrs['crs'],
        dst_transform=dst_grid.attrs['transform'],
        dst_crs=dst_grid.attrs['crs'],
        resampling=Resampling[RESAMPLE_TYPE])[0].copy()

    # Get the aggregated values on WRF grid - nominator
    ucp_2_wrf_denom = reproject(
        lcz_data_da_denom,
        dst_grid,
        src_transform=lcz_data_da_denom.attrs['transform'],
        src_crs=lcz_data_da_denom.attrs['crs'],
        dst_transform=dst_grid.attrs['transform'],
        dst_crs=dst_grid.attrs['crs'],
        resampling=Resampling[RESAMPLE_TYPE])[0].copy()

    hgt_urb2d = ucp_2_wrf_nom / ucp_2_wrf_denom

    ## In case nans occur, set to zero
    hgt_urb2d.values[0,np.isnan(hgt_urb2d[0,:,:])] = 0

    return hgt_urb2d

def _scale_hi(
        array,
):

    ''' Helper function to scale HI_URB2D to 100%'''

    return [(float(i) / sum(array) * 100.0) for i in array]

def _get_truncated_normal(
        mean,
        sd,
        low,
        upp,
):

    ''' Helper function to return bounded normal distribution'''

    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def _compute_hi_distribution(
        info,
        SAMPLE_SIZE=5000000,
        DIST_MARGIN=0.15,
        # HI_THRES_MIN=5,
):

    ''' Helper function to compute building height distribution'''

    # Read the look-up table
    ucp_table = pd.read_csv(
        importlib_resources.open_text('w2w.resources', 'LCZ_UCP_lookup.csv'),
        sep=',', index_col=0
    ).iloc[:17, :]

    # Initialize dataframe that stores building height distributions
    df_hi = pd.DataFrame(
        index = range(1,18,1),
        columns = ['0 - <5m', '5 - <10m', '10 - <15m', '15 - <20m',
                   '20 - <25m', '25 - <30m', '30 - <35m', '35 - <40m',
                   '40 - <45m', '45 - <50m', '50 - <55m', '55 - <60m',
                   '60 - <65m', '65 - <70m', '70 - <75m']
    )

    for i in info['BUILT_LCZ']:

        # LCZ 15 = paved, and considered to have no buildings (values = 0%)
        if not i == 15:

            # Create instance of a truncated normal distribution
            hi_inst = _get_truncated_normal(
                mean=ucp_table['MH_URB2D'].loc[i],
                sd=(ucp_table['MH_URB2D_MAX'].loc[i]-
                    ucp_table['MH_URB2D_MIN'].loc[i])/4,
                low=ucp_table['MH_URB2D_MIN'].loc[i],
                upp=ucp_table['MH_URB2D_MAX'].loc[i]
            )

            # populate with large enough sample for accuracy
            hi_sample = hi_inst.rvs(SAMPLE_SIZE)

            # Produce warning if approximated HI_URB2D distribution metrics
            # are not as expected: using a DIST_MARGIN % marging here.
            hi_metric = 'MH_URB2D_MIN'
            if not ucp_table[hi_metric].loc[i] * (1-DIST_MARGIN) < \
                hi_sample.min() < \
                ucp_table[hi_metric].loc[i] * (1+DIST_MARGIN):
                print("WARNING: MIN of HI_URB2D distribution not in "
                      f"expected range ({DIST_MARGIN}% marging) for LCZ class {i}: "
                      f"modelled: {np.round(hi_sample.min(),2)} | "
                      f"expected: [{(ucp_table[hi_metric].loc[i] * (1-DIST_MARGIN)).round(2)} - "
                      f"{(ucp_table[hi_metric].loc[i] * (1-DIST_MARGIN)).round(2)}]")

            hi_metric = 'MH_URB2D_MAX'
            if not ucp_table[hi_metric].loc[i] * (1-DIST_MARGIN) < \
                hi_sample.max() < \
                ucp_table[hi_metric].loc[i] * (1+DIST_MARGIN):
                print("WARNING: MAX of HI_URB2D distribution not in "
                      f"expected range ({DIST_MARGIN}% marging) for LCZ class {i}: "
                      f"modelled: {np.round(hi_sample.max(),2)} | "
                      f"expected: [{(ucp_table[hi_metric].loc[i] * (1-DIST_MARGIN)).round(2)} - "
                      f"{(ucp_table[hi_metric].loc[i] * (1+DIST_MARGIN)).round(2)}]")

            hi_metric = 'MH_URB2D'
            if not ucp_table[hi_metric].loc[i] * (1-DIST_MARGIN) < \
                hi_sample.mean() < \
                ucp_table[hi_metric].loc[i] * (1+DIST_MARGIN):
                print("WARNING: MEAN of HI_URB2D distribution not in "
                      f"expected range ({DIST_MARGIN}% marging) for LCZ class {i}: "
                      f"modelled: {np.round(hi_sample.mean(),2)} | "
                      f"expected: [{(ucp_table[hi_metric].loc[i] * (1-DIST_MARGIN)).round(2)} - "
                      f"{(ucp_table[hi_metric].loc[i] * (1+DIST_MARGIN)).round(2)}]")

            # Count the values within pre-set bins
            cnt = np.histogram(hi_sample, bins=np.arange(0,76,5))[0]
            cnt = cnt/(SAMPLE_SIZE/100) # Convert to %

            # Add to dataframe
            df_hi.loc[i,:] = cnt

        # Set nans to zero
        df_hi = df_hi.fillna(0)

    return df_hi

def _hi_resampler(
        info,
        RESAMPLE_TYPE,
        LCZ_BAND,
        HI_THRES_MIN=5,
):

    '''Helper function to resample ucp HI_URB2D_URB2D data to WRF grid'''

    # Read gridded data: LCZ and WRF grid
    src_data = xr.open_rasterio(info['src_file'])[LCZ_BAND, :, :]
    dst_grid = xr.open_rasterio(info['dst_gridinfo'])

    # Get mask of selected built LCZs
    lcz_urb_mask = xr.DataArray(
        np.in1d(src_data, info['BUILT_LCZ']).reshape(src_data.shape),
        dims=src_data.dims, coords=src_data.coords
    )

    # Get LCZ class values only.
    lcz_arr = src_data.values

    # Set LCZ classes not in BUILT_LCZ to 0
    lcz_arr[~lcz_urb_mask] = 0

    # Compute the building height densities.
    df_hi = _compute_hi_distribution(info)

    # Initialize array to store temp values
    hi_arr = np.zeros((15,dst_grid.shape[1],dst_grid.shape[2]))

    # Loop over the 15 height density classes.
    for hi_i in range(df_hi.shape[1]):

        print(f"Working on height interval {df_hi.columns[hi_i]} ...")
        lookup = df_hi.iloc[:, hi_i].loc[info['BUILT_LCZ']]

        # Make replacer object to map UCP values on LCZ class values
        replacer = np.zeros((max(info['BUILT_LCZ']) + 1,), object)
        replacer[lookup.index.values] = lookup
        lcz_data = np.array(replacer[lcz_arr], dtype='float')

        # Store into dataarray for resampling
        lcz_data_da = xr.Dataset(
            {'band': (['y', 'x'], lcz_data)},
            coords={'y': src_data.y.values, 'x': src_data.x.values},
            attrs={'transform': src_data.transform, 'crs': src_data.crs}
        ).to_array()

        # Get the aggregated values on WRF grid
        ucp_2_wrf = reproject(
            lcz_data_da,
            dst_grid,
            src_transform=lcz_data_da.attrs['transform'],
            src_crs=lcz_data_da.attrs['crs'],
            dst_transform=dst_grid.attrs['transform'],
            dst_crs=dst_grid.attrs['crs'],
            resampling=Resampling[RESAMPLE_TYPE])[0]

        ## In case nans occur, set to zero
        ucp_2_wrf.values[0, np.isnan(ucp_2_wrf[0, :, :])] = 0

        # Store UCPs in tmp hi_arr
        hi_arr[hi_i, :, :] = ucp_2_wrf[0, :, :]

        # For computational efficiency/storage, set values lower than
        # 5% (HI_THRES_MIN) to 0
        hi_arr[hi_arr < HI_THRES_MIN] = 0

    # re-scale HI_URB2D to 100% when summed over 118-132 indices!
    hi_arr_scaled = np.apply_along_axis(_scale_hi, 0, hi_arr)
    hi_arr_scaled[np.isnan(hi_arr_scaled)] = 0

    # Count max number of HI intervals over all grid cells.
    nbui_max = np.where(hi_arr_scaled, 1, 0).sum(axis=0).max()

    return hi_arr_scaled, nbui_max


def _lcz_resampler(
        info,
        frc_urb2d,
        LCZ_NAT_MASK,
        LCZ_BAND,
):

    '''Helper function to resample lcz classes to WRF grid'''

    # Read required gridded data, LCZ, WRF grid, and
    # original WRF (for original MODIS urban mask)
    src_data = xr.open_rasterio(info['src_file'])[LCZ_BAND, :, :]
    dst_grid = xr.open_rasterio(info['dst_gridinfo'])

    # Mask natural LCZs before majority filtering.
    if LCZ_NAT_MASK:
        src_data = src_data.where(
            src_data.isin(info['BUILT_LCZ'])
        ).copy()

    lcz_2_wrf = reproject(
        src_data,
        dst_grid,
        src_transform=src_data.attrs['transform'],
        src_crs=src_data.attrs['crs'],
        dst_transform=dst_grid.attrs['transform'],
        dst_crs=dst_grid.attrs['crs'],
        resampling=Resampling['mode'])[0].values

    # if LCZ 15 selected in 'BUILT_LCZ', rename to 11
    if 15 in info['BUILT_LCZ']:
        lcz_2_wrf[lcz_2_wrf == 15] = 11

    # Only keep LCZ pixels where FRC_URB2D > 0, for concistency
    frc_mask = frc_urb2d.values[0,:,:] != 0

    # Final LU_INDEX = 31 to 41 (included), as LCZ classes.
    lcz_resampled = lcz_2_wrf[0,frc_mask] + 30

    return frc_mask, lcz_resampled


def _adjust_greenfrac_landusef(
        info,
        dst_data,
        frc_mask,
):

    dst_data_orig = xr.open_dataset(info['dst_file'])

    # Adjust GREENFRAC and LANDUSEF
    # GREENFRAC is set as average / month from GREENFRAC
    # of original MODIS urban pixels
    wrf_urb = xr.DataArray(
        np.in1d(dst_data_orig['LU_INDEX'][0, :, :].values, [13])\
            .reshape(dst_data_orig['LU_INDEX'][0, :, :].shape),
        dims=dst_data_orig['LU_INDEX'][0, :, :].dims,
        coords=dst_data_orig['LU_INDEX'][0, :, :].coords
    )
    greenfrac_per_month = [
        dst_data_orig['GREENFRAC'].values[0, mm, wrf_urb].mean()
        for mm in range(12)
    ]

    # Loop over months and set average values
    for mm in range(12):
        dst_data['GREENFRAC'].values[0, mm, frc_mask] = \
            greenfrac_per_month[mm]

    # TODO: For lower resolution domains, this might not be valid?
    # Create new LANDUSEF with 41 levels instead of 21
    landusef_new = np.zeros(
        (41, dst_data.LANDUSEF.shape[2], dst_data.LANDUSEF.shape[3])
    )

    # Copy values from original file
    landusef_new[:21,:,:] = dst_data['LANDUSEF'][0, :21, :, :]

    # First set all values to zero for urban mask
    landusef_new[:, frc_mask] = 0 # First all to 0, so sum remains 1 in the end

    # LOOP over LCZ LU_INDEX values, and set to 1 there
    # So e.g. LANDUSE[0,31-1,:,1] = 1, where LU_INDEX = 31 (=LCZ 1)
    for lu_i in np.arange(31,42,1):
        lu_mask = dst_data.LU_INDEX == int(lu_i)
        landusef_new[int(lu_i)-1, lu_mask[0, :, :]] = 1
        del lu_mask

    # First store orginal attributes, then drop variable
    luf_attrs = dst_data.LANDUSEF.attrs
    dst_data = dst_data.drop_vars('LANDUSEF')

    # Expand axis to take shape (1,41,x,y)
    landusef_new = np.expand_dims(landusef_new, axis=0)

    # Add back to data-array, including (altered) attributes
    dst_data['LANDUSEF'] =     (
        ('Time', 'land_cat', 'south_north', 'west_east'),
        landusef_new
    )
    dst_data['LANDUSEF'] = dst_data.LANDUSEF.astype('float32')

    luf_attrs['description'] = 'Noah-modified 41-category IGBP-MODIS landuse'
    for key in luf_attrs.keys():
        dst_data['LANDUSEF'].attrs[key] = luf_attrs[key]

    return dst_data


def add_frc_lu_index_2_wrf(
        info,
        LCZ_BAND,
        FRC_THRESHOLD,
        LCZ_NAT_MASK,
):

    '''
    Add FRC_URB2D and adjusted LCZ-based LU_INDEX to WRF file
    Also alters LANDUSEF and GREENFRAC in line with LU_INDEX
    '''

    # Integrate FRC_URB2D
    ucp_key = 'FRC_URB2D'

    # Get the aggrated frc_urb values
    frc_urb = _ucp_resampler(
        info=info,
        ucp_key=ucp_key,
        RESAMPLE_TYPE='average',
        LCZ_BAND=LCZ_BAND,
        FRC_THRESHOLD = FRC_THRESHOLD,
    )

    # Add to geo_em* that that has no MODIS urban
    dst_data = xr.open_dataset(info['dst_nu_file'])

    # Make a FRC_URB field and store aggregated data.
    dst_data[ucp_key] = dst_data['LU_INDEX'].copy()
    dst_data[ucp_key] = (
                    ('Time', 'south_north', 'west_east'),
                    frc_urb.data
                )

    # Add proper attributes to the FRC_URB2D field
    dst_data[ucp_key].attrs['FieldType'] = np.intc(104)
    dst_data[ucp_key].attrs['MemoryOrder'] = "XY"
    dst_data[ucp_key].attrs['units'] = "-"
    dst_data[ucp_key].attrs['description'] = "ufrac"
    dst_data[ucp_key].attrs['stagger'] = "M"
    dst_data[ucp_key].attrs['sr_x'] = np.intc(1)
    dst_data[ucp_key].attrs['sr_y'] = np.intc(1)

    # Integrate LU_INDEX, also adjusts GREENFRAC and LANDUSEF
    frc_mask, lcz_resampled = _lcz_resampler(
        info=info,
        frc_urb2d=dst_data['FRC_URB2D'],
        LCZ_NAT_MASK=LCZ_NAT_MASK,
        LCZ_BAND=LCZ_BAND,
    )

    # 2) as LU_INDEX = 30 to 41, as LCZ classes.
    dst_data['LU_INDEX'].values[0,frc_mask] = lcz_resampled

    # Also adjust GREENFRAC and LANDUSEF
    dst_data = _adjust_greenfrac_landusef(info, dst_data, frc_mask)

    # Save to final _lcz_params file
    if os.path.exists(info['dst_lcz_params_file']):
        os.remove(info['dst_lcz_params_file'])
    dst_data.to_netcdf(info['dst_lcz_params_file'])

    return frc_mask


def _initialize_urb_param(
        info,
):

    ''' Helper function to initialize URB_PARAM in WRF geo_em file'''

    dst_data = xr.open_dataset(info['dst_lcz_params_file'])
    URB_PARAM = np.zeros([1, 132,
                          len(dst_data.south_north),
                          len(dst_data.west_east)])

    # Add to destination WRF file, with attributes
    dst_data['URB_PARAM'] = \
        (('Time', 'num_urb_params', 'south_north', 'west_east'), URB_PARAM)
    dst_data['URB_PARAM'].attrs['FieldType'] = np.intc(104)
    dst_data['URB_PARAM'].attrs['MemoryOrder'] = "XYZ"
    dst_data['URB_PARAM'].attrs['units'] = "dimensionless"
    dst_data['URB_PARAM'].attrs['description'] = "all urban parameters"
    dst_data['URB_PARAM'].attrs['stagger'] = "M"
    dst_data['URB_PARAM'].attrs['sr_x'] = np.intc(1)
    dst_data['URB_PARAM'].attrs['sr_y'] = np.intc(1)

    return dst_data

def add_urb_params_to_wrf(
        info,
        LCZ_BAND,
):

    ''' Map, aggregate and add lcz-based UCP values to WRF'''

    # Initialize empty URB_PARAM in final wrf file,
    # with all zeros and proper attributes
    dst_final = _initialize_urb_param(info)

    # get frc_mask, to only set values where FRC_URB2D > 0.
    frc_mask = dst_final.FRC_URB2D.values[0,:,:] != 0

    # Define the UCPs that need to be integrated,
    # together with their positions (index starts at 1) in URB_PARAMS
    # HGT_URB2D and HI_URB2D follow a different approach, see further.
    ucp_dict = {
        'LP_URB2D'  : 91,
        'MH_URB2D'  : 92,
        'STDH_URB2D': 93,
        'HGT_URB2D' : 94,
        'LB_URB2D'  : 95,
        'LF_URB2D'  : 96,   # 97, 98, 99, for all 4 directions
        'HI_URB2D'  : 118,  # Goes on until index 132
    }

    for ucp_key in ucp_dict.keys():

        print(f"Processing {ucp_key} ...")

        # Obtain aggregated LCZ-based UCP values
        if ucp_key in ['MH_URB2D', 'STDH_URB2D', 'LB_URB2D',
                       'LF_URB2D', 'LP_URB2D']:
            ucp_res = _ucp_resampler(
                info=info,
                ucp_key=ucp_key,
                RESAMPLE_TYPE='average',
                LCZ_BAND=LCZ_BAND,
            )
        elif ucp_key in ['HGT_URB2D']:
            ucp_res = _hgt_resampler(
                info=info,
                RESAMPLE_TYPE='average',
                LCZ_BAND=LCZ_BAND,
            )
        elif ucp_key in ['HI_URB2D']:
            ucp_res, nbui_max = _hi_resampler(
                info=info,
                RESAMPLE_TYPE='average',
                LCZ_BAND=LCZ_BAND
            )

        # Store UCPs in wrf destination file.
        if ucp_key == 'LF_URB2D':
            # Frontal area Index in N,E,S,W directions respectively
            # for WUDAPT LCZs they are considered all equal
            for i in range(4):
                ucp_res.values[:, frc_mask == 0] = 0
                dst_final['URB_PARAM'][:, ucp_dict[ucp_key]-1+i, :, :] = ucp_res
        if ucp_key == 'HI_URB2D':
            ucp_res[:,frc_mask==0] = 0
            dst_final['URB_PARAM'][0, (ucp_dict[ucp_key] - 1):, :, :] = ucp_res
        else:
            ucp_res.values[:, frc_mask == 0] = 0
            dst_final['URB_PARAM'].values[0, ucp_dict[ucp_key] - 1, :,:] = ucp_res

    # Make sure URB_PARAM is float32
    dst_final['URB_PARAM'] = dst_final.URB_PARAM.astype('float32')

    # Add/Change some additional global attributes,
    # including NBUI_MAX = max. nr. of HI intervals over the grid
    glob_attrs = {
        'NUM_LAND_CAT': 41,
        'FLAG_URB_PARAM': 1,
        'NBUI_MAX': np.intc(nbui_max),
    }
    for key in glob_attrs.keys():
        dst_final.attrs[key] = np.intc(glob_attrs[key])

    #TODO: add final repo link when done.
    #Add DESCRIPTION in attrs, referring to tool
    gh_repo = 'https://github.com/matthiasdemuzere/wrf-lcz-KL'
    dst_final.attrs['DESCRIPTION'] = \
        f"W2W.py tool used to create geo_em*.nc file: {gh_repo}"

    # Save back to file
    if os.path.exists(info['dst_lcz_params_file']):
        os.remove(info['dst_lcz_params_file'])
    dst_final.to_netcdf(info['dst_lcz_params_file'])

    return nbui_max

def create_extent_file(
        info,
        frc_mask,
):

    '''Create a domain file with an LCZ-based urban extent (excluding other LCZ-based info)'''

    dst_params = xr.open_dataset(info['dst_lcz_params_file'])
    dst_extent = dst_params.copy()

    lu_index = dst_extent.LU_INDEX.values
    lu_index[lu_index >= 31] = 13

    dst_extent.LU_INDEX.values = lu_index

    # Remove some unnecesary variables to reduce file size
    dst_extent = dst_extent.drop_vars(['FRC_URB2D','URB_PARAM'])

    # Reset LANDUSEF again to 21 classes.
    luf_attrs = dst_extent.LANDUSEF.attrs
    luf_values = dst_extent.LANDUSEF.values
    dst_extent = dst_extent.drop_vars('LANDUSEF')

    # Add back to data-array, including (altered) attributes
    dst_extent['LANDUSEF'] =     (
        ('Time', 'land_cat', 'south_north', 'west_east'),
        luf_values[:,:21,:,:]
    )
    dst_extent['LANDUSEF'].values[0, 12, frc_mask] = 1
    dst_extent['LANDUSEF'] = dst_extent.LANDUSEF.astype('float32')

    luf_attrs['description'] = 'Noah-modified 21-category IGBP-MODIS landuse'
    for key in luf_attrs.keys():
        dst_extent['LANDUSEF'].attrs[key] = luf_attrs[key]

    # Reset some other global attributes
    dst_extent.attrs['FLAG_URB_PARAM'] = np.intc(0)
    dst_extent.attrs['NUM_LAND_CAT'] = np.intc(21)

    # Save file.
    dst_extent.to_netcdf(info['dst_lcz_extent_file'])



def expand_land_cat_parents(
        info,
):

    # Get final domain number
    domain_nr = int(info['dst_file'][-5:-3])

    #list domain numbers to loop over
    domain_lst = list(np.arange(1,domain_nr,1))

    for i in domain_lst:

        ifile = f"{info['dst_file'][:-5]}{i:02d}.nc"

        try:
            da = xr.open_dataset(ifile)

        except Exception:
            print(f"WARNING: Parent domain {info['dst_file'][:-5]}{i:02d}.nc not found.\n"
                  f"Please make sure the parent domain files are in {info['io_dir']}\n"
                  f"Without this information, you will not be able to produce the boundary"
                  f"conditions with real.exe.")

        if int(da.attrs['NUM_LAND_CAT']) != 41:

            try:
                # Set number of land categories to 41
                da.attrs['NUM_LAND_CAT'] = np.intc(41)

                # Create new landusef array with expanded dimensions
                landusef_new = np.zeros(
                    (1, 41, da.LANDUSEF.shape[2], da.LANDUSEF.shape[3])
                )
                landusef_new[:, :21, :, :] = da['LANDUSEF'].values

                # First store orginal attributes, then drop variable
                luf_attrs = da.LANDUSEF.attrs
                da = da.drop_vars('LANDUSEF')

                # Add back to data-array, including (altered) attributes
                da['LANDUSEF'] = (
                    ('Time', 'land_cat', 'south_north', 'west_east'),
                    landusef_new
                )
                da['LANDUSEF'] = da.LANDUSEF.astype('float32')

                luf_attrs['description'] = 'Noah-modified 41-category IGBP-MODIS landuse'
                for key in luf_attrs.keys():
                    da['LANDUSEF'].attrs[key] = luf_attrs[key]

                ofile = ifile.replace('.nc', '_41.nc')
                da.to_netcdf(ofile)

            except Exception:
                err = traceback.format_exc()
                print(f'Cannot read change NUM_LAND_CAT and LANDUSEF dimensions\n{err}')

        else:
            print(f"Parent domain {info['dst_file'][:-5]}{i:02d}.nc "
                  f"already contains 41 LC classes")


def checks_and_cleaning(
        info,
):

    'Sanity checks and cleaning'

    print(f"Check 1: Urban class removed from "
          f"{info['dst_nu_file'].split('/')[-1]}?")
    ifile = info['dst_nu_file']
    da = xr.open_dataset(ifile)
    if 13 in da.LU_INDEX.values:
        print(f"WARNING: Urban land use still present")
    else:
        print(f"OK")
    print("")

    print(f"Check 2: LCZ Urban extent present in "
          f"{info['dst_lcz_extent_file'].split('/')[-1]}?")
    ifile = info['dst_lcz_extent_file']
    da = xr.open_dataset(ifile)
    if 13 in da.LU_INDEX.values:
        print(f"OK")
    else:
        print(f"WARNING: LCZ-based urban extent missing")
    print("")

    print(f"Check 3: Urban LCZ classes exists in "
          f"{info['dst_lcz_params_file'].split('/')[-1]}?")
    ifile = info['dst_lcz_params_file']
    da = xr.open_dataset(ifile)
    if 13 in da.LU_INDEX.values:
        print(f"WARNING: Urban extent still defined via LU_INDEX = 13?")
    else:
        LU_values = np.unique(da.LU_INDEX.values.flatten())
        LCZs = [int(i) for i in list(LU_values[LU_values >= 31] - 30)]
        print(f"OK: LCZ Classes ({LCZs}) present")
    print("")

    print(f"Check 4: URB_PARAMS matrix present in file "
          f"{info['dst_lcz_params_file'].split('/')[-1]}?")
    ifile = info['dst_lcz_params_file']
    da = xr.open_dataset(ifile)
    if not 'URB_PARAM' in list(da.keys()):
        print(f"WARNING: URB_PARAM matrix not present")
    else:
        print(f"OK")
    print("")

    print(f"Check 5: FRC_URB2D present in "
          f"{info['dst_lcz_params_file'].split('/')[-1]}?")
    ifile = info['dst_lcz_params_file']
    da = xr.open_dataset(ifile)
    if not 'FRC_URB2D' in list(da.keys()):
        print(f"WARNING: FRC_URB2D not present in {ifile}")
    else:
        FRC_URB2D = da.FRC_URB2D.values
        print(f"OK: FRC_URB2D values range between "
              f"{'{:0.2f}'.format(FRC_URB2D.min())} and "
              f"{'{:0.2f}'.format(FRC_URB2D.max())}")
    print("")

    print("Check 6: Do URB_PARAM variable values follow expected range in "
          f"{info['dst_lcz_params_file'].split('/')[-1]}?")
    ifile = info['dst_lcz_params_file']
    da = xr.open_dataset(ifile)

    # Take expected ranges from the look-up table,
    # add some margin for changes due to interpolation.
    ucp_table = pd.read_csv(
        importlib_resources.open_text('w2w.resources', 'LCZ_UCP_lookup.csv'),
        sep=',', index_col=0
    ).iloc[:17, :]

    ucp_dict = {
        'LP_URB2D'  : {
            'index': 91,
            'range': [0,1]
        },
        'MH_URB2D'  : {
            'index': 92,
            'range': [0, ucp_table['MH_URB2D'].max() + ucp_table['MH_URB2D'].std()]
        },
        'HGT_URB2D' : {
            'index': 94,
            'range': [0, ucp_table['MH_URB2D'].max() + ucp_table['MH_URB2D'].std()]
        },
        'LB_URB2D'  : {
            'index': 95,
            'range': [0, 5]
        },
        'LF_URB2D'  : {
            'index': 96,
            'range': [0, 5]
        },
        'LF_URB2D'  : {
            'index': 97,
            'range': [0, 5]
        },
        'LF_URB2D'  : {
            'index': 98,
            'range': [0, 5]
        },
        'LF_URB2D'  : {
            'index': 99,
            'range': [0, 5]
        },
    }

    def _check_range(darr, exp_range):

        total_len = len(darr)
        sel_len = ((darr >= exp_range[0]) & (darr <= exp_range[1])).sum(axis=0)

        if not (total_len - sel_len == 0):
            return -1
        else:
            return 0

    for ucp_key in ucp_dict.keys():

        darr = da.URB_PARAM[0,ucp_dict[ucp_key]['index']-1,:,:].values.flatten()
        exp_range = ucp_dict[ucp_key]['range']

        result = _check_range(darr, exp_range)

        if result == -1:
            print(f"WARNING: {ucp_key} exceeds expected value range")
        else:
            print(f"OK for {ucp_key}")
    print("")

    print("Check 7: Does HI_URB2D sum to 100% for urban pixels "
          f"in {info['dst_lcz_params_file'].split('/')[-1]}?")
    da = xr.open_dataset(info['dst_lcz_params_file'])
    hi_sum = da.URB_PARAM[0, 117:, :, :].sum(axis=0)
    hi_sum = hi_sum.where(hi_sum != 0, drop=True)

    if np.nanmax(np.abs((100 - hi_sum).values)) > 0.1:
        print(f"WARNING: Not all pixels have sum HI_URB2D == 100%")
    else:
        print(f"OK")
    print("")

    print("Check 8: Do FRC_URB and LCZs (from LU_INDEX) cover same extent "
          f"in {info['dst_lcz_params_file'].split('/')[-1]}?")
    frc_urb2d = xr.open_dataset(info['dst_lcz_params_file']).FRC_URB2D
    lu_index = xr.open_dataset(info['dst_lcz_params_file']).LU_INDEX
    frc_urb_res = xr.where(frc_urb2d != 0, 1, 0)
    lu_index_res = xr.where(lu_index >= 31, 1, 0)

    if int((frc_urb_res - lu_index_res).sum()) != 0:
        print(f"WARNING: FRC_URB and LCZs in LU_INDEX "
              f"do not cover same extent")
    else:
        print(f"OK")
    print("")

    print("Check 9: Extent and # urban pixels same for "
          "*_extent.nc and *_params.nc output file?")
    da_e = xr.open_dataset(info['dst_lcz_extent_file'])
    da_p = xr.open_dataset(info['dst_lcz_params_file'])
    da_e_res = xr.where(da_e.LU_INDEX == 13, 1,0)
    da_p_res = xr.where(da_p.LU_INDEX >=31, 1,0)

    if int((da_p_res - da_e_res).sum()) != 0:
        print(f"WARNING: Different # urban pixels (or extent) "
              f"according to LU_INDEX: "
              f" - extent: {int(da_e_res.sum().values)}"
              f" - params: {int(da_p_res.sum().values)}"
              )
    else:
        print(f"OK, urban extent the same. \n"
              f"Both files have {int(da_p_res.sum().values)} "
              f"urban pixels according to LU_INDEX")

    print("")
    print('Cleaning up ...')
    if os.path.exists(info['dst_gridinfo']):
        os.remove(info['dst_gridinfo'])

###############################################################################
##### __main__  scope
###############################################################################

if __name__ == "__main__":

    main()

###############################################################################
