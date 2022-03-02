import numpy as np

np.seterr(divide='ignore', invalid='ignore')
import pandas as pd
import rioxarray as rxr
import rasterio
import xarray as xr
from rasterio.warp import reproject, Resampling
from scipy.stats import mode, truncnorm
import os, sys
import argparse
from argparse import RawTextHelpFormatter
import traceback
from typing import Dict, Any
from typing import List
from typing import NamedTuple
from typing import Optional
from typing import Sequence
from typing import SupportsInt
from typing import Tuple
from typing import Union
from numpy.typing import NDArray
from scipy import stats
from pyproj import CRS

if sys.version_info >= (3, 9):  # pragma: >=3.9 cover
    import importlib.metadata as importlib_metadata
    import importlib.resources as importlib_resources
else:  # pragma: <3.9 cover
    import importlib_metadata
    import importlib_resources


def main(argv: Optional[Sequence[str]] = None) -> int:

    '''Add WUDAPT info to WRF's'''

    parser = argparse.ArgumentParser(
        description='PURPOSE: Add LCZ-based info to WRF geo_em.d0X.nc\n \n'
        'OUTPUT:\n'
        '- *_NoUrban.nc: MODIS Urban replaced by surrounding natural LC\n'
        '- *_LCZ_extent.nc: LCZ urban extent implemented, no LCZ UCPs yet\n'
        '- *_LCZ_params.nc: LCZ urban extent + UPC parameter values\n'
        '- *_d0X_41.nc: Parent domain files reflecting 41 Land categories',
        formatter_class=RawTextHelpFormatter,
    )

    # Required arguments
    parser.add_argument(
        type=str,
        dest='io_dir',
        help='Directory that contains geo_em.d0X.nc and LCZ.tif file',
    )
    parser.add_argument(
        type=str,
        dest='lcz_file',
        help='LCZ map file name',
    )
    parser.add_argument(
        type=str,
        dest='wrf_file',
        help='WRF geo_em* file name',
    )

    # Additional arguments
    parser.add_argument(
        '-V',
        '--version',
        action='version',
        version=f'%(prog)s {importlib_metadata.version("w2w")}',
    )
    parser.add_argument(
        '-b',
        '--built-lcz',
        nargs='+',
        metavar='',
        type=int,
        dest='built_lcz',
        help='LCZ classes considered as urban ' '(DEFAULT: 1 2 3 4 5 6 7 8 9 10)',
        default=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    )
    parser.add_argument(
        '-l',
        '--lcz-band',
        metavar='',
        type=int,
        dest='LCZ_BAND',
        help='Band to use from LCZ file (DEFAULT: 0). '
        'For maps produced with LCZ Generator, use 1',
        default=0,
    )
    parser.add_argument(
        '-f',
        '--frc-threshold',
        metavar='',
        type=float,
        dest='FRC_THRESHOLD',
        help='FRC_URB2D treshold value to assign pixel as urban ' '(DEFAULT: 0.2)',
        default=0.2,
    )
    parser.add_argument(
        '-n',
        '--npix-nlc',
        metavar='',
        type=int,
        dest='NPIX_NLC',
        help='Number of pixels to use for sampling neighbouring '
        'natural land cover (DEFAULT: 45)',
        default=45,
    )
    parser.add_argument(
        '--lcz-ucp',
        type=str,
        help='File with custom LCZ-based urban canopy parameters',
    )
    args = parser.parse_args(argv)

    # check if a custom LCZ UCP file was set and read it
    if args.lcz_ucp is not None:
        lookup_table = args.lcz_ucp
    else:
        lookup_table = importlib_resources.files(
            'w2w.resources',
        ).joinpath('LCZ_UCP_lookup.csv')

    # Execute the functions
    print('--> Set files structure')
    info = Info.from_argparse(args)
    ucp_table = pd.read_csv(lookup_table, index_col=0)

    print('--> Check LCZ integrity, in terms of ' 'class labels, projection and extent')
    check_lcz_integrity(
        info=info,
        LCZ_BAND=args.LCZ_BAND,
    )

    print('--> Replace WRF MODIS urban LC with surrounding natural LC')
    wrf_remove_urban(
        info=info,
        NPIX_NLC=args.NPIX_NLC,
    )

    print('--> Create temporary WRF grid .tif file for resampling')
    create_wrf_gridinfo(
        info=info,
    )

    print('--> Create LCZ-based geo_em file')
    nbui_max = create_lcz_params_file(
        info=info,
        FRC_THRESHOLD=args.FRC_THRESHOLD,
        LCZ_NAT_MASK=True,
        ucp_table=ucp_table,
    )

    print(
        '--> Create LCZ-based urban extent geo_em file '
        '(excluding other LCZ-based info)'
    )
    create_lcz_extent_file(
        info=info,
    )

    print('--> Expanding land categories of parent domain(s) to 41')
    expand_land_cat_parents(
        info=info,
    )

    print('\n--> Start sanity check and clean-up ...')
    checks_and_cleaning(
        info=info,
        ucp_table=ucp_table,
        nbui_max=nbui_max,
    )
    return 0


class Info(NamedTuple):
    """
    Immutable class representing the configuration with all files and directories
    """

    io_dir: str
    src_file: str
    src_file_clean: str
    dst_file: str
    dst_nu_file: str
    dst_gridinfo: str
    dst_lcz_extent_file: str
    dst_lcz_params_file: str
    BUILT_LCZ: List[int]

    @classmethod
    def from_argparse(cls, args: argparse.Namespace) -> 'Info':
        # Define output and tmp file(s), the latter is removed when done.
        return cls(
            io_dir=args.io_dir,
            src_file=os.path.join(args.io_dir, args.lcz_file),
            src_file_clean=os.path.join(
                args.io_dir, args.lcz_file.replace('.tif', '_clean.tif')
            ),
            dst_file=os.path.join(args.io_dir, args.wrf_file),
            dst_nu_file=os.path.join(
                args.io_dir, args.wrf_file.replace('.nc', '_NoUrban.nc')
            ),
            # TMP file, will be removed
            dst_gridinfo=os.path.join(
                args.io_dir, args.wrf_file.replace('.nc', '_gridinfo.tif')
            ),
            dst_lcz_extent_file=os.path.join(
                args.io_dir, args.wrf_file.replace('.nc', '_LCZ_extent.nc')
            ),
            dst_lcz_params_file=os.path.join(
                args.io_dir, args.wrf_file.replace('.nc', '_LCZ_params.nc')
            ),
            BUILT_LCZ=args.built_lcz,
        )


def _replace_lcz_number(
    lcz: xr.DataArray, lcz_to_change: NDArray[np.int_]
) -> xr.DataArray:
    lcz_expected = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])

    lcz_arr = lcz.data.flatten().astype(np.int32)
    df = pd.Series(lcz_arr, dtype=lcz_expected.dtype)

    d = dict(zip(lcz_to_change, lcz_expected))

    lcz_new = lcz.copy()
    lcz_new.data = df.map(d).values.reshape(lcz.shape).astype(np.int32)

    return lcz_new


def _check_lcz_wrf_extent(lcz: xr.DataArray, wrf: xr.Dataset) -> None:

    ERROR = '\033[0;31m'
    ENDC = '\033[0m'

    # Get bounding box coordinates
    lcz_xmin, lcz_ymin, lcz_xmax, lcz_ymax = lcz.rio.bounds()
    wrf_xmin, wrf_ymin, wrf_xmax, wrf_ymax = (
        float(wrf.XLONG_M.min()),
        float(wrf.XLAT_M.min()),
        float(wrf.XLONG_M.max()),
        float(wrf.XLAT_M.max()),
    )

    # Evaluate and throw error if WRF not within LCZ domain
    if (
        not (wrf_xmin > lcz_xmin)
        & (wrf_xmax < lcz_xmax)
        & (wrf_ymin > lcz_ymin)
        & (wrf_ymax < lcz_ymax)
    ):

        message = (
            f'{ERROR}ERROR: LCZ domain should be larger than '
            'WRF domain in all directions.\n'
            'LCZ bounds (xmin, ymin, xmax, ymax): '
            f'{lcz_xmin, lcz_ymin, lcz_xmax, lcz_ymax}\n'
            'WRF bounds (xmin, ymin, xmax, ymax): '
            f'{wrf_xmin, wrf_ymin, wrf_xmax, wrf_ymax}{ENDC}'
        )
        print(message)
        sys.exit(1)
    else:
        print('> LCZ domain is covering WRF domain')


def check_lcz_integrity(info: Info, LCZ_BAND: int) -> None:

    '''
    Check integrity of LCZ GeoTIFF file, which should have:
    - LCZ class labels between 1 and 17
    - WGS84 (EPSG:4326) projection
    - extent larger than WRF domain file, in all direction

    Note that this procedure directly select correct band in GeoTIFF file,
    so no need to select that specific band afterwards.

    Output:
        _clean.tif lcz file, used in the remainder of the tool.
    '''

    # Read the data
    lcz = rxr.open_rasterio(info.src_file)[LCZ_BAND, :, :]
    wrf = xr.open_dataset(info.dst_file)

    # If any of [101, 102, 103, 104, 105, 106, 107] is in the lcz tif file.
    lcz_100 = np.array(
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 101, 102, 103, 104, 105, 106, 107]
    )
    if any(x in lcz.data.flatten() for x in lcz_100[10:]):
        lcz = _replace_lcz_number(
            lcz=lcz,
            lcz_to_change=lcz_100,
        )
        print('> LCZ class labels renamed from 1 to 17')
    else:
        print('> LCZ labels as expected (1 to 17)')

    # Re-project when not WGS84 (EPSG:4326)
    if lcz.rio.crs != CRS.from_epsg(4326):
        lcz = lcz.rio.reproject('EPSG:4326')
        lcz.data = xr.where(lcz.data > 0, lcz.data, 0)
        lcz.data = lcz.data.astype(np.int32)
        print('> LCZ map reprojected to WGS84 (EPSG:4326)')
    else:
        print('> LCZ provided as WGS84 (EPSG:4326)')

    # Check if LCZ map exceeds domain of geo_em file in all directions
    _check_lcz_wrf_extent(lcz, wrf)

    # Write clean LCZ to file, used in all subsequent routines.
    lcz.rio.to_raster(info.src_file_clean, dtype=np.int32)


def _calc_distance_coord(
    lat1: float, lon1: float, lat2: float, lon2: float
) -> xr.DataArray:

    '''Calculate distance using coordinates
    This uses the spherical law of cosines
    '''

    earth_radius = 6371000  # Earth radius in m
    lat1r = lat1 * np.pi / 180.0
    lat2r = lat2 * np.pi / 180.0
    lon1r = lon1 * np.pi / 180.0
    lon2r = lon2 * np.pi / 180.0

    d = (
        np.arccos(
            np.sin(lat1r) * np.sin(lat2r)
            + np.cos(lat1r) * np.cos(lat2r) * np.cos(lon2r - lon1r)
        )
        * earth_radius
    )

    return d


def wrf_remove_urban(
    info: Info,
    NPIX_NLC: int,
) -> None:

    '''Remove MODIS urban extent from geo_em*.nc file'''

    # Make a copy of original dst file
    dst_data = xr.open_dataset(info.dst_file)

    # Read the relevant parameters
    luse = dst_data.LU_INDEX.squeeze()
    luf = dst_data.LANDUSEF.squeeze()
    greenf = dst_data.GREENFRAC.squeeze()
    lat = dst_data.XLAT_M.squeeze()
    lon = dst_data.XLONG_M.squeeze()
    newluse = luse.values.copy()
    newluf = luf.values.copy()
    newgreenf = greenf.values.copy()
    orig_num_land_cat = dst_data.NUM_LAND_CAT

    # Convert urban to surrounding natural characteristics
    for i in dst_data.south_north:
        for j in dst_data.west_east:
            if luse.isel(south_north=i, west_east=j) == 13:

                dis = _calc_distance_coord(
                    lat.where((luse != 13) & (luse != 17) & (luse != 21)),
                    lon.where((luse != 13) & (luse != 17) & (luse != 21)),
                    lat.isel(south_north=i, west_east=j),
                    lon.isel(south_north=i, west_east=j),
                )

                disflat = (
                    dis.stack(gridpoints=('south_north', 'west_east'))
                    .reset_index('gridpoints')
                    .drop_vars(['south_north', 'west_east'])
                )
                aux = luse.where(
                    dis <= disflat.sortby(disflat).isel(gridpoints=NPIX_NLC), drop=True
                )
                m = stats.mode(aux.values.flatten(), nan_policy='omit')[0]
                newluse[i, j] = int(m)

                auxg = (
                    greenf.where(
                        dis <= disflat.sortby(disflat).isel(gridpoints=NPIX_NLC),
                        drop=True,
                    )
                    .where(aux == newluse[i, j])
                    .mean(dim=['south_north', 'west_east'])
                )
                newgreenf[:, i, j] = auxg

            if luf.isel(south_north=i, west_east=j, land_cat=12) > 0.0:
                if orig_num_land_cat > 20:  # USING MODIS_LAKE
                    dis = _calc_distance_coord(
                        lat.where(
                            (luf.isel(land_cat=12) == 0.0)
                            & (luf.isel(land_cat=16) == 0.0)
                            & (luf.isel(land_cat=20) == 0.0)
                        ),
                        lon.where(
                            (luf.isel(land_cat=12) == 0.0)
                            & (luf.isel(land_cat=16) == 0.0)
                            & (luf.isel(land_cat=20) == 0.0)
                        ),
                        lat.isel(south_north=i, west_east=j),
                        lon.isel(south_north=i, west_east=j),
                    )
                else:  # USING MODIS (NO LAKES)
                    dis = _calc_distance_coord(
                        lat.where(
                            (luf.isel(land_cat=12) == 0.0)
                            & (luf.isel(land_cat=16) == 0.0)
                        ),
                        lon.where(
                            (luf.isel(land_cat=12) == 0.0)
                            & (luf.isel(land_cat=16) == 0.0)
                        ),
                        lat.isel(south_north=i, west_east=j),
                        lon.isel(south_north=i, west_east=j),
                    )

                disflat = (
                    dis.stack(gridpoints=('south_north', 'west_east'))
                    .reset_index('gridpoints')
                    .drop_vars(['south_north', 'west_east'])
                )
                aux = luse.where(
                    dis <= disflat.sortby(disflat).isel(gridpoints=NPIX_NLC), drop=True
                )
                m = stats.mode(aux.values.flatten(), nan_policy='omit')[0]
                newlu = int(m) - 1
                # newlu = int(mode(aux.values.flatten())[0])-1
                newluf[newlu, i, j] += luf.isel(
                    south_north=i, west_east=j, land_cat=12
                ).values
                newluf[12, i, j] = 0.0

    dst_data.LU_INDEX.values[0, :] = newluse[:]
    dst_data.LANDUSEF.values[0, :] = newluf[:]
    dst_data.GREENFRAC.values[0, :] = newgreenf[:]

    # Save to final _lcz_params file
    if os.path.exists(info.dst_nu_file):
        os.remove(info.dst_nu_file)
    dst_data.to_netcdf(info.dst_nu_file)


# Make WRF grid info available for Resampler (tmp file)
def create_wrf_gridinfo(info: Info) -> None:

    # Read  gridded WRF data
    dst_data = xr.open_dataset(info.dst_nu_file)

    # Create simpler WRF grid target.
    da_lu = xr.Dataset(
        {'LU_INDEX': (['y', 'x'], dst_data['LU_INDEX'][0, :, :].values)},
        coords={
            'y': dst_data.XLAT_M.values[0, :, 0],
            'x': dst_data.XLONG_M.values[0, 0, :],
        },
    )

    # Add projection information as attributes, save and read back in.
    da_lu.rio.write_crs('epsg:4326', inplace=True)
    da_lu.rio.to_raster(info.dst_gridinfo)


def _get_SW_BW(ucp_table: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:

    '''Get Street and Building Width'''

    # Street width extracted from S02012 Building heighht and H2W.
    SW = ucp_table['MH_URB2D'] / ucp_table['H2W']
    # Building Width according to bldfr_urb2d/(frc_urb2d-bldfr_urb2d)*sw
    BW = (
        ucp_table['BLDFR_URB2D'] / (ucp_table['FRC_URB2D'] - ucp_table['BLDFR_URB2D'])
    ) * SW

    return SW, BW


def _get_lcz_arr(src_data: xr.DataArray, info: Info) -> NDArray[np.int_]:

    '''Get LCZ data as array, setting non-built pixels to 0'''

    # Get mask of selected built LCZs
    lcz_urb_mask = xr.DataArray(
        np.in1d(src_data, info.BUILT_LCZ).reshape(src_data.shape),
        dims=src_data.dims,
        coords=src_data.coords,
    )

    # Get LCZ class values only.
    lcz_arr = src_data.data.astype(np.int32)

    # Set LCZ classes not in BUILT_LCZ to 0
    lcz_arr[~lcz_urb_mask] = 0

    return lcz_arr


def _ucp_resampler(
    info: Info,
    ucp_key: str,
    RESAMPLE_TYPE: str,
    ucp_table: pd.DataFrame,
    **kwargs: float,
) -> xr.DataArray:

    '''Helper function to resample lcz ucp data ('FRC_URB2D', 'MH_URB2D',
    'STDH_URB2D', 'LB_URB2D', 'LF_URB2D', 'LP_URB2D') to WRF grid'''

    # Read gridded data: LCZ and WRF grid
    src_data = rxr.open_rasterio(info.src_file_clean)[0, :, :]
    dst_grid = rxr.open_rasterio(info.dst_gridinfo)

    # Get Street and Building Width
    SW, BW = _get_SW_BW(ucp_table)

    # Get Look-up for FRC_values
    if ucp_key in ['LB_URB2D', 'LF_URB2D', 'LP_URB2D']:

        # Following Zonato et al (2020)
        LAMBDA_P = BW / (BW + SW)
        LAMBDA_F = 2 * ucp_table['MH_URB2D'] / (BW + SW)
        LAMBDA_B = LAMBDA_P + LAMBDA_F

        if ucp_key == 'LB_URB2D':
            lookup = LAMBDA_B.loc[info.BUILT_LCZ]
        elif ucp_key == 'LP_URB2D':
            lookup = LAMBDA_P.loc[info.BUILT_LCZ]
        elif ucp_key == 'LF_URB2D':
            lookup = LAMBDA_F.loc[info.BUILT_LCZ]

    elif ucp_key == 'STDH_URB2D':
        lookup = ((ucp_table['MH_URB2D_MAX'] - ucp_table['MH_URB2D_MIN']) / 4).loc[
            info.BUILT_LCZ
        ]
    else:
        lookup = ucp_table[ucp_key].loc[info.BUILT_LCZ]

    # Get mask of selected built LCZs
    lcz_arr = _get_lcz_arr(src_data, info)

    # Make replacer object to map UCP values on LCZ class values
    replacer = np.zeros((max(info.BUILT_LCZ) + 1,), object)
    replacer[lookup.index.values] = lookup
    lcz_data = np.array(replacer[lcz_arr], dtype='float')

    # Store into dataarray for resampling
    lcz_data_da = xr.Dataset(
        {'band': (['y', 'x'], lcz_data)},
        coords={'y': src_data.y.values, 'x': src_data.x.values},
        attrs={'transform': src_data.rio.transform(), 'crs': src_data.rio.crs},
    ).to_array()

    # Info: https://rasterio.readthedocs.io/en/latest/api/rasterio.warp.html?highlight=reproject(#rasterio.warp.reproject
    ucp_2_wrf = reproject(
        lcz_data_da,
        dst_grid,
        src_transform=lcz_data_da.rio.transform(),
        src_crs=lcz_data_da.rio.crs,
        dst_transform=dst_grid.rio.transform(),
        dst_crs=dst_grid.rio.crs,
        resampling=Resampling[RESAMPLE_TYPE],
    )[0]

    # In case of FRC_URB2D, filter for too low values
    if 'FRC_THRESHOLD' in kwargs.keys():
        ucp_2_wrf = ucp_2_wrf.where(ucp_2_wrf > kwargs['FRC_THRESHOLD'], 0)

    ## In case nans occur, set to zero
    ucp_2_wrf.values[0, np.isnan(ucp_2_wrf[0, :, :])] = 0

    return ucp_2_wrf


def _hgt_resampler(
    info: Info,
    RESAMPLE_TYPE: str,
    ucp_table: pd.DataFrame,
) -> xr.DataArray:

    '''Helper function to resample HGT_URB2D (=Area Weighted
    Mean Building Height ) data to WRF grid'''

    # Read gridded data: LCZ and WRF grid
    src_data = rxr.open_rasterio(info.src_file_clean)[0, :, :]
    dst_grid = rxr.open_rasterio(info.dst_gridinfo)

    # Street width extracted from S02012 Building heighht and H2W.
    SW, BW = _get_SW_BW(ucp_table)

    # Get Look-up for HGT values
    lookup_nom = BW.loc[info.BUILT_LCZ] ** 2 * ucp_table['MH_URB2D'].loc[info.BUILT_LCZ]
    lookup_denom = BW.loc[info.BUILT_LCZ] ** 2

    # Get mask of selected built LCZs
    lcz_arr = _get_lcz_arr(src_data, info)

    # Make replacer object for nominator
    replacer_nom = np.zeros((max(info.BUILT_LCZ) + 1,), object)
    replacer_nom[lookup_nom.index.values] = lookup_nom
    dataLcz_nom = np.array(replacer_nom[lcz_arr], dtype='float')

    # Make replacer object for denominator
    replacer_denom = np.zeros((max(info.BUILT_LCZ) + 1,), object)
    replacer_denom[lookup_denom.index.values] = lookup_denom
    dataLcz_denom = np.array(replacer_denom[lcz_arr], dtype='float')

    # Store into dataarray for resampling
    lcz_data_da_nom = xr.Dataset(
        {'band': (['y', 'x'], dataLcz_nom)},
        coords={'y': src_data.y.values, 'x': src_data.x.values},
        attrs={'transform': src_data.rio.transform(), 'crs': src_data.rio.crs},
    ).to_array()
    lcz_data_da_denom = xr.Dataset(
        {'band': (['y', 'x'], dataLcz_denom)},
        coords={'y': src_data.y.values, 'x': src_data.x.values},
        attrs={'transform': src_data.rio.transform(), 'crs': src_data.rio.crs},
    ).to_array()

    # Get the aggregated values on WRF grid - nominator
    ucp_2_wrf_nom = reproject(
        lcz_data_da_nom,
        dst_grid,
        src_transform=lcz_data_da_nom.rio.transform(),
        src_crs=lcz_data_da_nom.crs,
        dst_transform=dst_grid.rio.transform(),
        dst_crs=dst_grid.rio.crs,
        resampling=Resampling[RESAMPLE_TYPE],
    )[0].copy()

    # Get the aggregated values on WRF grid - nominator
    ucp_2_wrf_denom = reproject(
        lcz_data_da_denom,
        dst_grid,
        src_transform=lcz_data_da_denom.rio.transform(),
        src_crs=lcz_data_da_denom.crs,
        dst_transform=dst_grid.rio.transform(),
        dst_crs=dst_grid.rio.crs,
        resampling=Resampling[RESAMPLE_TYPE],
    )[0].copy()

    hgt_urb2d = ucp_2_wrf_nom / ucp_2_wrf_denom

    ## In case nans occur, set to zero
    hgt_urb2d.values[0, np.isnan(hgt_urb2d[0, :, :])] = 0

    return hgt_urb2d


def _get_truncated_normal_sample(
    lcz_i: int, ucp_table: pd.DataFrame, SAMPLE_SIZE: int = 100000
) -> NDArray[np.float_]:

    '''Helper function to return bounded normal distribution sample'''

    # Create instance of a truncated normal distribution
    low = ucp_table['MH_URB2D_MIN'].loc[lcz_i]
    mean = ucp_table['MH_URB2D'].loc[lcz_i]
    upp = ucp_table['MH_URB2D_MAX'].loc[lcz_i]
    sd = (
        ucp_table['MH_URB2D_MAX'].loc[lcz_i] - ucp_table['MH_URB2D_MIN'].loc[lcz_i]
    ) / 4

    hi_inst = truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

    # populate with large enough sample for accuracy
    hi_sample = hi_inst.rvs(SAMPLE_SIZE)

    return hi_sample


def _check_hi_values(
    lcz_i: int,
    hi_sample: NDArray[np.float_],
    ucp_table: pd.DataFrame,
    ERROR_MARGIN: float,
) -> None:
    hi_metrics = [
        'MH_URB2D_MIN',
        'MH_URB2D_MAX',
        'MH_URB2D',
    ]
    # Produce warning if approximated HI_URB2D distribution metrics
    # are not as expected: using a ERROR_MARGIN % marging here.
    for hi_metric in hi_metrics:

        if hi_metric == 'MH_URB2D_MIN':
            hi_sample_values = hi_sample.min()
        elif hi_metric == 'MH_URB2D_MAX':
            hi_sample_values = hi_sample.max()
        elif hi_metric == 'MH_URB2D':
            hi_sample_values = hi_sample.mean()

        if (
            not ucp_table[hi_metric].loc[lcz_i] * (1 - ERROR_MARGIN)
            < hi_sample_values
            < ucp_table[hi_metric].loc[lcz_i] * (1 + ERROR_MARGIN)
        ):
            print(
                f'WARNING: {hi_metric} distribution not in '
                f'expected range ({ERROR_MARGIN*100}% marging) for LCZ class {lcz_i}: '
                f'modelled: {np.round(hi_sample_values, 2)} | '
                f'expected: ['
                f'{np.round(ucp_table[hi_metric].loc[lcz_i] * (1 - ERROR_MARGIN),2)} - '
                f'{np.round(ucp_table[hi_metric].loc[lcz_i] * (1 - ERROR_MARGIN),2)}]'
            )


def _compute_hi_distribution(
    info: Info,
    ucp_table: pd.DataFrame,
    SAMPLE_SIZE: int = 100000,
    ERROR_MARGIN: float = 0.05,
) -> pd.DataFrame:

    '''Helper function to compute building height distribution'''

    # Initialize dataframe that stores building height distributions
    df_hi = pd.DataFrame(
        index=range(1, 18, 1),
        columns=[
            '0 - <5m',
            '5 - <10m',
            '10 - <15m',
            '15 - <20m',
            '20 - <25m',
            '25 - <30m',
            '30 - <35m',
            '35 - <40m',
            '40 - <45m',
            '45 - <50m',
            '50 - <55m',
            '55 - <60m',
            '60 - <65m',
            '65 - <70m',
            '70 - <75m',
        ],
    )

    for i in info.BUILT_LCZ:

        # LCZ 15 = paved, and considered to have no buildings (values = 0%)
        if not i == 15:

            # Make the sample
            hi_sample = _get_truncated_normal_sample(
                lcz_i=i,
                ucp_table=ucp_table,
                SAMPLE_SIZE=SAMPLE_SIZE,
            )

            # Check if values are within expected bounds.
            _check_hi_values(
                lcz_i=i,
                hi_sample=hi_sample,
                ucp_table=ucp_table,
                ERROR_MARGIN=ERROR_MARGIN,
            )

            # Count the values within pre-set bins
            count_bins = np.histogram(hi_sample, bins=np.arange(0, 76, 5))[0]
            count_bins = count_bins / (SAMPLE_SIZE / 100)  # Convert to %

            # Add to dataframe
            df_hi.loc[i, :] = count_bins

        # Set nans to zero
        df_hi = df_hi.fillna(0)

    return df_hi


def _scale_hi(array: NDArray[np.int_]) -> List[float]:

    '''Helper function to scale HI_URB2D to 100%'''

    scaled_array = [(float(i) / sum(array) * 100.0) for i in array]
    return scaled_array


def _hi_resampler(
    info: Info,
    RESAMPLE_TYPE: str,
    ucp_table: pd.DataFrame,
    HI_THRES_MIN: int = 5,
) -> Tuple[NDArray[np.float_], float]:

    '''Helper function to resample ucp HI_URB2D_URB2D data to WRF grid'''

    # Read gridded data: LCZ and WRF grid
    src_data = rxr.open_rasterio(info.src_file_clean)[0, :, :]
    dst_grid = rxr.open_rasterio(info.dst_gridinfo)

    # Get mask of selected built LCZs
    lcz_arr = _get_lcz_arr(src_data, info)

    # Compute the building height densities.
    df_hi = _compute_hi_distribution(info, ucp_table=ucp_table)

    # Initialize array to store temp values
    hi_arr = np.zeros((15, dst_grid.shape[1], dst_grid.shape[2]))

    # Loop over the 15 height density classes.
    for hi_i in range(df_hi.shape[1]):

        # print(f"Working on height interval {df_hi.columns[hi_i]} ...")
        lookup = df_hi.iloc[:, hi_i].loc[info.BUILT_LCZ]

        # Make replacer object to map UCP values on LCZ class values
        replacer = np.zeros((max(info.BUILT_LCZ) + 1,), object)
        replacer[lookup.index.values] = lookup
        lcz_data = np.array(replacer[lcz_arr], dtype='float')

        # Store into dataarray for resampling
        lcz_data_da = xr.Dataset(
            {'band': (['y', 'x'], lcz_data)},
            coords={'y': src_data.y.values, 'x': src_data.x.values},
            attrs={'transform': src_data.rio.transform(), 'crs': src_data.rio.crs},
        ).to_array()

        # Get the aggregated values on WRF grid
        ucp_2_wrf = reproject(
            lcz_data_da,
            dst_grid,
            src_transform=lcz_data_da.rio.transform(),
            src_crs=lcz_data_da.rio.crs,
            dst_transform=dst_grid.rio.transform(),
            dst_crs=dst_grid.rio.crs,
            resampling=Resampling[RESAMPLE_TYPE],
        )[0]

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
    info: Info,
    frc_urb2d: xr.DataArray,
    LCZ_NAT_MASK: bool,
) -> Tuple[NDArray[np.bool_], NDArray[np.float_]]:

    '''Helper function to resample lcz classes to WRF grid using majority'''

    # Read required gridded data, LCZ, WRF grid, and
    # original WRF (for original MODIS urban mask)
    src_data = rxr.open_rasterio(info.src_file_clean)[0, :, :]
    dst_grid = rxr.open_rasterio(info.dst_gridinfo)

    # Mask natural LCZs before majority filtering.
    if LCZ_NAT_MASK:
        src_data = src_data.where(src_data.isin(info.BUILT_LCZ)).copy()

    lcz_2_wrf = reproject(
        src_data,
        dst_grid,
        src_transform=src_data.rio.transform(),
        src_crs=src_data.rio.crs,
        dst_transform=dst_grid.rio.transform(),
        dst_crs=dst_grid.rio.crs,
        resampling=Resampling['mode'],
    )[0].values

    # if LCZ 15 selected in 'BUILT_LCZ', rename to 11
    if 15 in info.BUILT_LCZ:
        lcz_2_wrf[lcz_2_wrf == 15] = 11

    # Only keep LCZ pixels where FRC_URB2D > 0, for concistency
    frc_mask = frc_urb2d.values[0, :, :] != 0

    # Final LU_INDEX = 31 to 41 (included), as LCZ classes.
    lcz_resampled = lcz_2_wrf[0, frc_mask] + 30

    return frc_mask, lcz_resampled


def _adjust_greenfrac_landusef(
    info: Info,
    dst_data: xr.Dataset,
    frc_mask: NDArray[np.bool_],
) -> xr.Dataset:

    dst_data_orig = xr.open_dataset(info.dst_file)
    orig_num_land_cat = dst_data_orig.NUM_LAND_CAT

    # Adjust GREENFRAC and LANDUSEF
    # GREENFRAC is set as average / month from GREENFRAC
    # of original MODIS urban pixels
    wrf_urb = xr.DataArray(
        np.in1d(dst_data_orig['LU_INDEX'][0, :, :].values, [13]).reshape(
            dst_data_orig['LU_INDEX'][0, :, :].shape
        ),
        dims=dst_data_orig['LU_INDEX'][0, :, :].dims,
        coords=dst_data_orig['LU_INDEX'][0, :, :].coords,
    )
    greenfrac_per_month = [
        dst_data_orig['GREENFRAC'].values[0, mm, wrf_urb].mean() for mm in range(12)
    ]

    # Loop over months and set average values
    for mm in range(12):
        dst_data['GREENFRAC'].values[0, mm, frc_mask] = greenfrac_per_month[mm]

    # TODO: For lower resolution domains, this might not be valid?
    # Create new LANDUSEF with 41 levels instead of 21
    landusef_new = np.zeros(
        (41, dst_data.LANDUSEF.shape[2], dst_data.LANDUSEF.shape[3])
    )

    # Copy values from original file
    landusef_new[:orig_num_land_cat, :, :] = dst_data['LANDUSEF'][
        0, :orig_num_land_cat, :, :
    ]

    # First set all values to zero for urban mask
    landusef_new[:, frc_mask] = 0  # First all to 0, so sum remains 1 in the end

    # LOOP over LCZ LU_INDEX values, and set to 1 there
    # So e.g. LANDUSE[0,31-1,:,1] = 1, where LU_INDEX = 31 (=LCZ 1)
    for lu_i in np.arange(31, 42, 1):
        lu_mask = dst_data.LU_INDEX == int(lu_i)
        landusef_new[int(lu_i) - 1, lu_mask[0, :, :]] = 1
        del lu_mask

    # First store orginal attributes, then drop variable
    luf_attrs = dst_data.LANDUSEF.attrs
    dst_data = dst_data.drop_vars('LANDUSEF')

    # Expand axis to take shape (1,41,x,y)
    landusef_new = np.expand_dims(landusef_new, axis=0)

    # Add back to data-array, including (altered) attributes
    dst_data['LANDUSEF'] = (
        ('Time', 'land_cat', 'south_north', 'west_east'),
        landusef_new,
    )
    dst_data['LANDUSEF'] = dst_data.LANDUSEF.astype('float32')

    luf_attrs['description'] = 'Noah-modified 41-category IGBP-MODIS landuse'
    for key in luf_attrs.keys():
        dst_data['LANDUSEF'].attrs[key] = luf_attrs[key]

    return dst_data


def _add_frc_lu_index_2_wrf(
    info: Info,
    FRC_THRESHOLD: float,
    LCZ_NAT_MASK: bool,
    ucp_table: pd.DataFrame,
) -> xr.Dataset:

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
        ucp_table=ucp_table,
        FRC_THRESHOLD=FRC_THRESHOLD,
    )

    # Add to geo_em* that that has no MODIS urban
    dst_data = xr.open_dataset(info.dst_nu_file)

    # Make a FRC_URB field and store aggregated data.
    dst_data[ucp_key] = dst_data['LU_INDEX'].copy()
    dst_data[ucp_key] = (('Time', 'south_north', 'west_east'), frc_urb.data)

    # Add proper attributes to the FRC_URB2D field
    dst_data[ucp_key].attrs['FieldType'] = np.intc(104)
    dst_data[ucp_key].attrs['MemoryOrder'] = 'XY'
    dst_data[ucp_key].attrs['units'] = '-'
    dst_data[ucp_key].attrs['description'] = 'ufrac'
    dst_data[ucp_key].attrs['stagger'] = 'M'
    dst_data[ucp_key].attrs['sr_x'] = np.intc(1)
    dst_data[ucp_key].attrs['sr_y'] = np.intc(1)

    # Integrate LU_INDEX, also adjusts GREENFRAC and LANDUSEF
    frc_mask, lcz_resampled = _lcz_resampler(
        info=info,
        frc_urb2d=dst_data['FRC_URB2D'],
        LCZ_NAT_MASK=LCZ_NAT_MASK,
    )

    # 2) as LU_INDEX = 30 to 41, as LCZ classes.
    dst_data['LU_INDEX'].values[0, frc_mask] = lcz_resampled

    # Also adjust GREENFRAC and LANDUSEF
    dst_data = _adjust_greenfrac_landusef(info, dst_data, frc_mask)

    return dst_data


def _initialize_urb_param(dst_data: xr.Dataset) -> xr.Dataset:

    '''Helper function to initialize URB_PARAM in WRF geo_em file'''

    URB_PARAM = np.zeros([1, 132, len(dst_data.south_north), len(dst_data.west_east)])

    # Add to destination WRF file, with attributes
    dst_data['URB_PARAM'] = (
        ('Time', 'num_urb_params', 'south_north', 'west_east'),
        URB_PARAM,
    )
    dst_data['URB_PARAM'].attrs['FieldType'] = np.intc(104)
    dst_data['URB_PARAM'].attrs['MemoryOrder'] = 'XYZ'
    dst_data['URB_PARAM'].attrs['units'] = 'dimensionless'
    dst_data['URB_PARAM'].attrs['description'] = 'all urban parameters'
    dst_data['URB_PARAM'].attrs['stagger'] = 'M'
    dst_data['URB_PARAM'].attrs['sr_x'] = np.intc(1)
    dst_data['URB_PARAM'].attrs['sr_y'] = np.intc(1)

    return dst_data


def create_lcz_params_file(
    info: Info,
    FRC_THRESHOLD: float,
    LCZ_NAT_MASK: bool,
    ucp_table: pd.DataFrame,
) -> float:

    '''
    Create a domain file with all LCZ-based information:
    Map, aggregate and add lcz-based UCP values to the inner geo_em file.
    '''

    dst_data = _add_frc_lu_index_2_wrf(
        info=info,
        FRC_THRESHOLD=FRC_THRESHOLD,
        LCZ_NAT_MASK=LCZ_NAT_MASK,
        ucp_table=ucp_table,
    )

    # Initialize empty URB_PARAM in final wrf file,
    # with all zeros and proper attributes
    dst_final = _initialize_urb_param(dst_data=dst_data)

    # get frc_mask, to only set values where FRC_URB2D > 0.
    frc_mask = dst_final.FRC_URB2D.values[0, :, :] != 0

    # Define the UCPs that need to be integrated,
    # together with their positions (index starts at 1) in URB_PARAMS
    # HGT_URB2D and HI_URB2D follow a different approach, see further.
    ucp_dict = {
        'LP_URB2D': 91,
        'MH_URB2D': 92,
        'STDH_URB2D': 93,
        'HGT_URB2D': 94,
        'LB_URB2D': 95,
        'LF_URB2D': 96,  # 97, 98, 99, for all 4 directions
        'HI_URB2D': 118,  # Goes on until index 132
    }

    for ucp_key in ucp_dict.keys():

        print(f'> Processing {ucp_key} ...')

        # Obtain aggregated LCZ-based UCP values
        if ucp_key in ['MH_URB2D', 'STDH_URB2D', 'LB_URB2D', 'LF_URB2D', 'LP_URB2D']:
            ucp_res = _ucp_resampler(
                info=info,
                ucp_key=ucp_key,
                RESAMPLE_TYPE='average',
                ucp_table=ucp_table,
            )
        elif ucp_key in ['HGT_URB2D']:
            ucp_res = _hgt_resampler(
                info=info,
                RESAMPLE_TYPE='average',
                ucp_table=ucp_table,
            )
        elif ucp_key in ['HI_URB2D']:
            ucp_res, nbui_max = _hi_resampler(
                info=info,
                RESAMPLE_TYPE='average',
                ucp_table=ucp_table,
            )

        # Store UCPs in wrf destination file.
        if ucp_key == 'LF_URB2D':
            # Frontal area Index in N,E,S,W directions respectively
            # for WUDAPT LCZs they are considered all equal
            for i in range(4):
                ucp_res.values[:, frc_mask == 0] = 0
                dst_final['URB_PARAM'][:, (ucp_dict[ucp_key] - 1 + i), :, :] = ucp_res
        if ucp_key == 'HI_URB2D':
            ucp_res[:, frc_mask == 0] = 0
            dst_final['URB_PARAM'][0, (ucp_dict[ucp_key] - 1) :, :, :] = ucp_res
        else:
            ucp_res.values[:, frc_mask == 0] = 0
            dst_final['URB_PARAM'].values[0, (ucp_dict[ucp_key] - 1), :, :] = ucp_res

    # Make sure URB_PARAM is float32
    dst_final['URB_PARAM'] = dst_final.URB_PARAM.astype('float32')

    # Expand attribute title
    att_title = dst_final.attrs['TITLE']
    dst_final.attrs['TITLE'] = f'{att_title}, perturbed by W2W'

    # Add/Change some additional global attributes,
    # including NBUI_MAX = max. nr. of HI intervals over the grid
    glob_attrs: Dict[str, Union[int, SupportsInt]] = {
        'NUM_LAND_CAT': 41,
        'FLAG_URB_PARAM': 1,
        'NBUI_MAX': np.intc(nbui_max),
    }
    for key in glob_attrs.keys():
        dst_final.attrs[key] = np.intc(glob_attrs[key])

    # Add DESCRIPTION in attrs, referring to tool
    gh_ref = (
        'Demuzere, M., Argüeso, D., Zonato, A., & Kittner, J. (2021). \n'
        "W2W: A Python package that injects WUDAPT's Local Climate Zone \n"
        'information in WRF [Computer software]. \n'
        'https://github.com/matthiasdemuzere/w2w'
    )
    dst_final.attrs[
        'DESCRIPTION'
    ] = f'W2W.py tool used to create geo_em*.nc file:\n {gh_ref}'

    # Save back to file
    if os.path.exists(info.dst_lcz_params_file):
        os.remove(info.dst_lcz_params_file)
    dst_final.to_netcdf(info.dst_lcz_params_file)

    return nbui_max


def create_lcz_extent_file(info: Info) -> None:

    '''
    Create a domain file with an LCZ-based urban extent
    (excluding other LCZ-based info)
    '''

    dst_params = xr.open_dataset(info.dst_lcz_params_file)
    frc_mask = dst_params.FRC_URB2D.values[0, :, :] != 0

    dst_extent = dst_params.copy()

    dst_data_orig = xr.open_dataset(info.dst_file)
    orig_num_land_cat = dst_data_orig.NUM_LAND_CAT
    orig_luf_description = dst_data_orig.LANDUSEF.description

    lu_index = dst_extent.LU_INDEX.values
    lu_index[lu_index >= 31] = 13

    dst_extent.LU_INDEX.values = lu_index

    # Remove some unnecesary variables to reduce file size
    dst_extent = dst_extent.drop_vars(['FRC_URB2D', 'URB_PARAM'])

    # Reset LANDUSEF again to 21 classes.
    luf_attrs = dst_extent.LANDUSEF.attrs
    luf_values = dst_extent.LANDUSEF.values
    dst_extent = dst_extent.drop_vars('LANDUSEF')

    # Add back to data-array, including (altered) attributes
    dst_extent['LANDUSEF'] = (
        ('Time', 'land_cat', 'south_north', 'west_east'),
        luf_values[:, :orig_num_land_cat, :, :],
    )
    dst_extent['LANDUSEF'].values[0, 12, frc_mask] = 1
    dst_extent['LANDUSEF'] = dst_extent.LANDUSEF.astype('float32')

    luf_attrs['description'] = orig_luf_description
    for key in luf_attrs.keys():
        dst_extent['LANDUSEF'].attrs[key] = luf_attrs[key]

    # Reset some other global attributes
    dst_extent.attrs['FLAG_URB_PARAM'] = np.intc(0)
    dst_extent.attrs['NUM_LAND_CAT'] = np.intc(orig_num_land_cat)

    # Save file.
    dst_extent.to_netcdf(info.dst_lcz_extent_file)


def expand_land_cat_parents(info: Info) -> None:

    # Get final domain number
    domain_nr = int(info.dst_file[-5:-3])

    # list domain numbers to loop over
    domain_lst = list(np.arange(1, domain_nr, 1))

    for i in domain_lst:

        ifile = f'{info.dst_file[:-5]}{i:02d}.nc'

        if os.path.exists(ifile):

            da = xr.open_dataset(ifile)

            try:
                if int(da.attrs['NUM_LAND_CAT']) != 41:

                    orig_num_land_cat = da.attrs['NUM_LAND_CAT']
                    # Set number of land categories to 41
                    da.attrs['NUM_LAND_CAT'] = np.intc(41)

                    # Create new landusef array with expanded dimensions
                    landusef_new = np.zeros(
                        (1, 41, da.LANDUSEF.shape[2], da.LANDUSEF.shape[3])
                    )
                    landusef_new[:, :orig_num_land_cat, :, :] = da['LANDUSEF'].values

                    # First store orginal attributes, then drop variable
                    luf_attrs = da.LANDUSEF.attrs
                    da = da.drop_vars('LANDUSEF')

                    # Add back to data-array, including (altered) attributes
                    da['LANDUSEF'] = (
                        ('Time', 'land_cat', 'south_north', 'west_east'),
                        landusef_new,
                    )
                    da['LANDUSEF'] = da.LANDUSEF.astype('float32')

                    luf_attrs['description'] = (
                        'Noah-modified 41-category ' 'IGBP-MODIS landuse'
                    )
                    for key in luf_attrs.keys():
                        da['LANDUSEF'].attrs[key] = luf_attrs[key]

                    ofile = ifile.replace('.nc', '_41.nc')
                    da.to_netcdf(ofile)

                else:
                    print(f'> Parent domain d{i:02d}.nc already contains 41 LC classes')
            except Exception:
                err = traceback.format_exc()
                print(f'Cannot read NUM_LAND_CAT and LANDUSEF dimensions\n{err}')

        else:
            print(
                f'WARNING: Parent domain {info.dst_file[:-5]}{i:02d}.nc'
                f' not found.\n'
                f'Please make sure the parent domain files are '
                f'in {info.io_dir}\n'
                f'Without this information, you will not be able to produce '
                f'the boundary conditions with real.exe.'
            )


def checks_and_cleaning(info: Info, ucp_table: pd.DataFrame, nbui_max: float) -> None:

    'Sanity checks and cleaning'

    OKGREEN = '\033[0;32m'
    WARNING = '\033[0;35m'
    ENDC = '\033[0m'

    base_text = (
        f'> Check 1: Urban class removed from ' f"{info.dst_nu_file.split('/')[-1]}?"
    )
    ifile = info.dst_nu_file
    da = xr.open_dataset(ifile)
    if 13 in da.LU_INDEX.values:
        print(
            f'{base_text}\n{WARNING} WARNING: Urban land use ' f'still present {ENDC}'
        )
    else:
        print(f'{base_text}{OKGREEN} OK {ENDC}')

    base_text = (
        f'> Check 2: LCZ Urban extent present in '
        f"{info.dst_lcz_extent_file.split('/')[-1]}?"
    )
    ifile = info.dst_lcz_extent_file
    da = xr.open_dataset(ifile)
    if not 13 in da.LU_INDEX.values:
        print(
            f'{base_text}\n{WARNING} WARNING: LCZ-based urban ' f'extent missing {ENDC}'
        )
    else:
        print(f'{base_text}{OKGREEN} OK {ENDC}')

    base_text = (
        f'> Check 3: Urban LCZ classes exists in '
        f"{info.dst_lcz_params_file.split('/')[-1]}?"
    )
    ifile = info.dst_lcz_params_file
    da = xr.open_dataset(ifile)
    if 13 in da.LU_INDEX.values:
        print(
            f'{base_text}\n{WARNING} WARNING: Urban extent still '
            f'defined via LU_INDEX = 13? {ENDC}'
        )
    else:
        LU_values = np.unique(da.LU_INDEX.values.flatten())
        LCZs = [int(i) for i in list(LU_values[LU_values >= 31] - 30)]
        print(f'{base_text}{OKGREEN} OK: LCZ Classes ({LCZs}) ' f'present {ENDC}')

    base_text = (
        f'> Check 4: FRC_URB2D present in '
        f"{info.dst_lcz_params_file.split('/')[-1]}?"
    )
    ifile = info.dst_lcz_params_file
    da = xr.open_dataset(ifile)
    if 'FRC_URB2D' not in list(da.keys()):
        print(
            f'{base_text}\n{WARNING} WARNING: FRC_URB2D not '
            f'present in {ifile} {ENDC}'
        )
        frc_urb2d_present = 'NO'
    else:
        FRC_URB2D = da.FRC_URB2D.values
        print(
            f'{base_text}{OKGREEN} OK: FRC_URB2D values '
            f"range between {'{:0.2f}'.format(FRC_URB2D.min())} and "
            f"{'{:0.2f}'.format(FRC_URB2D.max())} {ENDC}"
        )
        frc_urb2d_present = 'YES'

    base_text = (
        f'> Check 5: URB_PARAMS matrix present in file '
        f"{info.dst_lcz_params_file.split('/')[-1]}?"
    )
    ifile = info.dst_lcz_params_file
    da = xr.open_dataset(ifile)
    if 'URB_PARAM' not in list(da.keys()):
        print(
            f'{base_text}\n{WARNING} WARNING: URB_PARAM matrix not ' f'present {ENDC}'
        )
        urb_param_present = 'NO'
    else:
        print(f'{base_text}{OKGREEN} OK {ENDC}')
        urb_param_present = 'YES'

    if urb_param_present == 'YES':
        base_text = (
            '> Check 6: Do URB_PARAM variable values follow expected '
            f"range in {info.dst_lcz_params_file.split('/')[-1]}?"
        )
        ifile = info.dst_lcz_params_file
        da = xr.open_dataset(ifile)

        # URB PAR Indices: https://ral.ucar.edu/sites/default/files/public/product-tool/NUDAPT_44_Documentation.pdf
        ucp_dict: Dict[str, Dict[str, Any]] = {
            'LP_URB2D': {'index': 91, 'range': [0, 1]},
            'MH_URB2D': {
                'index': 92,
                'range': [0, ucp_table['MH_URB2D'].max() + ucp_table['MH_URB2D'].std()],
            },
            'STDH_URB2D': {
                'index': 93,
                'range': [
                    0,
                    ucp_table['MH_URB2D_MAX'].max() + ucp_table['MH_URB2D'].std(),
                ],
            },
            'HGT_URB2D': {
                'index': 94,
                'range': [0, ucp_table['MH_URB2D'].max() + ucp_table['MH_URB2D'].std()],
            },
            'LB_URB2D': {'index': 95, 'range': [0, 5]},
            'LF_URB2D': {'index': 96, 'range': [0, 5]},
        }

        def _check_range(darr: NDArray[np.float_], exp_range: List[int]) -> int:
            total_len = len(darr)
            sel_len = ((darr >= exp_range[0]) & (darr <= exp_range[1])).sum(axis=0)

            if not (total_len - sel_len == 0):
                return -1
            else:
                return 0

        print(base_text)
        for ucp_key in ucp_dict.keys():

            darr = da.URB_PARAM[
                0, ucp_dict[ucp_key]['index'] - 1, :, :
            ].values.flatten()
            exp_range = ucp_dict[ucp_key]['range']

            result = _check_range(darr, exp_range)

            if result == -1:
                print(
                    f'{WARNING} WARNING: {ucp_key} exceeds '
                    f'expected value range {ENDC}'
                )
            else:
                print(f'{OKGREEN}   + OK for {ucp_key} {ENDC}')

        base_text = (
            '> Check 7: Does HI_URB2D sum to 100% for urban pixels '
            f"in {info.dst_lcz_params_file.split('/')[-1]}?"
        )
        da = xr.open_dataset(info.dst_lcz_params_file)
        hi_sum = da.URB_PARAM[0, 117:, :, :].sum(axis=0)
        hi_sum = hi_sum.where(hi_sum != 0, drop=True)

        if np.nanmax(np.abs((100 - hi_sum).values)) > 0.1:
            print(
                f'{base_text}\n{WARNING} WARNING: Not all pixels '
                f'have sum HI_URB2D == 100% {ENDC}'
            )
        else:
            print(f'{base_text}{OKGREEN} OK {ENDC}')

    if frc_urb2d_present == 'YES':
        base_text = (
            '> Check 8: Do FRC_URB and LCZs (from LU_INDEX) cover same extent '
            f"in {info.dst_lcz_params_file.split('/')[-1]}?"
        )
        frc_urb2d = xr.open_dataset(info.dst_lcz_params_file).FRC_URB2D
        lu_index = xr.open_dataset(info.dst_lcz_params_file).LU_INDEX
        frc_urb_res = xr.where(frc_urb2d != 0, 1, 0)
        lu_index_res = xr.where(lu_index >= 31, 1, 0)

        if int((frc_urb_res - lu_index_res).sum()) != 0:
            print(
                f'{base_text}\n{WARNING} WARNING: FRC_URB and LCZs in '
                f'LU_INDEX do not cover same extent {ENDC}'
            )
        else:
            print(f'{base_text}{OKGREEN} OK {ENDC}')

    base_text = (
        '> Check 9: Extent and # urban pixels same for '
        '*_extent.nc and *_params.nc output file?'
    )
    da_e = xr.open_dataset(info.dst_lcz_extent_file)
    da_p = xr.open_dataset(info.dst_lcz_params_file)
    da_e_res = xr.where(da_e.LU_INDEX == 13, 1, 0)
    da_p_res = xr.where(da_p.LU_INDEX >= 31, 1, 0)

    if int((da_p_res - da_e_res).sum()) != 0:
        print(
            f'{base_text}\n {WARNING} WARNING: Different '
            f'# urban pixels (or extent) '
            f'according to LU_INDEX: '
            f' - extent: {int(da_e_res.sum().values)}'
            f' - params: {int(da_p_res.sum().values)} {ENDC}'
        )
    else:
        print(
            f'{base_text}{OKGREEN} OK, urban extent the same '
            f'({int(da_p_res.sum().values)}) {ENDC}'
        )

    print('> Cleaning up ... all done!')
    if os.path.exists(info.dst_gridinfo):
        os.remove(info.dst_gridinfo)
    if os.path.exists(info.src_file_clean):
        os.remove(info.src_file_clean)

    print(
        f'\n\n ----------- !! NOTE !! --------- \n'
        f' Set nbui_max to {nbui_max} during compilation, '
        f'in order to optimize memory storage.\n\n'
    )


###############################################################################
##### __main__  scope
###############################################################################

if __name__ == '__main__':

    raise SystemExit(main())

###############################################################################
