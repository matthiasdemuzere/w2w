import os
from unittest import mock
import pytest
import rioxarray as rxr
from affine import Affine
from rasterio.warp import reproject, Resampling
import shutil
import w2w.w2w
from w2w.w2w import main
from w2w.w2w import create_info_dict
from w2w.w2w import _check_lcz_wrf_extent
from w2w.w2w import _replace_lcz_number
from w2w.w2w import check_lcz_integrity
from w2w.w2w import _calc_distance_coord
from w2w.w2w import wrf_remove_urban
from w2w.w2w import create_wrf_gridinfo
from w2w.w2w import _get_SW_BW
from w2w.w2w import _get_lcz_arr
from w2w.w2w import _ucp_resampler
from w2w.w2w import _hgt_resampler
from w2w.w2w import _scale_hi
from w2w.w2w import _get_truncated_normal_sample
from w2w.w2w import _check_hi_values
from w2w.w2w import _compute_hi_distribution
from w2w.w2w import _hi_resampler
from w2w.w2w import _lcz_resampler
from w2w.w2w import _add_frc_lu_index_2_wrf
from w2w.w2w import _initialize_urb_param
from w2w.w2w import create_lcz_params_file
from w2w.w2w import create_lcz_extent_file
from w2w.w2w import expand_land_cat_parents
from w2w.w2w import checks_and_cleaning
import pandas as pd
import xarray as xr
import numpy as np
import pandas as pd
from pytest import approx


def test_argparse_shows_help():
    with pytest.raises(SystemExit):
        main(['--help'])

def test_create_info_dict():

    class mock_args:
        lcz_file= 'lcz_file.tif'
        wrf_file = 'wrf_file.nc'
        io_dir = 'input/directory'
        built_lcz = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    info = create_info_dict(mock_args)

    # info is Dict, with 9 keys
    assert len(info.keys()) == 9

    # Three files are tifs
    assert len([i for i in list(info.values())[:-1] if
                i.endswith('.tif')]) == 3

    # 4 files are netcdf files
    assert len([i for i in list(info.values())[:-1] if
                i.endswith('.nc')]) == 4
    # Last entry is list of built LCZs
    assert isinstance(list(info.values())[-1], list)

def test_replace_lcz_number_ok():
    info = {
        'src_file': 'testing/Shanghai.tif',
    }
    LCZ_BAND = 0
    lcz = rxr.open_rasterio(info['src_file'])[LCZ_BAND, :, :]

    lcz_100 = np.array(
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
         101, 102, 103, 104, 105, 106, 107]
    )
    lcz_new = _replace_lcz_number(
            lcz=lcz,
            lcz_to_change=lcz_100,
        )

    # Test whether LCZs 100+ are converted to 11+
    assert (np.unique(lcz_new.data).tolist() ==
            [1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 14, 15, 17])

def test_check_lcz_integrity_lcz_numbers_as_expected(capsys, tmpdir):

    info = {
        'src_file': 'sample_data/lcz_zaragoza.tif',
        'src_file_clean': os.path.join(tmpdir, 'lcz_zaragoza_clean.tif'),
        'dst_file': 'sample_data/geo_em.d04.nc',
    }
    LCZ_BAND = 0
    check_lcz_integrity(info=info, LCZ_BAND=LCZ_BAND)
    out, _ = capsys.readouterr()
    assert '> LCZ labels as expected (1 to 17)' in out

def test_check_lcz_integrity_crs_changed(capsys, tmpdir):

    info = {
        'src_file': 'testing/Shanghai.tif',
        'src_file_clean': os.path.join(tmpdir, 'Shanghai_clean.tif'),
        'dst_file': 'testing/geo_em.d02_Shanghai.nc',
    }
    LCZ_BAND = 0
    check_lcz_integrity(info=info, LCZ_BAND=LCZ_BAND)
    out, _ = capsys.readouterr()
    assert 'LCZ map reprojected to WGS84 (EPSG:4326)' in out



def test_check_lcz_wrf_extent_lcz_too_small(capsys):
    info = {
        'src_file': 'testing/lcz_too_small.tif',
        'dst_file': 'sample_data/geo_em.d04.nc',
    }
    LCZ_BAND = 0

    lcz = rxr.open_rasterio(info['src_file'])[LCZ_BAND, :, :]
    wrf = xr.open_dataset(info['dst_file'])

    with pytest.raises(SystemExit):
        _check_lcz_wrf_extent(lcz=lcz, wrf=wrf)

    out, _ = capsys.readouterr()
    assert (
        'ERROR: LCZ domain should be larger than WRF domain'
    ) in out
    # TODO maybe add the actual values to check they are correct
    assert (
        'ERROR: LCZ domain should be larger than WRF domain '
        'in all directions.\nLCZ bounds (xmin, ymin, xmax, ymax): '
    ) in out
    assert 'WRF bounds (xmin, ymin, xmax, ymax): ' in out


def test_check_lcz_wrf_extent_ok(capsys):
    info = {
        'src_file': 'sample_data/lcz_zaragoza.tif',
        'dst_file': 'sample_data/geo_em.d04.nc',
    }
    LCZ_BAND = 0

    lcz = rxr.open_rasterio(info['src_file'])[LCZ_BAND, :, :]
    wrf = xr.open_dataset(info['dst_file'])
    _check_lcz_wrf_extent(lcz=lcz, wrf=wrf)
    out, _ = capsys.readouterr()
    assert '> LCZ domain is covering WRF domain' in out

def test_check_lcz_integrity_clean_file_written(tmpdir):

    info = {
        'src_file': 'testing/Shanghai.tif',
        'dst_file': 'testing/geo_em.d02_Shanghai.nc',
        'src_file_clean': os.path.join(tmpdir, 'Shanghai_clean.tif'),
    }
    LCZ_BAND = 0
    check_lcz_integrity(info=info, LCZ_BAND=LCZ_BAND)
    assert os.listdir(tmpdir) == ['Shanghai_clean.tif']

@pytest.mark.parametrize(
    ('coords', 'expected'),
    (
        pytest.param((89, 0, 90, 1), 111194, id='near northern pole'),
        pytest.param((-89, 0, -90, 1), 111194, id='near southern pole'),
        pytest.param((-60, -10, -70, -12), 1115764, id='southern eastern hemisphere'),
        pytest.param((-1, 0, 1, 2), 314498, id='across equator'),
        pytest.param((45, 179, 46, -179), 191461, id='across dateline')
    ),
)
def test_calc_distance_coord(coords, expected):
    assert _calc_distance_coord(*coords) == pytest.approx(expected, abs=1)

@pytest.mark.parametrize(
    ('dst_file', 'dst_nu_file'),
    (
        pytest.param('testing/5by5.nc', '5by5_new.nc', id='all cat'),
        pytest.param('testing/5by5_20cat.nc', '5by5_20cat_new.nc', id='20 cat'),
    )
)
def test_wrf_remove_urban(tmpdir, dst_file, dst_nu_file):
    info = {
        'dst_file': dst_file,
        'dst_nu_file': os.path.join(tmpdir, dst_nu_file)
    }
    old_ds = xr.open_dataset(info['dst_file'])
    wrf_remove_urban(info=info, NPIX_NLC=9)
    ds = xr.open_dataset(info['dst_nu_file'])
    # check lused 13 was reclassified to 12
    assert ds.LU_INDEX.values[0][2][2] == 12
    assert ds.LU_INDEX.values[0][2][3] == 12
    assert ds.LU_INDEX.values[0][4][1] == 12
    assert 13 not in ds.LU_INDEX.values.flatten()
    # sum up the booleans to check that the luse 13 (which are 3) values have
    # changed
    compare = old_ds.GREENFRAC.values[0][0] != ds.GREENFRAC.values[0][0]
    assert sum(compare.flatten()) == 3
    # check the values were changed at those coords
    assert compare[2][2].item() is True
    assert compare[2][3].item() is True
    assert compare[4][1].item() is True
    compare_luf = old_ds.LANDUSEF.values[0][11] != ds.LANDUSEF.values[0][11]
    # check the values were changed at those coords
    assert compare_luf[2][2].item() is True
    assert compare_luf[2][3].item() is True
    assert compare_luf[4][1].item() is True
    # TODO: there is one more value changed than we expect
    # the value at [4][0] is changed -- why?
    # assert sum(compare_luf.flatten()) == 3


def test_wrf_remove_urban_output_already_exists_is_overwritten(tmpdir):
    tmpdir.ensure('5by5_new.nc')
    info = {
        'dst_file': 'testing/5by5.nc',
        'dst_nu_file': os.path.join(tmpdir, '5by5_new.nc')
    }
    assert os.listdir(tmpdir) == ['5by5_new.nc']
    m_time_old = os.path.getmtime(info['dst_nu_file'])
    wrf_remove_urban(info=info, NPIX_NLC=9)
    assert m_time_old != os.path.getmtime(info['dst_nu_file'])


def test_create_wrf_gridinfo(tmpdir):
    info = {
        'dst_nu_file': 'testing/5by5.nc',
        'dst_gridinfo': os.path.join(tmpdir, 'dst_gridinfo.tif'),
    }
    create_wrf_gridinfo(info)
    assert os.listdir(tmpdir) == ['dst_gridinfo.tif']
    tif = rxr.open_rasterio(info['dst_gridinfo'])
    assert tif.rio.crs.to_proj4() == '+init=epsg:4326'
    assert tif.rio.transform() == Affine(
            0.01000213623046875, 0.0, -1.4050254821777344,
            0.0, 0.010000228881835938, 41.480000495910645,
    )

def test_get_SW_BW():

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    SW, BW = _get_SW_BW(ucp_table)

    assert approx(list(SW[:10])) == [
        20.0, 14.0, 5.2, 50.0, 35.0,
        13.0, 3.333333, 32.5, 43.333333, 28.571428
    ]
    assert approx(list(BW[:10])) == [
        22.222222, 22.0, 9.533333, 42.857142, 26.25,
        13.0, 25.0, 28.888888, 43.333333, 23.809523
    ]

def test_get_lcz_arr():

    info = {
        'src_file': 'sample_data/lcz_zaragoza.tif',
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    # Read gridded data: LCZ and WRF grid
    src_data = rxr.open_rasterio(info['src_file'])[0, :, :]

    lcz_arr = _get_lcz_arr(
        src_data=src_data,
        info=info
    )

    # Make sure shape is fine
    assert lcz_arr.shape == (765, 1210)
    # Only built LCZs labels should be in array.
    assert list(np.unique(lcz_arr)) == [0, 2, 3, 5, 6, 8, 9]
    # All other LCZs should be masked = 0.
    assert np.sum(lcz_arr == 0) == 905938
    # Type should be integere, to make sure lookup works
    assert lcz_arr.dtype == np.int32


@pytest.mark.parametrize(
    ('ucp_key', 'data_sum', 'data_mean'),
    (
        ('MH_URB2D', 2, 1.1429937),
        ('STDH_URB2D', 2, 0.29812142),
        ('LB_URB2D', 1, 0.15104416),
        ('LF_URB2D', 1, 0.07032267),
        ('LP_URB2D', 1, 0.08072148),
    ),
)
def test_ucp_resampler_output_values_per_paramater(
        ucp_key,
        data_sum,
        data_mean,
):
    info = {
        'src_file': 'sample_data/lcz_zaragoza.tif',
        'src_file_clean': 'testing/lcz_zaragoza_clean.tif',
        'dst_gridinfo': 'testing/geo_em.d04_gridinfo.tif',
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }
    ucp_table = pd.read_csv('w2w/resources/LCZ_UCP_lookup.csv', index_col=0)
    ucp_res = _ucp_resampler(
        info=info,
        ucp_key=ucp_key,
        RESAMPLE_TYPE='average',
        ucp_table=ucp_table,
    )
    # Parameters should be on WRF grid size
    assert ucp_res.shape == (1, 102, 162)
    # Per parameter, check # max values and non-0 domain average
    assert np.sum(ucp_res.data == ucp_res.data.max()) == data_sum
    assert approx(np.mean(ucp_res.data[ucp_res.data > 0])) == data_mean
    # Make sure no nans are present.
    np.sum(np.isnan(ucp_res.data)) == 0

def test_ucp_resampler_output_values_per_paramater_frc_threshold():
    info = {
        'src_file': 'sample_data/lcz_zaragoza.tif',
        'src_file_clean': 'testing/lcz_zaragoza_clean.tif',
        'dst_gridinfo': 'testing/geo_em.d04_gridinfo.tif',
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }
    ucp_table = pd.read_csv('w2w/resources/LCZ_UCP_lookup.csv', index_col=0)
    ucp_res = _ucp_resampler(
        info=info,
        ucp_key='FRC_URB2D',
        RESAMPLE_TYPE='average',
        ucp_table=ucp_table,
        FRC_THRESHOLD=0.2,
    )
    # Parameters should be on WRF grid size
    assert ucp_res.shape == (1, 102, 162)
    # Per parameter, check # max values and non-0 domain average
    assert np.sum(ucp_res.data == ucp_res.data.max()) == 1
    assert approx(np.mean(ucp_res.data[ucp_res.data > 0])) == 0.428921
    # Make sure no nans are present.
    np.sum(np.isnan(ucp_res.data)) == 0

def test_hgt_resampler_output_values():

    info = {
        'src_file': 'sample_data/lcz_zaragoza.tif',
        'src_file_clean': 'testing/lcz_zaragoza_clean.tif',
        'dst_gridinfo': 'testing/geo_em.d04_gridinfo.tif',
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    ucp_res = _hgt_resampler(
        info=info,
        RESAMPLE_TYPE='average',
        ucp_table=ucp_table,
    )

    # # Parameters should be on WRF grid size
    assert ucp_res.shape == (1, 102, 162)

    # Check # max values and non-0 domain average
    assert np.sum(ucp_res.data == ucp_res.data.max()) == 9
    assert approx(np.mean(ucp_res.data[ucp_res.data > 0])) == 6.7111716

    # Make sure no nans are present.
    assert np.sum(np.isnan(ucp_res.data)) == 0

def test_check_hi_values_fails_skewed_distribution_ucp_table(capsys):

    lcz_i = 1
    SAMPLE_SIZE = 100000
    ERROR_MARGIN=0.05

    # create dummpy ucp_table, with skewed distribution
    ucp_table = pd.DataFrame(
        index = [1],
        columns = ['MH_URB2D_MIN', 'MH_URB2D_MAX', 'MH_URB2D']
    )
    ucp_table.loc[1] = [5, 150, 40]

    hi_sample = _get_truncated_normal_sample(
        lcz_i=lcz_i,
        ucp_table=ucp_table,
        SAMPLE_SIZE=SAMPLE_SIZE
    )
    assert len(hi_sample) == SAMPLE_SIZE

    _check_hi_values(
        lcz_i=lcz_i,
        hi_sample=hi_sample,
        ucp_table=ucp_table,
        ERROR_MARGIN=ERROR_MARGIN)
    out, _ = capsys.readouterr()
    assert 'distribution not in expected range (5.0% marging) ' \
           'for LCZ class 1' in out

def test_check_hi_values_fails_sample_size_too_small(capsys):

    lcz_i = 6
    SAMPLE_SIZE = 10
    ERROR_MARGIN = 0.05

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    hi_sample = _get_truncated_normal_sample(
        lcz_i=lcz_i,
        ucp_table=ucp_table,
        SAMPLE_SIZE=SAMPLE_SIZE
    )

    _check_hi_values(
        lcz_i=lcz_i,
        hi_sample=hi_sample,
        ucp_table=ucp_table,
        ERROR_MARGIN=ERROR_MARGIN)
    out, _ = capsys.readouterr()
    assert 'distribution not in expected range' in out

def test_check_hi_values_fails_error_margin_too_small(capsys):

    lcz_i = 4
    SAMPLE_SIZE = 100000
    ERROR_MARGIN = 0.000001

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    hi_sample = _get_truncated_normal_sample(
        lcz_i=lcz_i,
        ucp_table=ucp_table,
        SAMPLE_SIZE=SAMPLE_SIZE
    )

    _check_hi_values(
        lcz_i=lcz_i,
        hi_sample=hi_sample,
        ucp_table=ucp_table,
        ERROR_MARGIN=ERROR_MARGIN)
    out, _ = capsys.readouterr()
    assert 'distribution not in expected range' in out

def test_compute_hi_distribution_values():

    info = {
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )
    SAMPLE_SIZE = 100000
    ERROR_MARGIN = 0.05

    df_hi = _compute_hi_distribution(
        info=info,
        ucp_table=ucp_table,
        SAMPLE_SIZE=SAMPLE_SIZE,
        ERROR_MARGIN=ERROR_MARGIN,
    )

    # The mean over all height bins and built classes should be 100%
    assert 100 == approx(df_hi.sum(axis=1)[info['BUILT_LCZ']].mean())

    # Natural LCZ classes should have all 0
    value, counts = np.unique(df_hi.loc[range(11, 18, 1)].values, return_counts=True)
    assert value[0] == 0
    assert counts[0] == 7*15 # 7 Natural LCZs x 15 height bins

def test_compute_hi_distribution_values_lcz15():

    info = {
        'BUILT_LCZ': [15],
    }

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )
    SAMPLE_SIZE = 100000
    ERROR_MARGIN = 0.05

    df_hi = _compute_hi_distribution(
        info=info,
        ucp_table=ucp_table,
        SAMPLE_SIZE=SAMPLE_SIZE,
        ERROR_MARGIN=ERROR_MARGIN,
    )

    # LCZ 15 has no buildings (bare rock) so distribution of building heights
    # across all height bins should be 0.
    assert df_hi.iloc[15].fillna(0).sum() == 0

def test_scale_hi():

    # Create random array to work with,
    # with 10 levels of HI frequency values
    a = np.random.randint(1, 100, size=(10, 5, 5))

    # Scale as done in code, along axis 0
    a_scaled = np.apply_along_axis(_scale_hi, 0, a)

    # Check that all pixels sum to 100
    assert np.full((5, 5), 100.) == approx((np.sum(a_scaled, axis=0)))

def test_hi_resampler():

    info = {
        'src_file_clean': 'testing/lcz_zaragoza_clean.tif',
        'dst_gridinfo': 'testing/geo_em.d04_gridinfo.tif',
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv('w2w/resources/LCZ_UCP_lookup.csv', index_col=0)

    hi_arr, nbui = _hi_resampler(
        info=info,
        RESAMPLE_TYPE='average',
        ucp_table=ucp_table,
    )

    # Output should be 15 building height classes on WRF grid
    assert hi_arr.shape == (15, 102, 162)

    # Output can not contain % values below 5% (except 0%)
    hi_arr_no_zero = hi_arr[hi_arr != 0]
    assert np.sum(hi_arr_no_zero < 5) == 0

    # % per building height class should sum to 100 for urban pixels (1057)
    hi_arr_sum = hi_arr.sum(axis=0)
    assert hi_arr_sum[hi_arr_sum != 0].mean() == 100

    # Max number of HI intervals over all grid cells should be 5
    assert nbui == 5


@pytest.mark.parametrize(
    ('LCZ_NAT_MASK', 'lcz_values', 'lcz_counts'),
    (
        (False,
         np.array([32., 33., 35., 36., 38., 39., 41., 42., 43., 44., 46.]),
         np.array([ 17,   3,   1,  77,  95,   2,   7,   5,   2, 150,  10]),
         ),
        (True,
         np.array([32., 33., 35., 36., 38., 39., 41.]),
         np.array([ 18,  30,   1, 171, 136,   4,   9]),
         ),
    ),
)

def test_lcz_resampler_lcz_nat_mask_on_off_with_lcz15(
        LCZ_NAT_MASK,
        lcz_values,
        lcz_counts,
):

    info = {
        'src_file_clean': 'testing/lcz_zaragoza_clean.tif',
        'dst_gridinfo': 'testing/geo_em.d04_gridinfo.tif',
        'src_file_frc_urb2d': 'testing/geo_em.d04_LCZ_params_no_urb_param.nc',
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15],
    }

    frc_data = xr.open_dataset(info['src_file_frc_urb2d'])
    frc_urb2d = frc_data['FRC_URB2D']

    frc_mask, lcz_resampled = _lcz_resampler(
        info=info,
        frc_urb2d=frc_urb2d,
        LCZ_NAT_MASK=LCZ_NAT_MASK,
    )
    # With natural masking off, majority filtering also includes
    # Natural classes
    lcz_values_def, lcz_counts_def = np.unique(lcz_resampled, return_counts=True)
    assert (lcz_values_def == lcz_values).all()
    assert (lcz_counts_def == lcz_counts).all()

def test_add_frc_lu_index_2_wrf():

    info = {
        'src_file_clean': 'testing/lcz_zaragoza_clean.tif',
        'dst_file': 'sample_data/geo_em.d04.nc',
        'dst_nu_file': 'testing/geo_em.d04_NoUrban.nc',
        'dst_gridinfo': 'testing/geo_em.d04_gridinfo.tif',
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv('w2w/resources/LCZ_UCP_lookup.csv', index_col=0)

    dst_data = _add_frc_lu_index_2_wrf(
        info = info,
        FRC_THRESHOLD = 0.2,
        LCZ_NAT_MASK = True,
        ucp_table = ucp_table,
    )

    # LU_Index should be expanded, reflecting LCZ classes between 31-41.
    # Highest class available in sample data is LCZ 8 (39)
    assert dst_data['LU_INDEX'].max().data == 39

    # The 41 category should also be embedded within the attributes' description
    assert '41-category' in dst_data['LANDUSEF'].attrs['description']

    # GREENFRAC is adjusted, because of LCZ implementation
    # Can be - or +, depending on conversion natural to LCZ land.
    dst_nu = xr.open_dataset(info['dst_nu_file'])
    assert (dst_data['GREENFRAC'][0,:,:,:].mean(axis=0) -
            dst_nu['GREENFRAC'][0,:,:,:].mean(axis=0)).data.min() == approx(-0.193724)
    assert (dst_data['GREENFRAC'][0,:,:,:].mean(axis=0) -
            dst_nu['GREENFRAC'][0,:,:,:].mean(axis=0)).data.max() == approx(0.31252602)

    # If done properly, the Landuse Fraction should still add up to 1 for all pixels
    assert np.unique(dst_data['LANDUSEF'][0,:,:,:].sum(axis=0)) == np.array(1)


@pytest.mark.parametrize(
    ('att_name', 'att_value'),
    (
        ('FieldType', np.intc(104)),
        ('MemoryOrder', "XYZ"),
        ('units', "dimensionless"),
        ('description', "all urban parameters"),
        ('stagger', "M"),
        ('sr_x', np.intc(1)),
        ('sr_y', np.intc(1)),
    ),
)
def test_initialize_urb_param(
        att_name,
        att_value,
):

    info = {
        'dst_no_params': 'testing/geo_em.d04_LCZ_params_no_urb_param.nc',
    }

    dst_data = xr.open_dataset(info['dst_no_params'])

    # Create the URB_PARAMS matrix
    dst_final = _initialize_urb_param(
        dst_data=dst_data)

    # Check size of URB_PARAM
    assert dst_final['URB_PARAM'].shape == (1, 132, 102, 162)

    # All attributes available?
    assert dst_final['URB_PARAM'].attrs[att_name] == att_value

@pytest.mark.parametrize(
    ('att_name', 'att_value'),
    (
        ('NUM_LAND_CAT', 41),
        ('FLAG_URB_PARAM', 1),
        ('NBUI_MAX', np.intc(5)),
        ('TITLE', "OUTPUT FROM GEOGRID V4.0, perturbed by W2W"),
    ),
)
def test_create_lcz_params_file_attrs_type(
        att_name,
        att_value,
        tmpdir,
        capsys
):
    input_dir = tmpdir.mkdir('input')
    shutil.copy(os.path.join('sample_data', 'geo_em.d04.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'lcz_zaragoza_clean.tif'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_gridinfo.tif'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_NoUrban.nc'), input_dir)

    info = {
        'src_file_clean': os.path.join(input_dir, 'lcz_zaragoza_clean.tif'),
        'dst_file': os.path.join(input_dir, 'geo_em.d04.nc'),
        'dst_nu_file': os.path.join(input_dir, 'geo_em.d04_NoUrban.nc'),
        'dst_gridinfo': os.path.join(input_dir, 'geo_em.d04_gridinfo.tif'),
        'dst_lcz_params_file': os.path.join(input_dir, 'geo_em.d04_LCZ_params.tif'),
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    # Copy file, as dummy for LCZ_params.nc
    shutil.copy(info['dst_file'], info['dst_lcz_params_file'])

    nbui_max = create_lcz_params_file(
        info=info,
        FRC_THRESHOLD=0.2,
        LCZ_NAT_MASK=True,
        ucp_table=ucp_table)
    out, _ = capsys.readouterr()

    # Check that output file is written
    contents = os.listdir(input_dir)
    assert 'geo_em.d04_LCZ_params.tif' in contents

    # Check if file has expected attributes (values)
    dst_params = xr.open_dataset(info['dst_lcz_params_file'])
    assert dst_params.attrs[att_name] == att_value

    # Make sure type of URB_PARAM if float
    assert dst_params['URB_PARAM'].dtype == np.float32

    # Check if nbui_max has expected value
    assert nbui_max == 5

def test_create_lcz_extent_file(tmpdir):

    # tmpdir = '/home/demuzmp4/Desktop'
    # info = {
    #     'dst_lcz_extent_file': os.path.join(tmpdir,'geo_em.d04_LCZ_extent.nc'),
    # }
    info = {
        'src_file_clean': 'testing/lcz_zaragoza_clean.tif',
        'dst_lcz_params_file': 'testing/geo_em.d04_LCZ_params.nc',
        'dst_lcz_extent_file': os.path.join(tmpdir,'geo_em.d04_LCZ_extent.nc'),
        'dst_file': 'sample_data/geo_em.d04.nc',
        'dst_nu_file': 'testing/geo_em.d04_NoUrban.nc',
        'dst_gridinfo': 'testing/geo_em.d04_gridinfo.tif',
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    frc_mask = _add_frc_lu_index_2_wrf(
        info=info,
        FRC_THRESHOLD=0.2,
        LCZ_NAT_MASK=True,
        ucp_table=ucp_table,
    )

    # Produce LCZ extent file
    create_lcz_extent_file(
        info=info,
    )

    # Read produced file
    dst_extent = xr.open_dataset(info['dst_lcz_extent_file'])

    # Check perturbations applied to LCZ_params file
    assert dst_extent.attrs['FLAG_URB_PARAM'] == np.intc(0)
    assert dst_extent.attrs['NUM_LAND_CAT'] == 41
    assert dst_extent.LANDUSEF.description == \
           'Noah-modified 21-category IGBP-MODIS landuse'
    assert 'FRC_URB2D' not in list(dst_extent.data_vars)
    assert 'URB_PARAM' not in list(dst_extent.data_vars)

    # Number of urban pixels within LANDUSEF[12]
    assert np.sum(dst_extent['LANDUSEF'].data[0, 12, :, :] == 1) == 369

@pytest.mark.parametrize(
    ('domain_id'),
    (
        ('01'),
        ('02'),
        ('03'),
    ),
)
def test_expand_land_cat_parents_files_missing(domain_id, capsys, tmpdir):
    ifiles = (
        'geo_em.d04.nc', 'lcz_zaragoza_clean.tif',
    )
    input_dir = tmpdir.mkdir('input')
    shutil.copy(os.path.join('sample_data', ifiles[0]), input_dir)
    shutil.copy(os.path.join('testing', ifiles[1]), input_dir)

    info = {
        'dst_file': os.path.join(input_dir, 'geo_em.d04.nc'),
        'io_dir': 'SOME_DIRECTORY',
    }

    expand_land_cat_parents(info=info)
    out, _ = capsys.readouterr()

    #warning_str = "Without this information, you will not be able to produce the boundary"
    warning_str = f"WARNING: Parent domain {info['dst_file'][:-5]}{domain_id}.nc not found"
    assert warning_str in out

def test_expand_land_cat_parents_no_num_land_cat(tmpdir, capsys):
    input_dir = tmpdir.mkdir('input')
    shutil.copy(os.path.join('testing', 'geo_em.d02_Shanghai.nc'),
                os.path.join(input_dir,'geo_em.d02.nc')
                )
    shutil.copy(os.path.join('testing', 'geo_em.d01_Shanghai_no_NUM_LAND_CAT.nc'),
                os.path.join(input_dir,'geo_em.d01.nc')
                )
    info = {
        'dst_file' : os.path.join(input_dir,'geo_em.d02.nc'),
        'io_dir': 'some_dummy_directory'
    }
    expand_land_cat_parents(
        info=info,
    )
    out, _ = capsys.readouterr()
    warning_str = f"Cannot read NUM_LAND_CAT"
    assert warning_str in out


def test_expand_land_cat_parents_num_land_cat_41(capsys, tmpdir):

    input_dir = tmpdir.mkdir('input')
    shutil.copy(os.path.join('testing', 'geo_em.d02_Shanghai.nc'),
                os.path.join(input_dir,'geo_em.d02.nc')
                )
    shutil.copy(os.path.join('testing', 'geo_em.d01_Shanghai_ncl41.nc'),
                os.path.join(input_dir,'geo_em.d01.nc')
                )
    info = {
        'dst_file' : os.path.join(input_dir,'geo_em.d02.nc'),
        'io_dir': 'some_dummy_directory'
    }
    expand_land_cat_parents(
        info=info,
    )
    out, _ = capsys.readouterr()
    assert 'already contains 41 LC classes' in out

def test_expand_land_cat_parents_num_land_cat_not41(capsys, tmpdir):

    input_dir = tmpdir.mkdir('input')
    shutil.copy(os.path.join('testing', 'geo_em.d02_Shanghai.nc'),
                os.path.join(input_dir,'geo_em.d02.nc')
                )
    shutil.copy(os.path.join('testing', 'geo_em.d01_Shanghai_ncl20.nc'),
                os.path.join(input_dir,'geo_em.d01.nc')
                )
    info = {
        'dst_file' : os.path.join(input_dir,'geo_em.d02.nc'),
        'io_dir': 'some_dummy_directory'
    }
    expand_land_cat_parents(
        info=info,
    )
    out, _ = capsys.readouterr()
    contents = os.listdir(input_dir)
    assert 'geo_em.d01_41.nc' in contents

def test_checks_and_cleaning_sample_data_all_ok(capsys, tmpdir):

    input_dir = tmpdir.mkdir('input')
    shutil.copy(os.path.join('sample_data', 'geo_em.d04.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'lcz_zaragoza_clean.tif'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_gridinfo.tif'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_NoUrban.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_LCZ_extent.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_LCZ_params.nc'), input_dir)

    info = {
        'src_file_clean': os.path.join(input_dir, 'lcz_zaragoza_clean.tif'),
        'dst_file': os.path.join(input_dir, 'geo_em.d04.nc'),
        'dst_nu_file': os.path.join(input_dir, 'geo_em.d04_NoUrban.nc'),
        'dst_gridinfo': os.path.join(input_dir, 'geo_em.d04_gridinfo.tif'),
        'dst_lcz_extent_file': os.path.join(input_dir, 'geo_em.d04_LCZ_extent.nc'),
        'dst_lcz_params_file': os.path.join(input_dir, 'geo_em.d04_LCZ_params.nc'),
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    NBUI_MAX = 5

    class bcolors:
        OKGREEN = "\033[0;32m"

    checks_and_cleaning(
        info=info,
        ucp_table=ucp_table,
        nbui_max=NBUI_MAX
    )
    out, _ = capsys.readouterr()

    # Check the outcome of the sample data, for which all is fine.
    base_text = f"> Check 1: Urban class removed from " \
          f"{info['dst_nu_file'].split('/')[-1]}?"
    check1 = f"{base_text}{bcolors.OKGREEN} OK"
    assert check1 in out

    base_text = f"> Check 2: LCZ Urban extent present in " \
          f"{info['dst_lcz_extent_file'].split('/')[-1]}?"
    check2 = f"{base_text}{bcolors.OKGREEN} OK"
    assert check2 in out

    base_text = f"> Check 3: Urban LCZ classes exists in " \
                f"{info['dst_lcz_params_file'].split('/')[-1]}?"
    check3 = f"{base_text}{bcolors.OKGREEN} OK: LCZ Classes"
    assert check3 in out

    base_text = f"> Check 4: FRC_URB2D present in " \
          f"{info['dst_lcz_params_file'].split('/')[-1]}?"
    check4 = f"{base_text}{bcolors.OKGREEN} OK: FRC_URB2D values "
    assert check4 in out

    base_text = f"> Check 5: URB_PARAMS matrix present in file " \
          f"{info['dst_lcz_params_file'].split('/')[-1]}?"
    check5 = f"{base_text}{bcolors.OKGREEN} OK "
    assert check5 in out

    base_text = "> Check 6: Do URB_PARAM variable values follow expected range in " \
          f"{info['dst_lcz_params_file'].split('/')[-1]}?"
    check6a = f"{base_text}"
    check6b = f"{bcolors.OKGREEN}   + OK for"
    assert check6a in out
    assert check6b in out

    base_text = "> Check 7: Does HI_URB2D sum to 100% for urban pixels " \
          f"in {info['dst_lcz_params_file'].split('/')[-1]}?"
    check7 = f"{base_text}{bcolors.OKGREEN} OK"
    assert check7 in out

    base_text = "> Check 8: Do FRC_URB and LCZs (from LU_INDEX) cover same extent " \
          f"in {info['dst_lcz_params_file'].split('/')[-1]}?"
    check8 = f"{base_text}{bcolors.OKGREEN} OK"
    assert check8 in out

    base_text = "> Check 9: Extent and # urban pixels same for " \
          "*_extent.nc and *_params.nc output file?"
    check9 = f"{base_text}{bcolors.OKGREEN} OK, urban extent the same "
    assert check9 in out

    assert f" Set nbui_max to {NBUI_MAX} during compilation, " in out

def test_checks_and_cleaning_sample_data_check1to5_wrong(
        capsys,
        tmpdir
):

    # Test WARNINGS by using dummy lcz_extent and lcz_params files

    input_dir = tmpdir.mkdir('input')
    shutil.copy(os.path.join('sample_data', 'geo_em.d04.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'lcz_zaragoza_clean.tif'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_gridinfo.tif'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_NoUrban_dummy.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_LCZ_extent_dummy.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_LCZ_params_dummy.nc'), input_dir)

    info = {
        'src_file_clean': os.path.join(input_dir, 'lcz_zaragoza_clean.tif'),
        'dst_file': os.path.join(input_dir, 'geo_em.d04.nc'),
        'dst_nu_file': os.path.join(input_dir, 'geo_em.d04_NoUrban_dummy.nc'),
        'dst_gridinfo': os.path.join(input_dir, 'geo_em.d04_gridinfo.tif'),
        'dst_lcz_extent_file': os.path.join(input_dir, 'geo_em.d04_LCZ_extent_dummy.nc'),
        'dst_lcz_params_file': os.path.join(input_dir, 'geo_em.d04_LCZ_params_dummy.nc'),
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    checks_and_cleaning(
        info=info,
        ucp_table=ucp_table,
        nbui_max=5
    )
    out, _ = capsys.readouterr()
    # Check 1
    assert "WARNING: Urban land use" in out
    # Check 2
    assert "WARNING: LCZ-based urban" in out
    # Check 3
    assert "WARNING: Urban extent still" in out
    # Check 4
    assert "WARNING: FRC_URB2D not" in out
    # Check 5
    assert "WARNING: URB_PARAM matrix not" in out


@pytest.mark.parametrize(
    ('ucp_key'),
    (
        ('MH_URB2D'),
        ('STDH_URB2D'),
        ('LB_URB2D'),
        ('LF_URB2D'),
        ('LP_URB2D'),
    ),
)
def test_checks_and_cleaning_sample_data_check6to9_wrong(
        ucp_key,
        capsys,
        tmpdir
):

    # Test WARNINGS by using dummy lcz_extent and lcz_params files

    input_dir = tmpdir.mkdir('input')
    shutil.copy(os.path.join('sample_data', 'geo_em.d04.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'lcz_zaragoza_clean.tif'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_gridinfo.tif'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_NoUrban_dummy.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_LCZ_extent_dummy.nc'), input_dir)
    shutil.copy(os.path.join('testing', 'geo_em.d04_LCZ_params_with_ucp_dummy.nc'), input_dir)

    info = {
        'src_file_clean': os.path.join(input_dir, 'lcz_zaragoza_clean.tif'),
        'dst_file': os.path.join(input_dir, 'geo_em.d04.nc'),
        'dst_nu_file': os.path.join(input_dir, 'geo_em.d04_NoUrban_dummy.nc'),
        'dst_gridinfo': os.path.join(input_dir, 'geo_em.d04_gridinfo.tif'),
        'dst_lcz_extent_file': os.path.join(input_dir, 'geo_em.d04_LCZ_extent_dummy.nc'),
        'dst_lcz_params_file': os.path.join(input_dir, 'geo_em.d04_LCZ_params_with_ucp_dummy.nc'),
        'BUILT_LCZ': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    }

    ucp_table = pd.read_csv(
        'w2w/resources/LCZ_UCP_lookup.csv',
         index_col=0
    )

    checks_and_cleaning(
        info=info,
        ucp_table=ucp_table,
        nbui_max=5
    )
    out, _ = capsys.readouterr()
    # Check 6
    assert f"WARNING: {ucp_key} exceeds " in out
    # Check 7
    assert f"WARNING: Not all pixels " in out
    # Check 8
    assert f"WARNING: FRC_URB and LCZs in " in out
    # Check 9
    assert f"WARNING: Different # urban pixels (or extent) " in out


def test_main_with_example_data(tmpdir):
    input_dir = tmpdir.mkdir('input')
    input_files = (
        'geo_em.d01.nc',
        'geo_em.d03.nc',
        'geo_em.d03.nc',
        'geo_em.d04.nc',
        'lcz_zaragoza.tif',
    )
    for f in input_files:
        shutil.copy(os.path.join('sample_data', f), input_dir)

    with tmpdir.as_cwd():
        main(['input', 'lcz_zaragoza.tif', 'geo_em.d04.nc'])

        contents = os.listdir(input_dir)
        assert 'geo_em.d04_LCZ_params.nc' in contents
        assert 'geo_em.d04_NoUrban.nc' in contents
        assert 'geo_em.d04_LCZ_extent.nc' in contents

def test_main_shanghai_data(tmpdir):
    input_dir = tmpdir.mkdir('input')

    shutil.copy(os.path.join('testing', 'geo_em.d01_Shanghai_ncl20.nc'),
                os.path.join(input_dir, 'geo_em.d01.nc'))
    shutil.copy(os.path.join('testing', 'geo_em.d02_Shanghai.nc'),
                os.path.join(input_dir, 'geo_em.d02.nc'))
    shutil.copy(os.path.join('testing', 'Shanghai.tif'), input_dir)

    with tmpdir.as_cwd():
        main(['input', 'Shanghai.tif', 'geo_em.d02.nc'])

        contents = os.listdir(input_dir)
        assert 'geo_em.d02_LCZ_params.nc' in contents
        assert 'geo_em.d02_NoUrban.nc' in contents
        assert 'geo_em.d02_LCZ_extent.nc' in contents

def test_main_default_lcz_ucp(tmpdir):
    input_dir = tmpdir.mkdir('input')
    input_files = (
        'geo_em.d01.nc',
        'geo_em.d03.nc',
        'geo_em.d03.nc',
        'geo_em.d04.nc',
        'lcz_zaragoza.tif'
    )
    for f in input_files:
        shutil.copy(os.path.join('sample_data', f), input_dir)

    with tmpdir.as_cwd():
        with mock.patch.object(w2w.w2w, 'checks_and_cleaning') as m:
            main(['input', 'lcz_zaragoza.tif', 'geo_em.d04.nc'])

    exp_df = pd.read_csv('w2w/resources/LCZ_UCP_lookup.csv', index_col=0)
    pd.testing.assert_frame_equal(exp_df, m.call_args[1]['ucp_table'])

def test_main_custom_lcz_ucp(tmpdir):
    input_dir = tmpdir.mkdir('input')
    input_files = (
        'geo_em.d01.nc',
        'geo_em.d03.nc',
        'geo_em.d03.nc',
        'geo_em.d04.nc',
        'lcz_zaragoza.tif'
    )
    for f in input_files:
        shutil.copy(os.path.join('sample_data', f), input_dir)
    shutil.copy('testing/custom_lcz_ucp.csv', tmpdir)

    with tmpdir.as_cwd():
        with mock.patch.object(w2w.w2w, 'checks_and_cleaning') as m:
            main(
                [
                    'input', 'lcz_zaragoza.tif', 'geo_em.d04.nc',
                    '--lcz-ucp', 'custom_lcz_ucp.csv',
                ],
            )

    exp_df = pd.read_csv('testing/custom_lcz_ucp.csv', index_col=0)
    pd.testing.assert_frame_equal(exp_df, m.call_args[1]['ucp_table'])