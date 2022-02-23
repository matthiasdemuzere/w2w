import os
from unittest import mock
import pytest
import rioxarray as rxr
from affine import Affine
import shutil
import w2w.w2w
from w2w.w2w import main
from w2w.w2w import _check_lcz_wrf_extent
from w2w.w2w import wrf_remove_urban
from w2w.w2w import create_wrf_gridinfo
from w2w.w2w import calc_distance_coord
import pandas as pd
import xarray as xr


def test_argparse_shows_help():
    with pytest.raises(SystemExit):
        main(['--help'])


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
    assert 'OK - LCZ domain is covering WRF domain' in out


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


def test_full_run_with_example_data(tmpdir):
    input_files = (
        'geo_em.d01.nc', 'geo_em.d03.nc', 'geo_em.d03.nc','geo_em.d04.nc',
        'lcz_zaragoza.tif',
    )
    input_dir = tmpdir.mkdir('input')
    for f in input_files:
        shutil.copy(os.path.join('sample_data', f), input_dir)

    with tmpdir.as_cwd():
        main(['input', 'lcz_zaragoza.tif', 'geo_em.d04.nc'])

        contents = os.listdir(input_dir)
        assert 'geo_em.d04_LCZ_params.nc' in contents
        assert 'geo_em.d04_NoUrban.nc' in contents
        assert 'geo_em.d04_LCZ_extent.nc' in contents


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
    assert calc_distance_coord(*coords) == pytest.approx(expected, abs=1)


def test_main_with_custom_lcz_ucp(tmpdir):
    input_files = (
        'geo_em.d01.nc', 'geo_em.d03.nc', 'geo_em.d03.nc','geo_em.d04.nc',
        'lcz_zaragoza.tif',
    )
    input_dir = tmpdir.mkdir('input')
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


def test_main_default_lcz_ucp(tmpdir):
    input_files = (
        'geo_em.d01.nc', 'geo_em.d03.nc', 'geo_em.d03.nc','geo_em.d04.nc',
        'lcz_zaragoza.tif',
    )
    input_dir = tmpdir.mkdir('input')
    for f in input_files:
        shutil.copy(os.path.join('sample_data', f), input_dir)

    with tmpdir.as_cwd():
        with mock.patch.object(w2w.w2w, 'checks_and_cleaning') as m:
            main(['input', 'lcz_zaragoza.tif', 'geo_em.d04.nc'])

    exp_df = pd.read_csv('w2w/resources/LCZ_UCP_lookup.csv', index_col=0)
    pd.testing.assert_frame_equal(exp_df, m.call_args[1]['ucp_table'])
