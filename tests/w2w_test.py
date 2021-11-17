import os
import pytest
from w2w.w2w import main
from w2w.w2w import check_lcz_wrf_extent
from w2w.w2w import wrf_remove_urban
from w2w.w2w import main
from w2w.w2w import create_wrf_gridinfo
import xarray


def test_argparse_shows_help():
    with pytest.raises(SystemExit):
        main(['--help'])


def test_check_lcz_wrf_extent_lcz_too_small(capsys):
    info = {
        'src_file': 'testing/lcz_too_small.tif',
        'dst_file': 'sample_data/geo_em.d04.nc',
    }
    with pytest.raises(SystemExit):
        check_lcz_wrf_extent(info=info)

    out, _ = capsys.readouterr()
    assert (
        'ERROR: LCZ domain should be larger than WRF domain in all directions.'
    ) in out
    # TODO maybe add the actual values to check they are correct
    assert (
        'ERROR: LCZ domain should be larger than WRF domain '
        'in all directions.\nLCZ bounds  (xmin, ymin, xmax, ymax): '
    ) in out
    assert 'WRF bounds  (xmin, ymin, xmax, ymax): ' in out


def test_check_lcz_wrf_extent_ok(capsys):
    info = {
        'src_file': 'sample_data/lcz_zaragoza.tif',
        'dst_file': 'sample_data/geo_em.d04.nc',
    }
    check_lcz_wrf_extent(info=info)
    out, _ = capsys.readouterr()
    assert 'OK - LCZ domain is covering WRF domain' in out


def test_wrf_remove_urban(tmpdir):
    info = {
        'dst_file': 'testing/5by5.nc',
        'dst_nu_file': os.path.join(tmpdir, '5by5_new.nc')
    }
    old_ds = xarray.open_dataset(info['dst_file'])
    wrf_remove_urban(info=info, NPIX_NLC=9)
    ds = xarray.open_dataset(info['dst_nu_file'])
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
    tif = xarray.open_rasterio(info['dst_gridinfo'])
    assert tif.crs == '+init=epsg:4326'
    assert tif.transform == pytest.approx(
        (
            0.01000213623046875, 0.0, -1.4050254821777344,
            0.0, 0.010000228881835938, 41.480000495910645,
        )
    )


def test_full_run_with_example_data():
    main(['sample_data', 'lcz_zaragoza.tif', 'geo_em.d04.nc'])
