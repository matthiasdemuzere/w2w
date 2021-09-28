import pytest
from w2w.w2w import main
from w2w.w2w import check_lcz_wrf_extent


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
