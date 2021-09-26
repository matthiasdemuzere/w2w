import pytest
from w2w.w2w import main


def test_argparse_shows_help():
    with pytest.raises(SystemExit):
        main(['--help'])
