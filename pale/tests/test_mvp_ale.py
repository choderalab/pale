"""
Unit and regression test for the mvp_ale package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import mvp_ale


def test_mvp_ale_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "mvp_ale" in sys.modules
