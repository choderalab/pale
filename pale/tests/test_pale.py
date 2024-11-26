"""
Unit and regression test for the pale package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import pale


def test_pale_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "pale" in sys.modules
