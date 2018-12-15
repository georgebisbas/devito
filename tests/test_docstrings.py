# Python provides the doctest module, which "searches for pieces of text that look
# like interactive Python sessions, and then executes those sessions to verify that
# they work exactly as shown".
#
# pytest integrates doctest: https://docs.pytest.org/en/latest/doctest.html
#
# Some good reasons to have a ``test_docstrings.py`` rather than calling
# ``pytest --doctest-modules`` directly:
# * inclusion in code coverage
# * skipping tests when using a devito backend (where they would fail, for
#   the most disparate reasons)

from importlib import import_module

import pytest
import doctest

from conftest import skipif

pytestmark = skipif(['yask', 'ops'])


@pytest.mark.parametrize('modname', [
    'dimension', 'equation', 'function', 'grid', 'operator',
    'data.decomposition', 'finite_differences.finite_difference',
    'ir.support.space'
])
def test_docstrings(modname):
    module = import_module('devito.%s' % modname)
    assert doctest.testmod(module).failed == 0
