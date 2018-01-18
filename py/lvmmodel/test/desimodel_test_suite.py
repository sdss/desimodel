# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
lvmmodel.test.lvmmodel_test_suite
===================================

Used to initialize the unit test framework via ``python setup.py test``.
"""
import unittest


# This is factored out separately from runtests() so that it can be used by
# python setup.py test.
def lvmmodel_test_suite():
    """Returns unittest.TestSuite of lvmmodel tests"""
    from os.path import dirname
    lvmmodel_dir = dirname(dirname(__file__))
    # print(lvmmodel_dir)
    return unittest.defaultTestLoader.discover(lvmmodel_dir,
        top_level_dir=dirname(lvmmodel_dir))


def runtests():
    """Run all tests in lvmmodel.test.test_*.
    """
    # Load all TestCase classes from lvmmodel/test/test_*.py
    tests = lvmmodel_test_suite()
    # Run them
    unittest.TextTestRunner(verbosity=2).run(tests)


if __name__ == "__main__":
    runtests()
