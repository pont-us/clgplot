#!/usr/bin/env python3

import unittest
from clgplot import gradient
from clgplot import x_for_half_max_y


class TestClgPlot(unittest.TestCase):

    def test_gradient(self):
        result = gradient([1,2],[5,7])
        self.assertEqual(2, result[0])
        self.assertEqual(2, result[1])

    def test_x_for_half_max_y(self):
        self.assertEqual(0, x_for_half_max_y([0, 0, 0], [0, 0, 0]))
        self.assertEqual(5, x_for_half_max_y([1, 2, 5, 16, 19, 99],
                                             [1, 1, 9, 11, 18, 14]))
        self.assertAlmostEqual(2.5, x_for_half_max_y([0, 2, 3, 8],
                                                     [1, 2, 4, 6]),
                               places=10)
        self.assertAlmostEqual(2.75, x_for_half_max_y([0, 2, 3, 8],
                                                      [0, 1, 5, 8]),
                               places=10)


if __name__ == "__main__":
    unittest.main()
