# Copyright 2023 by Filip Hajdyła.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

""" Tests for FindInSeq module """

import unittest

from Bio.SeqUtils.FindInSeq import FindRepeats

s = "ATGCATGCATGC"


class test_FindRepeats(unittest.TestCase):

    def test_find_all_subseqs(self):
        # no optional arguments
        r1 = FindRepeats(s)
        self.assertEqual(r1.subseqs, ["ATG",
                                      "CAT",
                                      "GCA",
                                      "TGC"])
        # custom length
        r2 = FindRepeats(s, length=4)
        self.assertEqual(r2.subseqs, ["ATGC",
                                      "CATG",
                                      "GCAT",
                                      "TGCA"])

    def test_count_subseqs(self):
        # no optional arguments
        r1 = FindRepeats(s)
        self.assertEqual(r1.counts, {"ATG": 3,
                                     "CAT": 2,
                                     "GCA": 2,
                                     "TGC": 3})

        # descending sort
        r2 = FindRepeats(s, sort="d")
        self.assertEqual(r2.counts, {"ATG": 3,
                                     "TGC": 3,
                                     "CAT": 2,
                                     "GCA": 2})

        # ascending sort
        r3 = FindRepeats(s, sort="a")
        self.assertEqual(r3.counts, {"CAT": 2,
                                     "GCA": 2,
                                     "ATG": 3,
                                     "TGC": 3})
        # threshold = 3
        r4 = FindRepeats(s, threshold=3)
        self.assertEqual(r4.counts, {"ATG": 3,
                                     "TGC": 3})


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
