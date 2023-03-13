# Copyright 2023 by Filip Hajdyła.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

class FindRepeats:
    """
    A class that is used to find and count all subsequences within the
    sequence. It creates the set of all possible subsequences and then
    uses it to count non-overlapping repeats of each subsequence in the set.
    Length of subsequences can be adjusted using `length` argument as
    described below. The output dictionary can be filtered by

    Args:
        - sequence (Seq object): The input sequence to search in.
        - length (int): The length of subsequences that are to be found. (default 3)
        - sort (str): Determines if and how (descending/ascending) `self.counts`
                      is to be sorted. It is sorted primarily by value and secondarily
                      by key in alphabetical order.
        - threshold (int): The minimum number of subsequence counts in
                           the sequence. Subsequences that have number
                           of counts below threshold will be filtered out
                           from the output dictionary. (default 0 -> no filtering)

    Attributes:
        - subseqs (set): The set of all possible subsequences of
                         `length == length`. Subsequences are searched for
                         in all reading frames which means they do overlap.
        - counts (dict): The dictionary containing the subsequences (keys) and
                         their counts (values). This dictionary can be filtered
                         using `threshold` argument or sorted using `sort` argument.

    ! ONLY FOR INTERNAL USE !

    TODO

    Methods:
        `_find_all_subseqs(sequence: str, length: int) -> set:`
            Finds all possible subsequences of
            `length == length` in all reading frames.

        `_count_subseqs(sequence: str, subseqs: set, length: int, threshold: int) -> dict:`
            Finds number of appearances of all `subseqs` in `sequence`
            Creates a sorted (descending) and filtered
            (`value >= threshold`) `counts` dictionary.
    """

    def __init__(self, sequence, length=3, threshold=0, sort=None):

        self._sequence = sequence.upper()
        self._length = length
        self._threshold = threshold
        self._sort = sort
        self.subseqs = self._find_all_subseqs(self._sequence,
                                              self._length)
        self.counts = self._count_subseqs(self._sequence,
                                          self.subseqs,
                                          self._sort,
                                          self._threshold)

    def _find_all_subseqs(self, sequence: object, length: int) -> set:

        subseqs = []
        for j in range(length):
            seq = sequence[j:]
            subseqs += [seq[i : i + length]
                        for i in range(0, len(seq), length)
                        if len(seq[i : i + length]) == length]
        return sorted(list(set(subseqs)))

    def _count_subseqs(self, sequence: object, subseqs: set,
                       sort: str, threshold: int) -> dict:

        counts = {subseq: sequence.count(subseq) for subseq in subseqs}

        if sort in ["d", "descending"]:
            sort_mode = True
        elif sort in ["a", "ascending"]:
            sort_mode = False
        if sort:
            counts = dict(sorted(counts.items(),
                                 key=lambda x: x[1],
                                 reverse=sort_mode))

        def _filter(pair) -> bool:
            key, value = pair
            if value >= threshold:
                return True
            else:
                return False
        counts = dict(filter(_filter, counts.items()))
        return counts

# TODO:
# class FindORFs
# class FindSNPs
# class FindSubseqs (maybe some other name?)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
