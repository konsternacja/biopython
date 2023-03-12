class FindRepeats:
    """
    A class that finds all possible subsequences of a given length in all
    reading frames and counts the number of appearances of each subsequence
    in a given sequence for values equal or above a given threshold.

    Args:
        sequence (Seq object): The input sequence to search.
        rep_len (int): The length of the desired substring.
        threshold (int): The minimum number of appearances of a substring
                         to be counted. Threshold affects `self.counts`
                         output but does not trim off `self.subseqs` set.

    Attributes:
        sequence (Seq object): The input sequence to search.
        rep_len (int): The length of the desired substring.
        threshold (int): The minimum number of appearances of a substring
                         to be counted. Threshold affects `self.counts`
                         output but does not trim off `self.subseqs` set.
        subseqs (set): The set of all possible subsequences of
                       `length == rep_len` in all reading frames.
        counts (dict): The dictionary containing the subsequence (key) and
                       their number of appearance (value) in sequence after
                       being filtered by threshold.

    Methods:
        `_find_all_subseqs(sequence: str, rep_len: int) -> set:`
            Finds all possible subsequences of
            `length == rep_len` in all reading frames.

        `_count_subseqs(sequence: str, subseqs: set, rep_len: int, threshold: int) -> dict:`
            Finds number of appearances of all `subseqs` in `sequence`
            Creates a sorted (descending) and filtered
            (`value >= threshold`) `counts` dictionary.


    """
    def __init__(self, sequence, rep_len=3, threshold=0):
        """
        Initializes a FindRepeats object.

        Args:
        - sequence (Seq object): DNA/RNA/Protein input sequence object or
                                 string.
        - rep_len (int): Desired substring length. Default is 3.
        - threshold (int): Minimum number of appearances of subsequence in
                           sequence. Default is 0.
        """
        self.sequence = sequence.upper()
        self.rep_len = rep_len
        self.threshold = threshold
        self.subseqs = self._find_all_subseqs(self.sequence,
                                              self.rep_len)
        self.counts = self._count_subseqs(self.sequence,
                                          self.subseqs,
                                          self.rep_len,
                                          self.threshold)

    def _find_all_subseqs(self, sequence: object, rep_len: int) -> set:
        """Finds all possible subsequences of `length == rep_len`
        in all reading frames.

        Args:
            sequence (Seq object): DNA/RNA/Protein input sequence object or string.
            rep_len (int): Desired substring length.

        Returns:
            set: Set of all possible subsequences of `length == rep_len`
            in all reading frames.
        """
        for j in range(rep_len):
            seq = sequence[j:]  # shift the reading frame
            subseqs = set(
                seq[i : i + rep_len]
                for i in range(0, len(seq), rep_len)
                if len(seq[i : i + rep_len]) == rep_len
            )
        return subseqs

    def _count_subseqs(self, sequence: object, subseqs: set, rep_len: int,
                       threshold: int) -> dict:
        """Finds number of appearances of all `subseqs` in `sequence`.
        Creates a sorted (descending) and filtered (value >= `threshold`)
        `counts` dictionary.

        Args:
            sequence (Seq object): DNA/RNA/Protein input sequence object or string.
            subseqs (set): Subsequence set created with `_find_all_subseqs()`
            rep_len (int): Desired substring length.
            threshold (int): Minimum number of appearances of subsequence in sequence.

        Returns:
            dict: sorted and filtered dictionary containing
            subsequence (key) and number of their appearance
            (value) in sequence.
        """

        # TODO repeats are now overlapping
        # implement overlap=True/False modes

        counts = {}
        for j in range(rep_len):
            seq = sequence[j:]
            counts.update({subseq: seq.count(subseq) for subseq in subseqs})

        # sort by value (descending)
        counts = dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))

        # filter out reps below threshold
        def _filter(pair) -> bool:
            key, value = pair
            if value >= threshold:
                return True
            else:
                return False
        counts = dict(filter(_filter, counts.items()))
        return counts


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
