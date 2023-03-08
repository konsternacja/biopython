def repeats(sequence, n):
    # TODO refactor it to class
    counts = {}
    for j in range(n):
        seq = sequence[j:]  # shift the reading frame
        subseqs = set(
            seq[i : i + n] for i in range(0, len(seq), n)
            if len(seq[i : i + n]) == n
        )
        counts.update({subseq: seq.count(subseq) for subseq in subseqs})
    return dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))


r = repeats("ATGATATATTTTTATTATTGTAGTTAGTAGTAGTA", 3)
print(r)