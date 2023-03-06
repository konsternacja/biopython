def repeats(sequence, n):
    counts = {}
    for j in range(n):
        # shift reading frame
        seq = sequence[j:]
        subseqs = [
            seq[i : i + n]
            for i in range(0, len(seq), n)
            if len(seq[i : i + n]) == n
        ]
        print(subseqs)
        subseqs = set(subseqs)
        for s in subseqs:
            counts[s] = seq.count(s)
    counts = dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))
