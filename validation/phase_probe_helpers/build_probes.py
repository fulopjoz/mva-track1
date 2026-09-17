#!/usr/bin/env python3
"""Build a non-redundant 31-mer probe set from a reference interval (FASTA).

Output: one probe per line to stdout (used with `seqkit grep -f <file>`).
Only kmers occurring exactly once within the interval are emitted (avoids
low-complexity duplication); the reference interval is expected to be unique
in the genome for the reads we target (BUB1B / MVA regulator region).
"""
import sys

K = 31

def main(region_fa: str, out: str):
    seq = []
    with open(region_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    seq = "".join(seq)

    counts: dict[str, int] = {}
    for i in range(len(seq) - K + 1):
        kmer = seq[i:i + K]
        counts[kmer] = counts.get(kmer, 0) + 1

    seen = set()
    with open(out, "w") as fh:
        for i in range(len(seq) - K + 1):
            kmer = seq[i:i + K]
            if counts[kmer] == 1 and kmer not in seen:
                seen.add(kmer)
                fh.write(kmer + "\n")
    sys.stderr.write(f"interval len={len(seq)} unique-kmers-written={len(seen)}\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
