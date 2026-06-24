#!/usr/bin/env python3
"""
Reference-free alignment-quality scorer.

Takes a set of *rotated* sequences (MARS output), aligns them with an
independent aligner (Clustal Omega), and scores the resulting MSA with a
sum-of-pairs consistency metric. This measures whether a given set of
rotations is 'as alignable' as another -- the relevant quality criterion,
per Blackshields (2008) / Edgar (2014) showing MSA quality is guide-tree
robust even when the specific alignment differs.

Usage:
  python3 score_quality.py clustalo out1.fasta [out2.fasta ...]
"""
import subprocess
import sys
import os
import numpy as np


def read_fasta(path):
    seqs = []
    sid = None
    buf = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None:
                seqs.append("".join(buf))
            sid = line[1:]
            buf = []
        else:
            buf.append(line.upper())
    if sid is not None:
        seqs.append("".join(buf))
    return seqs


def clustalo_align(path, clustalo, tmpdir):
    aln = os.path.join(tmpdir, "aln.fasta")
    subprocess.run([clustalo, "--infile", path, "--outfile", aln,
                    "--force", "--outfmt", "fasta"], check=True,
                   capture_output=True)
    return read_fasta(aln)


def consistency(aln):
    """Sum-of-pairs identity over all columns (gap-gap excluded)."""
    n = len(aln)
    if n < 2:
        return 1.0
    L = len(aln[0])
    total = 0.0
    pairs = 0
    for col in range(L):
        chars = [a[col] for a in aln]
        for i in range(n):
            for j in range(i + 1, n):
                ci, cj = chars[i], chars[j]
                if ci == "-" and cj == "-":
                    continue
                pairs += 1
                if ci == cj and ci != "-":
                    total += 1.0
    return total / pairs if pairs else 0.0


def main():
    clustalo = sys.argv[1]
    tmpdir = "/tmp/opencode/vrf/qtmp"
    os.makedirs(tmpdir, exist_ok=True)
    print("%-40s %8s %8s" % ("file", "consistency", "rel"))
    base = None
    for path in sys.argv[2:]:
        aln = clustalo_align(path, clustalo, tmpdir)
        c = consistency(aln)
        if base is None:
            base = c
        print("%-40s %10.5f %8.4f" % (os.path.basename(path), c, c / base))


if __name__ == "__main__":
    main()
