#!/usr/bin/env python3
"""Compute GC content for each sequence in a FASTA file.

Usage:
    python gc_content.py example.fasta
"""
import sys


def read_fasta(path):
    """Yield (name, sequence) pairs from a FASTA file."""
    name, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.upper())
    if name is not None:
        yield name, "".join(parts)


def gc_content(seq):
    """Return GC fraction of a sequence (0-1). Ignores N and other ambiguous bases."""
    counted = [b for b in seq if b in "ACGT"]
    if not counted:
        return 0.0
    return sum(b in "GC" for b in counted) / len(counted)


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    gc_values = []
    for name, seq in read_fasta(sys.argv[1]):
        gc = gc_content(seq)
        gc_values.append(gc)
        at = gc_content(seq.translate(str.maketrans("ACGT", "CATG")))
        print(f"{name}\tLength = {len(seq)}\tGC = {100 * gc:.1f}%\tAT = {100 * at:.1f}%")
    mean_gc = sum(gc_values) / len(gc_values) if gc_values else 0.0
    print(f"Mean GC = {100 * mean_gc:.1f}%")


if __name__ == "__main__":
    main()
