#!/usr/bin/env python3
"""
Compute trivially derived OEIS sequences from existing b-files.

A007869 = (A000088 + A000171) / 2  — graphs with even number of edges
A054928 = (A000273 + A003086) / 2  — digraphs with even number of arcs

Author: opus-2026-03-08-S48
"""

import sys
sys.path.insert(0, '.')

def load_bfile(fname):
    """Load b-file into dict {n: value}."""
    d = {}
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                d[int(parts[0])] = int(parts[1])
    return d


def save_bfile(d, fname, start_n=None):
    """Save dict to b-file."""
    keys = sorted(d.keys())
    if start_n is not None:
        keys = [k for k in keys if k >= start_n]
    with open(fname, 'w') as f:
        for k in keys:
            f.write(f"{k} {d[k]}\n")
    print(f"Written {fname} with {len(keys)} entries (n={keys[0]}..{keys[-1]})")


def derive_average(fname_a, fname_b, fname_out, name, start_n=1):
    """Compute (A + B) / 2."""
    a = load_bfile(fname_a)
    b = load_bfile(fname_b)

    # Find common range
    min_n = max(min(a.keys()), min(b.keys()), start_n)
    max_n = min(max(a.keys()), max(b.keys()))

    result = {}
    for n in range(min_n, max_n + 1):
        if n in a and n in b:
            s = a[n] + b[n]
            assert s % 2 == 0, f"n={n}: sum {s} is odd!"
            result[n] = s // 2

    print(f"\n=== {name} ===")
    print(f"Range: n={min_n}..{max_n} ({len(result)} terms)")
    for n in range(min_n, min(min_n + 10, max_n + 1)):
        print(f"  a({n}) = {result[n]}")

    save_bfile(result, fname_out, start_n=min_n)
    return result


if __name__ == "__main__":
    # A007869 = (A000088 + A000171) / 2
    # A000088 starts at n=0, A000171 starts at n=1
    # Need to align: A000171(n) = self-comp graphs on n nodes, A000088(n) = all graphs on n nodes
    derive_average('b000088.txt', 'b000171.txt', 'b007869.txt',
                   'A007869: graphs with even number of edges', start_n=1)

    # A054928 = (A000273 + A003086) / 2
    # A000273 starts at n=0, A003086 starts at n=1
    derive_average('b000273.txt', 'b003086.txt', 'b054928.txt',
                   'A054928: digraphs with even number of arcs', start_n=1)
