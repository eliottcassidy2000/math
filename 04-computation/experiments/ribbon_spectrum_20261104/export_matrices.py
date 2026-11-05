#!/usr/bin/env python3
"""export_matrices.py -- write the span-w profile transfer system to a text file for the C enumerator.
Format: d nletters nfinal / final indices / per letter: mask(bits s-1 for s in letter) nnz then triples j i m."""
import sys
from ribbon import transfer_system, letters
w = int(sys.argv[1]) if len(sys.argv) > 1 else 2
order, states, mats, final = transfer_system(w)
L = letters(w)
with open(f"span{w}_matrices.txt", "w") as f:
    f.write(f"{len(order)} {len(L)} {len(final)}\n")
    f.write(" ".join(map(str, final)) + "\n")
    for a in L:
        mask = sum(1 << (s - 1) for s in a)
        items = sorted(mats[a].items())
        f.write(f"{mask} {len(items)}\n")
        for (j, i), m in items:
            f.write(f"{j} {i} {m}\n")
print("written", f"span{w}_matrices.txt", len(order), "states")
