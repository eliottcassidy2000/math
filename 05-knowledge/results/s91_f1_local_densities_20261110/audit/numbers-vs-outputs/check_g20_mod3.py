#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: class factors modulo 3 (all three classes) and 9 of H over all 4^n words of G(2,0),
from the local histograms g20_local/full_hist_n{n}.json (not committed; read in place)."""
import json
import sys

import numpy as np

MOD = 720720
d0 = sys.argv[1]
for n in range(8, 17):
    try:
        d = json.load(open(f"{d0}/full_hist_n{n}.json"))
    except OSError:
        print(n, "missing")
        continue
    h = np.zeros(MOD)
    h[np.array(d["residues"], dtype=np.int64)] = np.array(d["counts"], dtype=np.float64)
    tot = h.sum()
    out = [f"n = {n:2d} ({int(tot)} words)"]
    for p in (3, 9, 5, 7):
        c = np.zeros(p)
        np.add.at(c, np.arange(MOD) % p, h)
        out.append(f"mod {p}: " + " ".join(f"{p * x / tot:.4f}" for x in c))
    print(" | ".join(out))
