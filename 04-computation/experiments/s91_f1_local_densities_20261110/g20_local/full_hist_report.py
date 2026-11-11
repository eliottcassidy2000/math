#!/usr/bin/env python3
"""S91 (HYP-9660 (D)): residue classes of H over ALL 4^n words of the two-apex grammar G(2,0), n <= 16.

Input: full_hist_n{n}.json written by full_hist.py (the helper's exhaustive histograms of H modulo 720720 =
16 * 9 * 5 * 7 * 11 * 13; n = 16 is 4,294,967,296 words).  The histograms themselves are not kept in the repository
(regenerate with `python3 full_hist.py 16`, about 5 minutes); this script's output is.

Printed for each n: the class factors  p * share(H = c mod p)  at c = 0 and the extreme factors among c != 0, for
p = 3, 5, 7, 11, 13 and for 9; and the extreme factors (M/2) * share over the odd classes modulo M = 8 and 16.
A factor 1 means no bias.  For comparison the one-apex grammar F1 inside G(2,0) has factors 1 - 1/p^2 + 1/p at its
distinguished class in the limit (THM-5128).

Usage: python3 full_hist_report.py [NMAX]
"""
import json
import os
import sys

import numpy as np

MOD = 720720


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    here = os.path.dirname(os.path.abspath(__file__))
    rows = []
    for n in range(4, nmax + 1):
        path = os.path.join(here, f"full_hist_n{n}.json")
        if not os.path.exists(path):
            continue
        d = json.load(open(path))
        h = np.zeros(MOD, dtype=np.float64)
        h[np.array(d["residues"], dtype=np.int64)] = np.array(d["counts"], dtype=np.float64)
        tot = h.sum()
        assert abs(tot - 4.0 ** n) < 0.5 and h[0::2].sum() == 0
        line = [f"{n:2d}"]
        for p in (3, 5, 7, 11, 13, 9):
            c = np.zeros(p)
            np.add.at(c, np.arange(MOD) % p, h)
            f = p * c / tot
            line.append(f"{f[0]:.4f} [{f[1:].min():.4f}, {f[1:].max():.4f}]")
        for m in (8, 16):
            c = np.zeros(m)
            np.add.at(c, np.arange(MOD) % m, h)
            f = (m / 2) * c[1::2] / tot
            line.append(f"[{f.min():.4f}, {f.max():.4f}]")
        rows.append(line)
    hdr = ["n"] + [f"p={p}: L(0) [min, max of L(c != 0)]" for p in (3, 5, 7, 11, 13, 9)] + ["odd classes mod 8", "mod 16"]
    print(" | ".join(hdr))
    for line in rows:
        print(" | ".join(line))


if __name__ == "__main__":
    main()
