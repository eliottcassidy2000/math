#!/usr/bin/env python3
"""
regular7_doubling_Htest_kps1.py — kind-pasteur-2026-06-09-S1
BRANCH E follow-up (HYP-2337): does the H-coincidence H(T[K2]) = H(SCblow) seen at
the regular classes (C_3 at n=3 with 45=45, regular T_5 with 15565=15565, plus iso
verified there) extend to the THREE regular n=7 classes (doubles on 14 vertices)?
Context: among n=5 SC classes ONLY the regular one had H(K2)=H(SC) (and was iso);
all non-regular SC classes split. Hypothesis under test:
  HYP-E: T regular => H(T[K2](T)) = H(SCblow(T))  (and conjecturally iso).

The 3 regular n=7 classes (S18h/THM-027): alpha_2=7 (BIBD/Paley, H=189),
alpha_2=10 (H=171), alpha_2=14 (H=175). Found by rejection sampling + canon.

Output: 05-knowledge/results/regular7_doubling_Htest_kps1.out
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from skew_doubling_core_kps1 import (canon, H_count, M_of, A_of, scores,
                                     D_skew, D_blowup, D_scblowup, is_DRT)

OUTPATH = os.path.join(ROOT, "05-knowledge", "results",
                       "regular7_doubling_Htest_kps1.out")


def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.integers(2):
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A


def main():
    out = open(OUTPATH, "w", encoding="utf-8")

    def w(line=""):
        out.write(line + "\n")
        out.flush()
        print(line)

    w("=== regular7_doubling_Htest_kps1 — H(K2) vs H(SCblow) vs H(Dskew), regular n=7 ===")
    rng = np.random.default_rng(20260609)
    found = {}
    tries = 0
    t0 = time.time()
    while len(found) < 3 and tries < 2_000_000:
        tries += 1
        A = random_tournament(7, rng)
        if scores(A) != (3, 3, 3, 3, 3, 3, 3):
            continue
        C = canon(A)
        key = C.tobytes()
        if key not in found:
            found[key] = C
            w("found regular class %d after %d tries (H=%d, DRT=%s)" %
              (len(found), tries, H_count(C), is_DRT(C)))
    w("sampling time %.1fs; %d distinct regular n=7 classes" %
      (time.time() - t0, len(found)))
    w("")
    w("%6s %4s | %12s %12s %12s | K2==SC?" % ("H(T)", "DRT", "H(Dskew)", "H(TK2)", "H(SCbl)"))
    for C in sorted(found.values(), key=lambda c: H_count(c)):
        hT = H_count(C)
        t1 = time.time()
        hD = H_count(D_skew(C)[0])
        hK = H_count(D_blowup(C)[0])
        hS = H_count(D_scblowup(C)[0])
        w("%6d %4s | %12d %12d %12d | %s   (%.1fs)" %
          (hT, is_DRT(C), hD, hK, hS, hK == hS, time.time() - t1))
    w("")
    w("=== done ===")
    out.close()


if __name__ == "__main__":
    main()
