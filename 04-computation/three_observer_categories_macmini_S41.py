"""The three observer-categories under the affine group (mac-mini-2026-06-22-S41).
A quantity on speed sets splits by what it ignores under translation (=move observer/anchor) and
scaling (=change units):
  1. OBSERVER-RELATIVE  (scale-inv, translation-SENSITIVE): meas(safe)/coverage   -> TILING model (anchored), FINEST
  2. METRIC-DIFFERENCE  (translation-inv, scale-SENSITIVE) : the gap multiset      -> metric winding, the THIRD category
  3. OBSERVER-BLIND/AFFINE (both invariant)               : additive energy, H, order -> TOURNAMENT model, COARSEST
This explains S39 (coverage != additive energy: categories 1 vs 3) and kps S31m (the H-level analogy
breaks: H is category 3, the LRC coverage is category 1 -- two levels finer).
See 07-reflections/three-observer-categories-tiling-is-relative-tournament-is-blind.md."""
from fractions import Fraction as Fr
from collections import Counter

def meas_safe(S):
    iv = []
    for s in S:
        if s == 0: continue
        w = Fr(1, 7) / abs(s)
        for k in range(abs(s)):
            c = Fr(k, abs(s)); lo = (c - w/2) % 1; hi = lo + w
            if hi <= 1: iv.append((lo, hi))
            else: iv.append((lo, Fr(1))); iv.append((Fr(0), hi - 1))
    iv.sort(); tot = Fr(0); clo = chi = None
    for lo, hi in iv:
        if chi is None: clo, chi = lo, hi
        elif lo <= chi: chi = max(chi, hi)
        else: tot += chi - clo; clo, chi = lo, hi
    if chi is not None: tot += chi - clo
    return 1 - tot

def A(S):
    c = Counter()
    for a in S:
        for b in S: c[a + b] += 1
    return sum(v * v for v in c.values())

def diffset(S): return tuple(sorted(a - b for a in S for b in S if a != b))

if __name__ == "__main__":
    base = list(range(1, 14)); tr = [s + 5 for s in base]; sc = [3 * s for s in base]
    print("category 1 OBSERVER-RELATIVE  meas(safe): base=%.4f +5=%.4f x3=%.4f  (scale-inv, transl-SENS)"
          % (float(meas_safe(base)), float(meas_safe(tr)), float(meas_safe(sc))))
    print("category 3 AFFINE             A(E)      : base=%d +5=%d x3=%d  (transl-inv AND scale-inv)"
          % (A(base), A(tr), A(sc)))
    print("category 2 METRIC-DIFFERENCE  diffset== : +5:%s x3:%s  (transl-inv, scale-SENS)"
          % (diffset(base) == diffset(tr), diffset(base) == diffset(sc)))
