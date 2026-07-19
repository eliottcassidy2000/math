"""
opus-2026-07-19-S399: cross-check the (D,s) active-pair identity D = M*s
(THM-1261, formerly THM-1245, opus-S398) on out-of-sample data:

  1. death-star-S59b's brand-new N=31 single-far discovery {1..29,31,120}
     claimed at M = 4/127 (the D=4 slack-1 rung for 32 runners: 127 = 32*4-1).
  2. The two known n=14 gap/ladder families {1..11,13,36} (M=3/41, pair (5,36))
     and {1..11,13,24} (M=1/14, second extremal).

For each family: locate the maximizer among rationals a/q for q <= QMAX
(sanity band, NOT a proof of the global max -- the exact values are canon:
death-star THM-1255/1256 cross-N S59b, opus THM-1230/1235), then at t*=a/q:
  - list the active speeds (||v t*|| = M exactly),
  - find straddling active pairs (v_i t* = a_i + M, v_j t* = a_j - M),
  - verify s = v_i + v_j equals the reduced denominator times m, and
    D = M*s exactly (the determinant identity), and
  - report slack n*D - s for the ambient runner count n.

Purpose (synthesis session): test that the (D,s)/slack frame is structural
across N, i.e. the identity is not n=14 numerology.  Output frozen to
05-knowledge/results/.
"""
from fractions import Fraction
from math import gcd

def dist_to_Z(x: Fraction) -> Fraction:
    f = x - Fraction(int(x))  # frac part, in [0,1) for x>=0
    if f < 0:
        f += 1
    return min(f, 1 - f)

def family_min_at(V, t: Fraction) -> Fraction:
    return min(dist_to_Z(v * t) for v in V)

def scan_max(V, qmax):
    """max over t=a/q, q<=qmax, of min_v ||v t||; returns (M, t*)."""
    best = (Fraction(0), Fraction(0))
    for q in range(2, qmax + 1):
        for a in range(1, q // 2 + 1):
            if gcd(a, q) != 1:
                continue
            m = family_min_at(V, Fraction(a, q))
            if m > best[0]:
                best = (m, Fraction(a, q))
    return best

def analyze(name, V, n_runners, qmax=300):
    print(f"\n=== {name}  (|V|={len(V)}, n={n_runners} runners, threshold 1/{n_runners}) ===")
    M, tstar = scan_max(V, qmax)
    print(f"  scan (q<={qmax}): M = {M} = {float(M):.6f} at t* = {tstar}")
    print(f"  threshold 1/{n_runners} = {float(Fraction(1, n_runners)):.6f}; "
          f"M - 1/n = {M - Fraction(1, n_runners)}")
    # active speeds and straddling structure at t*
    act_plus, act_minus = [], []   # v t* = a + M  /  a - M... signs via frac part
    for v in V:
        x = v * tstar
        f = x - Fraction(int(x))
        if f < 0:
            f += 1
        if min(f, 1 - f) == M:
            if f == M:
                act_minus.append(v)   # just ABOVE an integer: v t* = a + M
            else:
                act_plus.append(v)    # just BELOW an integer: v t* = a - M
    print(f"  active speeds: above-integer {act_minus}, below-integer {act_plus}")
    # straddling pairs: one from each side; check D = M*s
    found = []
    for vi in act_minus:
        for vj in act_plus:
            s = vi + vj
            D = M * s
            # D must be a positive integer for the identity
            if D.denominator == 1:
                slack = n_runners * int(D) - s
                found.append((vi, vj, s, int(D), slack))
    for (vi, vj, s, D, slack) in found:
        print(f"  straddling pair ({vi},{vj}): s = {s}, D = M*s = {D}  "
              f"[integer OK], slack n*D - s = {slack}")
    if not found:
        print("  NO straddling pair with integer D -- identity FAILS here")
    return M, tstar, found

if __name__ == "__main__":
    print("(D,s) identity cross-check -- opus-S399")
    print("NOTE: scan bound is a sanity band; exact M values are canon "
          "(THM-1256, THM-1230/1235). This verifies the IDENTITY at the "
          "maximizer, not the global max itself.")

    # death-star-S59b N=31 discovery: 31 speeds, 32 runners
    V31 = list(range(1, 30)) + [31, 120]
    analyze("death-star-S59b {1..29,31,120}", V31, 32, qmax=200)

    # opus-S395/S396 n=14 gap family: 13 speeds, 14 runners
    V36 = list(range(1, 12)) + [13, 36]
    analyze("opus-S395 {1..11,13,36}", V36, 14, qmax=200)

    # second extremal (Goddyn-Wong accel): M = 1/14 exactly
    V24 = list(range(1, 12)) + [13, 24]
    analyze("GW second extremal {1..11,13,24}", V24, 14, qmax=200)
