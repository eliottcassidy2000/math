#!/usr/bin/env python3
"""
monad-explorer-2026-06-06-S701
LRC: is the self-conjugate runner n/2 LOAD-BEARING in AP_n (even n)?

Builds DIRECTLY on S700 (HYP-2259): for even n the negation involution
sigma:a->-a on Z/n has fixed set {0, n/2}; the second apex n/2 is the AP's own
self-conjugate runner, slack (distance 1/2) at every tight DIVISION-GRID time.
S700's handoff conjecture was: n/2 stays slack OFF the grid, so it can be deleted,
reducing even-n LRC to the odd-n single-apex residual.

This script tests that OFF the grid using exact arithmetic (M_exact scans the full
finite candidate set of breakpoints: half-integer-over-v local maxima + pairwise
pinch times m/(v_i +- v_j)).  Definitions:

  A runner r in S is LOAD-BEARING iff  M(S) < M(S \ {r})  (deleting it strictly
  raises the gap, i.e. r is forced below threshold at the reduced optimum).

  At a given argmax t* of S, r is a BINDER iff ||r t*|| = M(S) (it attains the min).
  r is a UNIQUE binder iff it is the only runner attaining the min at t*.

Questions answered exactly for even n = 4,6,...,18:
 Q1. M(AP_n) (expect 1/n) and M(AP_n \ {n/2}).  Is n/2 load-bearing?
 Q2. At EVERY argmax t* of AP_n, is n/2 ever a (unique) binder?  [grid claim S700:
     at the tight grid times it is slack; here we sweep the FULL argmax set.]
 Q3. The reduced set R_n = {1,..,n-1}\{n/2} is negation-symmetric with NO fixed
     point.  What standard fraction is M(R_n)?  Compare 1/n, 1/(n-1), 2/(2n-1),...
 Q4. Where is the argmax of R_n, and what are its binders?  Does n/2 sit slack there
     (consistent with deletability) or below threshold (load-bearing)?
 Q5. Cross-check on the near-floor families V* (AP with 12->24 at n=14) and the
     doubled-AP atoms: is n/2 ever the UNIQUE binder at their optima?
"""
from fractions import Fraction
from math import gcd


def frac(x):
    r = x % 1
    return min(r, 1 - r)


def candidate_times(V):
    """Full finite breakpoint set for max_t min_i ||v_i t|| (exact)."""
    c = set()
    for v in V:
        v = abs(v)
        for k in range(2 * v):                 # half-integer/v: local maxima of ||vt||
            c.add(Fraction(2 * k + 1, 2 * v) % 1)
    for i in range(len(V)):
        for j in range(len(V)):
            for s in (1, -1):
                d = V[i] + s * V[j]            # pinch: two runners equidistant
                if d:
                    d = abs(d)
                    for k in range(d + 1):
                        c.add(Fraction(k, d) % 1)
    c.discard(Fraction(0))
    return c


def M_argmax_binders(V):
    """Return (M, argmax_times sorted, {t: sorted binder speeds}) exactly."""
    V = list(V)
    best = Fraction(0)
    vals = {}
    for t in candidate_times(V):
        mn = min(frac(v * t) for v in V)
        vals[t] = mn
        if mn > best:
            best = mn
    argmax = sorted(t for t, m in vals.items() if m == best)
    binders = {}
    for t in argmax:
        binders[t] = sorted(v for v in V if frac(v * t) == best)
    return best, argmax, binders


def analyze(n):
    AP = list(range(1, n))
    half = n // 2
    R = [v for v in AP if v != half]           # AP_n minus the self-conjugate runner

    M_ap, arg_ap, bind_ap = M_argmax_binders(AP)
    M_r, arg_r, bind_r = M_argmax_binders(R)

    # Is n/2 a binder at ANY argmax of AP_n?  (S700 grid claim says slack there.)
    half_ever_binds = any(half in bind_ap[t] for t in arg_ap)
    half_unique_binder = any(bind_ap[t] == [half] for t in arg_ap)

    # Load-bearing test
    load_bearing = M_ap < M_r

    # At the reduced optimum, where does the n/2 runner sit?
    half_dists_at_R_argmax = sorted({frac(half * t) for t in arg_r})

    return {
        "n": n, "M_ap": M_ap, "argmax_ap": arg_ap, "bind_ap": bind_ap,
        "M_r": M_r, "argmax_r": arg_r, "bind_r": bind_r,
        "half": half, "half_ever_binds_AP": half_ever_binds,
        "half_unique_binder_AP": half_unique_binder,
        "load_bearing": load_bearing,
        "half_dists_at_R_argmax": half_dists_at_R_argmax,
    }


def fmt_t(ts):
    return ", ".join(str(t) for t in ts)


if __name__ == "__main__":
    print("=" * 78)
    print("S701: Is the self-conjugate runner n/2 load-bearing in AP_n? (even n)")
    print("=" * 78)
    print("Convention: AP_n = {1,...,n-1}, n-1 runners, conjectured floor M=1/n.")
    print("R_n = AP_n \\ {n/2} = negation-symmetric, n-2 runners, floor >= 1/(n-1).")
    print()
    rows = []
    for n in range(4, 20, 2):
        d = analyze(n)
        rows.append(d)
        print(f"--- n = {n}  (n/2 = {d['half']}) ---")
        print(f"  M(AP_n)            = {d['M_ap']}   (1/n = {Fraction(1,n)})  "
              f"argmax t = {{{fmt_t(d['argmax_ap'])}}}")
        # show binders at each AP argmax
        for t in d['argmax_ap']:
            mark = "  <== n/2 BINDS" if d['half'] in d['bind_ap'][t] else ""
            print(f"      t={t}: binders={d['bind_ap'][t]}{mark}")
        print(f"  n/2 ever a binder of AP_n? {d['half_ever_binds_AP']}    "
              f"unique binder? {d['half_unique_binder_AP']}")
        print(f"  M(R_n)=M(AP_n\\{{{d['half']}}}) = {d['M_r']} = {float(d['M_r']):.6f}"
              f"   (1/(n-1)={Fraction(1,n-1)}={float(Fraction(1,n-1)):.6f})")
        print(f"      argmax t = {{{fmt_t(d['argmax_r'])}}}")
        for t in d['argmax_r'][:6]:
            print(f"        t={t}: binders={d['bind_r'][t]}   "
                  f"||(n/2)t||={frac(d['half']*t)}")
        print(f"  n/2 LOAD-BEARING (M(AP)<M(R))? {d['load_bearing']}    "
              f"dist of n/2 at R-optima: {d['half_dists_at_R_argmax']}")
        print()

    print("=" * 78)
    print("SUMMARY TABLE")
    print("=" * 78)
    print(f"{'n':>3} {'M(AP_n)':>10} {'M(R_n)':>12} {'M(R_n)~':>10} "
          f"{'1/(n-1)':>9} {'load?':>6} {'n/2 binds AP?':>13}")
    for d in rows:
        print(f"{d['n']:>3} {str(d['M_ap']):>10} {str(d['M_r']):>12} "
              f"{float(d['M_r']):>10.5f} {str(Fraction(1,d['n']-1)):>9} "
              f"{str(d['load_bearing']):>6} {str(d['half_ever_binds_AP']):>13}")

    # ---- Cross-check on near-floor families at n=14 ----
    print()
    print("=" * 78)
    print("CROSS-CHECK at n=14: is runner 7 = n/2 ever the UNIQUE binder?")
    print("=" * 78)
    fams = {
        "AP_14            ": list(range(1, 14)),
        "V* (AP, 12->24)  ": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24],
        "2AP (doubled)    ": [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26],
    }
    for tag, V in fams.items():
        if 7 not in V and 14 not in V:
            seven = None
        M, arg, bind = M_argmax_binders(V)
        # which runner plays the self-conjugate role: 7 in AP/V*, 14 in 2AP (=2*7)
        sc = 7 if 7 in V else (14 if 14 in V else None)
        uniq = any(bind[t] == [sc] for t in arg) if sc else None
        ever = any(sc in bind[t] for t in arg) if sc else None
        print(f"  {tag} M={M}={float(M):.5f}  sc-runner={sc}  "
              f"sc ever binds? {ever}  sc UNIQUE binder? {uniq}")
        for t in arg[:4]:
            print(f"      t={t}: binders={bind[t]}")
