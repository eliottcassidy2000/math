#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spectrum_complement_even_kpswf13.py   (kind-pasteur 2026-06-22-S36, THREAD 2)

THE COMPLEMENT-EVEN REFRAME of the Node-3 spectrum sum (HYP-2867).

SETUP (exact rationals; reuses the kpswf12 arc machinery):
  G_P    = {x : ||p x|| >= 1/14  for all p in P}          (witness/lonely arc set)
  GOOD   = covc(E) = cover_set(E)^c = {x : >=1 of the 7 sectors [r/7,(r+1)/7) is MISSED
                                       by {frac(e_i x): e_i in E}}    (1/n-decay, 7-vanishing)
  chat(n)= Fourier coeffs of 1_{G_P}   (REAL: G_P = -G_P)
  ghat(n)= Fourier coeffs of 1_{GOOD}  (REAL: GOOD = -GOOD)
  SPEC   = sum_{n!=0} chat(n) conj(ghat(n)),   R' = 1 + SPEC/baseline,  baseline=meas(G_P)meas(GOOD)

THE COMPLEMENT INVOLUTION (the tournament T^op = v->N-v).  Two distinct circle-level actions:
  (A) x -> -x  on the CIRCLE.  Both G_P, GOOD are invariant => chat, ghat REAL.  (already known)
  (B) e -> N-e on the CLUSTER SPEEDS (the v->N-v of HYP-2867).  This is the new axis.

KEY EXACT FACT we exploit:
  frac((N-e) x) at x = j/N (the N-point Bohr grid) equals frac(-e j/N) = 1 - frac(e j/N).
  So on the grid (1/N)Z, the speed map e->N-e induces the phase reflection ph -> 1-ph = (-ph),
  which is a SECTOR REFLECTION r -> 6-r (sectors of [0,1) into 7 parts, r=0..6 -> 6-r, with the
  boundary sector handled by the half-open convention).  Hence the sector-COVER predicate
  "all 7 sectors hit" is INVARIANT under e->N-e ON THE GRID.  Off the grid the two clusters
  E and N-E give literally the SAME cover set as a SET (because frac((N-e)x)=frac(-ex) only needs
  Nx in Z to be exact, but the cover predicate over ALL i is reflection-symmetric: if {frac(e_i x)}
  covers all 7 sectors then {1-frac(e_i x)} = {frac(-e_i x)} = {frac((N-e_i)x)} also covers all 7,
  exactly, for EVERY x, since r->6-r is a bijection of the 7 sectors -- see verify_cover_reflection).

  CONSEQUENCE we TEST: cover_set(N-E) == cover_set(-E) == cover_set(E) for EVERY x (set equality),
  hence ghat is UNCHANGED by e->N-e.  This makes the *complement-even symmetrization*
  E_sym = E (multiset) "u" (N-E) a candidate that does not change GOOD at all.  We test that and
  the finer question: which Z_2-graded PART of E actually moves ghat / SPEC.

THE SPLIT WE COMPUTE:
  (1) cluster reflection invariance: meas / ghat of cover_set(E) vs cover_set(N-E) vs cover_set(-E).
  (2) complement-even part E_even and complement-odd part E_odd of the cluster (w.r.t. a chosen
      center c=N/2): pair e with N-e; e is its own partner iff 2e=N. Compute SPEC for:
        - the full E,
        - E_evensym = the symmetric closure E u (N-E) (even under e<->N-e),
        - E restricted to fixed points / to one representative per pair,
      and see whether SPEC(E_evensym) == SPEC(E) (the conjecture: even part determines SPEC).
  (3) frequency Z_2 grading: split SPEC = SPEC_even_n + SPEC_odd_n by parity of n, and by n mod
      gcd-structure; report which part carries the (negative) deviation and which is ~mean-zero.

OUTPUT (the deliverable): raw numbers for the even/odd split of SPEC + the conjecture verdict.
All quantities EXACT (Fraction) where the arc algebra allows; spectral sums high-precision and
cross-checked against the exact real-space SPEC.
"""
import sys, itertools, cmath
from fractions import Fraction as F
from math import gcd, pi
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, danger_arcs, safe_set,
    cover_set, sector_of, fourier_num_of_arcs,
)

def chat(arcs, n):
    if n == 0:
        return complex(float(meas(arcs)), 0.0)
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)

def lattice_gcd(P):
    return reduce(gcd, [abs(x) for x in P if x != 0]) if any(P) else 0

# ------------------------------------------------------------------ exact SPEC (real space)
def spec_exact(P, E):
    gp = safe_set(P); cov = cover_set(E); covc = complement(cov)
    mGP = meas(gp); mC = meas(covc); base = mGP * mC
    if base == 0:
        return None
    inter = meas(intersect(gp, covc))
    SPEC = inter - base
    Rp = inter / base
    return dict(mGP=mGP, mC=mC, base=base, inter=inter, SPEC=SPEC, Rp=Rp, gp=gp, covc=covc)

# ------------------------------------------------------------------ (1) cover reflection test
def verify_cover_reflection(E, Ns):
    """For each candidate N, test cover_set(E) == cover_set([N-e]) == cover_set([-e]) as EXACT sets.
       Returns dict with set-equality booleans and measures."""
    E = sorted(set(int(e) for e in E))
    covE = cover_set(E)
    mE = meas(covE)
    negE = sorted(set(-e for e in E))          # may contain negatives; cover_set uses frac so ok? it uses 7*e grid
    out = {"E": E, "mE": mE, "rows": []}
    # cover_set assumes nonneg e for its breakpoint loop range(7*e+1); handle -e by using |e| equivalently:
    # frac(-e x) sector = 6 - sector(frac(e x)) (boundary aside) -> cover predicate identical.
    # We verify directly via a fine exact check on the union of ALL breakpoints of E and reflections.
    for N in Ns:
        NE = sorted(set(N - e for e in E))
        # build a common breakpoint set from E and NE (positive parts); use abs for the grid
        same_negE = _cover_equal(E, [-e for e in E])
        same_NE = _cover_equal(E, NE)
        mNE = _cover_meas(NE)
        out["rows"].append((N, NE, bool(same_NE), bool(same_negE), mNE, mNE == mE))
    return out

def _phases_cover(E, x):
    """does {frac(e x): e in E} hit all 7 sectors? exact."""
    hit = set()
    for e in E:
        hit.add(int((F(e) * x) % 1 * 7))
    return len(hit) == 7

def _all_breakpoints(E, extra=()):
    bps = set([F(0), F(1)])
    for e in E:
        ae = abs(int(e))
        if ae == 0:
            continue
        for m in range(7 * ae + 1):
            v = F(m, 7 * ae)
            if 0 < v < 1:
                bps.add(v)
    for v in extra:
        bps.add(v)
    return sorted(bps)

def _cover_equal(E1, E2):
    """exact set-equality of cover_set(E1) and cover_set(E2) by testing midpoints of the
       common refinement of all breakpoints."""
    bps = _all_breakpoints(list(E1))
    bps2 = _all_breakpoints(list(E2))
    allb = sorted(set(bps) | set(bps2))
    for a, b in zip(allb, allb[1:]):
        mid = (a + b) / 2
        if _phases_cover(E1, mid) != _phases_cover(E2, mid):
            return False
    return True

def _cover_meas(E):
    bps = _all_breakpoints(list(E))
    tot = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        if _phases_cover(E, mid):
            tot += (b - a)
    return tot

# ------------------------------------------------------------------ (2) complement-even/odd of cluster
def cluster_even_odd(E, N):
    """Pair e<->N-e (center N/2). Return (pairs, fixed, evensym multiset, one-rep set)."""
    E = sorted(set(int(e) for e in E))
    Eset = set(E)
    pairs = []; fixed = []; seen = set()
    for e in E:
        if e in seen:
            continue
        partner = N - e
        if partner == e:
            fixed.append(e); seen.add(e)
        elif partner in Eset:
            pairs.append((e, partner)); seen.add(e); seen.add(partner)
        else:
            # e has no partner inside E -> it is "complement-asymmetric"
            pass
    asym = [e for e in E if (N - e) not in Eset]
    # evensym closure: add the missing reflections
    evensym = sorted(set(E) | set(N - e for e in E))
    return dict(pairs=pairs, fixed=fixed, asym=asym, evensym=evensym, N=N)

# ------------------------------------------------------------------ (3) frequency Z_2 grading of SPEC
def spec_freq_grading(gp, covc, g, N=4000):
    """Split SPEC = sum 2 Re[chat(n)conj(ghat(n))] (n=1..N) by:
        - parity of n (even n vs odd n)
        - n in 7Z vs n not in 7Z
        - n in gZ (P-lattice) -- though chat=0 off gZ so this is automatic
       Returns the spectral SPEC plus each graded subsum, with per-n top terms."""
    tot = 0.0; even_n = 0.0; odd_n = 0.0; sev = 0.0; nonsev = 0.0
    terms = []
    for n in range(1, N + 1):
        cn = chat(gp, n); gn = chat(covc, n)
        term = 2.0 * (cn * gn.conjugate()).real
        tot += term
        if n % 2 == 0: even_n += term
        else: odd_n += term
        if n % 7 == 0: sev += term
        else: nonsev += term
        if abs(term) > 1e-5:
            terms.append((n, term))
    terms.sort(key=lambda t: -abs(t[1]))
    return dict(tot=tot, even_n=even_n, odd_n=odd_n, sev=sev, nonsev=nonsev, top=terms[:16])

# ================================================================== main
def main():
    print("#" * 96)
    print("# THREAD 2: COMPLEMENT-EVEN REFRAME of the Node-3 spectrum sum (kind-pasteur kpswf13)")
    print("#   e<->N-e (cluster T^op) vs x->-x (circle).  Does the complement-EVEN cluster part")
    print("#   determine SPEC?  Plus the frequency-n Z_2 grading of SPEC.")
    print("#" * 96)

    # The representative (P,E) bank (same as kpswf12 + the {2,3,4,5,6} killer + a deliberately
    # complement-asymmetric cluster to expose the odd part).
    bank = [
        ([1, 2, 3, 12, 13], list(range(8)),                         "k=8 consec (AP-like, R'~0.73)"),
        ([1, 2, 3, 4, 5],   list(range(8)),                         "k=8 consec small P (R'~0.63)"),
        ([1, 2, 3],         list(range(10)),                        "k=10 consec |P|=3 (R'~0.59)"),
        ([2, 3, 4, 5, 6],   list(range(8)),                         "P={2,3,4,5,6} killer"),
        ([5, 7, 11],        list(range(10)),                        "k=10 coprime P (R'~1.21)"),
        ([1, 2, 6],         [0, 4, 6, 8, 10, 12, 14, 15, 16, 17],   "wide d>=2 (R'~0.84)"),
        ([1, 2, 3, 12, 13], [0, 1, 3, 7, 12, 20, 30, 31],          "complement-ASYMMETRIC cluster"),
    ]

    # ---------------- (1) CLUSTER REFLECTION INVARIANCE of GOOD ----------------
    print("\n" + "=" * 96)
    print("(1) CLUSTER REFLECTION INVARIANCE:  cover_set(E) =?= cover_set(N-E) =?= cover_set(-E)")
    print("    (if TRUE for all x, then ghat is UNCHANGED by e->N-e: GOOD only sees the cluster")
    print("     up to speed-reflection.  Tested EXACTLY on the common breakpoint refinement.)")
    print("=" * 96)
    for P, E, lab in bank:
        Es = sorted(set(E))
        N_natural = max(Es) + min(Es)               # reflection center fixing the set's span
        N_max = max(Es)
        res = verify_cover_reflection(Es, [N_natural, N_max, 2 * max(Es)])
        print(f"\n  {lab}:  E={Es}   meas(cover_set(E))={float(res['mE']):.6f}")
        print(f"    cover_set(-E) == cover_set(E)? {_cover_equal(Es, [-e for e in Es])}"
              f"   (the x->-x / speed-negation invariance)")
        for N, NE, sNE, snE, mNE, meq in res["rows"]:
            print(f"    N={N:3d}: N-E={NE}  cover(N-E)==cover(E)? {sNE}  meas={float(mNE):.6f}  (meas eq? {meq})")

    # ---------------- (2) DOES THE COMPLEMENT-EVEN PART DETERMINE SPEC? ----------------
    print("\n" + "=" * 96)
    print("(2) COMPLEMENT-EVEN SYMMETRIZATION:  SPEC(E) vs SPEC(E u (N-E))  [the conjecture]")
    print("    center N = max(E)+min(E).  pairs e<->N-e; fixed: 2e=N; asym: N-e not in E.")
    print("    CONJECTURE (HYP-2867): R'/SPEC controlled by the complement-EVEN part; odd cancels.")
    print("=" * 96)
    for P, E, lab in bank:
        Es = sorted(set(E))
        N = max(Es) + min(Es)
        eo = cluster_even_odd(Es, N)
        base = spec_exact(P, Es)
        sym = spec_exact(P, eo["evensym"])
        print(f"\n  {lab}:  P={P}")
        print(f"    E          = {Es}")
        print(f"    N={N}: pairs={eo['pairs']}  fixed={eo['fixed']}  asym(odd)={eo['asym']}")
        print(f"    E_evensym  = {eo['evensym']}")
        if base:
            print(f"    SPEC(E)        = {float(base['SPEC']):+.6f}   R'(E)        = {float(base['Rp']):.6f}"
                  f"   base={float(base['base']):.6f}")
        if sym:
            print(f"    SPEC(E_sym)    = {float(sym['SPEC']):+.6f}   R'(E_sym)    = {float(sym['Rp']):.6f}"
                  f"   base={float(sym['base']):.6f}")
        if base and sym:
            dS = float(sym['SPEC'] - base['SPEC']); dR = float(sym['Rp'] - base['Rp'])
            asym_empty = (len(eo['asym']) == 0)
            print(f"    Delta SPEC (sym-E) = {dS:+.6f}   Delta R' = {dR:+.6f}"
                  f"   {'[E already even: must be 0]' if asym_empty else '[E has odd part]'}")

    # ---------------- (3) FREQUENCY Z_2 GRADING OF SPEC ----------------
    print("\n" + "=" * 96)
    print("(3) FREQUENCY GRADING:  SPEC = SPEC[even n] + SPEC[odd n] = SPEC[7|n] + SPEC[7 nmid n]")
    print("    Which graded part carries the (negative) deviation; which is ~mean-zero?")
    print("=" * 96)
    for P, E, lab in bank:
        Es = sorted(set(E)); g = lattice_gcd(P)
        b = spec_exact(P, Es)
        if not b:
            continue
        fr = spec_freq_grading(b['gp'], b['covc'], g, N=4000)
        print(f"\n  {lab}:  P={P} gcd={g}")
        print(f"    SPEC(spectral,|n|<=4000) = {fr['tot']:+.6f}   (exact {float(b['SPEC']):+.6f},"
              f" |diff|={abs(fr['tot']-float(b['SPEC'])):.1e})")
        print(f"    by parity:   even-n = {fr['even_n']:+.6f}   odd-n = {fr['odd_n']:+.6f}")
        print(f"    by 7-divis:  7|n    = {fr['sev']:+.6f}   7 nmid n = {fr['nonsev']:+.6f}")
        tops = ", ".join(f"n={n}:{t:+.4f}" for n, t in fr['top'][:8])
        print(f"    top terms:   {tops}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
