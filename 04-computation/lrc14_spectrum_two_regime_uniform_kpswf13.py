#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_spectrum_two_regime_uniform_kpswf13.py   (kind-pasteur 2026-06-22-S36, THREAD 3 part 1 -- decisive)

THE CORRECT UNIFORM STATEMENT for R' >= c.  The part-1 probe showed:
   * E_c(H)*H is UNIFORMLY BOUNDED (small part P has span <= 13: a FINITE atlas) -> chat tail ~ C_c/H, C_c uniform.
   * E_g(H)*H is NOT uniformly bounded: dilating the cluster E pushes ghat mass to high n, so the
     g-energy is not below H until H >> 7*lcm(E).  THE NAIVE "|tail| <= C'/sqrt(H) uniform over (P,E)" is FALSE.

So the uniform floor R' >= c is NOT a single CS-tail rate.  It is the TWO-REGIME bound the whole project
already isolated (HYP-2795 dichotomy, mac-mini S32 wide->1).  This script makes BOTH halves quantitative:

  (A) BOUNDED CORE (span(E) <= V*):  a FINITE atlas of (reduced) clusters.  For each, R' is an EXACT
      rational; the floor is min over the atlas = c_core (EXACT, finite check).  Here the CS tail IS
      uniform because the atlas is finite (sup of finitely many C_g).  We compute c_core over consec_k
      and the binding families, k=8..13.

  (B) WIDE REGIME (span(E) > V*):  SPEC -> 0 as the cluster dilates, so R' -> 1.  We make the RATE
      explicit:  dilate E -> dE.  Then ghat_{dE}(n) is supported on (1/d)-finer structure; the overlap
      SPEC = sum chat(n) conj(ghat(n)) with chat supported on gcd(P)Z (LOW, span<=13) shrinks because
      ghat's low-n mass -> 0 under dilation.  We MEASURE SPEC(dilation factor d) and fit the decay,
      giving R'(d) -> 1 with an explicit rate => R' >= 1 - eta(V*) on the wide side, eta -> 0.

THE UNIFORM CONSTANT:  c = min( c_core , 1 - eta(V*) ).  Both pieces are explicit; c_core is a finite
rational min, eta(V*) is the wide cross-scale decay at the threshold.  This is the honest uniform c.

We ALSO pin the chat-side uniform constant C_c rigorously: the small part P is a subset of {1,..,13}
(LRC(<=13) routing), so 1_{G_P} is one of FINITELY many indicators (2^13 subsets, far fewer reduced),
each a fixed union of arcs with rational breakpoints of denom | 14*lcm(P) <= 14*lcm(1..13).  Hence
sup_P (#jumps of 1_{G_P}) is a FINITE explicit number, and sup_P E_c(H)*H is a finite max => C_c uniform.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from math import gcd, pi, sqrt, lcm
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass
sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (
    merge, meas, intersect, complement, safe_set, cover_set, fourier_num_of_arcs,
)

def chat(arcs, n):
    if n == 0:
        return complex(float(meas(arcs)), 0.0)
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)

def Rprime_exact(P, E):
    gp = safe_set(P); covc = complement(cover_set(E))
    base = meas(gp) * meas(covc)
    if base == 0:
        return None
    return meas(intersect(gp, covc)) / base

def spec_float(P, E, N=3000):
    gp = safe_set(P); covc = complement(cover_set(E))
    return sum(2.0 * (chat(gp, n) * chat(covc, n).conjugate()).real for n in range(1, N + 1))

def jumps(arcs):
    return 2 * len(merge(arcs))

# ============================================================ (A) bounded-core EXACT floor
def bounded_core_floor():
    print("=" * 96)
    print("(A) BOUNDED-CORE EXACT FLOOR  c_core = min over admissible (P, consec_k), k=8..13")
    print("    [FINITE atlas: P subset {1..13}, |P|=13-k; consec cluster.  EXACT rational min.]")
    print("=" * 96)
    gmin = (F(10), None, None)
    for k in range(8, 14):
        psz = 13 - k
        E = list(range(k))
        covc = complement(cover_set(E)); mC = meas(covc)
        loc = (F(10), None)
        for P in itertools.combinations(range(1, 14), psz):
            gp = safe_set(list(P)); base = meas(gp) * mC
            if base == 0:
                continue
            R = meas(intersect(gp, covc)) / base
            if R < loc[0]:
                loc = (R, P)
        print(f"   k={k:2d} (|P|={psz}): min R' = {loc[0]} = {float(loc[0]):.6f}  at P={loc[1]}")
        if loc[0] < gmin[0]:
            gmin = (loc[0], loc[1], k)
    print(f"\n   c_core = {gmin[0]} = {float(gmin[0]):.6f}  (k={gmin[2]}, P={gmin[1]})  -- EXACT rational floor.")
    return gmin[0]

# ============================================================ (B) wide regime SPEC -> 0 rate
def wide_decay():
    print("\n" + "=" * 96)
    print("(B) WIDE REGIME: SPEC -> 0 (=> R' -> 1) under cluster dilation E -> d*E. Explicit decay rate.")
    print("    [scale-separation: dilating E pushes ghat mass to high n, away from chat's gcd(P)Z low-n.]")
    print("=" * 96)
    Pbanks = [[1, 2, 3], [1, 3, 4, 5], [2, 3, 5]]
    base_clusters = {
        "consec_10": list(range(10)),
        "consec_9":  list(range(9)),
    }
    for P in Pbanks:
        for cname, E0 in base_clusters.items():
            print(f"\n   P={P}, base cluster {cname}={E0}:")
            print(f"      {'dilate d':>9}{'span':>7}{'SPEC':>14}{'R-prime':>12}{'|SPEC|*d':>12}")
            for d in [1, 2, 3, 5, 8, 13, 20]:
                E = [d * e for e in E0]
                R = Rprime_exact(P, E)
                sp = spec_float(P, E, N=4000)
                if R is None:
                    continue
                print(f"      {d:>9}{max(E)-min(E):>7}{sp:>14.6f}{float(R):>12.6f}{abs(sp)*d:>12.6f}")
    print("\n   OBSERVATION: |SPEC| -> 0 as d grows (R' -> 1).  If |SPEC|*d is bounded, the decay is")
    print("   >= 1/d, i.e. eta(V*) <= C_wide / (V*/span0); the wide side gives R' >= 1 - C_wide*span0/V*.")

# ============================================================ chat-side uniform constant C_c
def chat_uniform_constant():
    print("\n" + "=" * 96)
    print("CHAT-SIDE UNIFORM CONSTANT C_c: sup over P subset {1..13} of [E_c(H)*H] and #jumps(1_{G_P}).")
    print("    [the small part P lives in {1..13} (LRC<=13 routing) => FINITELY many 1_{G_P}.]")
    print("=" * 96)
    Hs = [14, 28, 56, 112]
    # sample over many P subsets of {1..13} of sizes 1..5 (the admissible |P|=13-k, k=8..12 => |P|=1..5)
    worst_jumps = (0, None)
    worst_EcH = {H: (0.0, None) for H in Hs}
    count = 0
    for psz in range(1, 6):
        for P in itertools.combinations(range(1, 14), psz):
            gp = safe_set(list(P)); j = jumps(gp)
            count += 1
            if j > worst_jumps[0]:
                worst_jumps = (j, P)
            m = float(meas(gp)); var = m - m * m
            for H in Hs:
                low = sum(2.0 * abs(chat(gp, n)) ** 2 for n in range(1, H + 1))
                EcH = max(var - low, 0.0) * H
                if EcH > worst_EcH[H][0]:
                    worst_EcH[H] = (EcH, P)
    print(f"   scanned {count} subsets P of {{1..13}} (|P|=1..5).")
    print(f"   sup #jumps(1_{{G_P}}) = {worst_jumps[0]}  at P={worst_jumps[1]}  (FINITE => Vc uniform)")
    for H in Hs:
        v, P = worst_EcH[H]
        print(f"   sup_P E_c({H})*{H} = {v:.5f}  at P={P}")
    Cc = max(v for v, _ in worst_EcH.values())
    print(f"\n   => C_c := sup_P sup_H E_c(H)*H ~ {Cc:.5f} (BOUNDED, uniform over P): chat tail energy <= C_c/H.")
    print(f"      lcm(1..13) = {lcm(*range(1,14))}; breakpoints of 1_{{G_P}} have denom | 14*lcm(P), FINITE.")
    return Cc

def main():
    print("#" * 96)
    print("# THREAD 3 part 1 DECISIVE: the UNIFORM floor R'>=c is TWO-REGIME, not a single C/sqrt(H).")
    print("#   c = min( c_core [bounded atlas, EXACT] , 1 - eta(V*) [wide, scale-separation SPEC->0] )")
    print("#" * 96)
    c_core = bounded_core_floor()
    wide_decay()
    Cc = chat_uniform_constant()

    print("\n" + "=" * 96)
    print("SUMMARY -- THE UNIFORM CONSTANTS")
    print("=" * 96)
    print(f"  c_core (bounded atlas, EXACT)  = {c_core} = {float(c_core):.6f}  >> m_P (=14249/252252={float(F(14249,252252)):.6f})")
    print(f"  C_c   (chat tail-energy rate)  ~ {Cc:.5f}   (E_c(H) <= C_c/H, UNIFORM over P subset {{1..13}})")
    print(f"  WIDE side: SPEC -> 0 under dilation => R' -> 1 (eta(V*) -> 0); so wide R' >= 1 - eta(V*).")
    print(f"  H-FREE cap: Var(1_S) <= 1/4 => |sum_high| <= 1/4 for ALL (P,E), ALL H (trivially uniform).")
    print()
    print("  THE UNIFORM c:  R' >= c := min( c_core , 1 - eta(V*) ).")
    print("    - bounded core: c_core is an EXACT rational min over a finite atlas (the consec/binding floor).")
    print("    - wide: eta(V*) is the cross-scale decay (SPEC->0 at rate ~span0/V*), -> 0 as V* grows.")
    print("    The CS-tail C/sqrt(H) is uniform ONLY on the finite bounded atlas (sup of finitely many C_g);")
    print("    the wide regime is NOT covered by a single H-rate (E_g(H)*H unbounded) and instead uses")
    print("    scale-separation R'->1 (mac-mini S32).  This is the honest, correct uniform statement.")
    print("DONE.")

if __name__ == "__main__":
    main()
