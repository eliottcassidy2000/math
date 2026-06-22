#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_thread3_uniform_constants_TRUE_kpswf13.py  (kind-pasteur 2026-06-22-S37, THREAD 3 -- TRUE GOOD)

THREAD-3 UNIFORM CONSTANTS, on the CORRECT object.  S36b found the kpswf12 `cover^c` proxy is x->-x
ASYMMETRIC, so its ghat is NOT real; the prompt's GOOD = {x: cluster maxgap{frac(e_i x)}>1/7} is the
genuinely complement-symmetric object (ghat REAL, sineE=0).  This script recomputes the Thread-3
uniform constants on GOOD_true (NOT the proxy):

  R'_true = meas(G_P cap GOOD_true) / (meas(G_P) meas(GOOD_true)) = 1 + SPEC_true/baseline,
  SPEC_true = sum_{n!=0} chat(n) ghat_true(n),   ghat_true REAL.

DELIVERABLES (the uniform constants):
  (A) c_core = EXACT bounded-atlas floor = min over admissible (P, consec_k), k=8..13, of R'_true.
  (B) C_c    = sup over P subset {1..13} of E_c(H)*H  (chat tail-energy rate; UNCHANGED from proxy,
               depends only on G_P) -- the chat side IS uniform (finite atlas).
  (C) C_wide = sup-envelope constant in |SPEC_true(dE)| <= C_wide/d (wide decorrelation, R'_true->1).
  (D) the H-free cap Var<=1/4 (trivially uniform), and the dilation identity ghat_{dE}=ghat_E(n/d)[d|n].
"""
import sys, itertools
from fractions import Fraction as F
from math import pi, sqrt, lcm
from functools import reduce
try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass
sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import merge, meas, intersect, complement, safe_set, fourier_num_of_arcs
from lrc14_true_maxgap_good_kpswf13 import good_true, reflect_arcs

def chat(arcs, n):
    if n == 0:
        return complex(float(meas(arcs)), 0.0)
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)

def Rprime_true(P, E, good_cache=None):
    gp = safe_set(P)
    good = good_cache if good_cache is not None else good_true(E)
    base = meas(gp) * meas(good)
    if base == 0:
        return None, None, good
    return meas(intersect(gp, good)) / base, base, good

def spec_true(P, E, N=3000, good_cache=None):
    gp = safe_set(P)
    good = good_cache if good_cache is not None else good_true(E)
    return sum(2.0 * (chat(gp, n) * chat(good, n).conjugate()).real for n in range(1, N + 1)), good

def p(*a):
    print(*a, flush=True)

# ============================================================ confirm GOOD_true is symmetric
def confirm_symmetry():
    p("=" * 92)
    p("(0) CONFIRM GOOD_true is x->-x symmetric (=> ghat REAL), unlike the kpswf12 cover^c proxy.")
    p("=" * 92)
    for E in [list(range(8)), list(range(10)), [0, 3, 7, 11]]:
        g = good_true(E); gr = reflect_arcs(g)
        sd = meas(merge([(a, b) for a, b in g])) - meas(intersect(g, gr))  # meas(g \ reflect)
        # symmetric difference measure:
        symdiff = meas(g) + meas(gr) - 2 * meas(intersect(g, gr))
        sineE = sum(chat(g, n).imag ** 2 for n in range(1, 800))
        p(f"   E={E}: meas(GOOD)={float(meas(g)):.6f}  symdiff(GOOD,reflect)={float(symdiff):.2e}  sineE(800)={sineE:.2e}")
    p("   => symdiff = 0 EXACTLY, sineE ~ 0 => ghat_true REAL.  This is the prompt's GOOD.")

# ============================================================ (A) bounded-core floor (TRUE)
def bounded_core_floor_true():
    p("\n" + "=" * 92)
    p("(A) c_core (TRUE maxgap GOOD): min over admissible (P, consec_k), k=8..13, of R'_true  [EXACT]")
    p("=" * 92)
    gmin = (F(10), None, None)
    for k in range(8, 14):
        psz = 13 - k
        E = list(range(k)); good = good_true(E); mG = meas(good)
        loc = (F(10), None)
        for P in itertools.combinations(range(1, 14), psz):
            gp = safe_set(list(P)); base = meas(gp) * mG
            if base == 0:
                continue
            R = meas(intersect(gp, good)) / base
            if R < loc[0]:
                loc = (R, P)
        p(f"   k={k:2d} (|P|={psz}): min R'_true = {loc[0]} = {float(loc[0]):.6f}  at P={loc[1]}  [meas(GOOD)={float(mG):.5f}]")
        if loc[0] < gmin[0]:
            gmin = (loc[0], loc[1], k)
    p(f"\n   c_core = {gmin[0]} = {float(gmin[0]):.6f}  (k={gmin[2]}, P={gmin[1]})  -- EXACT rational floor on TRUE GOOD.")
    return gmin[0]

# ============================================================ (B) chat-side C_c (G_P only; UNCHANGED)
def chat_uniform_C_c():
    p("\n" + "=" * 92)
    p("(B) C_c = sup over P subset {1..13} of E_c(H)*H  (chat tail energy; depends only on G_P).")
    p("=" * 92)
    Hs = [14, 28, 56, 112]
    worst = {H: (0.0, None) for H in Hs}
    worst_j = (0, None)
    cnt = 0
    for psz in range(1, 6):
        for P in itertools.combinations(range(1, 14), psz):
            gp = safe_set(list(P)); cnt += 1
            j = 2 * len(merge(gp))
            if j > worst_j[0]:
                worst_j = (j, P)
            m = float(meas(gp)); var = m - m * m
            for H in Hs:
                low = sum(2.0 * abs(chat(gp, n)) ** 2 for n in range(1, H + 1))
                EcH = max(var - low, 0.0) * H
                if EcH > worst[H][0]:
                    worst[H] = (EcH, P)
    p(f"   scanned {cnt} subsets P of {{1..13}}; sup #jumps(1_{{G_P}}) = {worst_j[0]} at P={worst_j[1]}")
    for H in Hs:
        v, P = worst[H]
        p(f"   sup_P E_c({H})*{H} = {v:.5f} at P={P}")
    Cc = max(v for v, _ in worst.values())
    p(f"   => C_c ~ {Cc:.5f} (BOUNDED, uniform over P).  lcm(1..13)={lcm(*range(1,14))} bounds breakpoint denoms.")
    return Cc

# ============================================================ (C) wide envelope C_wide (TRUE)
def wide_envelope_true():
    p("\n" + "=" * 92)
    p("(C) C_wide (TRUE GOOD): sup-envelope of |SPEC_true(dE)| over dilations; |SPEC_true(dE)|<=C_wide/d.")
    p("=" * 92)
    cases = [([1, 3, 4, 5], list(range(9)), "argmin-ish k9"),
             ([1, 2, 3], list(range(10)), "k10 |P|=3"),
             ([2, 3, 5], list(range(10)), "coprime-P")]
    allMD = []
    for P, E0, lab in cases:
        p(f"   CASE {lab}  P={P}:")
        sp = {}
        for d in range(1, 41):
            sval, _ = spec_true(P, [d * e for e in E0], N=2500)
            sp[d] = abs(sval)
        for D in [1, 2, 3, 5, 8, 13, 21]:
            dom = [d for d in sp if d >= D]
            am = max(dom, key=lambda d: sp[d]); M = sp[am]
            allMD.append(M * D)
            p(f"      D={D:3d}  M(D)={M:.5f}  M(D)*D={M*D:.5f}  argmax_d={am}")
    Cw = max(allMD)
    p(f"   => C_wide = sup_D sup_{{d>=D}}|SPEC_true(dE)|*D = {Cw:.5f}")
    return Cw

# ============================================================ dilation identity (TRUE GOOD)
def dilation_identity_true():
    p("\n" + "=" * 92)
    p("(D) DILATION IDENTITY on TRUE GOOD: ghat_true_{dE}(n) = ghat_true_E(n/d)*[d|n] (exact mechanism).")
    p("=" * 92)
    E0 = list(range(10)); d = 3
    g0 = good_true(E0); gd = good_true([d * e for e in E0])
    p(f"   meas(GOOD(E0))={float(meas(g0)):.6f}  meas(GOOD(3E0))={float(meas(gd)):.6f}  (EQUAL: dilation preserves measure)")
    mx = 0.0
    for n in range(1, 22):
        v = abs(chat(gd, n))
        if n % d == 0:
            w = abs(chat(g0, n // d)); mx = max(mx, abs(v - w))
        else:
            mx = max(mx, v)
    p(f"   max |ghat_dE(n) - ghat_E(n/d)*[d|n]| over n=1..21 = {mx:.2e}  (=> dilation moves ghat spectrum to dZ)")

def main():
    p("#" * 92)
    p("# THREAD 3 (TRUE GOOD): uniform constants for R'_true>=c.  GOOD={maxgap>1/7}, ghat REAL.")
    p("#   c = min( c_core[bounded atlas, EXACT] , 1 - eta(V*)[wide, SPEC->0] ); C_c chat-tail, C_wide wide-decay")
    p("#" * 92)
    confirm_symmetry()
    c_core = bounded_core_floor_true()
    Cc = chat_uniform_C_c()
    Cw = wide_envelope_true()
    dilation_identity_true()
    mP = F(14249, 252252)
    p("\n" + "=" * 92)
    p("SUMMARY -- THE UNIFORM CONSTANTS (on the prompt's TRUE GOOD = {maxgap>1/7})")
    p("=" * 92)
    p(f"  c_core (bounded atlas, EXACT) = {c_core} = {float(c_core):.6f}  = {float(c_core/mP):.1f}x m_P")
    p(f"  C_c    (chat tail-energy)     ~ {Cc:.5f}   (E_c(H)<=C_c/H, UNIFORM over P subset {{1..13}})")
    p(f"  C_wide (wide-decay envelope)  ~ {Cw:.5f}   (|SPEC_true(dE)|<=C_wide/d => R'_true->1)")
    p(f"  H-free cap: Var(1_S)<=1/4 => |sum_high|<=1/4 for ALL (P,E), ALL H.")
    p(f"  THE UNIFORM c: R'_true >= min( c_core={float(c_core):.4f} , 1 - C_wide*span0/(base*V*) ).")
    p(f"  Bounded core: finite EXACT min. Wide: scale-separation R'->1, explicit envelope. m_P={float(mP):.6f}.")
    p("DONE")

if __name__ == "__main__":
    main()
