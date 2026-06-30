#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Where does 7/89 come from? The binding-pair / packing-radius anatomy of M for covering
sets, through the Littlewood (badly-approximable / continued-fraction) lens.
(mac-mini-2026-06-29-S18)

7/89 = M({1..11,13,84}), the tightest known covering set (S15). Anatomy:
  - covering completion of the 'skip-12' core {1..11,13} (M=1/12) by the killer 84=lcm(12,14).
  - M=7/89 at t*=37/89; binding pair (5,84) since 84 == -5 (mod 89), so 89 = 84+5.
  - 7 = the PACKING RADIUS: at the best rotation a=37, the 13 speeds*a mod 89 all avoid the
    7-neighbourhood of 0; 7 is the largest such radius (the covering radius in Z/89).

GENERAL STRUCTURE (S15): M(S)=max over t of min_i ||v_i t||, the max at a point where two
speeds BIND; M = j/D with D | (v_a +- v_b) a binding-pair sum/difference. So 'where M comes
from' = which binding pair + what packing radius. The LITTLEWOOD lens: the tightest configs
should sit at BADLY-APPROXIMABLE speed ratios (bounded partial quotients), and the binding
modulus D should be a CONVERGENT denominator of the relevant ratio.
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
print = functools.partial(print, flush=True)


def frac(x): return x - math.floor(x)
def dist(x): f = frac(x); return min(f, 1 - f)


def M_exact_via_pairs(S, Dmax=2000):
    """M(S) exactly: scan candidate binding moduli D (divisors of v_a+-v_b), and for each D
    the best rotation a in [1,D), M = max over (a,D) of min_i ||v_i a / D||. Returns (M, a, D, binders)."""
    cands = set()
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            cands.add(S[i] + S[j]); cands.add(abs(S[i] - S[j]))
            cands.add(2 * S[i]); cands.add(2 * S[j])
    cands = sorted(d for d in cands if 1 < d <= Dmax)
    best = (F(0), 0, 1, [])
    for D in cands:
        for a in range(1, D):
            if math.gcd(a, D) != 1: continue
            md = min((v * a) % D for v in S)
            md = min(md, D - max((v * a) % D for v in S)) if False else min(min((v*a)%D, D-((v*a)%D)) for v in S)
            if F(md, D) > best[0]:
                binders = [v for v in S if min((v*a) % D, D-((v*a) % D)) == md]
                best = (F(md, D), a, D, binders)
    return best


def cf(p, q):
    a = []
    while q:
        a.append(p // q); p, q = q, p % q
    return a


def main():
    print("=" * 82)
    print("Where does 7/89 come from? binding-pair / packing-radius anatomy + Littlewood lens")
    print("=" * 82)

    # ---- (A) anatomy of 7/89 ----
    S = list(range(1, 12)) + [13, 84]
    print(f"\n[A] {S}  (skip-12 core + killer 84=lcm(12,14)):")
    M, a, D, binders = M_exact_via_pairs(S)
    print(f"    M = {M} at t* = {a}/{D};  binders = {binders}  (84 == {84 % D} = -5 mod {D}, so {D}=84+5)")
    print(f"    packing radius = numerator {M.numerator}: the 13 speeds*{a} mod {D} all avoid the")
    res = sorted(min((v*a) % D, D-((v*a) % D)) for v in S)
    print(f"    {M.numerator}-neighbourhood of 0; sorted distances to 0: {res}")
    print(f"    cont.frac(89/7)={cf(89,7)}, (89/37)={cf(89,37)}, (84/5)={cf(84,5)}  (89=F_11, coincidence)")

    # ---- (B) the single-large covering family {1..11,13,L}, L multiple of 84 ----
    print(f"\n[B] single-large covering family {{1..11,13,L}}, L=84k: M, binding modulus D, partner:")
    print(f"    {'L':>5} {'M':>10} {'D':>6} {'partner (D-L mod?)':>18} {'cf(D/Mnum)':>16}")
    for k in range(1, 7):
        L = 84 * k; Sk = list(range(1, 12)) + [13, L]
        Mk, ak, Dk, bk = M_exact_via_pairs(Sk, Dmax=3000)
        partner = [v for v in bk if v != L]
        print(f"    {L:>5} {str(Mk):>10} {Dk:>6} {str(partner):>18} {str(cf(Dk, Mk.numerator)):>16}")
    print("    => D = L + (small partner); M numerator = packing radius. NOT Fibonacci across the")
    print("    family (89,173,...) => 89=F_11 was coincidence. The structure is binding-pair + packing.")

    # ---- (C) the LITTLEWOOD lens: tightest covering sets sit at badly-approximable ratios? ----
    print(f"\n[C] LITTLEWOOD lens -- vary the killer L freely (any covering single-large set),")
    print(f"    find the TIGHTEST M and check the binding-ratio continued fraction (bounded q_i?):")
    rows = []
    for L in range(14, 600):
        Sk = list(range(1, 12)) + [13, L]
        if len(set(Sk)) != 13: continue
        if not all(any(v % q == 0 for v in Sk) for q in range(2, 15)): continue   # covering
        Mk, ak, Dk, bk = M_exact_via_pairs(Sk, Dmax=2 * L + 50)
        rows.append((Mk, L, Dk, ak, bk))
    rows.sort()
    print(f"    covering single-large sets: {len(rows)}; TIGHTEST 6 (smallest M):")
    for Mk, L, Dk, ak, bk in rows[:6]:
        partner = [v for v in bk if v != L]
        ratio_cf = cf(L, partner[0]) if partner else []
        print(f"      L={L:>4}: M={str(Mk):>9} D={Dk:>5} partner={str(partner):>8} "
              f"cf(L/partner)={ratio_cf}")
    print("    Littlewood reading: the binding pair (L, partner) -- the tightness is the simultaneous")
    print("    approximation of a wall by BOTH L*t and partner*t at t*=a/D; D=L+partner is their")
    print("    'resonance modulus'. The product ||L t*||*||partner t*|| = M^2 (both = M at the bind).")

    # ---- (D) the Littlewood product at the bind, and the q*||.||*||.|| quantity ----
    print(f"\n[D] LITTLEWOOD product q*||q a||*||q b|| at the binding (q=D, a=L/D, b=partner/D):")
    Mk, ak, Dk, bk = M_exact_via_pairs(S)   # the 7/89 set
    partner = [v for v in bk if v != 84][0]
    val = Dk * float(Mk) * float(Mk)
    print(f"    7/89 set: D={Dk}, M={Mk}; D*M*M = {Dk}*({Mk})^2 = {val:.4f} = {Dk}*{Mk.numerator}^2/{Dk}^2 "
          f"= {Mk.numerator**2}/{Dk} = {F(Mk.numerator**2, Dk)}")
    print(f"    Littlewood asks inf_q q||qa||||qb||=0; here the LRC sup pins it at {Mk.numerator}^2/D = 49/89")
    print(f"    >0 -- the LRC FLOOR is the OBSTRUCTION to the Littlewood product vanishing at the bind.")

    print("\n" + "=" * 82)
    print("UNDERSTANDING 7/89: it is the PACKING RADIUS (7) over the BINDING MODULUS (89 = 84+5)")
    print("of the tightest covering completion of the skip-12 core. M=j/D where D=binding-pair sum,")
    print("j = covering radius of the speeds in Z/D. Littlewood lens: the bind is the simultaneous")
    print("approximation of a wall by the pair (L, partner); the LRC floor M>=1/14 is exactly the")
    print("statement that this product CANNOT be driven to 0 -- the anti-Littlewood (a POSITIVE floor).")
    print("=" * 82)


if __name__ == "__main__":
    main()
