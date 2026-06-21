#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_doubled_prime_generalization_macmini.py  (mac-mini 2026-06-21, THREAD D NEW LEAD (i))

DOES THE LRC(14) PROOF METHOD EXTEND TO LRC(n)=LRC(2p) FOR OTHER PRIMES p?

The LRC(14) S3-sector route reduces to: measSp(E) = P(all p sectors hit) <= cap_k for
primitive k-sets, where sector(y)=floor(p*frac(y)) for the apex prime p=7 (n=2p=14).

The method's pillars (all p-keyed):
  - covering reduction (THM-527/530/531): tournament-side, structural, p-agnostic
  - sector cover = coupon-cover over p sectors  (p-keyed)
  - Delsarte LP over the binary Krawtchouk K_j(.;p-1)  (p-keyed)
  - L7 torus-discrepancy Koksma bound D <= 2(p-1)/p / something  (p-keyed; 14/p was p=7)

QUESTION: for p=5 (n=10), p=11 (n=22), p=13 (n=26):
  (a) is consec [1..k] STILL the measSp-maximizer over primitive k-sets?
  (b) does the apex D=0 law generalize:  D_{a,b}=0 <=> p|ab ?
  (c) what is the analog Koksma constant (was 14/p = 2(p-1)/p * p/... )?

This is the first systematic test of whether the whole machine is p-portable.
We use SMALL p (5) where exhaustive checks are cheap, and confirm the same phenomena.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd

def make_sector(p):
    def sector(yf): return int(p * yf)
    return sector

def breakpoints(E, p):
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for t in range(0, p * e):
            bp.add(Fr(t, p * e))
    return sorted(bp)

def measSp(E, p):
    """P_v[ all p sectors hit by {e v mod 1 : e in E} ]."""
    sector = make_sector(p)
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E, p); tot = Fr(0)
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        if len(set(sector((e * mid) % 1) for e in E)) == p:
            tot += (b - a)
    return tot

def mu_2far(a, b, p):
    """joint sector occupancy of (sector(a v), sector(b v)) -- the apex D law object."""
    sector = make_sector(p)
    bp = {Fr(0), Fr(1)}
    for f in (a, b):
        for t in range(0, p * f): bp.add(Fr(t, p * f))
    vs = sorted(bp); cell = {}
    for x, y in zip(vs, vs[1:]):
        mid = (x + y) / 2
        key = (sector((a*mid) % 1), sector((b*mid) % 1))
        cell[key] = cell.get(key, Fr(0)) + (y - x)
    return cell

def D_apex(a, b, p):
    cell = mu_2far(a, b, p)
    return sum(abs(cell.get((i, j), Fr(0)) - Fr(1, p*p)) for i in range(p) for j in range(p))

def main():
    print("#"*80)
    print("# DOUBLED-PRIME GENERALIZATION of the LRC(14) method")
    print("#"*80)

    for p in [3, 5, 7, 11]:
        n = 2 * p
        print(f"\n{'='*70}\np = {p}  (n = {n})\n{'='*70}")

        # (a) consec-max test: smallest interesting k for the cover is k>=p (need p sectors)
        #     test k = p, p+1 : is consec [1..k] the maximizer over primitive k-sets?
        for k in [p, p+1]:
            cons = list(range(1, k+1)); pc = measSp(cons, p)
            # exhaustive over primitive k-sets, min=1, max<=k+4 (cheap window)
            beat = []; total = 0
            Nmax = k + 4
            for combo in itertools.combinations(range(2, Nmax+1), k-1):
                E = (1,) + combo
                g = 0
                for e in E: g = gcd(g, e)
                if g != 1: continue
                total += 1
                v = measSp(E, p)
                if v > pc + Fr(1, 10**12):
                    beat.append((E, float(v)))
            flag = "consec is MAX" if not beat else f"BEATEN by {len(beat)}: {sorted(beat,key=lambda t:-t[1])[:2]}"
            print(f"  k={k}: measS{p}(consec)={float(pc):.5f}  searched {total} prim sets (max<={Nmax})  -> {flag}")

        # (b) apex D=0 law: D_{a,b}=0 <=> p|ab ?
        viol = 0; checked = 0
        for q in range(1, 22):
            for r in range(1, 22):
                if gcd(q, r) != 1: continue
                checked += 1
                d = D_apex(q, r, p)
                law = (q % p == 0) or (r % p == 0)
                if (d == 0) != law:
                    viol += 1
        print(f"  apex law D=0 <=> {p}|ab : {viol} violations over {checked} coprime pairs (a,b<=21)")

        # (c) Koksma constant: sup D*min(a,b) over non-apex coprime pairs; compare to 2(p-1)
        sup_Dp = Fr(0); arg = None
        for q in range(1, 30):
            if q % p == 0: continue
            for r in range(1, 30):
                if r % p == 0: continue
                if gcd(q, r) != 1: continue
                d = D_apex(q, r, p)
                m = min(q, r)
                if d * m > sup_Dp: sup_Dp = d * m; arg = (q, r)
        # the proved bound for p=7 was D <= 14/p_low = 2(p-1)... actually 14/p where 14=2p.
        # general Koksma: D <= 2*p*(p-1)/p / p?  let's just report sup D*min and the ratio to (n/min)
        print(f"  Koksma: sup D*min(a,b) over non-apex = {sup_Dp} = {float(sup_Dp):.4f} at {arg}; "
              f"compare n=2p={n}, p^2*2/p={2*p}")

    print("\n=== INTERPRETATION ===")
    print("If consec-max holds, apex law holds, and Koksma constant scales for all p,")
    print("the whole S3-sector route is p-PORTABLE: LRC(2p) reduces the same way.")
    print("\nDONE.")

if __name__ == "__main__":
    main()
