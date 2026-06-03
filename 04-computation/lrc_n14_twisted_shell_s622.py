#!/usr/bin/env python3
"""
S622 — the TWISTED-SHELL DODGE for C'(14).
A time t=a/m (gcd(a,m)=1) gives  M(S) >= min_i ||a s_i / m|| = (1/m) min_i dist(a s_i mod m, 0).
This BEATS 1/14 iff every residue a s_i mod m lies in the central band-complement
    r in ( m/14 , 13m/14 ),   i.e.  min(r, m-r) > m/14.
 - m <= 13: band excludes only 0 -> works iff no s_i ≡ 0 (mod m) (twist a irrelevant): the 1-clock.
 - 15 <= m <= 27: the multiplier a is ESSENTIAL — it rotates all runners off the danger band.
   These live on the shells (21=3·7, 25=5^2, 27=3^3).
CLAIM to test: every multiple-of-14 config has a twisted-shell escape with SMALL m
(=> a constructive proof of C'(14)).  We record the minimal escaping m over hard configs.
"""
from math import gcd
import itertools, random

def escapes(S, m):
    """is there a multiplier a coprime to m with all a*s_i mod m in the safe band (dist> m/14)?
       returns a witnessing a, or None.  (m/14 compared via 14*min(r,m-r) > m, strict)."""
    units = [a for a in range(1, m) if gcd(a, m) == 1]
    for a in units:
        ok = True
        for s in S:
            r = (a*s) % m
            if 14*min(r, m-r) <= m:    # in danger band (dist <= 1/14)
                ok = False; break
        if ok:
            return a
    return None

def min_escaping_m(S, mmax=60):
    for m in range(2, mmax+1):
        a = escapes(S, m)
        if a is not None:
            twist = (m > 13)            # for m<=13 any unit works; twist essential only for m>=15
            return m, a, twist
    return None, None, None

def hard_configs():
    """multiple-of-14 configs engineered to block many small clocks (contain multiples of
       2,3,...,13) so the escape is forced onto a shell."""
    cfgs = []
    # 1) AP-with-swap family (near tight)
    base = list(range(1, 14))
    for w in (1, 2, 3):
        for i in range(13):
            S = base[:]; S[i] = 14*w
            if len(set(S)) == 13: cfgs.append(tuple(sorted(S)))
    # 2) "block every small clock": include 14 and one multiple of each of 2..13 (highly composite picks)
    random.seed(1)
    for _ in range(4000):
        S = {14}
        # force multiples of several small primes/clocks
        for m in random.sample([3,4,5,6,7,8,9,10,11,12,13], random.randint(6,11)):
            S.add(m*random.randint(1,3))
        while len(S) < 13:
            S.add(random.randint(1,40))
        S = tuple(sorted(S))
        if len(S) == 13:
            g = 0
            for x in S: g = gcd(g, x)
            if g == 1: cfgs.append(S)
    return list(dict.fromkeys(cfgs))

if __name__ == "__main__":
    cfgs = hard_configs()
    print(f"Testing twisted-shell escape on {len(cfgs)} multiple-of-14 configs.\n")
    fails = []; mdist = {}; twist_needed = 0; maxm = 0; ex_twist = []
    for S in cfgs:
        m, a, tw = min_escaping_m(S)
        if m is None:
            fails.append(S); continue
        mdist[m] = mdist.get(m, 0) + 1
        maxm = max(maxm, m)
        if tw:
            twist_needed += 1
            if len(ex_twist) < 8: ex_twist.append((S, m, a))
    print(f"configs with NO escape up to m=60: {len(fails)}  {fails[:3]}")
    print(f"max minimal-escaping-m over all configs: {maxm}")
    print(f"configs needing a genuine twist (escape m>=15): {twist_needed}/{len(cfgs)}")
    print("distribution of minimal escaping m:", dict(sorted(mdist.items())))
    print("\nexamples needing a twisted-shell escape (S, m, multiplier a):")
    for S, m, a in ex_twist:
        print(f"  m={m}  a={a}  S={S}")
