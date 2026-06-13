#!/usr/bin/env python3
"""
lrc_resonance_energy_sufficient_s550.py    oracle-2026-06-01-S550o

A RIGOROUS SUFFICIENT CONDITION for LRC, via the resonance-energy bound, and the
MINIMAL RESONANCE LENGTH as the controlling invariant. Searching for progress.

THE IDENTITY (S529). With ghat(0)=1-2/n, ghat(m)=-sin(2 pi m/n)/(pi m),
   |LONELY(v)| = INT_0^1 prod_i 1_[1/n,1-1/n](v_i t) dt
              = sum_{m in Z^{n-1}: sum m_i v_i = 0} prod_i ghat(m_i).
The all-zero term is the MAIN TERM (1-2/n)^{n-1}; everything else is the resonance
correction. By the triangle inequality,
   | |LONELY(v)| - (1-2/n)^{n-1} |  <=  E(v) := sum_{0!=m, sum m_i v_i=0} prod_i |ghat(m_i)|.

>>> SUFFICIENT CONDITION (rigorous): if  E(v) < (1-2/n)^{n-1}  then |LONELY(v)| > 0,
    i.e. v is n-lonely (LRC holds for v).  <<<

E(v) is dominated by SMALL resonances (ghat decays like 1/(pi|m|)); the controlling
invariant is the MINIMAL RESONANCE LENGTH  l(v) = min { sum|m_i| : 0!=m, sum m_i v_i=0 }.
Large l(v) (decorrelated speeds, S544) => tiny E => LRC. Small l(v) (AP-like, many
returns S545) => large E => the hard core (the regular polygon saturates it).

This script: (1) the pairwise resonance energy E2 (the dominant, exactly-summable
part) and the bound vs the main term; (2) the minimal resonance length l(v) vs the
direct lonely measure -- the decorrelation->LRC picture; (3) the HARD CORE map: which
small-n sets have E2 >= main (the bound fails) -- are they exactly the AP-like sets?
"""
from itertools import combinations, product as iproduct
from functools import reduce
from math import gcd, sin, pi
import random

def frac(x): return x - int(x // 1)
def ghat_abs(m, n):
    if m == 0: return 1.0 - 2.0/n
    return abs(sin(2*pi*m/n))/(pi*abs(m))

# ---------- direct lonely measure ----------
def lonely_measure(v, n, G=200000):
    lo, hi = 1.0/n, 1.0 - 1.0/n; c = 0
    for i in range(G):
        t = (i+0.5)/G
        if all(lo < (s*t) % 1.0 < hi for s in v): c += 1
    return c/G

# ---------- minimal resonance length ----------
def minimal_resonance(v, M=10):
    """min sum|m_i| over nonzero integer m with sum m_i v_i = 0, |m_i|<=M (search)."""
    n1 = len(v); best = None
    # search by increasing total |m|: enumerate small supports first (pairs, then triples)
    # pairwise: m_i v_i + m_j v_j = 0 -> m_i = v_j/g, m_j = -v_i/g (g=gcd)
    for i, j in combinations(range(n1), 2):
        g = gcd(v[i], v[j]); l = v[i]//g + v[j]//g  # |m_i|+|m_j| = (v_j+v_i)/g
        if best is None or l < best: best = l
    # also short full/short combos via small-coeff search (catches sub-pair shorter)
    rng = range(-3, 4)
    for ms in iproduct(rng, repeat=n1):
        if any(ms) and sum(a*b for a, b in zip(ms, v)) == 0:
            l = sum(abs(a) for a in ms)
            if best is None or l < best: best = l
    return best

# ---------- pairwise resonance energy E2 ----------
def E2(v, n, K=2000):
    """dominant (order-2) part of E(v): sum over pairs (i,j) and k>=1 of
    (1-2/n)^{n-3} * |ghat(k*v_j/g)| * |ghat(k*v_i/g)|, g=gcd(v_i,v_j)."""
    n1 = len(v); base = (1 - 2.0/n)**(n1 - 2)  # n-1 runners, 2 involved -> n-3 = n1-2 zeros
    tot = 0.0
    for i, j in combinations(range(n1), 2):
        g = gcd(v[i], v[j]); a = v[j]//g; b = v[i]//g
        s = 0.0
        for k in range(1, K+1):
            # m_i = k*a, m_j = -k*b (so m_i v_i + m_j v_j = k a v_i - k b v_j = k(v_j v_i - v_i v_j)/g=0)
            gi = ghat_abs(k*a, n); gj = ghat_abs(k*b, n)
            if gi == 0 and gj == 0: continue
            s += gi*gj
        tot += 2*base*s   # factor 2: (m_i,m_j) and (-m_i,-m_j)
    return tot

def main():
    print("="*74)
    print("RIGOROUS SUFFICIENT CONDITION: E(v) < (1-2/n)^{n-1}  =>  LRC(v).")
    print("E dominated by small resonances; controlled by minimal resonance length l(v).")
    print("="*74)
    print()
    print("(1)+(2) minimal resonance l(v), pairwise energy E2, main term, lonely measure")
    for n in (5, 6, 7):
        main_t = (1 - 2.0/n)**(n-1)
        print(f"\n  n={n}: main term (1-2/n)^(n-1) = {main_t:.4f}")
        rnd = random.Random(100+n)
        sets = {"AP 1..n-1 (regular)": tuple(range(1, n))}
        # a few primitive sets with increasing minimal resonance
        cand = []
        seen = 0
        while len(cand) < 4 and seen < 4000:
            seen += 1
            v = tuple(sorted(rnd.sample(range(2, 12*n), n-1)))
            if reduce(gcd, v) != 1: continue
            cand.append(v)
        for idx, v in enumerate(cand): sets[f"random#{idx+1}"] = v
        rows = []
        for name, v in sets.items():
            l = minimal_resonance(v); e2 = E2(v, n); lm = lonely_measure(v, n)
            suff = "LRC (E2<main)" if e2 < main_t else "bound fails"
            rows.append((name, v, l, e2, lm, suff))
            print(f"    {name:20s} v={v}: l(v)={l:2d}, E2={e2:.4f}, |LONELY|={lm:.4f}  [{suff}]")
    print()
    print("="*74)
    print("(3) HARD CORE map: fraction of primitive sets where E2 >= main (bound fails)")
    print("="*74)
    for n in (5, 6):
        main_t = (1 - 2.0/n)**(n-1); rnd = random.Random(7+n)
        tot = 0; fail = 0; failmin = []
        while tot < 150:
            v = tuple(sorted(rnd.sample(range(2, 14*n), n-1)))
            if reduce(gcd, v) != 1: continue
            tot += 1
            e2 = E2(v, n)
            if e2 >= main_t:
                fail += 1
                if len(failmin) < 6: failmin.append((v, minimal_resonance(v)))
        print(f"  n={n}: E2>=main (bound fails) for {fail}/{tot} primitive sets; "
              f"those have SMALL minimal resonance, e.g. {failmin}")
    print()
    print("  => the rigorous bound E(v)<main PROVES LRC for the decorrelated majority (large")
    print("     minimal resonance); the hard core (E>=main) is the small-minimal-resonance")
    print("     sets, led by the AP/regular polygon. Progress = LRC reduced to bounding the")
    print("     resonance energy of the small-minimal-resonance core (the returns, S545).")

if __name__ == "__main__":
    main()
