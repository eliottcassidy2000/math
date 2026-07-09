#!/usr/bin/env python3
r"""
lrc14_detuned_dispatch_monad_S3.py   (monad-explorer-2026-07-09-S3, HYP-5727/THM-668)

THE DETUNED-HARMONIC DISPATCH (d = 1) -- explicit witnesses + the d >= 2 frontier map.

THM-668: v = g*H u {delta}, g >= 2, |H| = 12, g not| delta  ==>  M(v) >= 1/13 > 1/14.
Witness: tau* = (u* + j*)/g with u* an LRC(13) time for H and j* the branch aligning
the detuned runner (coset spacing gcd(delta,g)/g <= 1/2 => clearance >= 1/4).

PART 1: explicit exact witnesses on the S2 residual instances + ratio>13 + g=7,2 cases.
        (u* found on the pair-sum rulers of H -- mac-mini THM-666; clearances exact.)
PART 2: the UNION-COVERAGE HUNT at d = 2 (the open frontier): primitive covering
        ratio>13 families g*H11 u {d1, d2}; which are caught by [phi-interval
        composition] / [direct branch search] / neither.
"""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import numpy as np

_s = open('/home/bigo/math/04-computation/lrc14_clamp_port_composition_monad_S2.py').read()
exec(_s[:_s.rfind('if __name__')])
_s2 = open('/home/bigo/math/04-computation/lrc14_phi_interval_composition_monad_S2.py').read()
_b = _s2[:_s2.rfind('if __name__')]
_b = _b[_b.index('def is_primitive'):]
exec(_b)

def nearint(x):
    y = x % 1
    return min(y, 1 - y)

def exact_witness_H(H, target=F(1, 14)):
    """Find u* on the pair-sum rulers of H with min_m ||m u*|| as large as possible;
    returns (clearance, u*). THM-666 (mac-mini): the maximizer lives at p/(h_i+h_j)."""
    H = sorted(set(H))
    best = (F(0), None)
    dens = sorted(set(a + b for i, a in enumerate(H) for b in H[i:]))
    for den in dens:
        for p in range(1, den):
            if gcd(p, den) != 1 and False:
                pass
            u = F(p, den)
            c = min(nearint(m * u) for m in H)
            if c > best[0]:
                best = (c, u)
    return best

def dispatch_witness(g, H, delta):
    """THM-668 witness: exact tau*, clearance table."""
    cH, u = exact_witness_H(H)
    # branch pigeonhole
    bestj = (F(0), None)
    for j in range(g):
        c = nearint(F(delta) * (u + j) / g)
        if c > bestj[0]:
            bestj = (c, j)
    cD, j = bestj
    tau = (u + j) / g
    v = sorted([g * m for m in H] + [delta])
    clear = min(nearint(vi * tau) for vi in v)
    return dict(u=u, cH=cH, j=j, cD=cD, tau=tau, minclear=clear, v=v)

if __name__ == "__main__":
    print("=" * 100)
    print("PART 1 -- THM-668 explicit witnesses (exact rationals)")
    print("=" * 100)
    cases = [
        ("S2 fail-1: 14*(AP\\6) + 83", 14, [1,2,3,4,5,7,8,9,10,11,12,13], 83),
        ("S2 fail-2: + 85",            14, [1,2,3,4,5,7,8,9,10,11,12,13], 85),
        ("S2 2pert: 14*(AP\\5) + 69",  14, [1,2,3,4,6,7,8,9,10,11,12,13], 69),
        ("ratio>13: 14*{1..11,20}+167",14, [1,2,3,4,5,6,7,8,9,10,11,20], 167),
        ("ratio>13, gcd7: +161",       14, [1,2,3,4,5,6,7,8,9,10,11,20], 161),
        ("g=7: 7*{1..12} + 90",        7,  [1,2,3,4,5,6,7,8,9,10,11,12], 90),
        ("g=2: 2*{1..12} + 13",        2,  [1,2,3,4,5,6,7,8,9,10,11,12], 13),
    ]
    for name, g, H, delta in cases:
        d = dispatch_witness(g, H, delta)
        ok = d['minclear'] >= F(1, 14)
        print(f"  [{name}]")
        print(f"    u* = {d['u']} (H-clearance {d['cH']} = {float(d['cH']):.4f}), "
              f"branch j* = {d['j']} (detuned clearance {float(d['cD']):.4f})")
        print(f"    tau* = {d['tau']}, min clearance = {d['minclear']} = "
              f"{float(d['minclear']):.5f} {'>= 1/14 LONELY' if ok else '< 1/14 FAIL'}")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 2 -- THE d = 2 FRONTIER (union-coverage hunt, primitive covering, ratio > 13)")
    print("=" * 100)
    rng = random.Random(57270709)
    # constructions: g = 14, H11 subset of {1..11, 14..20} (11 multipliers), D = {d1, d2}
    # covering: harmonics 14*H cover q iff exists m in H with q | 14m; D must carry the rest.
    def covers_q(g, H, D, q):
        return any((g * m) % q == 0 for m in H) or any(dd % q == 0 for dd in D)
    tested = comp_fire = branch_fire = neither = 0
    neither_list = []
    for trial in range(120):
        Hbase = list(range(1, 12))
        # ensure ratio > 13: swap one multiplier up
        i = rng.randrange(11)
        Hbase[i] = rng.randint(14, 22)
        H = sorted(set(Hbase))
        if len(H) != 11:
            continue
        g = 14
        # missing q's after harmonics
        missing = [q for q in range(2, 15) if not any((g*m) % q == 0 for m in H)]
        # d1 carries them (multiple of lcm(missing) not divisible by 14), d2 adversarial
        from math import lcm
        L = 1
        for q in missing:
            L = lcm(L, q)
        if L == 1:
            L = 13
        k1 = rng.randint(1, max(1, (20 * g) // L))
        d1 = L * k1
        if d1 % 14 == 0:
            d1 += L
        d2 = 14 * rng.choice([m for m in range(1, 22) if m not in H]) + rng.choice([-1, 1])
        v = sorted(set([g * m for m in H] + [d1, d2]))
        if len(v) != 13 or not is_primitive(v) or not is_covering(v):
            continue
        Vmax = v[-1]
        if Vmax > 320 or max(v) <= 13 * min(v):   # need ratio > 13 (outside spread13)
            continue
        tested += 1
        got = False
        for S_L in sorted(set([Vmax // 3, Vmax // 2, 2 * Vmax // 3, int(0.8 * Vmax), Vmax - 1])):
            fires, _, _, _ = phi_composition(v, S_L)
            if fires:
                got = True
                break
        if got:
            comp_fire += 1
            continue
        # direct branch search: u* from the 11 harmonics, then best j over g branches
        cH, u = exact_witness_H(H)
        found = False
        for j in range(g):
            tau = (u + j) / g
            if min(nearint(vi * tau) for vi in v) >= F(1, 14):
                found = True
                break
        if found:
            branch_fire += 1
        else:
            neither += 1
            neither_list.append((v, [d1, d2]))
    print(f"  tested {tested} primitive covering ratio>13 d=2 constructions:")
    print(f"    phi-composition fires: {comp_fire}")
    print(f"    branch search (u* + j) fires: {branch_fire}")
    print(f"    NEITHER (the true d=2 frontier): {neither}")
    for v, D in neither_list[:6]:
        print(f"      open: {v}  (detuned {D})")
    sys.stdout.flush()
