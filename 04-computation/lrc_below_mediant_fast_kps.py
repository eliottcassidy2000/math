#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FAST integer-arithmetic filtered scan for below-mediant LRC configs (a_max growth).
kind-pasteur-2026-06-19-S9.  (v2 of lrc_below_mediant_filtered_kps.py; integer probes.)

M(S) < mediant=2/(2k+1)  <=>  no time t=m/q has min_v ||vt|| >= mediant.
Integer test at t=m/q:  min_v min(vm mod q, q - vm mod q) >= ceil_threshold,
where ||vt||>=2/(2k+1) <=> (2k+1)*min(vm%q, q-vm%q) >= 2q.
Probes ordered with the strong killers (q=k+1,k,k-1,...) first.
Survivors get an EXACT Fraction maxmin (rare).  SOUND: never discards a below-mediant config.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None

def fn(x):
    r = x - (x.numerator // x.denominator); return min(r, 1 - r)

def maxmin_exact(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(i, len(S)):
            for d in (S[i] + S[j], S[i] - S[j]):
                if d == 0: continue
                d = abs(d)
                for m in range(1, d): cand.add(Fraction(m, d))
    best = Fraction(0); bt = Fraction(0)
    for t in cand:
        mv = min(fn(v * t) for v in S)
        if mv > best: best, bt = mv, t
    return best, bt

def primitive(S): return reduce(gcd, S) == 1

def ordered_probe_denoms(k, Qmax):
    """denominators ordered so the strongest 'killer' clocks come first."""
    front = [k + 1, k, k - 1, k + 2, k - 2, 2 * k + 1, 2 * k - 1, 2 * k, 2 * k + 3,
             2 * k + 2, k + 3, k - 3, 2 * k - 3, 3 * k + 2, 3 * k - 1]
    seen = set(); order = []
    for q in front + list(range(2, Qmax + 1)):
        if 2 <= q <= Qmax and q not in seen:
            seen.add(q); order.append(q)
    return order

def make_probes(k, Qmax):
    """list of (m,q) reduced, ordered."""
    pr = []
    for q in ordered_probe_denoms(k, Qmax):
        for m in range(1, q):
            if gcd(m, q) == 1:
                pr.append((m, q))
    return pr

def scan(k, box, Qprobe=None):
    twok1 = 2 * k + 1                  # mediant = 2/(2k+1)
    floor = Fraction(1, k + 1)
    mediant = Fraction(2, twok1)
    if Qprobe is None: Qprobe = 3 * k + 3
    probes = make_probes(k, Qprobe)
    survivors = []; ns = 0; npass = 0
    for tail in itertools.combinations(range(2, box + 1), k - 1):
        S = (1,) + tail
        if reduce(gcd, S) != 1: continue
        ns += 1
        skip = False
        for (m, q) in probes:
            # min over v of min(vm%q, q-vm%q); test >= 2q/(2k+1)  <=> (2k+1)*mn >= 2q
            mn = q
            for v in S:
                r = (v * m) % q
                d = r if r < q - r else q - r
                if d < mn:
                    mn = d
                    if twok1 * mn < 2 * q:   # this v already below mediant at this t
                        break
            if twok1 * mn >= 2 * q:
                skip = True; break
        if skip: continue
        npass += 1
        M, ts = maxmin_exact(S)
        if floor < M < mediant:
            survivors.append((M, S, ts))
    survivors.sort(key=lambda x: x[0])
    return survivors, ns, npass

def express(M, k):
    p, q = M.numerator, M.denominator
    e = p * (k + 1) - q
    return (p if e == 1 else None), e, p, q

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks", type=str, default="9,10,11,12,13")
    ap.add_argument("--boxes", type=str, default="")  # optional explicit "k:box,..."
    args = ap.parse_args()
    ks = [int(x) for x in args.ks.split(",")]
    boxmap = {}
    if args.boxes:
        for kv in args.boxes.split(","):
            a, b = kv.split(":"); boxmap[int(a)] = int(b)
    print("=== FAST below-mediant scan: a_max(k) growth test ===")
    for k in ks:
        box = boxmap.get(k, int(2.9 * k) + 2)
        max_box_a = (2 * box) // (k + 1)
        survivors, ns, npass = scan(k, box)
        floor = Fraction(1, k + 1)
        if not survivors:
            print(f"  k={k:2d} box={box} (catch a<={max_box_a}): NO below-mediant. "
                  f"scanned={ns} passed={npass}", flush=True)
            continue
        s2 = survivors[0][0]; a, e, p, q = express(s2, k)
        avals = sorted(set(express(M, k)[0] for M, _, _ in survivors if express(M, k)[0]))
        g = s2 - floor
        print(f"  k={k:2d} box={box} (catch a<={max_box_a}): sigma_2={s2} (a={a},e={e}) "
              f"g*k^2={float(g)*k*k:.4f} #below={len(survivors)} mediant-a={avals} "
              f"maxS={max(survivors[0][1])} wit={survivors[0][1]}", flush=True)
