#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c95 -- HYP-8040: WORKING THE SMALL-WITNESS LAW
(mac-mini-S125 handoff (ii-a); the Z-face of the HYP-8030 fused theorem).

LAW (target): every duty-satisfying two-far family W = S u {x,y}
(S an 11-subset of {1..13}, 14 <= x < y) is beaten to >= 3/41 by a witness
of denominator q <= 41 -- except the divisibility (lcm) channel, which is
then closed separately by the drift/resonance argument at huge fars.

PART (a) -- THE 78 FINITE BAD_q ANALYSES, EXACT.
For each complement pair {a,b} (S = {1..13} \\ {a,b}) and q in [14, 41]:
  t_q = ceil(3q/41)   (witness needs dist_q(v*j) >= t_q for ALL v in W)
  A_q(S) = {j in [1,q-1] : all s in S clear}   (S-clear rotations)
  E_q = {(x mod q, y mod q) : no j in A_q clears both x and y}  (escapes)
  structure: divisibility part D_q = {(u,v) : u = 0 or v = 0}; CLEAN modulus
  iff E_q subset of D_q -- then any all-heights escape needs q | x or q | y.
Deliverable: per complement, the list of clean moduli; if every complement
has >= 2 coprime-ish clean moduli, escapes force divisibility chains and
the law reduces to fars >= X0 (explicit).

PART (b) -- THE DRIFT CONSTANTS.
Per complement: the 11-part's exact 3/41-safe set (breakpoint enumeration
over denominators <= 2*13; components + widest width w0(S)).  The drift
lemma needs x >= 2/w0 and y >= ~5x; print X0-candidates per complement and
the global max.  (The near-equal resonance check: (1+m)*6/41 < 1 for
m <= 5 -- arithmetic printed for the record.)

GATES: mac-mini's resistant families ({1..9,11,13,20,108} = 9/113 exit at
q <= 41 etc.) + F3(13) = 3/41 at t = 17/41 + the sweep's verdict shape
(every family exits at q <= 41) reproduced on a random sample.
"""
import sys, time, random
from fractions import Fraction as F
from math import gcd
from itertools import combinations
sys.path.insert(0, '04-computation')
from lrc14_ladder_realization_crossN_kps_S128c86 import M_exact_arg

random.seed(95)
TH = F(3, 41)

def distq(v, j, q):
    r = (v*j) % q
    return min(r, q - r)

def exit_modulus(W, qmax=41):
    """smallest q in [14, qmax] with a witness j: all dist_q >= ceil(3q/41)"""
    for q in range(14, qmax+1):
        t = -(-3*q // 41)
        for j in range(1, q//2 + 1):
            if all(distq(v, j, q) >= t for v in W):
                return q, j
    return None, None

def gates():
    print("== gates ==", flush=True)
    ok = True
    for W, wantM in ((list(range(1,10))+[11,13,20,108], F(9,113)),
                     (list(range(1,12))+[13,36], F(3,41))):
        M, arg = M_exact_arg(W)
        q, j = exit_modulus(W)
        good = (M == wantM) and (q is not None)
        ok &= good
        print(f"  {W[:5]}...: M={M} (want {wantM}) exit q={q} j={j} -> {'OK' if good else 'FAIL'}", flush=True)
    # random duty-satisfying two-far sample must all exit at q <= 41
    fails = 0; n = 0
    for _ in range(300):
        ab = random.sample(range(1, 14), 2)
        S = [v for v in range(1, 14) if v not in ab]
        x = random.randint(14, 300); y = random.randint(x+1, 301)
        W = sorted(set(S + [x, y]))
        if len(W) != 13: continue
        # duty: multiples of removed {a,b} must exist among fars (else family
        # is not duty-satisfying; skip -- the law's hypothesis)
        if not all(any(v % d == 0 for v in W) for d in range(2, 14)): continue
        if any(v % 14 == 0 for v in W): continue
        n += 1
        q, j = exit_modulus(W)
        if q is None: fails += 1
    print(f"  random duty two-far sample: {n} families, exit-fails at q<=41: {fails} "
          f"{'OK' if fails == 0 else '(!)'}  [mac-mini sweep: 0/804,119]", flush=True)
    return ok

def partA():
    print("\n== (a) the 78 finite BAD_q analyses ==", flush=True)
    summary = {}
    for ab in combinations(range(1, 14), 2):
        S = [v for v in range(1, 14) if v not in ab]
        clean = []; useless = []; rows = []
        for q in range(14, 42):
            t = -(-3*q // 41)
            A = [j for j in range(1, q) if all(distq(s, j, q) >= t for s in S)]
            if not A:
                useless.append(q); continue
            # escapes: (u,v) residue pairs no j in A clears
            esc = []
            for u in range(q):
                # precompute u-clear js
                uj = [j for j in A if distq(u, j, q) >= t]
                for v in range(q):
                    if not any(distq(v, j, q) >= t for j in uj):
                        esc.append((u, v))
            ndiv = sum(1 for (u, v) in esc if u % q and v % q)  # non-divisibility escapes
            # divisibility part: u==0 or v==0 always escape iff 0 in danger (t>0: yes)
            is_clean = (ndiv == 0)
            if is_clean: clean.append(q)
            rows.append((q, len(A), len(esc), ndiv))
        summary[ab] = (clean, useless, rows)
        ex = [r for r in rows if r[3] > 0]
        print(f"  S = 13\\{set(ab)}: clean moduli {clean if clean else 'NONE'}; "
              f"useless {useless if useless else '-'}; "
              f"moduli with NON-divisibility escapes: {[(q, nd) for q, _, _, nd in ex][:6]}", flush=True)
    ncleanall = sum(1 for ab in summary if summary[ab][0])
    two_coprime = 0
    for ab, (clean, _, _) in summary.items():
        found = any(gcd(q1, q2) <= 2 for i, q1 in enumerate(clean) for q2 in clean[i+1:])
        if found: two_coprime += 1
    print(f"  => complements with >= 1 clean modulus: {ncleanall}/78; "
          f"with a near-coprime clean pair: {two_coprime}/78", flush=True)
    return summary

def partB():
    print("\n== (b) drift constants: the 11-part's exact 3/41-safe widths ==", flush=True)
    worst = None
    for ab in combinations(range(1, 14), 2):
        S = [v for v in range(1, 14) if v not in ab]
        # exact safe set {t in [0,1]: all ||s t|| >= 3/41}: breakpoints at
        # (k +- 3/41)/s; build intervals exactly
        pts = set()
        for s in S:
            for k in range(s+1):
                for sgn in (-1, 1):
                    t = F(k, s) + sgn*TH/s
                    if 0 <= t <= 1: pts.add(t)
        pts = sorted(pts)
        best = F(0)
        for i in range(len(pts)-1):
            lo, hi = pts[i], pts[i+1]
            mid = (lo + hi)/2
            if all(min(mid*s - int(mid*s), 1 - (mid*s - int(mid*s))) >= TH for s in S):
                if hi - lo > best: best = hi - lo
        if worst is None or best < worst[0]: worst = (best, ab)
        # print only extremes to keep output light
    print(f"  minimum widest-safe-interval over all 78 complements: w0 = {worst[0]} ~ {float(worst[0]):.5f} at S = 13\\{set(worst[1])}", flush=True)
    x0 = 2/worst[0]
    print(f"  drift lemma constants: x >= 2/w0 = {float(x0):.0f}, y >= ~5x  => witness exists (elementary two-step drift)", flush=True)
    print(f"  near-equal resonance check: (1+m)*6/41 < 1 for m <= 5: "
          f"{[(m, float((1+m)*F(6,41))) for m in range(1, 6)]} -- all < 1 OK", flush=True)

if __name__ == "__main__":
    if gates():
        partA()
        partB()
    else:
        print("!! GATES FAILED", flush=True)
