#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
GROUNDING THE MULTI-TIGHTENER CONFINEMENT (THM-612 open gap), opus-2026-07-03-S61.

Confinement conjecture: a PRIMITIVE tight (M(S)=1/14) 13-family has q*=14 (NOT 14m, m>=2).
Proved for one tightener (Lemma C). Open for >=2 tighteners. Here: exact M(S) and q* (denominator of the
argmax), to (a) verify the even-block reduction (2*{1..13}: M=1/14, q*=28, imprimitive), (b) independently
probe structured f=2,3,4 primitivizations for a PRIMITIVE tight q*>14 (expect none), (c) verify the
covering/extremity facts my proof attempt rests on.

exact_M(S): F(t)=min_v ||vt|| is piecewise-linear; a local max where runners i,j bind with opposite
slopes sits at t=(k_i+k_j)/(v_i+v_j), common value (v_i k_j - v_j k_i)/(v_i+v_j). Enumerate pairs & k,
evaluate F exactly (Fraction), take the max => (M, argmax t*, q*=den(t*)).
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def norm(x):  # ||x|| = dist to nearest integer, exact
    x = x - int(x)
    if x < 0: x += 1
    return min(x, 1 - x)
def F(S, t):
    return min(norm(v * t) for v in S)
def exact_M(S):
    # local maxima of F=min_v||vt|| are at pair-crossings t=s/(vi+vj) or s/|vi-vj|, or single peaks
    # (2k+1)/(2v). Enumerate all candidate t in [0,1), evaluate F exactly, take the max.
    S = sorted(set(S)); cands = set()
    for v in S:
        for k in range(v): cands.add(Fr(2*k+1, 2*v))
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for den in (S[i]+S[j], abs(S[i]-S[j])):
                if den == 0: continue
                for s in range(den): cands.add(Fr(s, den))
    best = (Fr(0), Fr(0))
    for t in cands:
        val = F(S, t)
        if val > best[0]: best = (val, t)
    return best  # (M, argmax t*)

def prim(S): return reduce(gcd, S) == 1
def cov(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))

print("="*96); print(" (a) even-block reduction: M and q* (=denominator of argmax t*)"); print("="*96)
for name, S in [("AP {1..13}", list(range(1,14))),
                ("even block 2*{1..13}", [2*j for j in range(1,14)]),
                ("2*({1..13}\\{6}) u {5}", sorted([2*j for j in range(1,14) if j!=6]+[5])),
                ("deep well {1..12,182}", sorted(list(range(1,13))+[182]))]:
    M, t = exact_M(S)
    print(f"  {name:<26} M={str(M):>8}  =1/14? {M==Fr(1,14)}  t*={str(t):>8} q*={t.denominator:>4} "
          f"primitive={prim(S)} covering={cov(S)}")

print("\n"+"="*96); print(" (b) search structured f=2,3,4 primitivizations of even blocks for PRIMITIVE tight q*>14"); print("="*96)
# even block on U0 = a subset of {1..13} of size e; F = f odd tighteners; e+f=13.
# take U0 = {1..e} (or {1..13} minus some), E=2*U0, add f odd tighteners <= 27, require primitive tight q*=28.
found = []
count = 0
for e in [12, 11, 10]:
    f = 13 - e
    U0s = [list(range(1, e+1))]
    if e < 13:
        U0s.append([x for x in range(1, e+2) if x != 6][:e])  # skip 6 variant
    for U0 in U0s:
        E = [2*u for u in U0]
        odds = [w for w in range(1, 28) if w % 2 == 1 and w not in E]
        for F_ in itertools.combinations(odds, f):
            S = sorted(E + list(F_)); count += 1
            if len(set(S)) != 13 or not prim(S): continue
            M, t = exact_M(S)
            if M == Fr(1,14) and t.denominator > 14:
                found.append((S, t.denominator))
print(f"  searched {count} structured even-block+odd-tightener families (e=10,11,12).")
print(f"  PRIMITIVE tight (M=1/14) with q*>14 found: {len(found)}")
for S, q in found[:8]: print(f"     q*={q}  S={S}")
if not found:
    print("  => 0 hits (independent confirmation of mac-mini's confinement search on this slice).")

print("\n"+"="*96); print(" (c) covering/extremity facts the proof rests on"); print("="*96)
# fact: for odd w, ||w(t+1/2)|| = 1/2 - ||wt||  (so D_w disjoint from D_w+1/2)
import random
random.seed(0); ok=True
for _ in range(100000):
    w = 2*random.randint(0,50)+1; t=Fr(random.randint(0,997),998)
    if norm(w*(t+Fr(1,2))) != Fr(1,2)-norm(w*t): ok=False; break
print(f"  odd w: ||w(t+1/2)|| = 1/2 - ||wt||  (D_w ∩ (D_w+1/2)=∅):  verified {ok}")
print("  => a single odd tightener covers <=1 of each 2-orbit {t,t+1/2}; f>=2 needed (Lemma C).")
print("DONE.")
