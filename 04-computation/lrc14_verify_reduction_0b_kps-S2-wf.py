#!/usr/bin/env python3
"""
PART 5/6 (fast): decoupling-bound check + adversarial meas(G_C) minimization.
Imports the interval machinery from the part-0 script.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random, sys, importlib.util

# load helpers from the part-0 file
spec = importlib.util.spec_from_file_location(
    "p0", "04-computation/lrc14_verify_reduction_0_kps-S2-wf.py")
# The part-0 file runs prints on import; suppress by redirecting stdout during import.
import io, contextlib
p0 = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(p0)

lonely_set = p0.lonely_set
lonely_measure = p0.lonely_measure
danger = p0.danger
complement = p0.complement
intersect = p0.intersect
measure = p0.measure
M = p0.M

C_drop6 = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13]
TARGET = F(7, 858)

def L_exact(C, w, h=F(1, 14)):
    G = lonely_set(C, h)
    Dw = danger(w, h)
    return measure(intersect(G, complement(Dw)))

def components(C, h=F(1, 14)):
    return len(lonely_set(C, h))

print("="*78)
print("PART 5: decoupling bound  L(S) >= (6/7)meas(G_C) - r/(7w)  -- exact check")
print("="*78)
mGc = lonely_measure(C_drop6); r = components(C_drop6)
print(f"drop-6 core: meas(G_C)={mGc}={float(mGc):.8f}, r={r} components")
viol = 0; minL = None
for m in range(1, 60):
    w = 14*m
    L = L_exact(C_drop6, w)
    bound = F(6,7)*mGc - F(r, 7*w)
    if L < bound: viol += 1
    if minL is None or L < minL[0]: minL = (L, w)
    if m <= 8:
        print(f"   w={w:4d}: L={float(L):.8f}  (6/7)meas-r/(7w)={float(bound):.8f}  L>=bound? {L>=bound}")
print(f"   violations of decoupling bound over m=1..59: {viol}")
print(f"   min L over m=1..59: {minL[0]}={float(minL[0]):.8f} at w={minL[1]}")
thresh = F(r,1)/(F(6)*mGc)
print(f"   bound > 0 requires w >= r/(6 meas) = {thresh} = {float(thresh):.2f}  (i.e. m>={float(thresh)/14:.2f})")
print(f"   So for the drop-6 core (smallest meas), the single-w decoupling bound is")
print(f"   POSITIVE only for w >= {float(thresh):.1f}; smaller w must be exact-checked")
print(f"   (and they are: min L over all w=14m is {float(minL[0]):.8f} > 0).")
print()

print("="*78)
print("PART 6: ADVERSARIAL meas(G_C) minimization -- can it beat 7/858?")
print("="*78)

# 6a coordinated multiples of 14
print("6a. {1..k} + coordinated multiples of 14:")
best = (F(1), None); tested=0
for ncore in [9,10,11]:
    core=list(range(1,ncore+1)); nl=12-ncore
    for combo in itertools.combinations([14*j for j in range(1,18)], nl):
        C=sorted(core+list(combo))
        if len(set(C))!=12: continue
        tested+=1; m=lonely_measure(C)
        if m<best[0]: best=(m,C)
print(f"   tested {tested}; MIN meas(G_C)={best[0]}={float(best[0]):.8f} at {best[1]}; <7/858? {best[0]<TARGET}")

# 6b coordinated multiples of 7,12,13,84
print("6b. {1..k} + coordinated multiples of 7,12,13,84,91,168:")
best2=(F(1),None); tested2=0
mults=sorted(set(b*j for b in [7,12,13,84,91,168] for j in range(1,8) if b*j<=600))
for ncore in [8,9,10]:
    core=list(range(1,ncore+1)); nl=12-ncore; cnt=0
    for combo in itertools.combinations(mults, nl):
        C=sorted(set(core+list(combo)))
        if len(C)!=12: continue
        cnt+=1
        if cnt>40000: break
        tested2+=1; m=lonely_measure(C)
        if m<best2[0]: best2=(m,C)
print(f"   tested {tested2}; MIN meas(G_C)={best2[0]}={float(best2[0]):.8f} at {best2[1]}; <7/858? {best2[0]<TARGET}")

# 6c smart greedy descent: candidate replacements limited, accept best move per round
print("6c. greedy descent (best-improvement, candidates 1..150) from drop-6:")
cur=list(C_drop6); cur_m=lonely_measure(cur)
for rnd in range(40):
    bestmove=None
    for idx in range(12):
        for newv in range(1,151):
            if newv in cur: continue
            cs=sorted(cur[:idx]+[newv]+cur[idx+1:])
            if len(set(cs))!=12: continue
            m=lonely_measure(cs)
            if m<cur_m and (bestmove is None or m<bestmove[0]):
                bestmove=(m,cs)
    if bestmove is None: break
    cur_m,cur=bestmove[0],bestmove[1]
print(f"   final meas(G_C)={cur_m}={float(cur_m):.8f} at {cur}; <7/858? {cur_m<TARGET}")

# 6d random primitive 12-sets, wider net incl very spread speeds
print("6d. random primitive 12-sets, speeds up to 300:")
random.seed(424242)
br=(F(1),None); nb=0; nt=0
for V in [40,80,150,300]:
    for _ in range(6000):
        C=sorted(random.sample(range(1,V+1),12))
        if reduce(gcd,C)!=1: continue
        nt+=1; m=lonely_measure(C)
        if m<br[0]: br=(m,C)
        if m<TARGET: nb+=1
print(f"   tested {nt}; #(<7/858)={nb}; MIN={br[0]}={float(br[0]):.8f} at {br[1]}")
print()
print("="*78)
print("VERDICT DATA")
print("="*78)
overall_min = min([best[0],best2[0],cur_m,br[0]])
print(f"min meas(G_C) found anywhere this script = {overall_min} = {float(overall_min):.8f}")
print(f"target 7/858 = {float(TARGET):.8f}")
print(f"any 12-set with meas(G_C) < 7/858 found? {overall_min < TARGET}")
