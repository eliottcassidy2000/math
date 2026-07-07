#!/usr/bin/env python3
"""
death-star-2026-07-07-S1 -- mu_1/7 AP-minimality across ALL cluster sizes k=8..13.

Route-1's floor = min over k of mu_1/7(cluster of size k); binding at k=13 (=477/1078).
opus-S130 checked k=8..13 by DESCENT only.  I showed descent MISSES structured minimizers
(it did for E[maxgap]).  So re-check each k with STRUCTURED adversaries + descent, exact-confirm.
mu_1/7 is dilation/translation-invariant, so AP_k = {1..k}.
"""
from fractions import Fraction as F
import random
from math import gcd
from functools import reduce

THR = F(1, 7)

def regions(E):
    E = list(E)
    kd = max(max(abs(j) for j in E), max(abs(a-b) for a in E for b in E))
    bps = set([F(0), F(1)])
    for d in range(1, kd+1):
        for m in range(1, d): bps.add(F(m, d))
    bps = sorted(bps); out = []
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2; fl = {j: (j*mid).__floor__() for j in E}
        order = sorted(E, key=lambda j: (j*mid - fl[j])); n = len(order); gaps = []
        for s in range(n):
            j1 = order[s]; j2 = order[(s+1) % n]
            if s < n-1: c = F(j2-j1); b0 = F(-(fl[j2]-fl[j1]))
            else: c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]])+1)
            gaps.append((c, b0))
        sub = set([a, b])
        for i in range(n):
            ci, bi = gaps[i]
            for jx in range(i+1, n):
                cj, bj = gaps[jx]
                if ci != cj:
                    xc = (bj-bi)/(ci-cj)
                    if a < xc < b: sub.add(xc)
        sub = sorted(sub)
        for u, v in zip(sub, sub[1:]):
            m2 = (u+v)/2; cb, bb = max(gaps, key=lambda cb: cb[0]*m2+cb[1]); out.append((u, v, cb, bb))
    return out

def mu17_exact(E):
    tot = F(0)
    for u, v, c, b0 in regions(E):
        if c == 0:
            if b0 > THR: tot += (v-u)
        elif c > 0:
            xc = (THR-b0)/c; lo = xc if xc > u else u
            if lo < v: tot += (v-lo)
        else:
            xc = (THR-b0)/c; hi = xc if xc < v else v
            if hi > u: tot += (hi-u)
    return tot

def mu17_num(E, res):
    E = list(E); thr = 1.0/7; cnt = 0
    for r in range(res):
        x = (r+0.5)/res
        ph = sorted((e*x) % 1.0 for e in E)
        mg = ph[0]+1.0-ph[-1]; prev = ph[0]
        for p in ph[1:]:
            g = p-prev
            if g > mg: mg = g
            prev = p
        if mg > thr: cnt += 1
    return cnt/res

def primitive(E):
    g = reduce(gcd, E); return tuple(sorted(e//g for e in E))

def struct_adversaries(k):
    """structured k-families to test against the AP."""
    outs = []
    outs.append([2*i for i in range(1, k)] + [k-1 if (k-1) % 2 else k])  # ~prim-sat analog
    outs.append([2*i for i in range(1, k)] + [1])           # 2*AP + unit
    outs.append(list(range(1, k)) + [k+1])                  # {1..k-1, k+1}
    outs.append(list(range(1, k)) + [2*k])                  # {1..k-1, far}
    outs.append([1] + [2*i for i in range(1, k)])           # unit + evens
    outs.append(list(range(1, k+1))[:-1] + [k+3])
    # clean to distinct k-sets
    res = []
    for E in outs:
        E = sorted(set(E))
        if len(E) == k: res.append(E)
    return res

if __name__ == "__main__":
    print("=== mu_1/7 AP-minimality across k=8..13 (structured + descent adversaries) ===\n")
    print(f"{'k':>3} {'mu_1/7(AP_k)':>18} {'~':>8} {'min structured':>15} {'descent min~':>13} {'AP-min?':>8}")
    random.seed(41)
    for k in range(8, 14):
        AP = list(range(1, k+1))
        muAP = mu17_exact(AP); f = float(muAP)
        # structured
        best_s = (2.0, None)
        for E in struct_adversaries(k):
            m = float(mu17_exact(E))
            if m < best_s[0]: best_s = (m, tuple(E))
        # descent (numeric)
        best_d = (2.0, None)
        for trial in range(30):
            H = random.choice([k+2, k+6, k+14, 2*k+8])
            if H < k: H = k+2
            E = sorted(random.sample(range(1, H+1), k)); cur = mu17_num(E, max(2500, 150*max(E)))
            for step in range(45):
                i = random.randrange(k); new = random.randint(1, random.choice([k+6, 2*k+6, 3*k]))
                if new in E: continue
                cand = sorted(set(E[:i]+E[i+1:]+[new]))
                if len(cand) != k: continue
                c = mu17_num(cand, max(2500, 150*max(cand)))
                if c < cur - 2e-3: E, cur = cand, c
            if cur < best_d[0]: best_d = (cur, primitive(E))
        d_exact = float(mu17_exact(list(best_d[1])))
        apmin = (best_s[0] >= f - 1e-9) and (d_exact >= f - 1e-4)
        print(f"{k:>3} {str(muAP):>18} {f:>8.5f} {best_s[0]:>15.5f} {d_exact:>13.5f} {'YES' if apmin else 'NO!':>8}")
    print("\n  If all YES: mu_1/7 AP-minimality holds at every relevant k (floor 477/1078 solid at k=13).")
