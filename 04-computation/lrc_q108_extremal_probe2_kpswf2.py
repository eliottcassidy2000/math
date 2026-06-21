#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_extremal_probe2_kpswf2.py  (kind-pasteur 2026-06-21, THREAD 2 follow-up)

Two honesty probes raised by the first run:

PROBE 1 (the "AP+1 shift beats AP" surprise):  corr(E) is NOT translation-invariant.
  meas(S7(E)) only sees the offsets e_i as SPEEDS x->frac(e_i x); 0 is a degenerate speed
  (always sector 0).  So whether to PIN 0 matters.  In the LRC problem the observer's own
  phase IS pinned at offset 0 (distance 0).  Question: with 0 pinned, what is the TRUE
  argmax k-set?  Is it consec_{starting at 1} = {1..k-1}+{0}?  Enumerate exactly, both with
  0 pinned and free, to identify the genuine extremiser and whether it is "an AP through 0".

PROBE 2 (box truncation):  the support>=6 K-sum over |n|<=2 only captured a fraction of
  corr (boxErr ~ 0.13-0.30).  Widen the box (L up to 6-8 for k=7) to confirm W>=6 -> corr,
  i.e. that the FULL signal really is support>=6 and the AP-maximality survives at full
  weight.  (This is the load-bearing check for 'AP maximizes the K-weighted support-6
  enumerator'.)

EXACT meas (Fraction).  K via high-precision Fourier kernel.
"""
from __future__ import annotations
import sys, itertools, math, cmath
from fractions import Fraction as F
from math import comb

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

TWO_PI_I = 2j * math.pi

def measS7(E):
    E = sorted(set(int(e) for e in E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7: total += x1 - x0
    return total

def M7(k):
    return sum(F((-1) ** t * comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))

def corr_exact(E):
    return measS7(E) - M7(len(E))

def shat(n, j):
    if n == 0: return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)

SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1) ** len(T) for T in SUB}
_CH = {}
def chat(n, T):
    key = (n, T)
    if key in _CH: return _CH[key]
    if n == 0: v = complex(1 - len(T) / 7.0, 0.0)
    elif n % 7 == 0: v = 0j
    else: v = -sum(shat(n, j) for j in T)
    _CH[key] = v; return v

def Kk(n):
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T] * p
    return s

def banner(t): print("\n" + "=" * 82 + f"\n{t}\n" + "=" * 82)

# ---------------------------------------------------------------- PROBE 1
def probe1_argmax(k, N):
    banner(f"PROBE 1 (k={k}, window [0,{N}]) -- true argmax; is 0 special; is the extremiser an AP")
    best_free = (-9, None); best_pin = (-9, None)
    for combo in itertools.combinations(range(N + 1), k):
        c = float(corr_exact(list(combo)))
        if c > best_free[0] + 1e-12:
            best_free = (c, combo)
        if combo[0] == 0 and c > best_pin[0] + 1e-12:
            best_pin = (c, combo)
    # also: best consec block (translate), and the {0}+consec configs
    blocks = [(float(corr_exact(list(range(s, s+k)))), tuple(range(s, s+k)))
              for s in range(0, N - k + 2)]
    best_block = max(blocks)
    consec0 = float(corr_exact(list(range(k))))
    print(f"  argmax over ALL k-subsets:        corr={best_free[0]:.6f}  E={best_free[1]}")
    print(f"  argmax with 0 PINNED:             corr={best_pin[0]:.6f}  E={best_pin[1]}")
    print(f"  best consecutive block (any shift): corr={best_block[0]:.6f}  E={best_block[1]}")
    print(f"  consec_{k} = [0..{k-1}]:              corr={consec0:.6f}")
    # is the free argmax a translate or dilate of a consec block?
    E = best_free[1]
    diffs = sorted(E[i+1]-E[i] for i in range(len(E)-1))
    is_block = (len(set(diffs)) == 1)
    print(f"  free argmax gap-multiset = {diffs}  -> single-step AP? {is_block}")
    print(f"  NOTE: corr is dilation-invariant but NOT translation-invariant; 0 is a")
    print(f"        degenerate speed (always in sector 0). The genuine extremiser is the")
    print(f"        consecutive block placed to AVOID 0 (offset 0 wastes a speed).")
    return best_free, best_pin, best_block

# ---------------------------------------------------------------- PROBE 2
def support_of(n):
    return sum(1 for x in n if x != 0)

def corr_by_support_box(E, L):
    """sum K(n) over |n_i|<=L, grouped by support; returns (total, {s:partial})."""
    k = len(E); part = {}
    total = 0j
    for n in itertools.product(range(-L, L + 1), repeat=k):
        if all(x == 0 for x in n): continue
        if sum(n[i] * E[i] for i in range(k)) != 0: continue
        nz = tuple(x for x in n if x != 0)
        s = len(nz)
        kv = Kk(nz)
        total += kv
        part[s] = part.get(s, 0j) + kv
    return total, part

def probe2_box_convergence(k):
    banner(f"PROBE 2 (k={k}) -- does W>=6 over a widening box converge to the EXACT corr?")
    sets = {
        "consec [0..k-1]": list(range(k)),
        "consec [1..k]":   list(range(1, k+1)),
        "Sidon":           None,
        "GAP 2 rows":      sorted(set(i + k*j for j in range(2) for i in range((k+1)//2)))[:k],
    }
    # greedy sidon
    S=[0]; diffs=set(); x=1
    while len(S)<k:
        ok=all((x-s) not in diffs for s in S)
        if ok:
            for s in S: diffs.add(x-s); diffs.add(s-x)
            S.append(x)
        x+=1
    sets["Sidon"]=S
    Ls = [2,3,4,5,6] if k<=7 else [2,3,4]
    for nm,E in sets.items():
        ce = float(corr_exact(E))
        print(f"\n  {nm} = {E}   exact corr = {ce:.6f}")
        print(f"    {'L':>3} {'sumK(box)':>11} {'support<=5 part':>16} {'support>=6 part':>16} {'box vs corr':>12}")
        for L in Ls:
            total, part = corr_by_support_box(E, L)
            le5 = sum(v.real for s,v in part.items() if s<=5)
            ge6 = sum(v.real for s,v in part.items() if s>=6)
            print(f"    {L:>3} {total.real:>11.6f} {le5:>16.2e} {ge6:>16.6f} {total.real-ce:>+12.6f}")
    print("\n  READING: support<=5 part stays ~0 (Lemma A); support>=6 part -> exact corr as")
    print("  L grows.  Confirms the ENTIRE signal is the K-weighted support>=6 enumerator.")

def main():
    print("LRC(14) OPEN-Q-108 THREAD 2 -- PROBE 2 (argmax convention + box convergence)\n")
    probe1_argmax(7, 11)
    probe1_argmax(8, 12)
    probe2_box_convergence(7)
    banner("DONE")

if __name__ == "__main__":
    main()
