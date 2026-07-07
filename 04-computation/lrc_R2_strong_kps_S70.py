#!/usr/bin/env python3
r"""
lrc_R2_strong_kps_S70.py   (kps-S70, HYP-5077)

STRONG R2 test: my first adversary's DESCENT was weak (stuck at 0.857 for k=8, never near
the AP min 0.720) -- the MISTAKE-102 trap (isolated structured minimizers). Proper tests:
 (A) LOCAL-MINIMUM test: from the exact AP winner, exhaustively try all 1- and 2-element
     moves; can anything go below the AP? (tests whether the spread AP is even a local min)
 (B) STRUCTURED large-diameter adversary aimed at the resonant regime: spread APs with a
     defect/outlier, dilated near-APs, AP-with-one-doubled (GW-type), interleaved APs --
     the families most likely to beat a pure spread AP if R2 is false.
 (C) exact-rational confirmation of the AP winner + best challenger.
Question restated for (A'): does anything go below T_k (route fails) or below the spread-AP
min (R2 false)?
"""
import random, math
from fractions import Fraction as F

def PA2(E, res=20000):
    E = sorted(set(E)); n = len(E); c = 0
    for r in range(res):
        x = (r + .5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        best0 = ph[0] + 1 - ph[-1]; besth = None
        # gap@0
        g0 = ph[0] + 1 - ph[-1]
        for i in range(n-1):
            if ph[i] <= 0.0 < ph[i+1]: g0 = ph[i+1]-ph[i]
        # since 0 < ph[0] generically, gap@0 = wrap
        g0 = ph[0] + 1 - ph[-1]
        gh = ph[0] + 1 - ph[-1]
        for i in range(n-1):
            if ph[i] <= 0.5 < ph[i+1]: gh = ph[i+1]-ph[i]; break
        else:
            gh = ph[0] + 1 - ph[-1]
        if max(g0, gh) > 1/7: c += 1
    return c / res

def PA2_exact(E):
    E = sorted(set(E)); THETA = F(1,7)
    bps = {F(0), F(1)}
    for i in range(len(E)):
        e = E[i]
        for m in range(1, e): bps.add(F(m, e))
        for m in range(0, e): bps.add(F(2*m+1, 2*e))
        for j in range(i+1, len(E)):
            d = E[j]-E[i]
            for m in range(1, d): bps.add(F(m, d))
    bps = sorted(bps); total = F(0)
    for lo, hi in zip(bps[:-1], bps[1:]):
        mid = (lo+hi)/2; fl = {e:(e*mid).__floor__() for e in E}
        order = sorted(E, key=lambda e: e*mid-fl[e]); p = order
        aff = {e:(e, F(-fl[e])) for e in E}
        s0 = aff[p[0]][0]-aff[p[-1]][0]; b0 = aff[p[0]][1]+1-aff[p[-1]][1]; g0=(s0,b0)
        pv = [e*mid-fl[e] for e in p]; idx=None
        for i in range(len(p)-1):
            if pv[i] <= F(1,2) < pv[i+1]: idx=i; break
        if idx is None: sh,bh = s0,b0
        else: sh = aff[p[idx+1]][0]-aff[p[idx]][0]; bh = aff[p[idx+1]][1]-aff[p[idx]][1]
        gh=(sh,bh); sub={lo,hi}
        for (s,b) in (g0,gh):
            if s!=0:
                xc=(THETA-b)/s
                if lo<xc<hi: sub.add(xc)
        if g0[0]!=gh[0]:
            xc=(gh[1]-g0[1])/(g0[0]-gh[0])
            if lo<xc<hi: sub.add(xc)
        sub=sorted(sub)
        for u,v in zip(sub[:-1],sub[1:]):
            m2=(u+v)/2
            if max(g0[0]*m2+g0[1], gh[0]*m2+gh[1])>THETA: total+=v-u
    return total

Tk = {8: 0.6185, 13: 0.0565}
winner = {8: [5,7,9,11,13,15,17,19], 13: [5+7*j for j in range(13)]}

for k in (8, 13):
    T = Tk[k]; W = winner[k]; wv = PA2(W, 24000)
    print("=" * 88)
    print(f"k={k}: AP winner {W}, PA_2={wv:.4f} (T_k={T:.4f})")
    print("=" * 88)
    # (A) LOCAL-MIN test: all 1- and 2-element moves
    localbest = (wv, "winner", W)
    steps = [-4,-3,-2,-1,1,2,3,4]
    # 1-element
    for i in range(k):
        for s in steps:
            E = W[:]; E[i] += s
            E = sorted(set(e for e in E if e > 0))
            if len(E) != k: continue
            v = PA2(E, 16000)
            if v < localbest[0]: localbest = (v, f"move v{i}{s:+d}", E)
    # 2-element (sample)
    rng = random.Random(k)
    for _ in range(600):
        E = W[:]; i,j = rng.sample(range(k), 2)
        E[i] += rng.choice(steps); E[j] += rng.choice(steps)
        E = sorted(set(e for e in E if e > 0))
        if len(E) != k: continue
        v = PA2(E, 12000)
        if v < localbest[0]: localbest = (v, "2-move", E)
    print(f"  (A) local-min: best neighbor PA_2 = {localbest[0]:.4f} ({localbest[1]}: {localbest[2]})")
    print(f"      => winner is a LOCAL MIN: {localbest[0] >= wv - 0.002}")
    # (B) structured large-diameter adversary
    struct = {}
    for d in (2,3,5,7,11):
        struct[f"spreadAP d={d}"] = [1+d*j for j in range(k)]
        struct[f"spreadAP d={d} outlier"] = [1+d*j for j in range(k-1)] + [1+d*(k-1)+d*10]
        struct[f"spreadAP d={d} defect"] = [1+d*j for j in range(k) if j != k//2] + [1+d*(k//2)+1]
    struct["GW (AP one doubled)"] = list(range(1,k)) + [2*(k-1)]
    struct["interleave 2AP"] = sorted([1+4*j for j in range(k//2+1)] + [3+4*j for j in range(k-k//2-1)])[:k]
    struct["dilated {1..k} x large"] = [97*j for j in range(1,k+1)]
    bstruct = (9.9, None, None)
    below = []
    for nm, E in struct.items():
        E = sorted(set(E))
        if len(E) != k: continue
        v = PA2(E, 16000)
        if v < bstruct[0]: bstruct = (v, nm, E)
        if v < T - 0.003: below.append((v,nm,E))
        if v < wv - 0.004: print(f"      *** {nm} PA_2={v:.4f} BEATS winner {wv:.4f}: {E}")
    print(f"  (B) structured large-diam min = {bstruct[0]:.4f} ({bstruct[1]}: {bstruct[2]})")
    print(f"      below T_k? {'YES ***' if below else 'NO'}; beats winner? {'YES (R2 false)' if bstruct[0] < wv-0.004 else 'NO'}")
    # (C) exact confirm winner
    if k == 8:
        ex = PA2_exact(W)
        print(f"  (C) EXACT PA_2(winner) = {ex} = {float(ex):.5f}")
    print()
print("VERDICT below.")
