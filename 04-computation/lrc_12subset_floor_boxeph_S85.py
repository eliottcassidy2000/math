#!/usr/bin/env python3
"""
Is the 12-subset floor TRUE for primitive covering families?  (boxeph-2026-07-18-S85)
====================================================================================

THM-1010 reduces the compact core to:  M(V\{v_max}) >= (1/14)(1+1/rho),  rho=v_max/v_2nd.
Before proving it, TEST it -- skeptically, hunting adversarial (near-dilated-AP) cases.
For each primitive covering 13-family with rho<13 we record which elementary tool closes it:
  (S) SIEVE on V'=V\{v_max}: V' misses some q with 1/q >= (1/14)(1+1/rho)  => M(V')>=1/q, descent closes;
  (D) descent-from-v_max floor rho*M(V')/(rho+1) >= 1/14  (== the 12-subset floor via v_max);
  (R) best-removal recursion  max_r rho_r*M(V\{r})/(rho_r+1) >= 1/14;
  (G) direct: V' good-set (thr 1/14) meets v_max's safe set  (mu0(V')>1/7 suffices).
Goal: find families where the v_max floor (D) FAILS, and see if (S)/(R)/(G) still close them,
and whether the (D)-failures are exactly near-dilated-AP (=> rigidity, not descent).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, itertools

def norm(x):
    r = x % 1
    return min(r, 1 - r)

def exact_M(V):
    if len(V) == 1: return F(1,2)
    best = F(0); qs = set([14])
    for i in range(len(V)):
        for j in range(i, len(V)):
            s = V[i]+V[j]
            for d in range(1, s+1):
                if s % d == 0: qs.add(d)
    for q in qs:
        for a in range(1, q):
            if gcd(a,q)==1:
                m = min(min((v*a)%q, q-(v*a)%q) for v in V)
                c = F(m,q)
                if c>best: best=c
    return best

def is_covering(V, n=14):
    return all(any(v % b == 0 for v in V) for b in range(2, n+1))

def gcd_all(V):
    return reduce(gcd, V)

def missed_qs(V, n=13):
    return [q for q in range(2, n+1) if all(v % q != 0 for v in V)]

def mu0(V, thr=F(1,14)):
    """Exact measure of {t in [0,1): min_v ||v t|| >= thr}."""
    # danger(v) = union over k of (k/v - thr/v, k/v + thr/v).  Collect covered mass.
    events = []  # (position, +1 start / -1 end) of danger intervals, on [0,1) with wrap
    for v in V:
        h = thr / v  # half-width
        for k in range(v):
            c = F(k, v)
            lo, hi = c - h, c + h
            # split wrap
            lo %= 1; hi_ = c + h
            # represent interval [c-h, c+h] mod 1 as up to two [0,1) pieces
            a = (c - h) % 1
            b = (c + h) % 1
            if (c - h) < 0 or (c + h) > 1 or a > b:
                # wraps: [a,1) and [0,b]
                events.append((a, 1)); events.append((F(1), -1))
                events.append((F(0), 1)); events.append((b, -1))
            else:
                events.append((a, 1)); events.append((b, -1))
    events.sort(key=lambda e: (e[0], -e[1]))
    covered = F(0); depth = 0; prev = F(0)
    for pos, d in events:
        if depth > 0:
            covered += pos - prev
        depth += d
        prev = pos
    good = F(1) - covered
    return good

def vmax_floor_ok(V):
    Vs = sorted(V); vmax, v2 = Vs[-1], Vs[-2]
    rho = F(vmax, v2)
    Vp = Vs[:-1]
    return exact_M(Vp) >= (F(1,14))*(1 + 1/rho), rho, exact_M(Vp), (F(1,14))*(1+1/rho)

def best_removal_ok(V):
    Vs = sorted(V); best = F(0)
    for idx, r in enumerate(Vs):
        kept = Vs[:idx]+Vs[idx+1:]
        Bk = max(kept); rho = F(r, Bk)
        fl = rho*exact_M(kept)/(rho+1)
        if fl > best: best = fl
    return best >= F(1,14), best

def sieve_descent_ok(V):
    """Removing v_max uncovers q with 1/q >= (1/14)(1+1/rho)?"""
    Vs = sorted(V); vmax, v2 = Vs[-1], Vs[-2]
    rho = F(vmax, v2); need = F(1,14)*(1+1/rho)
    Vp = Vs[:-1]
    for q in missed_qs(Vp, 13):
        if F(1,q) >= need:
            return True, q
    return False, None

# ---- build a stress pool: random primitive covering + adversarial dilated-AP+killer ----
def pool():
    fams = []
    random.seed(11)
    # random primitive covering, rho<13
    cnt=0; att=0
    while cnt<40 and att<60000:
        att+=1
        V=sorted(random.sample(range(1,60),13))
        if gcd_all(V)!=1: continue
        if not is_covering(V): continue
        if F(V[-1],V[-2])>=13: continue
        fams.append(("rand",V)); cnt+=1
    # adversarial: dilated AP {d,..,12d} + a covering killer (near-tight kept)
    for d in [1,2,3,4,5]:
        ap = [d*i for i in range(1,13)]
        # killer w must make it covering + primitive + rho<13
        for w in range(1, 400):
            V = sorted(ap+[w])
            if len(set(V))!=13: continue
            if gcd_all(V)!=1: continue
            if not is_covering(V): continue
            if F(max(V),sorted(V)[-2])>=13: continue
            fams.append((f"AP{d}+{w}",V))
            break
    # more adversarial: near-dilated-AP (perturb one AP element) + killer
    for d in [2,3]:
        for drop in range(1,13):
            ap=[d*i for i in range(1,13) if i!=drop]
            for extra in [d*13, d*14, d*15]:
                base=sorted(ap+[extra])
                if len(set(base))!=12: continue
                for w in range(1,400):
                    V=sorted(base+[w])
                    if len(set(V))!=13: continue
                    if gcd_all(V)!=1: continue
                    if not is_covering(V): continue
                    if F(max(V),sorted(V)[-2])>=13: continue
                    fams.append((f"nAP{d}",V)); break
                break
    return fams

if __name__ == "__main__":
    print("="*100)
    print("IS THE 12-SUBSET FLOOR TRUE?  primitive covering, rho<13.  (D)=v_max floor")
    print("="*100)
    fams = pool()
    nD=nS=nR=nG=0; N=0; lrc_ok=True
    D_fail=[]
    for name,V in fams:
        N+=1
        M=exact_M(V)
        if M < F(1,14): lrc_ok=False; print(f"  !!! LRC(14) VIOLATION {V} M={M}")
        okD,rho,Mp,need = vmax_floor_ok(V)
        okS,q = sieve_descent_ok(V)
        okR,rf = best_removal_ok(V)
        m0 = mu0(sorted(V)[:-1])
        okG = m0 > F(1,7)
        nD+=okD; nS+=okS; nR+=okR; nG+=okG
        if not okD:
            D_fail.append((name,V,rho,Mp,need,okS,okR,okG,m0,M))
    print(f"\nsampled {N} primitive covering rho<13 families;  LRC(14) holds on all: {lrc_ok}")
    print(f"  (D) v_max 12-subset floor holds:        {nD}/{N}")
    print(f"  (S) sieve-on-V' closes (missed small q): {nS}/{N}")
    print(f"  (R) best-removal recursion reaches 1/14: {nR}/{N}")
    print(f"  (G) mu0(V') > 1/7 (measure closes):      {nG}/{N}")
    print(f"\n--- families where the v_max floor (D) FAILS ({len(D_fail)}): are they near-dilated-AP? ---")
    for name,V,rho,Mp,need,okS,okR,okG,m0,M in D_fail[:20]:
        print(f"  {name:9s} rho={float(rho):.2f} M(V')={Mp}({float(Mp):.4f}) need={float(need):.4f} "
              f"M(V)={M}  [S={okS} R={okR} G={okG} mu0(V')={float(m0):.3f}]  V={V}")
    print("\nREADING: if (D) fails only on near-dilated-AP (M(V')~1/13) AND (S or R or G) still")
    print("closes them, the compact residual = {descent OR sieve OR measure}, with the")
    print("dilated-AP stratum handled by rigidity. If some family fails ALL of D,S,R,G -> genuine crux.")
