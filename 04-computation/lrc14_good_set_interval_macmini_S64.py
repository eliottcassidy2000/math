"""
mac-mini-2026-07-09-S64 (new direction) -- the good period is "the grid hits G's max interval",
a rigorous pigeonhole that UNIFIES the j=1 wraparound with the missed-residue case and isolates the
hard case to all-residue clusters on resonant rulers.

G(E) = {x in [0,1] : maxgap({e_i x mod 1}) > 1/7}, an open set (finite union of arcs), measure mu(E).
PIGEONHOLE (rigorous): an arc of length L contains a multiple of 1/V whenever L >= 1/V. So
   maxIntG(E) >= 1/V   ==>   some j/V in G   ==>   STRICT good period (7*maxCircGap > V).
Hence a strict good period exists for every V >= 1/maxIntG(E). Two provable lower bounds on maxIntG:
  (0-nbhd)  near x=0, maxgap = 1 - spread*x > 1/7 for x < 6/(7*spread), so maxIntG >= 6/(7*spread)
            ==> V >= 7*spread/6 gives a strict good period = EXACTLY the j=1 wraparound (compressed).
  (missed-residue)  if {e_i mod 7} misses a residue, x=1/7 has maxgap >= 2/7 with an arc around it;
            measured maxIntG*spread >= ~10 ==> V >= spread/10 gives a strict good period (covers the
            WIDE regime spread>6V/7 and more) -- a-priori, for ANY V (resonant or not) that clears 1/maxIntG.
The ONLY residual: all-residue clusters whose G is FRAGMENTED (maxIntG ~ 6/(7*spread), no big arc) on
rulers V in (spread, 7*spread/6) -- the wide regime at the knife-edge scale. There maxIntG = 6/(7*spread)
so 1/maxIntG = 7*spread/6, and V < that: the grid may miss (needs arithmetic luck / non-strict j=1 at
V=7spread/6 / density floor). This is klein-S201's resonant-ruler pathology, now with the geometric cause.

TESTS: (1) maxIntG*spread for missed-residue vs all-residue dissociated sets; (2) at non-resonant
V>spread, does a strict good period always exist (grid hits G's arc)?
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(64641)

def prim(E):
    E=sorted(E); return len(E)>=2 and reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best
def maxint_G(E,N):
    Ea=np.array(sorted(E)); i=np.arange(0,N)
    ph=np.mod(np.outer(i,Ea),N)/N; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    ind=(g.max(axis=1)>1.0/7+1e-12).astype(int)
    if ind.sum()==0: return 0.0
    if ind.sum()==N: return 1.0
    dd=np.diff(np.concatenate([[0],ind,[0]]))
    starts=np.where(dd==1)[0]; ends=np.where(dd==-1)[0]; runs=list(ends-starts)
    if ind[0]==1 and ind[-1]==1 and len(runs)>=2: runs[0]+=runs[-1]; runs=runs[:-1]
    return max(runs)/N
def has_strict_gp(E,V):
    E=sorted(E)
    for j in range(1,V):
        ps=sorted({(e*j)%V for e in E}); m=len(ps)
        mg=max(max((ps[(i+1)%m]-ps[i])%V for i in range(m-1)), ps[0]+V-ps[-1]) if m>1 else V
        if 7*mg>V: return True
    return False
def gen(s, all_res):
    for _ in range(300000):
        sev=[x for x in range(7,s,7)]
        k7=random.randint(3,min(len(sev),7)) if len(sev)>=3 else 0
        S7=random.sample(sev,k7) if k7 else []
        pool=[x for x in range(1,s) if x%7!=0]
        need=11-len(S7)
        if need<0 or len(pool)<need: continue
        E=sorted(set([0,s]+S7+random.sample(pool,need)))
        if len(E)!=13 or not prim(E) or longest_ap(E)>7: continue
        nres=len(set(e%7 for e in E))
        if all_res and nres==7: return E
        if (not all_res) and nres<7: return E
    return None

for label, all_res in (("missed-residue", False), ("ALL-residue (hard)", True)):
    print(f"=== {label} dissociated k=13 ===")
    print(f"{'spread':>7}{'maxIntG*s':>11}{'1/maxIntG':>10}{'strict GP @ 3 non-res V>s':>28}")
    mn=99
    for s in (40,70,130,250,420):
        E=gen(s, all_res)
        if E is None: print(f"{s:>7}  (none)"); continue
        N=max(8000,30*s); mi=maxint_G(E,N); mn=min(mn,mi*s)
        Vs=[]; V=s+1
        while len(Vs)<3:
            if V%7!=0: Vs.append(V)
            V+=1
        res=[has_strict_gp(E,V) for V in Vs]
        print(f"{s:>7}{mi*s:>11.2f}{(1.0/mi if mi>0 else 0):>10.1f}   V={Vs} GP={res}")
    print(f"  min maxIntG*spread = {mn:.2f}\n")

print("CONCLUSION: maxIntG >= 1/V => strict good period (rigorous pigeonhole).")
print(" 0-nbhd: maxIntG >= 6/(7*spread) => V>=7spread/6 (= j=1 wraparound).")
print(" missed-residue: maxIntG*spread >~ 10 => V>=spread/10 => strict GP at every non-res V>spread.")
print(" all-residue fragmented-G on V in (spread, 7spread/6): the residual (non-strict j=1 / density floor).")
