"""
S73b (v2): Ferromagnetic merge + creative niche sets, with SET SIZE controlled.
Frame (merging kps's two-maps + the antiferromagnetic-tournament + my Joukowski bridge):
  TOURNAMENT I(Omega,x): real-rooted = hard-core REPULSIVE / ANTI-ferromagnet (frustrated).
  LRC coverage Q(z): Lee-Yang CIRCLE = FERROMAGNETIC Ising; cap = the 7-cyclotomic FM ideal.
  Joukowski w=z+R^2/z bridges circle->real; DIP = off-circle = FRUSTRATION; AP = unfrustrated FM ground state.
PART A: FIXED SIZE (13 runners = LRC(14)) -- fair structural FM ranking. Is the AP the FM ground state?
PART B: COOLING -- adding runners (k=7..13) = lowering temperature -> spread shrinks toward the FM circle.
Signal: circle spread max|r|/min|r| (1=perfect FM circle); Joukowski max|Im w| (frustration / real-rootedness defect).
"""
from fractions import Fraction as F
from math import comb, cos, pi
import numpy as np

DEMOIVRE=sorted(2*cos(2*pi*j/7) for j in (1,2,3))
def sector_of(p): return int((p%1)*7)
def missdist(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        t=7-len(set(sector_of(e*((x0+x1)/2)) for e in E))
        if 0<=t<=6: q[t]+=x1-x0
    return q
def analyze(E):
    q=missdist(E)
    if any(q[t]==0 for t in range(7)): return None
    c=[float(q[t]) for t in range(6,-1,-1)]; z=np.roots(c)
    R=(float(q[0]/q[6]))**(1/6)
    rad=sorted(abs(zi) for zi in z); spread=rad[-1]/rad[0]
    w=[zi+R*R/zi for zi in z]; maxIm=max(abs(wi.imag) for wi in w)
    return dict(R=R,p0=float(q[0]),spread=spread,maxIm=maxIm)

def first(seq,n,start=1):
    out=[]; x=start
    while len(out)<n:
        v=seq(x)
        if v not in out: out.append(v)
        x+=1
    return out
primes=[2,3,5,7,11,13,17,19,23,29,31,37,41]
fib=[1,2,3,5,8,13,21,34,55,89,144,233,377]
dyad=[2**i for i in range(13)]
squares=[i*i for i in range(1,14)]
# Mian-Chowla greedy Sidon (OEIS A005282)
mc=[1,2,4,8,13,21,31,45,66,81,97,123,148]
sets13={
 "AP {1..13} (Paley/regular ground)": list(range(1,14)),
 "2*AP {2..26} (dilation)":           list(range(2,27,2)),
 "GW {1..11,13,24}":                  list(range(1,12))+[13,24],
 "primes":                            primes,
 "Fibonacci":                         fib,
 "dyadic {1,2,..,4096}":              dyad,
 "squares {1,4,..,169}":              squares,
 "Sidon (Mian-Chowla)":               mc[:13],
 "random spread":                     [1,4,9,15,23,30,41,52,66,77,88,99,113],
}
print("="*92)
print(" PART A -- FIXED SIZE (13 runners = LRC(14)): fair FM ranking (smallest spread = most ferromagnetic)")
print("="*92)
print(f"{'13-set (+observer 0)':<36}{'R':>8}{'p0(cov)':>9}{'circle spread':>14}{'frustration|Im|':>16}")
rows=[]
for name,sp in sets13.items():
    r=analyze([0]+list(sp))
    if r is None: print(f"{name:<36}{'-- degenerate --':>40}"); continue
    rows.append((name,r))
for name,r in sorted(rows,key=lambda kv: kv[1]['spread']):
    star=" <== FM GROUND STATE" if r['spread']==min(x[1]['spread'] for x in rows) else ""
    print(f"{name:<36}{r['R']:>8.3f}{r['p0']:>9.4f}{r['spread']:>14.3f}{r['maxIm']:>16.5f}{star}")
print("-"*92)
print(f" de Moivre/cyclotomic FM ideal: 2cos(2pi j/7)={['%.3f'%x for x in DEMOIVRE]}; perfect circle spread=1.0")

print("\n"+"="*92)
print(" PART B -- COOLING: AP/consec at growing size k (= lowering temperature) -> spread -> 1 (FM circle)")
print("="*92)
print(f"{'consec {0..k-1}':<20}{'k':>4}{'R':>8}{'circle spread':>14}{'frustration|Im|':>16}")
for k in range(7,14):
    r=analyze(tuple(range(k)))
    if r is None: print(f"{'consec':<20}{k:>4}{'  -- degenerate (q_t=0) --':>30}"); continue
    print(f"{'consec':<20}{k:>4}{r['R']:>8.3f}{r['spread']:>14.3f}{r['maxIm']:>16.5f}")
print(" => more runners = colder = tighter Lee-Yang circle: the FM ground state is the cold/full-AP limit.")
