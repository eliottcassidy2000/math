"""
mac-mini-2026-07-07-S41 (HYP-4817, part 3) -- corrected 13-point descents + the
WHERE-DOES-E[U]-LIVE decomposition (mechanism for the floor).

FIXES from part 2: 'opus stretched' is {0,2,3,4,5,6,7,8,9,10,12,17,28} (13 points,
skips 11); part-2 records 0.0800/0.2441 were 14-point artifacts. All descents now
assert 13 points.

NEW QUESTION (the how/why): at the E[U_{1/7}]-minimizing families, where in x-space
does the surviving avoidance mass live?
  Decompose E[U] = int U(x) dx into Farey shells: x within delta_q = 1/(2 q Qmax D)
  ... simpler: shells around rationals p/q, q <= Qcut, of width w_q = c/q^2 (Farey
  cell scale), remainder = 'generic x'.
  If the mass concentrates near small-q rationals -> the provable route is the
  rational-neighborhood (three-gap) bound; if generic -> decorrelation.
Also: verify PZ >= (7/6)E[U] domination on the bank, and the Bonferroni-by-weight
sandwich main+w3 <= E[U] <= main+w3+w4 on corrected families.
"""
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations
import random as rnd
rnd.seed(414141)

THETA=1/7; MP=14249/252252; BAR_EU=(6/7)*MP
GRID_F=200_000; xs_f=(np.arange(GRID_F)+0.5)/GRID_F
GRID_C=30_000;  xs_c=(np.arange(GRID_C)+0.5)/GRID_C

def U_vec(E, xs):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return np.clip(g-THETA,0,None).sum(axis=1)

def U_stats(E, xs):
    U=U_vec(E,xs); return U.mean(), (U*U).mean(), (U>0).mean()

def norm(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]
    g=reduce(gcd,E[1:]) if len(E)>1 else 1
    return tuple(e//g for e in E) if g else tuple(E)

BANK={
 'AP {1..13}': list(range(1,14)),
 'GW {1..11,13,24}': list(range(1,12))+[13,24],
 'death-star 2*{1..12}u{13}': [2*i for i in range(1,13)]+[13],
 'monad record evens+{11,13}': [2,4,6,8,10,11,12,13,14,16,18,20,22],
 'opus stretched (13pt, no 11)': [0,2,3,4,5,6,7,8,9,10,12,17,28],
 'monad S1-min': [6,8,10,11,12,13,14,16,18,25,36,43,61],
 'S41p2 EU-min 13pt': [0,30,36,45,50,54,60,63,70,72,75,90,108],
 'S41p2 PZ-min 13pt': [0,2,4,5,6,7,8,9,10,11,12,14,16],
 'random big': sorted(rnd.sample(range(1,2000),13)),
}
for name,E in BANK.items(): assert len(set(E))==13, name

print("=== corrected bank at theta=1/7 (bars: E[U]>=%.5f, PZ>=%.5f) ===" % (BAR_EU,MP))
print(f"{'family':30s} {'E[U]':>7s} {'PZ':>7s} {'(7/6)EU':>8s} {'mu':>6s} {'PZ>=(7/6)EU?':>12s}")
for name,E in BANK.items():
    EU,EU2,mu=U_stats(E,xs_f); pz=EU*EU/EU2
    print(f"{name:30s} {EU:7.4f} {pz:7.4f} {7*EU/6:8.4f} {mu:6.3f} {str(pz>=7*EU/6-1e-9):>12s}")

def descend(obj, starts, iters=1200):
    results=[]
    for E0 in starts:
        E=list(E0); assert len(set(E))==13
        cur=obj(E)
        for it in range(iters):
            i=rnd.randrange(13); cand=E.copy()
            mv=rnd.random()
            if mv<0.5: cand[i]=max(0,cand[i]+rnd.choice([-2,-1,1,2]))
            elif mv<0.8: cand[i]=max(0,cand[i]+rnd.choice([-9,-6,-4,4,6,9]))
            else: cand[i]=rnd.randint(0,120)
            if len(set(cand))!=13: continue
            cv=obj(cand)
            if cv<cur: E,cur=cand,cv
        results.append((cur,norm(E)))
    results.sort(); return results

STARTS=[list(range(1,14)),
        [2,4,6,8,10,11,12,13,14,16,18,20,22],
        [0,2,3,4,5,6,7,8,9,10,12,17,28],
        [6,8,10,11,12,13,14,16,18,25,36,43,61],
        [0,30,36,45,50,54,60,63,70,72,75,90,108],
        sorted(rnd.sample(range(1,150),13))]

print("\n=== 13-point descents (enforced) ===")
resU=descend(lambda E:U_stats(E,xs_c)[0], STARTS)
print("min E[U] results:")
for cur,E in resU[:3]: print(f"  {cur:.5f}  {E}  (n={len(E)})")
EU,EU2,mu=U_stats(list(resU[0][1]),xs_f)
print(f"finalist fine: E[U]={EU:.5f} ({EU/BAR_EU:.2f}x bar); (7/6)E[U]={7*EU/6:.5f} vs m_P={MP:.5f}; mu={mu:.3f}")

def pzobj(E):
    a,b,_=U_stats(E,xs_c); return a*a/b if b>0 else 9.9
resPZ=descend(pzobj, STARTS)
print("min PZ results:")
for cur,E in resPZ[:3]: print(f"  {cur:.5f}  {E}  (n={len(E)})")
EU,EU2,mu=U_stats(list(resPZ[0][1]),xs_f)
print(f"finalist fine: PZ={EU*EU/EU2:.5f} ({(EU*EU/EU2)/MP:.2f}x m_P); E[U]={EU:.5f}; mu={mu:.3f}")

# ---- WHERE does E[U] live? Farey-shell decomposition at the minimizers ----
print("\n=== Farey-shell decomposition of E[U] (mechanism) ===")
def farey_decomp(E, Qcut=8, cwidth=0.5):
    U=U_vec(E,xs_f)
    total=U.mean()
    assigned=np.zeros(GRID_F,dtype=bool)
    rows=[]
    for q in range(1,Qcut+1):
        wq=cwidth/(q*q*14)   # Farey-cell-scale window; c chosen small vs 1/q^2
        mask=np.zeros(GRID_F,dtype=bool)
        for p in range(0,q+1):
            if gcd(p,q)!=1 and not (p in (0,) and q==1): continue
            lo,hi=p/q-wq,p/q+wq
            mask |= (xs_f>=lo)&(xs_f<=hi)
        mask &= ~assigned
        rows.append((q, mask.mean(), (U*mask).mean()))
        assigned |= mask
    rem=~assigned
    rows.append(('gen', rem.mean(), (U*rem).mean()))
    return total, rows

for name in ('AP {1..13}','S41p2 EU-min 13pt','monad record evens+{11,13}','random big'):
    E=BANK[name]; total,rows=farey_decomp(E)
    parts=' '.join(f"q{q}:{um/total*100:4.1f}%" for q,mm,um in rows)
    print(f"{name:30s} E[U]={total:.4f}  mass-> {parts}")

# ---- Bonferroni-by-weight sandwich on CORRECTED families ----
print("\n=== sandwich main+w3 <= E[U] <= main+w3+w4 (corrected 13pt families) ===")
def phi(m,u=THETA):
    if m==0: return complex(u,0)
    return (1-np.exp(-2j*np.pi*m*u))/(2j*np.pi*m)
def w3(E,u=THETA,TMAX=60):
    tot=0.0
    for (a,b,c) in combinations(E,3):
        w=(b-c,c-a,a-b); g=reduce(gcd,[abs(v) for v in w]); w=tuple(v//g for v in w)
        tot+=sum(2*(phi(t*w[0],u)*phi(t*w[1],u)*phi(t*w[2],u)).real for t in range(1,TMAX+1))
    return -tot
def w4(E,u=THETA,MMAX=4):
    tot=0.0; Earr=list(E)
    for idx in combinations(range(13),4):
        e=[Earr[i] for i in idx]
        for m1 in range(-MMAX,MMAX+1):
            if m1==0: continue
            for m2 in range(-MMAX,MMAX+1):
                if m2==0: continue
                for m3 in range(-MMAX,MMAX+1):
                    if m3==0: continue
                    m4v=-(m1+m2+m3)
                    if m4v==0 or abs(m4v)>MMAX: continue
                    if m1*e[0]+m2*e[1]+m3*e[2]+m4v*e[3]!=0: continue
                    tot+=(phi(m1,u)*phi(m2,u)*phi(m3,u)*phi(m4v,u)).real
    return tot
main=(1-THETA)**13
for name in ('AP {1..13}','monad record evens+{11,13}','opus stretched (13pt, no 11)','S41p2 EU-min 13pt'):
    E=BANK[name]; EU,_,_=U_stats(E,xs_f)
    a=main+w3(E); b=a+w4(E)
    ok='SANDWICH OK' if a<=EU<=b else 'VIOLATED'
    print(f"{name:30s} lower={a:8.4f}  E[U]={EU:7.4f}  upper={b:7.4f}   {ok}")
