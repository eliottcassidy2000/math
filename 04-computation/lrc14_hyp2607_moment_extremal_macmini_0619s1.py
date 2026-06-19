#!/usr/bin/env python3
"""
lrc14_hyp2607_moment_extremal_macmini_0619s1.py  (mac-mini-2026-06-19-S1)

Crack at HYP-2607 (the UNIFIED crux all 8 angles converge on):
   prove L_y(E) = E[g(N_E)] <= L_y(consec_k)  (consec maximizes the empty-sector moment fnl).
N_E(x) = # of sectors j in {1..6} EMPTY of the orbit {frac(e_i x)} (sector 0 always full).
g(t) integer-root duals (THM-534): k=8 g=(t-1)(t-2)(t-4)(t-5)/40; k=9,10 g=-(t-2)(t-3)(t-6)/36.
L_y = sum_r y_r S_r, S_r(E)=E[C(N,r)]=sum_{|A|=r} meas{sectors A all empty}.
GOAL: (1) confirm L_y closes (L_y(consec_k) <= cap_k, k=8,9,10) exactly;
      (2) decompose into component moment inequalities: which S_r does consec extremize, and
          with the dual signs, do they ALL point at consec? (the THM-534 mechanism)
      (3) test whether each component S_r(E) is itself extremized by consec (a cleaner target):
          consec MIN S_1 (=> max avg #sectors-hit), consec MAX S_2,S_4 — over a thorough bank.
      (4) the SCALE-INVARIANT reframe: S_r are scale-invariant; is there a majorization
          S_r(E) <= / >= S_r(consec) provable by an AP-orbit rearrangement?
"""
import sys, itertools
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)

def avoid_measure(E, A):
    """meas{x: orbit {frac(e x)} avoids the union of sectors in A (A subset {1..6})}.
       = meas{x: no e_i x mod 1 in any [j/7,(j+1)/7), j in A}."""
    E=sorted(set(E))
    forb=[(F(j,7),F(j+1,7)) for j in A]
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for (lo,hi) in forb:
            for t in (lo,hi):
                for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if 0<=z<=1); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; ok=True
        for e in E:
            p=(e*xm)%1
            for (lo,hi) in forb:
                if lo<=p<hi: ok=False; break
            if not ok: break
        if ok: tot+=x1-x0
    return tot

def S_r(E, r):
    return sum(avoid_measure(E,A) for A in itertools.combinations(range(1,7),r)) if r>0 else F(1)

# integer-root duals g(t): coefficients y_r in binomial basis g(t)=sum_r y_r C(t,r)
def g_func(k):
    if k==8: return lambda t:F((t-1)*(t-2)*(t-4)*(t-5),40)
    if k in (9,10): return lambda t:F(-(t-2)*(t-3)*(t-6),36)
    return lambda t:F((t-3)*(t-4),12)  # k=11,12,13
def Ly(E,k):
    g=g_func(k);
    # L_y = E[g(N)] = sum over x-cells; easier: g(t)=sum_r y_r C(t,r), and E[C(N,r)]=S_r.
    # invert: y_r from finite differences of g at 0..6
    yr=[]
    for r in range(7):
        s=F(0)
        for i in range(r+1): s+=F((-1)**(r-i)*comb(r,i))*g(i)
        yr.append(s/1)  # y_r = r-th forward difference of g at 0 / ... actually Newton: g(t)=sum binom(t,r)*Delta^r g(0)
    # g(t)=sum_r Delta^r g(0) C(t,r); Delta^r g(0)=sum_i (-1)^{r-i}C(r,i)g(i)
    val=F(0)
    for r in range(7):
        if yr[r]!=0: val+=yr[r]*S_r(E,r)
    return val, yr

H=F(1,14)
def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H/u)%1; b=(c+H/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mgm(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    dz=mgm([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
import functools
@functools.lru_cache(None)
def cap(k): return min(measGP(list(P)) for P in itertools.combinations(range(1,14),13-k))

print("(1,2) L_y closure + dual coefficients + component moments (consec), k=8,9,10")
for k in (8,9,10):
    E=list(range(k)); ly,yr=Ly(E,k); capk=cap(k)
    srs=[S_r(E,r) for r in range(5)]
    print(f"  k={k}: L_y(consec)={ly}={float(ly):.5f}  cap_k={float(capk):.5f}  {'CLOSES' if ly<=capk else 'NO'}")
    print(f"        dual y_r (binomial)={[str(y) for y in yr if y!=0]}")
    print(f"        S_r(consec) r=0..4: {[float(s) for s in srs]}")

print("\n(3) component extremality over a thorough bank (k=8, spread<=12): "
      "consec MIN S_1? MAX S_2,S_4?")
k=8; cE=list(range(k))
s1c,s2c,s4c=S_r(cE,1),S_r(cE,2),S_r(cE,4)
v1=v2=v4=0; checked=0
for rest in itertools.combinations(range(1,13),k-1):
    E=[0]+list(rest); checked+=1
    if S_r(E,1)<s1c-F(1,10**12): v1+=1   # someone has SMALLER S_1 than consec
    if S_r(E,2)>s2c+F(1,10**12): v2+=1   # someone has LARGER S_2
    if S_r(E,4)>s4c+F(1,10**12): v4+=1
print(f"  consec S_1={float(s1c):.4f} S_2={float(s2c):.4f} S_4={float(s4c):.4f}; checked {checked}")
print(f"  shapes with S_1<consec (consec NOT min): {v1}")
print(f"  shapes with S_2>consec (consec NOT max): {v2}")
print(f"  shapes with S_4>consec (consec NOT max): {v4}")
print("  => if v1=v2=v4=0, consec extremizes each component, so L_y(E)<=L_y(consec) follows")
print("     from the dual signs (y_1<0, y_2,y_4>0). That would PROVE HYP-2607 for k=8.")
print("\nDONE.")
