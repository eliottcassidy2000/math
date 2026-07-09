# Is the ABSOLUTE resonant sum Sum_{Vmax|n.e}|What(n)| < (6/7)^k for DISSOCIATED clusters? (kps-S98)
# If yes => |R| < (6/7)^k a-priori with NO cancellation (triangle ineq) => good-period existence,
# closing the dissociated branch cleanly. Uses FULL LEM-011 |What| (with sigma-factor). Truncate at
# height H, watch convergence in H (geometric tail 0.371/coord). Critical window V=floor(7s/6), NOT V=k.
import numpy as np
from math import floor, gcd, sin, pi
from functools import reduce
from itertools import combinations, product
import random
random.seed(98017)
K=13; lead=(6.0/7.0)**K
def babs(m):
    if m%7==0: return 0.0
    return abs(sin(pi*m/7))/(pi*abs(m))
def What_abs(nvec):
    act=[m for m in nvec if m!=0]; r=len(act)
    if any(m%7==0 for m in act): return 0.0
    sig=sum(nvec); pb=1.0
    for m in act: pb*=babs(m)
    if sig==0: sf=6.0/7.0
    else:
        if sig%7==0: return 0.0
        sf=babs(sig)
    return (6.0/7.0)**((K-1)-r)*pb*sf
def prim(E):
    E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
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
def abs_resonant_sum(E,V,H,supmax=3):
    E=sorted(E); coords=list(range(1,K)); rng=[m for m in range(-H,H+1) if m!=0 and m%7!=0]; tot=0.0
    for s in range(1,supmax+1):
        for combo in combinations(coords,s):
            for vals in product(rng,repeat=s):
                ne=sum(vals[t]*E[combo[t]] for t in range(s))
                if ne%V!=0: continue
                nvec=[0]*(K-1)
                for t,ci in enumerate(combo): nvec[ci-1]=vals[t]
                tot+=What_abs(nvec)
    return tot
print("ABSOLUTE resonant sum Sum_{V|n.e}|What(n)| / (6/7)^k for DISSOCIATED clusters (longest-AP<=6).")
print("If < 1 => |R| < (6/7)^k a-priori, NO cancellation. Convergence in height H (support<=3).")
print(f"{'family':>12}{'s':>5}{'V':>5}{'lAP':>4}{'H=3':>8}{'H=5':>8}{'H=7':>8}")
def emit(lab,E,s):
    E=sorted(set(E)); V=floor(7*s/6)
    if len(E)!=K: return
    b3=abs_resonant_sum(E,V,3)/lead; b5=abs_resonant_sum(E,V,5)/lead; b7=abs_resonant_sum(E,V,7)/lead
    print(f"{lab:>12}{s:>5}{V:>5}{longest_ap(E):>4}{b3:>8.4f}{b5:>8.4f}{b7:>8.4f}")
    return max(b3,b5,b7)
# hard 7-structured (mac-mini MISTAKE-128) + generic dissociated at several spreads
emit("7-struct",[0,7,14,21,26,29,37,44,51,58,67,75,82],84)
worst=0
for s in (60,90,120):
    for _ in range(5000):
        mid=sorted(random.sample(range(1,s),K-2)); E=[0]+mid+[s]
        if len(set(E))==K and prim(E) and longest_ap(E)<=6:
            v=emit("dissoc",E,s)
            if v and v>worst: worst=v
            break
print(f"\nworst dissociated abs-resonant-sum/lead (H<=7) = {worst:.4f}")
print("=> if <1 and converging in H, the dissociated |R|<(6/7)^k is A-PRIORI (absolute, no cancellation).")
print("   the tail beyond H is <= (0.371)^? geometric, addable explicitly.")
