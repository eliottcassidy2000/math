import numpy as np
from fractions import Fraction as F
from math import comb
GRID=400000; xs=(np.arange(GRID)+0.5)/GRID
def gaps_matrix(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    return np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
def mu_17(E):
    return float((gaps_matrix(E).max(axis=1)>1/7).mean())
def EU(E):  # E_x[ sum_gaps (gap-1/7)_+ ] = E_x[uncovered length];  mu_17 >= EU (Paley-Zygmund, U in [0,1])
    g=gaps_matrix(E); return float(np.clip(g-1/7,0,None).sum(axis=1).mean())

print("=== E[U] = E_x[uncovered length] as a PROVABLE lower bound on mu_17 (mu_17 >= E[U]), opus-S131 ===\n")
print(f"{'k':>3} {'E':>22} {'mu_17':>8} {'E[U]':>8}  (mu_17>=E[U])")
def show(k,name,E):
    print(f"{k:>3} {name:>22} {mu_17(E):>8.4f} {EU(E):>8.4f}")
for k in range(8,14):
    show(k,"AP {0..k-1}",list(range(k)))
# is E[U] minimized at AP too? check structured + random for k=13
print()
import random; random.seed(5)
ap=list(range(13)); euap=EU(ap); minEU=euap; worst=("AP",ap)
show(13,"AP",ap)
for name,E in [("12AP+out20",list(range(12))+[20]),("two6blocks",list(range(6))+list(range(50,57))),
               ("hole@6+13",[i for i in range(13) if i!=6]+[13]),("step3",list(range(0,39,3)))]:
    if len(set(E))==13:
        e=EU(E); 
        if e<minEU: minEU=e; worst=(name,E)
        show(13,name,E)
for _ in range(40):
    E=sorted(random.sample(range(1,80),13)); e=EU(E)
    if e<minEU: minEU=e; worst=(name,E)
print(f"\nk=13: E[U](AP)={euap:.4f}; min E[U] over structured+40 random = {minEU:.4f} (at {worst[0]})")
print(f"  => if inf_E E[U] > 0, then mu_17 >= E[U] gives a UNIFORM density floor (positivity).")

# EXACT E[U] via inclusion-exclusion first two terms (structure-INDEPENDENT) to see the analytic handle
print("\n=== analytic: E[U] = 1 - E[covered]; I-E = k/7 - C(k,2)*3/196 + (triples...) ===")
for k in range(8,14):
    ie2 = F(k,7) - comb(k,2)*F(3,196)   # covered ~ first two I-E terms (upper approx to E[U]=1-covered from below)
    print(f"  k={k}: 1 - [k/7 - C(k,2)*3/196] = {float(1-ie2):.4f}  (2-term I-E estimate of E[U], triples+ correct it)")
