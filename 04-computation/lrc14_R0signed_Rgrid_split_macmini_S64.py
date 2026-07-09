"""
mac-mini-2026-07-09-S64 (handoff) -- the R0-signed / R_grid-absolute split of the Weyl residual.

kps-S97: E_grid[W] = (6/7)^k + R, R = sum_{n!=0, V|n.e} What(n) (Poisson tail over L_V(e)).
Split L_V(e) = {n.e = 0} (exact relations, grid-INDEPENDENT) U {n.e = mV, m!=0} (wraparound shells):
  R0     = sum_{n.e=0,   n!=0} What(n)   -- SIGNED, keep whole:  E_x[W] := (6/7)^k + R0 = continuum mean int_0^1 W(x)dx  (>0, fixed).
  R_grid = sum_{n.e=mV, m!=0} What(n)     -- ABSOLUTE, bound:      recedes as V grows (needs |n| >~ V/spread, 0.371-decay).
Then  E_grid[W] = E_x[W] + R_grid, and existence (E_grid[W]>0)  <==  E_x[W] > |R_grid|.
The winning side is the CONTINUUM floor (density-floor object); only the decaying wraparound is bounded in abs.

TEST for 7-structured dissociated k=13 sets across V (7|V; wide spread>6V/7 AND compressed):
 (Q1) sign/size of R0 (dissociated => few exact relations => R0 ~ 0 => E_x ~ (6/7)^k).
 (Q2) is |R_grid| < |R| (split strictly better-conditioned than kps's full-R route)?
 (Q3) does E_x[W] > |R_grid| hold, with what margin, in the WIDE regime (the handoff item)?
 (Q4) the crossover: knife-edge (spread=6V/7) should give E_x = |R_grid| EXACTLY (E_grid=0). Confirm.
E_x[W] computed by a FINE coprime grid (V_fine=100003, prime, gcd(.,7)=1) where R_grid(V_fine)~0.
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(640064)
K=13; lead=(6.0/7.0)**K
Vf=100003   # fine grid ~ continuum (prime, coprime to 7)

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
def meanW(E,V):
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return np.maximum(g-1.0/7,0).sum(axis=1).mean()
def make_7struct(k,s,mn):
    sev=[x for x in range(0,s+1,7)]
    if len(sev)<mn: return None
    S7=sorted(random.sample(sev,random.randint(mn,min(len(sev),k-1))))
    pool=[x for x in range(1,s) if x%7!=0]; need=k-len(S7)
    if need<0 or len(pool)<need: return None
    E=sorted(set(S7+[0,s]+random.sample(pool,max(0,need))))
    return E if len(E)==k else None

print(f"K={K}  (6/7)^K={lead:.6f}   E_x[W] via fine grid Vf={Vf}\n")
print(f"{'regime':>12}{'V':>5}{'E_x[W]':>9}{'R0/lead':>9}{'|R|/lead':>10}{'|Rg|/lead':>11}{'|Rg|/|R|':>10}{'E_x>|Rg|?':>10}{'margin':>9}")
# collect sup |R_grid|/E_x per regime (the split's actual existence ratio; <1 => E_grid>0)
sup_ratio={'compressed':0.0,'wide':0.0}; fails={'compressed':0,'wide':0}; ntot={'compressed':0,'wide':0}
examples=[]
for _ in range(4000):
    s=random.randint(40,120); E=make_7struct(K,s,4)
    if E is None or not prim(E) or longest_ap(E)>K-6: continue
    mx=max(E)
    Ex=meanW(E,Vf); R0=Ex-lead
    for t in range(1,8):
        V=7*((mx)//7+t)
        if V<=mx: continue
        Eg=meanW(E,V); R=Eg-lead; Rg=Eg-Ex
        reg='wide' if 7*mx>6*V else ('compressed' if 7*mx<6*V else 'knife')
        if reg=='knife':
            examples.append(('KNIFE',V,Ex,R0/lead,abs(R)/lead,abs(Rg)/lead,abs(Rg)/max(abs(R),1e-12),Ex>abs(Rg),Ex-abs(Rg)))
            continue
        ratio=abs(Rg)/Ex; ntot[reg]+=1
        if Ex<=abs(Rg): fails[reg]+=1
        if ratio>sup_ratio[reg]: sup_ratio[reg]=ratio
        if len(examples)<8:
            examples.append((reg,V,Ex,R0/lead,abs(R)/lead,abs(Rg)/lead,abs(Rg)/max(abs(R),1e-12),Ex>abs(Rg),Ex-abs(Rg)))
for e in examples[:12]:
    print(f"{e[0]:>12}{e[1]:>5}{e[2]:>9.5f}{e[3]:>9.3f}{e[4]:>10.4f}{e[5]:>11.4f}{e[6]:>10.3f}{str(e[7]):>10}{e[8]:>9.5f}")
print()
for reg in ('compressed','wide'):
    if ntot[reg]:
        print(f"{reg:>12}: n={ntot[reg]:>5}  sup |R_grid|/E_x = {sup_ratio[reg]:.4f}  (<1 => E_grid[W]>0 => strict good period);  E_x<=|R_grid| FAILURES: {fails[reg]}")
print()
print("READING:")
print(" R0/lead ~ 0 (dissociated => few exact relations => E_x ~ (6/7)^k, a CLEAN explicit floor).")
print(" |R_grid|/|R| < 1 => the split strips the (grid-independent) continuum part off the residual.")
print(" sup|R_grid|/E_x < 1 in a regime => that regime has a STRICT good period via the split (E_grid>0).")
print(" KNIFE (spread=6V/7): E_x = |R_grid| exactly (E_grid=0) => the split's equality case = j=1's non-strict job.")
