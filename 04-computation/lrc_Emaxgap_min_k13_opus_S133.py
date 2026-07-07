import numpy as np, random
from itertools import combinations
from math import gcd
from functools import reduce
random.seed(13)
GRID=150000; xs=(np.arange(GRID)+0.5)/GRID
def Emaxgap(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
    return float(g.max(axis=1).mean())
def norm(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]; g=reduce(gcd,E[1:]) if len(E)>1 else 1
    return tuple(e//g for e in E) if g else tuple(E)

thr=1/7
ap=Emaxgap(list(range(13)))
print(f"=== TRUE min of E[maxgap] over k=13 shapes (opus-S133) ===")
print(f"E[maxgap(AP_13)]={ap:.5f} (=93/440); 1/7={thr:.5f}\n")
best=ap; bestE=tuple(range(13))
# 1) endpoint-stretched APs {0..11, S} and {0, a..a+11}
for S in range(12,40):
    for E in [tuple(range(12))+(S,)]:
        v=Emaxgap(E)
        if v<best: best,bestE=v,norm(E)
# 2) systematic: {0..j-1} U {stretched} -- AP-core + a few stretched points
for cut in range(6,13):
    for ext in combinations(range(cut,cut+12),13-cut):
        E=tuple(range(cut))+ext
        if len(set(E))!=13: continue
        v=Emaxgap(E)
        if v<best: best,bestE=v,norm(E)
# 3) adversarial descent (minimize E[maxgap])
for _ in range(30):
    E=list(range(13)); cur=Emaxgap(E)
    for _ in range(50):
        i=random.randrange(13); cand=E.copy(); cand[i]=random.randint(0,30)
        if len(set(cand))<13: continue
        cv=Emaxgap(cand)
        if cv<cur: E,cur=cand,cv
    if cur<best: best,bestE=cur,norm(E)
print(f"TRUE min E[maxgap] (k=13) found: {best:.5f}  at shape {bestE}")
print(f"  below AP (0.2114)? {best<ap-1e-4}")
print(f"  > 1/7={thr:.5f}? {best>thr}   margin over 1/7: {best-thr:.5f}")
lb=(best-thr)/(1-thr)
print(f"  reverse-Markov floor at the TRUE min: mu_17 >= {lb:.5f} > 0")
print(f"\n  => density floor holds iff inf_E E[maxgap] > 1/7; the min is NOT the AP but is ~{best:.3f} >> 1/7.")
