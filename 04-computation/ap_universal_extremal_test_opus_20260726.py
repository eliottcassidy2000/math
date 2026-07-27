import numpy as np
from itertools import combinations
from math import gcd
from functools import reduce
# TEST boxeph-S205's recurrence: is the AP the SIMULTANEOUS extremal across independent functionals?
# (1) additive-triple count E3 (Schur/LRC THM-730): #{(a,b,c) in S^3 : a+b=c}  -- AP should MAXIMIZE
# (2) mean max-gap of the multiplicative orbit {frac(s x)} (LRC reach)  -- AP should be near-extremal
# If AP is extremal for BOTH (independent functionals), that's the 'universal extremal' recurrence.
def e3(S):
    Sset=set(S); return sum(1 for a in S for b in S if a+b in Sset)
def norm(S):
    S=sorted(set(S)); S=[x-S[0] for x in S]; g=reduce(gcd,S[1:]) if len(S)>1 else 1
    return tuple(x//g for x in S) if g else tuple(S)
GRID=40000; xs=(np.arange(GRID)+0.5)/GRID
def Emaxgap(S):
    S=[s for s in S if s>0] or [1]
    ph=np.mod(np.outer(xs,np.array(S,float)),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
    return float(g.max(axis=1).mean())

print("=== boxeph-S205 recurrence TEST: is the AP the simultaneous extremal? (opus) ===")
for n in [6,7,8]:
    AP=tuple(range(1,n+1))
    e3_AP=e3(AP)
    # census all n-subsets of [1..n+6], normalized, find E3-maximizers
    seen={}
    for S in combinations(range(1,n+7),n):
        sh=norm(S)
        if len(sh)==n and sh not in seen: seen[sh]=(e3(list(sh)), Emaxgap(list(sh)))
    e3max=max(v[0] for v in seen.values())
    e3maxers=[s for s,v in seen.items() if v[0]==e3max]
    apsh=norm(AP)
    print(f"\n n={n}: AP={AP} (norm {apsh}) E3(AP)={e3_AP}")
    print(f"   E3 max over {len(seen)} shapes = {e3max}; #maximizers={len(e3maxers)}; AP among them? {apsh in e3maxers}")
    print(f"   E3-maximizers (sample): {e3maxers[:3]}")
    # is AP also the Emaxgap-min among the E3-maximizers?
    emg={s:seen[s][1] for s in e3maxers}
    print(f"   among E3-maximizers, Emaxgap: AP={seen[apsh][1]:.4f}, min={min(emg.values()):.4f} => AP min-maxgap among them? {abs(seen[apsh][1]-min(emg.values()))<1e-3}")
print("\n  => AP maximizing E3 (Schur) is THM-730; testing if it ALSO co-optimizes the reach functional (independent).")
