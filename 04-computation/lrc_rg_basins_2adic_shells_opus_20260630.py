"""
TASK 2 - basins of the RG flow over 13-sets. The 2-adic shell profile (sizes by valuation) determines the
attractor (deepest shell) and the floor (min apex gap). Classify: CUSP (odd part covers Z_7), DOUBLET-floor
(a distinct-residue 2-shell, g=0.198), MULTIPLICITY (g<0.198), SINGLETON/other (g>0.198, no doublet).
"""
import cmath, math, random
from collections import Counter
def gmult(res,p=7):
    if not res: return None
    w=cmath.exp(2j*cmath.pi/p)
    return min(abs(sum(w**((k*x)%p) for x in res))**2 for k in range(1,p))
def shells(S):
    sh={}
    for s in S:
        v=0; u=s
        while u%2==0: u//=2; v+=1
        sh.setdefault(v,[]).append(u)
    return sh  # valuation -> odd parts
def classify(S):
    sh=shells(S)
    levels=sorted(sh)
    gaps=[gmult([u%7 for u in sh[v]]) for v in levels]
    cusp = any(set(u%7 for u in sh[v])==set(range(7)) for v in levels)
    floor=min(g for g in gaps if g is not None and g>1e-9) if any(g and g>1e-9 for g in gaps) else None
    # multiplicity: some shell has repeated residue mod 7
    mult = any(len(set(u%7 for u in sh[v]))<len(sh[v]) for v in levels)
    deepest_size=len(sh[max(levels)])
    if cusp: cat="CUSP(meas0)"
    elif mult and floor is not None and floor<0.198-1e-3: cat="MULTIPLICITY(<0.198)"
    elif floor is not None and abs(floor-0.1981)<1e-3: cat="DOUBLET(0.198)"
    else: cat="OTHER(>0.198)"
    return cat, round(floor,4) if floor else None, deepest_size, [len(sh[v]) for v in levels]
random.seed(7)
for label,(lo,hi) in [("dense {1..20}",(1,20)),("mid {1..60}",(1,60)),("sparse {1..300}",(1,300))]:
    cats=Counter(); floors=Counter(); samples=2000
    for _ in range(samples):
        S=random.sample(range(lo,hi+1),13)
        c,fl,ds,prof=classify(S)
        cats[c]+=1; floors[fl]+=1
    print(f"\n{label}: basin distribution over {samples} random 13-sets:")
    for c,n in cats.most_common():
        print(f"   {c:>22}: {100*n/samples:5.1f}%")
print()
print("DETERMINANT: the basin = the 2-adic shell profile + residue structure (a STATIC decomposition).")
print("  * CUSP needs 7 odd elements covering Z_7 (rare unless dense & AP-like).")
print("  * DOUBLET-floor needs a 2-element shell with distinct residues mod 7 (common in dense sets).")
print("  * MULTIPLICITY needs a shell with repeated residue mod 7 (common in sparse sets: collisions).")
print("  * deepest shell = the attractor (singleton g=1 generic; doublet 0.198 when 2 max-valuation elts).")
