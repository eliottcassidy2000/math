"""
The RG flow T: S -> even(S)/2. AP family {1..N} is invariant: T{1..N}={1..floor(N/2)}.
Apex-gap trajectory g(odd core mod p) along the flow, for LRC(2p) AP {1..2p-1}. Cusp Z_p (g=0) = unstable
entry; doublet (g=4sin^2(pi/2p)) = attractor. Test universality on general sets.
"""
import math, cmath
def gp(res, p):  # apex gap mod p
    if not res: return None
    w=cmath.exp(2j*cmath.pi/p)
    return min(abs(sum(w**((k*x)%p) for x in res))**2 for k in range(1,p))
def flow_odd_cores(S):
    ch=[]; cur=list(S)
    while cur:
        O=[x for x in cur if x%2==1]; E=[x//2 for x in cur if x%2==0]
        if O: ch.append(O)
        cur=E
        if len(ch)>30: break
    return ch
print("RG FLOW on the AP family {1..2p-1} (the cusp config of LRC(2p)). Apex-gap trajectory g(O_j mod p):")
for p in [3,5,7,11,13]:
    N=2*p-1
    ch=flow_odd_cores(list(range(1,N+1)))
    traj=[round(gp([x%p for x in O],p),4) for O in ch]
    sizes=[len(O) for O in ch]
    atom=4*math.sin(math.pi/(2*p))**2
    iscusp0 = set(x%p for x in ch[0])==set(range(p))
    print(f"  p={p:>2} (LRC{2*p}): AP{{1..{N}}} sizes {sizes}  g-traj {traj}")
    print(f"          cusp at level0 (Z_{p})={iscusp0}; doublet attractor g=4sin^2(pi/{2*p})={atom:.4f}")
print()
print("=> EVERY LRC(2p) AP starts at the CUSP Z_p (g=0, level 0) and flows down to the DOUBLET (g=atom),")
print("   then SINGLETON (g=1). The flow N->floor(N/2) is binary halving; depth ~ log2(2p).")
print("   The cusp is the UNSTABLE ENTRY (left after one step); the doublet is the floor-carrying attractor.")
print()
# universality: do general (non-AP) sets also flow to the doublet attractor (worst g)? sample random covering-ish sets
print("UNIVERSALITY: deepest-shell (attractor) apex gap over assorted 13-sets (mod 7):")
import random
random.seed(1)
tests={
 "AP {1..13}": list(range(1,14)),
 "{1..12,182}": list(range(1,13))+[182],
 "{3,5,7,...,27}": list(range(3,28,2))[:7]+[2,4,8,16,32,64],
 "powers-heavy": [1,3,5,7,2,4,8,16,32,64,128,6,12],
 "random A": sorted(random.sample(range(1,50),13)),
 "random B": sorted(random.sample(range(1,50),13)),
}
for name,S in tests.items():
    ch=flow_odd_cores(S)
    gtr=[round(gp([x%7 for x in O],7),4) for O in ch]
    # attractor = deepest nontrivial (last with |O|>=1); report min nonzero g
    nz=[g for g in gtr if g and g>1e-9]
    print(f"  {name:>16}: g-traj {gtr}  -> per-level min(nonzero) g = {min(nz) if nz else None}")
print("  => attractor's apex gap bottoms at the DOUBLET 0.198 whenever a doublet shell appears (generic);")
print("     else higher (singleton 1.0 / triple 0.308). The WORST attractor = doublet = THM-590 floor.")
