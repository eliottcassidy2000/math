# opus-2026-07-20-S401 -- HYP-8220: TESTING kps's BORSUK-ULAM PREDICTION.
#
# kind-pasteur-S31av: 14 = |D_7|, the heptagon's dihedral group.  Because
# p = 7 = 3 mod 4, -1 is a NON-residue, so the complement/reflection is an
# ANTI-automorphism -- orientation reversing -- and the Z/2 acts FREELY.  A free
# Z/2 is the hypothesis of BORSUK-ULAM, not Brouwer, and the prediction is that
# witnesses carry no reflection-fixed representative but come in ANTIPODAL PAIRS
# (t*, -t*), with the certificate an odd-degree/parity invariant.
#
# THIS IS TESTABLE against my maximizer data, and it sits on the SURVIVING side
# of the THM-1185 triage (measure methods are blind to the extremal families;
# pointwise/topological ones are not).  It also bears on my THM-1200, where I
# called the tournament-7 and the LRC-7 "structurally unrelated" -- if the
# maximizer set is a free Z/2-set, that claim needs scoping.
from fractions import Fraction as F
from itertools import combinations
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
def gap_at(V,t): return min(fn(t*v) for v in V)
def maximizers(V):
    """ALL critical points attaining the max, exactly."""
    D=set()
    for a,b in combinations(V,2):
        D.add(a+b)
        if a!=b: D.add(abs(a-b))
    for v in V: D.add(2*v)
    best=F(0); args=[]
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=gap_at(V,F(p,q))
            if g>best: best,args=g,[F(p,q)]
            elif g==best and F(p,q) not in args: args.append(F(p,q))
    return best, sorted(set(args))
print("(1) ARE THE MAXIMIZERS A FREE Z/2-SET UNDER t -> -t (i.e. t -> 1-t)?")
print("    Borsuk-Ulam predicts: NO fixed point, and an even count in antipodal pairs.")
print()
for nm,V in [("{1,...,13}  (extremal)",list(range(1,14))),
             ("{1..11,13,24} (extremal)",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("{1..11,13,36} (3/41)",[1,2,3,4,5,6,7,8,9,10,11,13,36]),
             ("{1,...,12,14} (1/13)",list(range(1,13))+[14]),
             ("{2,...,14}",list(range(2,15)))]:
    M,args=maximizers(V)
    anti=[1-t for t in args]
    closed = sorted(set(anti))==args
    fixed=[t for t in args if 1-t==t]     # t = 1/2 only
    pairs=len(args)//2
    print(f"  {nm:26s} M={str(M):8s}  #maximizers={len(args)}")
    print(f"      maximizers: {[str(t) for t in args]}")
    print(f"      closed under t -> 1-t: {closed}   fixed points: {fixed}   "
          f"=> {pairs} antipodal pair(s)")
print()
print("(2) THE D_7 IRREP COUNT: 14 = 1^2 + 1^2 + 2^2 + 2^2 + 2^2 -- THREE 2-dim irreps.")
print("    kps predicts the free-reflection witnesses organise by those three.")
V=list(range(1,14)); M,args=maximizers(V)
print(f"    {{1,...,13}}: {len(args)} maximizers = {len(args)//2} antipodal pairs"
      f"   {'== 3, matching the three 2-dim irreps' if len(args)//2==3 else ''}")
print(f"    numerators over 14: {[t.numerator for t in args]}")
print(f"    these are the units mod 14? {sorted(t.numerator for t in args)} vs "
      f"{[u for u in range(1,14) if __import__('math').gcd(u,14)==1]}")
