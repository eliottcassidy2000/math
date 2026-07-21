from fractions import Fraction as F
from itertools import combinations
print("WOWII-103 MOTIF = INFLATION/DECOUPLING: a small CORE + pendant LEAVES that pump ONE")
print("invariant (independence 9) while a coupled invariant (ecc, bipartite) stays fixed,")
print("breaking a conjectured inequality alpha <= 8. Found by exhaustive 2^11 search + Lean.")
print()
print("THE REPO ALREADY USES THIS MOTIF. Check: the LRC extremal {1..11,13,24} is a")
print("LEAF-INFLATION of the base family {1..13} -- swap the 'core' speed 12 for a SACRIFICIAL")
print("pendant 24 = 2*12 that reproduces the same danger comb at t*=1/14 but is a 'leaf'.")
def fn(x):
    r=x-int(x); r+= 1 if r<0 else 0
    return min(r,1-r)
def M(V):
    D=set()
    for a,b in combinations(V,2): D.add(a+b); D.add(abs(a-b))
    for v in V: D.add(2*v)
    best=F(0)
    for q in sorted(D):
        if q<2: continue
        for p in range(1,q):
            g=min(fn(F(p,q)*v) for v in V)
            if g>best: best=g
    return best
base=list(range(1,14))
infl=[1,2,3,4,5,6,7,8,9,10,11,13,24]
gap =[1,2,3,4,5,6,7,8,9,10,11,13,36]
for nm,V in [("{1..13} base",base),("{1..11,13,24} inflated",infl),("{1..11,13,36} gap",gap)]:
    m=M(V); print(f"  {nm:24s} M={m}   ({'=1/14 tight' if m==F(1,14) else str(float(m))})")
print("  => {1..11,13,24}: replace 12 by 24 (a '2x leaf'); SAME M=1/14 at t*=1/14 because")
print("     24*(1/14)=12/7, ||12/7||=||-2/7||... the leaf reproduces the core's danger. This is")
print("     EXACTLY the WOWII pendant-inflation: a leaf that pumps the config without changing")
print("     the binding invariant. (The (1/14,3/41) gap family {1..11,13,36} is another leaf.)")
print()
print("THE DECOUPLING VIEW (THM-1820 is the repo's own WOWII). WOWII-103 breaks a conjectured")
print("COUPLING of independence to eccentricity; THM-1820 broke the conjectured coupling of")
print("H-maximal to 3-cycle-maximal (they are DECOUPLED, opposite strata). SAME phenomenon:")
print("a conjectured inequality between two invariants dies to a construction that decouples")
print("them. So the transferable rule: to break an invariant-inequality, DECOUPLE via inflation.")
print()
print("TRANSFERABLE METHOD (WOWII -> repo):")
print("  1. INFLATION: small core + sacrificial 'leaf' (LRC: 2x speed; JC: append identity")
print("     coordinate; tournament: dominated/dominating vertex) pumps one invariant.")
print("  2. EXHAUSTIVE small-case certification (repo already does; e.g. 59048 nullcone sweep).")
print("  3. LEAN certification of the finite invariants (repo has TournamentH7) -- the WOWII PR")
print("     Lean-certified alpha=9,b=10,ecc=30/11 by decide+native. The repo could Lean-certify")
print("     its extremal LRC/tournament invariants the same way.")
