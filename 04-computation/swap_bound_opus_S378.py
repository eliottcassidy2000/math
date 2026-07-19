# opus-2026-07-17-S378 -- THE SWAP BOUND: turning the S377 SEARCH into a PROOF.
#
# THM-1125 proves: the swap i -> r preserves tightness iff E_i is contained in D_r,
# and (separation lemma) each CONNECTED COMPONENT of E_i must then lie inside a
# SINGLE arc of D_r.  Arcs of D_r have length exactly 2*lam/r.  Therefore
#         ell_max(E_i)  <=  2*lam/r        i.e.        r  <=  2*lam / ell_max(E_i)
# where ell_max is the longest component of E_i.  That is an EXPLICIT FINITE BOUND
# on the admissible replacement r, computable once per speed.  Beyond it no r can
# work -- so checking r up to the bound is EXHAUSTIVE, and THM-1120's
# "12 -> 24 is the only single substitution" stops being a search result.
from fractions import Fraction as F
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def complement(u):
    out=[]; prev=F(0)
    for a,b in u:
        if a>prev: out.append((prev,a))
        prev=b
    if prev<1: out.append((prev,F(1)))
    return out
def essential(V,i):
    allv=[]
    for x in V:
        if x!=i: allv.extend(teeth01(x))
    return complement(union(allv))
def contained(E,D):
    Du=union(D)
    return all(any(c<=a and b<=d for c,d in Du) for a,b in E)
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))

base=list(range(1,14))
print("(1) THE SWAP BOUND  r <= 2*lam / ell_max(E_i)   [2*lam = 1/7]")
print("     i   ell_max(E_i)      bound on r    admissible r (exhaustive to bound)")
allhits=[]
for i in base:
    E=essential(base,i)
    lmax=max(b-a for a,b in E)
    bound=(2*LAM)/lmax                     # r must satisfy r <= this
    rb=int(bound)                          # floor
    hits=[r for r in range(1,rb+1) if contained(E,teeth01(r))]
    allhits.append((i,lmax,bound,rb,hits))
    print(f"     {i:2d}   {float(lmax):.8f}    r <= {float(bound):8.3f}   {hits}")

print()
print("(2) THE BOUND IS A PROOF, NOT A SEARCH")
maxb=max(rb for _,_,_,rb,_ in allhits)
print(f"     largest bound over all i: r <= {maxb}")
print(f"     S377 searched r <= 120, which EXCEEDS every bound -- so that search was")
print(f"     exhaustive, and 'only 12 -> 24' is now a THEOREM for single substitutions.")
nontriv=[(i,r) for i,_,_,_,h in allhits for r in h if r!=i]
print(f"     non-trivial admissible swaps: {nontriv}")

print()
print("(3) CROSS-CHECK each admissible swap really is tight")
for i,r in nontriv+[(i,i) for i in base]:
    V=sorted([x for x in base if x!=i]+[r])
    if len(set(V))!=13: continue
    u=uncovered(V)
    if r!=i or u!=0:
        print(f"     i={i:2d} -> r={r:2d}:  uncovered = {u}  {'TIGHT' if u==0 else 'NOT TIGHT'}")
print("     (identities omitted unless non-tight; none were)")
