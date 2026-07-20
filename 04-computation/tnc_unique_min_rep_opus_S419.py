import sympy as sp
from itertools import combinations_with_replacement as cwr
u=sp.symbols('u')
def CT(R,N,mm):
    return sp.Poly(sp.expand(sp.sympify(R)**mm),u).coeff_monomial(u**(N*mm))
def charges(R,N):
    P=sp.Poly(sp.sympify(R),u)
    return sorted(k-N for (k,),c in zip(P.monoms(),P.coeffs()) if c!=0)
def min_rep_m(chs):
    """smallest m>=1 such that 0 is a sum of m charges; and whether the rep is UNIQUE"""
    from itertools import product
    for m in range(1,25):
        reps=[c for c in cwr(chs,m) if sum(c)==0]
        if reps: return m,len(reps)
    return None,None
print("SUFFICIENT CONDITION: if 0 has a UNIQUE minimal representation as a sum of charges,")
print("then CT(m0) is a single nonzero product => Lambda is NOT a nullcone element => TNC.")
print()
print("   R                       N   charges         m0   #minimal-reps   CT(m0)   TNC-forced?")
for R,N in [("2+3*u**2",1),("2+3*u**4",2),("1+u**3+u**4",2),("1+u+u**4",2),
            ("1+u**3+u**5",2),("1+u**4+u**5",3),("1+u**2+u**5",3),
            ("1+u+u**2+u**3+u**4",2),("1+u**2+u**3+u**5",3)]:
    chs=charges(R,N); m0,nrep=min_rep_m(chs); ctv=CT(R,N,m0) if m0 else None
    forced = (nrep==1 and ctv!=0)
    print(f"   {R:24s} {N}   {str(chs):16s} {str(m0):4s}   {str(nrep):8s}     {str(ctv):6s}   {forced}")
print()
print("  => whenever #minimal-reps = 1, CT(m0) != 0 and TNC is FORCED (no all-order argument")
print("     needed).  This closes: ALL BINOMIALS, and every R with a unique minimal 0-rep.")
print()
print("THE RESIDUAL is exactly #minimal-reps >= 2 with the reps CANCELLING.  Search for a")
print("trinomial where the minimal representation is NON-UNIQUE (so cancellation is POSSIBLE):")
found=[]
for N in [2,3]:
    for j in range(1,7):
        for d in range(j+1,9):
            if j==N or d==N: continue      # keep r_N=0
            R=1+u**j+u**d
            chs=charges(R,N); m0,nrep=min_rep_m(chs)
            if nrep and nrep>=2:
                found.append((str(R),N,chs,m0,nrep,CT(f"1+u**{j}+u**{d}",N,m0)))
for R,N,chs,m0,nrep,ct in found[:8]:
    print(f"   {R:16s} N={N}: charges {chs}, m0={m0}, #reps={nrep}, CT(m0)={ct}")
print(f"   trinomials with non-unique minimal rep found: {len(found)}")
if found:
    print("   -> these are the ONLY place cancellation can even start; with unit coeffs CT!=0,")
    print("      but TUNED coefficients could zero CT(m0). That tuning = the collision locus.")
