"""
WITNESS-FIRST search for the odd-n covering-min. Fix witness (q,a) and gap m; the SAFE BAND is
B={r mod q: min(r,q-r)>=m}. Build a covering set by picking, for each q' in {2..n}, a multiple of q'
whose residue *a mod q lands in B. The band guides the construction. Compute true M, track the min.
"""
import math
from fractions import Fraction
def M_exact(S,Qmax):
    best=Fraction(0)
    for q in range(2,Qmax+1):
        bb=0
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
                if m<=bb:break
            if m>bb:bb=m
        v=Fraction(bb,q)
        if v>best:best=v
    return best
def is_cov(S,n): return all(any(x%q==0 for x in S) for q in range(2,n+1))
def smart_covering_min(n):
    Qmax=2*(n*n-n+1)
    best=Fraction(1); bestS=None; bestw=None
    for q in range(n+1, 2*n*n):
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            # for gap m near q/n..q/(n-1)
            for m in range(max(1,q//(n+1)), q//(n-1)+1):
                band=set(r for r in range(q) if min(r,q-r)>=m)
                S=[]; ok=True
                for qp in range(2,n+1):
                    found=None
                    for k in range(1, 3*q):
                        sp=k*qp
                        if (sp*a)%q in band: found=sp; break
                    if found is None: ok=False; break
                    S.append(found)
                if not ok: continue
                Sset=sorted(set(S))
                if len(Sset)<n-1:  # need exactly n-1 distinct; pad with band elements covering nothing new
                    extra=[x for x in range(1,3*q) if (x*a)%q in band and x not in Sset]
                    for e in extra:
                        if len(Sset)>=n-1: break
                        Sset=sorted(set(Sset+[e]))
                if len(Sset)!=n-1 or not is_cov(Sset,n): continue
                M=M_exact(Sset,Qmax)
                if M<best: best=M; bestS=Sset; bestw=(a,q,m)
        if q>3*n*n//2: break
    return best,bestS,bestw
for n in [7,9]:
    M,S,w=smart_covering_min(n)
    Phi6=n*n-n+1
    print(f"n={n}(ODD): smart covering-min={M}={float(M):.5f} at {S}  witness(a,q,m)={w}  (construction {n}/{Phi6}={n/Phi6:.5f}, 1/n={1/n:.5f})")
    if S:
        odds=[x for x in S if x%2]; print(f"   odds={odds} (cover odd q=n); descent odd-core mod {min(7,n)}: {sorted(set(x% (7 if n>=7 else n) for x in odds))}")
