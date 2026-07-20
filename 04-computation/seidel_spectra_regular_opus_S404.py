"""opus-2026-07-20-S404 part 3 -- SEIDEL SPECTRA, REGULAR TWO-GRAPHS, and the ORDER-14
conference matrix.  Extends the S404 census; complements kps THM-1415 (tournament side)."""
import itertools, numpy as np
from collections import defaultdict
def pairs(k): return list(itertools.combinations(range(k),2))
def seidel(g,k):
    """S = J - I - 2A  (Seidel matrix).  Switching = conjugation by diag(+-1),
       so the SPECTRUM is a switching invariant AND an iso invariant."""
    S=np.ones((k,k),dtype=int)-np.eye(k,dtype=int)
    P=pairs(k)
    for i,(a,b) in enumerate(P):
        if (g>>i)&1: S[a,b]=S[b,a]=-1
    return S
def canon_graph(g,k,_c={}):
    if k not in _c:
        P=pairs(k); idx={p:i for i,p in enumerate(P)}
        _c[k]=[[idx[tuple(sorted((s[a],s[b])))] for (a,b) in P] for s in itertools.permutations(range(k))]
    best=None
    for act in _c[k]:
        h=0
        for i in range(len(act)):
            if (g>>i)&1: h|=1<<act[i]
        if best is None or h<best: best=h
    return best
print("="*76)
print("(3) SEIDEL SPECTRA of the switching classes, and REGULAR TWO-GRAPHS")
print("="*76)
print("  A REGULAR two-graph <=> its Seidel matrix has exactly TWO distinct eigenvalues.")
print("  These are the equiangular-line / conference-matrix objects.")
print()
print("   n  #distinct Seidel spectra  #REGULAR two-graph spectra")
for n in range(3,7):
    specs=set()
    for g in range(1<<len(pairs(n))):
        ev=np.round(np.linalg.eigvalsh(seidel(g,n).astype(float)),6)
        specs.add(tuple(sorted(ev)))
    reg=[ev for ev in specs if len(set(ev))==2]
    print(f"  {n:2d}  {len(specs):23d}  {len(reg):24d}"
          f"   {'<- the unique regular two-graph on 6 points' if n==6 and reg else ''}")
    for ev in reg:
        vals=sorted(set(ev))
        print(f"        regular spectrum: {vals} multiplicities "
              f"{[list(ev).count(x) for x in vals]}")

print()
print("="*76)
print("(4) THE ORDER-14 PALEY CONFERENCE MATRIX  ->  regular two-graph on 14 points")
print("="*76)
print("  Symmetric conference matrix of order n <=> regular two-graph on n vertices.")
print("  Paley construction needs n-1 an odd prime power = 1 mod 4.  n=14: n-1=13, and")
print("  13 = 1 mod 4.  So C_14 EXISTS, built from the quadratic residues mod 13.")
p=13
chi=[0]*p
QR={ (x*x)%p for x in range(1,p) }
for x in range(1,p): chi[x]= 1 if x in QR else -1
C=np.zeros((p+1,p+1),dtype=int)
for i in range(1,p+1):
    C[0,i]=1; C[i,0]=1
for i in range(1,p+1):
    for j in range(1,p+1):
        if i!=j: C[i,j]=chi[(i-j)%p]
print(f"  QR mod 13 = {sorted(QR)}   (13 = number of LRC speeds; 14 = the repo apex)")
print(f"  C is {C.shape[0]}x{C.shape[1]}, symmetric: {np.array_equal(C,C.T)}")
print(f"  zero diagonal: {bool(np.all(np.diag(C)==0))};  entries in {{0,+-1}}: {set(np.unique(C))}")
CCt=C@C.T
print(f"  C C^T = {C.shape[0]-1} * I ?  {np.array_equal(CCt,(p)*np.eye(p+1,dtype=int))}"
      f"   (conference-matrix defining identity)")
ev=np.round(np.linalg.eigvalsh(C.astype(float)),6)
print(f"  Seidel spectrum: {sorted(set(ev))} with multiplicities "
      f"{[list(ev).count(x) for x in sorted(set(ev))]}  -> TWO eigenvalues = REGULAR")
print(f"  => a regular two-graph on 14 points exists, eigenvalues +-sqrt(13).")
print()
print("  HONEST NOTE: this is a real object at n=14 built from 13 residues, but the repo")
print("  rule from HYP-8230 applies -- coincidence of INDEX is not evidence of mechanism.")
print("  The LRC 14 is 2*7 (parity, lam=1/14); this 14 is 13+1 with 13 = 1 mod 4.")
print("  Different arithmetic. Recorded as a catalogued object, NOT as a bridge.")
