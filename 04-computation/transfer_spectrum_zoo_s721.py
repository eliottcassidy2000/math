#!/usr/bin/env python3
"""
S721 — Using the transfer-spectrum mindset to make progress: a ZOO of C-finite tournament families.

HYP-2381 predicts, for an order-r C-finite tournament invariant family, char poly x^r + ...:
  - e_last = product of eigenvalues = +-1   (PARITY / unimodular transfer = the Pfaffian face, S713);
  - e1 = sum of eigenvalues = the GEOMETRY (corners/branching), should be a small fixed integer per shape;
  - the middle symmetric functions = the TEMPERATURE (additive ~ all roots 1; multiplicative ~ root >1).

We TEST this across many "Toeplitz" staircase tournaments (tile (i,j), j>i+1, oriented by a rule f(j-i)),
indexed by n. For each family we compute H(n) (#Hamiltonian paths) and Pf(n) (Pfaffian of the skew-adj),
auto-discover the minimal linear recurrence, and read off (order, e1=sum roots, e_last=product, dominant
root=temperature). This (a) tests the unimodular-transfer conjecture (t-0107.2), (b) maps intermediate
temperatures (t-0107.1), (c) checks whether Pf is C-finite with related spectrum (t-0102/t-0105).

No numpy/sympy.
"""
from fractions import Fraction as Fr

# ---- recurrence finder (min order, integer/rational coeffs) ----
def gauss(A,b):
    A=[r[:] for r in A]; b=b[:]; n=len(A)
    for k in range(n):
        p=next((i for i in range(k,n) if A[i][k]!=0),None)
        if p is None: return None
        A[k],A[p]=A[p],A[k]; b[k],b[p]=b[p],b[k]
        inv=A[k][k]; A[k]=[x/inv for x in A[k]]; b[k]/=inv
        for i in range(n):
            if i!=k and A[i][k]!=0:
                f=A[i][k]; A[i]=[A[i][j]-f*A[k][j] for j in range(n)]; b[i]-=f*b[k]
    return b
def find_rec(seq, maxord=5):
    n=len(seq)
    for r in range(1,maxord+1):
        if n<2*r+1: break
        A=[[Fr(seq[i-1-j]) for j in range(r)] for i in range(r,2*r)]
        rhs=[Fr(seq[i]) for i in range(r,2*r)]
        c=gauss(A,rhs)
        if c is None: continue
        if all(sum(c[j]*seq[i-1-j] for j in range(r))==seq[i] for i in range(r,n)):
            return r,[int(x) if x.denominator==1 else x for x in c]
    return None,None
def charpoly_from_rec(c):
    # recurrence s_n = c0 s_{n-1}+...+c_{r-1} s_{n-r}  =>  char x^r - c0 x^{r-1} - ... - c_{r-1}
    r=len(c); poly=[1]+[-c[i] for i in range(r)]   # [1, -c0, -c1, ..., -c_{r-1}]
    e1=c[0]                      # sum of roots = -coeff(x^{r-1}) = c0
    e_last=((-1)**(r+1))*c[-1]   # product of roots = (-1)^r * coeff(x^0)/1 ; coeff x^0 = -c_{r-1}
    return poly,e1,e_last

# ---- tournaments ----
def ham_paths(A):
    n=len(A); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            d=dp[mask][v]
            if not d or not (mask>>v)&1: continue
            for w in range(n):
                if (mask>>w)&1: continue
                if A[v][w]: dp[mask|(1<<w)][w]+=d
    return sum(dp[(1<<n)-1][v] for v in range(n))
def pf(M):
    n=len(M)
    if n==0: return 1
    if n%2: return 0
    t=0;r0=M[0]
    for j in range(1,n):
        a=r0[j]
        if a==0: continue
        rest=[k for k in range(1,n) if k!=j]
        t+=((-1)**(j-1))*a*pf([[M[r][c] for c in rest] for r in rest])
    return t
def skew(A):
    n=len(A); return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]

def toeplitz_tournament(n, rule):
    """base path n-1->...->0 (consecutive higher->lower); tile (i,j), j>i+1, oriented by rule(j-i)."""
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if j==i+1: A[j][i]=1                    # base path arc
            else:
                if rule(j-i): A[i][j]=1             # 'up': lower i beats higher j
                else:        A[j][i]=1              # 'down': higher beats lower
    return A

RULES = {
  "up(all)":        lambda g: True,
  "down(all)":      lambda g: False,
  "gap even->up":   lambda g: g%2==0,
  "gap odd->up":    lambda g: g%2==1,
  "gap mod3==0->up":lambda g: g%3==0,
  "gap<=2->up":     lambda g: g<=2,
  "QR(gap mod7)":   lambda g: (g%7) in (1,2,4),   # quadratic residues mod 7
}

if __name__=="__main__":
    print("="*92)
    print("S721 — transfer-spectrum zoo: testing frozen geometry/parity + running temperature (HYP-2381)")
    print("="*92)
    NS=list(range(3,15))   # n=3..14
    print(f"\nFamilies indexed by n={NS[0]}..{NS[-1]}. For H(n) and Pf(n): minimal recurrence + spectrum.\n")
    print(f"{'family':16} | {'inv':3} | {'ord':3} | {'e1=sum':6} | {'e_last=prod':11} | unimodular? | dominant root (temp)")
    print("-"*92)
    for name,rule in RULES.items():
        Hs=[ham_paths(toeplitz_tournament(n,rule)) for n in NS]
        # Pf only defined for even n
        Pfs=[abs(pf(skew(toeplitz_tournament(n,rule)))) for n in NS if n%2==0]
        for inv,seq in (("H",Hs),("Pf",Pfs)):
            r,c=find_rec(seq)
            if r is None:
                print(f"{name:16} | {inv:3} | {'>5?':3} |   --   |     --      |     --      | seq={seq[:6]}...")
                continue
            poly,e1,el=charpoly_from_rec(c)
            uni = "YES" if el in (1,-1) else "no("+str(el)+")"
            # dominant root via crude bisection on |companion| -> use power iteration magnitude estimate
            # estimate growth ratio from last terms
            ratio = (seq[-1]/seq[-2]) if len(seq)>=2 and seq[-2]!=0 else float('nan')
            print(f"{name:16} | {inv:3} | {r:3} | {e1:6} | {el:11} | {uni:^11} | char {poly}  growth~{ratio:.3f}")
    print("-"*92)
    print("\nReading: e_last in {+1,-1} confirms UNIMODULAR transfer = parity (Pfaffian face, S713).")
    print("e1 (sum of roots) = the geometry; the dominant root = the temperature (1=additive, >1=hot).")
    print("Families with the SAME (order,e1,e_last) but different middle = different TEMPERATURE only.")
    print("="*92)
