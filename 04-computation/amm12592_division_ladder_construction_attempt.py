"""Parity-corrected division ladder.

e_{k,i} must be = C(a_k,i) (mod 2).  Write e = f + 2g with f_{k,i} =
C(a_k,i) mod 2, i.e. F_k(u) = sum_{i subset a_k} u^i.  By reduction A,
sum_k F_k(u)(1+u)^{L_k} = u^{m-1} (mod 2), so

    W(u) := [ eps u^{m-1} - sum_k F_k(u)(1+u)^{L_k} ] / 2

has INTEGER coefficients, and it remains to solve sum_k G_k(1+u)^{L_k} = W
by the division ladder, then set e = f + 2g and test |e| <= C(a_k,i).
"""
import sys
from fractions import Fraction
from math import comb
sys.path.insert(0,'/Users/e/Documents/GitHub/math/04-computation')
from amm12592_exact_block_profile_solver import profile

def poly_div_full(R, L):
    D=[comb(L,i) for i in range(L+1)]; R=R[:]
    q=[0]*max(1,len(R)-L)
    for j in range(len(R)-1, L-1, -1):
        c=R[j]
        if c:
            q[j-L]=c
            for i,dv in enumerate(D): R[j-L+i]-=c*dv
    return q, R[:L]

def ladder(target, L, a, m):
    R=target[:]; 
    while len(R)>1 and R[-1]==0: R.pop()
    G=[]
    for k in range(m):
        q,r=poly_div_full(R, L[k])
        while len(q)>1 and q[-1]==0: q.pop()
        if len(q)-1 > a[k]:
            return None, f"quotient degree {len(q)-1} > a_{k}={a[k]}"
        G.append(q); R=r[:]
        while len(R)>1 and R[-1]==0: R.pop()
        if all(c==0 for c in R):
            G.extend([[0]]*(m-1-k)); break
    if any(c!=0 for c in R): return None, "nonzero remainder"
    return G, None

def attempt(m, C, eps):
    a=profile(m,C)
    if a is None: return "no profile", None
    L=[m-1-k-a[k] for k in range(m)]
    F=[[1 if comb(a[k],i)%2 else 0 for i in range(a[k]+1)] for k in range(m)]
    # W = [eps u^{m-1} - sum F_k (1+u)^{L_k}]/2
    acc=[0]*(2*m+2)
    for k in range(m):
        for i,c in enumerate(F[k]):
            if c:
                for e in range(L[k]+1): acc[i+e]+=c*comb(L[k],e)
    W=[0]*(2*m+2); W[m-1]=eps
    W=[W[j]-acc[j] for j in range(len(W))]
    if any(w%2 for w in W): return "reduction A FAILED", None
    W=[w//2 for w in W]
    while len(W)>1 and W[-1]==0: W.pop()
    G,err=ladder(W,L,a,m)
    if G is None: return err, None
    worst=None
    for k in range(m):
        Gk=G[k] if k<len(G) else [0]
        for i in range(a[k]+1):
            g=Gk[i] if i<len(Gk) else 0
            e=F[k][i]+2*g
            box=comb(a[k],i)
            if abs(e)>box:
                rr=Fraction(abs(e),box)
                if worst is None or rr>worst[0]: worst=(rr,k,i,e,box)
    return ("OK" if worst is None else "BOX"), worst

print("parity-corrected division ladder")
for m in (4,8,16,32,64,128):
    row=[]
    for C in (Fraction(8,5),Fraction(5,3),Fraction(7,4),Fraction(15,8),Fraction(2)):
        best=None
        for eps in (1,-1):
            st,w=attempt(m,C,eps)
            if st=="OK": best="OK"; break
            if best is None: best=(st,w)
        if best=="OK": row.append(f"C={C}:OK")
        else:
            st,w=best
            extra=f"(x{float(w[0]):.2f} at k={w[1]},i={w[2]})" if w else ""
            row.append(f"C={C}:{st}{extra}")
    print(f"m={m:4d}  " + "  ".join(row))
