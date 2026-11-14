"""Wave 24, THM-5162 (C'): lifting strong approximation to prime powers.
u1 = JOOZ^-1 JOZZ = I + N (N^3 = 0, entries in Z[1/2]).  For odd p, u1^p = I + p N' with
N' = N + ((p-1)/2) N^2 (mod p), a nonzero nilpotent element of sl_3(F_p).  We compute the F_p-span of the
SL_3(F_p)-conjugates of N' (the adjoint submodule it generates) by closing under conjugation by generators of
SL_3(F_p) (elementary matrices).  If the span is all of sl_3(F_p) (dimension 8), the closure of Gamma meets the
first congruence subgroup K_1 in a subgroup mapping onto K_1/K_2, hence equals K_1 (Frattini: K_1 is uniform for
p odd), so the closure of Gamma cap SL_3 in SL_3(Z_p) is SL_3(Z_p)."""
import numpy as np
from fractions import Fraction as F
def modp(x,p): return (x.numerator*pow(x.denominator,-1,p))%p
N=[[F(0),F(0),F(0)],[F(-5,4),F(0),F(-2)],[F(-1,4),F(0),F(0)]]
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
N2=mm(N,N); assert mm(N2,N)==[[0]*3]*3
def span_dim(p):
    Np=np.array([[ (modp(N[i][j],p)+((p-1)//2)*modp(N2[i][j],p))%p for j in range(3)] for i in range(3)],dtype=np.int64)
    gens=[]
    for i in range(3):
        for j in range(3):
            if i!=j:
                E=np.eye(3,dtype=np.int64); E[i,j]=1; Ei=np.eye(3,dtype=np.int64); Ei[i,j]=p-1
                gens.append((E,Ei))
    basis=[]  # row-echelon over F_p of flattened matrices
    def reduce(v):
        v=v.copy()%p
        for b,piv in basis:
            if v[piv]: v=(v-v[piv]*b)%p
        return v
    def add(v):
        v=reduce(v)
        nz=np.nonzero(v)[0]
        if len(nz)==0: return False
        piv=nz[0]; v=(v*pow(int(v[piv]),-1,p))%p
        for k,(b,pv) in enumerate(basis):
            if b[piv]: basis[k]=((b-b[piv]*v)%p,pv)
        basis.append((v,piv)); return True
    queue=[Np]; add(Np.flatten())
    while queue:
        X=queue.pop()
        for E,Ei in gens:
            Y=(E@X@Ei)%p
            if add(Y.flatten()): queue.append(Y)
    return len(basis), Np
for p in [3,5,7,11,13]:
    d,Np=span_dim(p); print("p=%d: N' mod p = %s, trace %d; dim of the adjoint submodule it generates: %d (sl_3 has 8)"%(p,Np.tolist(),int(np.trace(Np))%p,d))
