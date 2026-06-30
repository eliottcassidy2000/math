"""Exact cusp orders in the cuspidal group of X_0(14)=Z/6, and complement R=[-1] (negation)."""
import math
from sympy import Matrix, ilcm, gcd
divs=[1,2,7,14]
M=Matrix([[ (14//c)*(math.gcd(c,d)**2)//d for d in divs] for c in divs])  # 24*ord_c(eta(d tau))
# principal divisors: D(r) = (1/24) M r, r in Z^4, sum r =0, AND 24 | (M r)_c. 
# Cuspidal group = Z^4_0 / L, L = { (1/24) M r integer }. Compute order of [c]-[oo](=e_c - e_{14}).
# Approach: a degree-0 divisor v (integer) is principal iff v = (1/24) M r for integer r,sum r=0.
# order of e_c - e_14 = min k>0 s.t. 24k(e_c - e_14) in column space of M over Z (with sum-0 r).
def in_image(v):  # is there integer r with M r = v? (v must be integer vector)
    # solve M r = v over Q, check integer solution exists (M rank 3, ones in kernel-cokernel)
    aug=M.row_join(Matrix(v))
    sol=M.gauss_jordan_solve(Matrix(v)) if False else None
    try:
        r=M.solve(Matrix(v))  # particular rational solution
    except Exception:
        return False
    # M singular (rank3); use least squares over Q via pinv won't give integrality. Instead: r exists integer iff v in Z-span of columns of M. Use Smith.
    return None
# Smith normal form of M: M = U D V, U,V unimodular. Then M r = v solvable in Z iff D | (U v).
from sympy.matrices.normalforms import smith_normal_form
from sympy import diag
# get U, V too via algorithm: use smith_normal_form for D; for U,V use the .T trick is unavailable -> compute manually
def smith_with_transforms(A):
    A=A.copy(); m,n=A.shape
    U=Matrix.eye(m); V=Matrix.eye(n)
    def swap_rows(i,j):
        A.row_swap(i,j); U.row_swap(i,j)
    def swap_cols(i,j):
        A.col_swap(i,j); V.col_swap(i,j)
    t=0
    import sympy
    while t<min(m,n):
        # find pivot nonzero
        piv=None
        for i in range(t,m):
            for j in range(t,n):
                if A[i,j]!=0: piv=(i,j);break
            if piv:break
        if not piv: break
        swap_rows(piv[0],t); swap_cols(piv[1],t)
        done=False
        while not done:
            done=True
            for i in range(t+1,m):
                if A[i,t]!=0:
                    q=A[i,t]//A[t,t]
                    A[i,:]=A[i,:]-q*A[t,:]; U[i,:]=U[i,:]-q*U[t,:]
                    if A[i,t]!=0:
                        swap_rows(i,t); done=False
            for j in range(t+1,n):
                if A[t,j]!=0:
                    q=A[t,j]//A[t,t]
                    A[:,j]=A[:,j]-q*A[:,t]; V[:,j]=V[:,j]-q*V[:,t]
                    if A[t,j]!=0:
                        swap_cols(j,t); done=False
        t+=1
    return U,A,V
U,D,V=smith_with_transforms(M)
# M r = v solvable in Z  <=>  (U v)_i divisible by D[i,i] for i<rank, and (U v)_i==0 for i>=rank
diagD=[D[i,i] for i in range(4)]
rank=sum(1 for d in diagD if d!=0)
def solvable(v):
    Uv=U*Matrix(v)
    for i in range(4):
        if i<rank:
            if diagD[i]==0: 
                if Uv[i]!=0: return False
            elif Uv[i]%diagD[i]!=0: return False
        else:
            if Uv[i]!=0: return False
    return True
print("Smith diag of M (24*ord matrix):",diagD)
names={1:'T (transitive, d=1)',2:'+ (dom-out, d=2)',7:'- (dom-in, d=7)',14:'S (regular, d=14)'}
print("\nOrder of each cusp difference [c]-[d=14] in the cuspidal group (=Z/6):")
e={d:[1 if divs[i]==d else 0 for i in range(4)] for d in divs}
base=14
for c in divs:
    if c==base: 
        print(f"   {names[c]:>24}: BASE (identity, order 1)"); continue
    v0=[24*(e[c][i]-e[base][i]) for i in range(4)]  # 24*(e_c - e_base); principal iff this = M r
    order=None
    for k in range(1,13):
        vk=[k*x for x in v0]
        if solvable(vk): order=k;break
    print(f"   {names[c]:>24}: order {order}")
print("\nComplement R (tournament) swaps + <-> - and fixes T,S. On Z/6 the fixed pts of [-1] are {0,3} (2-torsion).")
print("  => R = [-1] (NEGATION on 14a): T=0 (id), S=3 (2-torsion), {+,-} = a +/- pair (3-torsion {2,4}).")
print("  SC classes {T,S} = the 2-TORSION (Z/2, R-fixed); complement pair {+,-} = the 3-TORSION (Z/3, R-swapped).")
print("  Z/6 = Z/2 x Z/3 = (the SC mirror) x (the complement-pair) ; 6 = phi(14) = units = EXISTENCE column.")
