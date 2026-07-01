"""TOURNAMENT + DIHEDRAL synthesis: the 6 unit atoms {1,3,5,9,11,13}/14 + the added antipode 7/14=1/2 form the
7 vertices at odd/14 = the ROOTS OF z^7=-1 (a regular heptagon), with dihedral D_7 (order 14 = the LRC 14).
vertex k <-> (2k+1)/14. Build the natural rotational tournament; compute invariants."""
import numpy as np
from itertools import combinations, permutations
V=list(range(7)); pos=[(2*k+1)/14 for k in V]  # vertex k at (2k+1)/14
print("7 vertices (odd/14):", [f"{2*k+1}/14" for k in V], " antipode = vertex 3 = 7/14 = 1/2")
# roots of z^7 = -1
z=[np.exp(2j*np.pi*p) for p in pos]
print("z_k^7 all = -1?", np.allclose([zz**7 for zz in z], -1))
# D_7 action: rotation rho: k->k+1 (add 2/14=1/7 to position); reflection iota: x->-x  => k->6-k
rho=lambda k:(k+1)%7
iota=lambda k:(6-k)%7
print(f"iota (antipode x->-x): pairs {[(k,iota(k)) for k in range(3)]}, FIXED vertex {[k for k in V if iota(k)==k]} (=7/14)")
print(f"  => 6 UNITS = non-fixed vertices {[k for k in V if iota(k)!=k]} in 3 iota-pairs; antipode = the reflection center.")
# units check: (2k+1) coprime to 14  <=> 2k+1 != 7  <=> k != 3
print(f"  (2k+1) a UNIT mod 14 (coprime): {[2*k+1 for k in V if np.gcd(2*k+1,14)==1]} = the 6 units; non-unit: 7 (=vertex 3)")
# ---- natural ROTATIONAL tournament: i->j iff (j-i) mod 7 in {1,2,3} (beat the next three) ----
def rot_tour(conn):
    A=np.zeros((7,7),int)
    for i in V:
        for j in V:
            if i!=j and ((j-i)%7) in conn: A[i,j]=1
    return A
A=rot_tour({1,2,3})
scores=A.sum(1)
print(f"\nROTATIONAL tournament (beat next 3): scores={list(scores)} (regular, all 3)")
# 3-cycles and transitive triples
c3=sum(1 for t in combinations(V,3) if A[t[0],t[1]]+A[t[1],t[2]]+A[t[2],t[0]] in (0,3) and (A[t[0],t[1]]==A[t[1],t[2]]==A[t[2],t[0]]))
def is_cyclic(t):
    a,b,c=t; e=[(a,b),(b,c),(c,a)]; d=[A[x,y] for x,y in e]; 
    return sum(d)==2 and set([A[a,b],A[b,c],A[c,a]])=={0,1} and (A[a,b]+A[b,c]+A[c,a] in (1,2))
cyc=0
for t in combinations(V,3):
    a,b,c=t
    # cyclic iff each vertex has in and out within the triple: sum of the 3 directed indicators (a->b,b->c,c->a) is 1 or 2 but specifically a 3-cycle
    s=A[a,b]+A[b,c]+A[c,a]
    if s==1 or s==2:  # not all same direction => but need actual 3-cycle
        pass
# cleaner: count cyclic triangles = C(n,3) - sum C(score,2)
cyc=int(len(list(combinations(V,3))) - sum(s*(s-1)//2 for s in scores))
print(f"  cyclic 3-cycles = C(7,3) - sum C(score,2) = 35 - {sum(s*(s-1)//2 for s in scores)} = {cyc}   (= 14 = |D_7| = the LRC 14!)")
print(f"  transitive triples = 35 - {cyc} = {35-cyc} (= 21 = C(7,2) = #edges)")
# Hamiltonian paths (Redei: odd)
def ham_paths(A):
    cnt=0
    for p in permutations(V):
        if all(A[p[i],p[i+1]]==1 for i in range(6)): cnt+=1
    return cnt
hp=ham_paths(A)
print(f"  Hamiltonian paths = {hp} (Redei: ODD? {hp%2==1})")
# self-converse via iota: iota maps arc i->j to iota(i)->iota(j); check T-reversed ~ T under iota
Aop=1-A-np.eye(7,dtype=int); Aop=np.where(np.eye(7,dtype=int)==1,0,1-A)  # reverse
P=np.zeros((7,7),int)
for k in V: P[k,iota(k)]=1
Aiota=P.T@A@P
print(f"  iota is an ANTI-automorphism (T -> T^op)? {np.array_equal(Aiota, (1-A)*(1-np.eye(7,dtype=int)))}  => self-converse (repo SC-flavor)")
# skew-adjacency spectrum (Gauss sums)
S=A-A.T  # skew-symmetric +-1
ev=np.linalg.eigvals(S.astype(float))
print(f"  skew-adjacency eigenvalues (imag parts): {sorted(np.round(ev.imag,3))}")
print(f"    |nonzero eigenvalue| = {np.round(max(abs(ev)),4)}  (sqrt(7)={np.sqrt(7):.4f}? Gauss-sum structure)")
