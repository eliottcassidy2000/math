"""Wave 24 (P52): unipotents A^-1 B from same-length invertible blocks, and the primes where they fail to
generate the Heisenberg group U3(F_p).  In (H,E,lam) coordinates A^-1 B = [[1,0,0],[alpha,1,beta],[gamma,0,1]],
lower unitriangular in the basis order (H, lam, E); its image in U/[U,U] = F_p^2 is (gamma, beta)."""
import itertools, sympy as sp
from math import gcd
from fractions import Fraction
def step(s,m):
    H,E,L,h=s
    if m=='O': return (H,0,2*L,h)
    if m=='Z': return (H+L-E,h,2*L-h,h)
    return (H+L-E,H-E,H+L,H)
def blockmat(w):
    cols=[]
    for e in ((1,0,0),(0,1,0),(0,0,1)):
        s=step((e[0],e[1],e[2],0),'J')
        for ch in w: s=step(s,ch)
        cols.append(s[:3])
    return sp.Matrix([[cols[c][r] for c in range(3)] for r in range(3)])
def odd(n):
    n=abs(int(n))
    while n and n%2==0: n//=2
    return n
J=blockmat('')
assert J==sp.Matrix([[1,-1,1],[1,-1,0],[1,0,1]])
vecs=[]
for n in range(2,7):
    ws=['O'+''.join(v)+'Z' for v in itertools.product('OZ',repeat=n-2)]
    Ms={w:blockmat(w) for w in ws}
    for A,B in itertools.combinations(ws,2):
        U=Ms[A].inv()*Ms[B]
        assert U[0,:].tolist()==[[1,0,0]] and U[1,1]==1 and U[2,1]==0 and U[2,2]==1, (A,B,U)
        N=U-sp.eye(3); assert (N**3)==sp.zeros(3)
        gam,bet=sp.Rational(U[2,0]),sp.Rational(U[1,2])
        vecs.append((A,B,gam,bet))
print("checked: A^-1 B is unipotent of the stated shape for all", len(vecs), "same-length pairs, |w| = 2..6")
G=0; ex=None
for (A,B,g1,b1),(C,D,g2,b2) in itertools.combinations(vecs,2):
    d=g1*b2-g2*b1
    if d!=0:
        num=odd(sp.fraction(d)[0]); G=gcd(G,num)
        if ex is None and num%3!=0: ex=(A,B,C,D,d)
print("odd part of gcd of all 2x2 determinants det[(gamma,beta),(gamma',beta')]:", G, sp.factorint(G) if G>1 else "")
# per-prime witness and smallest
for p in [3,5,7,11,13]:
    w=next(((A,B,C,D,g1*b2-g2*b1) for (A,B,g1,b1),(C,D,g2,b2) in itertools.combinations(vecs,2)
            if (g1*b2-g2*b1)!=0 and odd(sp.fraction(g1*b2-g2*b1)[0])%p!=0),None)
    print("p=%d witness:"%p, w)
# the witness pair inside the seven blocks of wave 23
A,B,C,D=[blockmat(w) for w in ['OOZ','OZZ','OOOZ','OZZZ']]
u1=A.inv()*B; u2=C.inv()*D
print('u1 =',u1.tolist(),' u2 =',u2.tolist())
print('det of abelianized images =',u1[2,0]*u2[1,2]-u2[2,0]*u1[1,2],'  commutator [u1,u2] =',(u1*u2*u1.inv()*u2.inv()).tolist())
# opposite-flag condition for J: (J)_{H,E} and (J^-1)_{H,E}
Ji=J.inv(); print("J_{H,E} =",J[0,1]," (J^-1)_{H,E} =",Ji[0,1], " det J =",J.det())
