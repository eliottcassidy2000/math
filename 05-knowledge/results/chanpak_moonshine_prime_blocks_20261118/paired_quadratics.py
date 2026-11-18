"""Independent visit tensor replay for simultaneous paired quadratics."""
from math import comb,isqrt,gcd,prod
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'moonshine_prime_blocks'))
from clone_replay import inflate,vectors,TT,read
from replay import profile
from dihedral_replay import count
DATA=(
 ((34,32,11,19,39,12),(6,4,24,51,195,87,143,47),(25,12,12,6),(175,84,84,42),6,1),
 ((62,16,26,50,96,70,15),(6,4,24,51,195,343,143,303,95),(67,32,32,16),(469,224,224,112),16,3),
 ((66,64,3,23,39,79,28),(6,4,336,295,459,343,395,303,3),(53,30,30,15),(371,210,210,105),15,-7),
 ((190,16,2,22,100,142,175,30),(6,4,152,115,195,471,391,11,159),(61,41,23,16),(427,224,224,112),16,-3),
)
def run():
    for k,(A,B,C,D,a,b) in enumerate(DATA):
        X=tuple(vectors(T,(0,1),1,len(T)-1) for T in (A,B))
        for r,c in zip(X,(C,D)):
            assert r==dict(zip(((1,1),(1,2),(2,1),(2,2)),c))
        for U,V in [(TT(1),TT(1)),(TT(2),TT(2)),(TT(2),TT(3))]:
            h,F=profile(U);g,G=profile(V)
            p=a*F*G+b*h*g
            x,y=[count(inflate(T,{0:U,1:V})) for T in (A,B)]
            assert y==7*p
            assert x==p+(9*(h*G-F*g) if k==3 else 0)
        print(a,b,'diagonal TT2',a*16+b,(a*16+b)*7)
    for U,p in (((22,60,56,49,32,1),393577),((126,124,88,112,96,68,0),393241)):
        assert all(p%i for i in range(2,isqrt(p)+1))
        h,B=profile(U)
        assert h*h+6*B*B==p
    for k in range(1,5):
        U=TT(k);h,B=profile(U);C=read(U,3)
        T=inflate((2,4,1),{0:U,1:U})
        assert profile(T)==(2*h*h+(B-h)**2,2*(B*B+(C-B)**2+(B-h)**2))
    U=(126,252,56,112,96,64,132,61)
    h,B=profile(U); C=read(U,3)
    assert (h,B,C)==(71,1976,21591)
    p=4144*(B**2+(C-B)**2+(B-h)**2)-209*(2*h*h+(B-h)**2)
    F=(2,2,3,5,181,1447,103399)
    assert p==1624855313581 and 7*p==11373987195067
    assert prod(F)==p-1 and pow(6,p-1,p)==1
    assert all(all(q%d for d in range(2,isqrt(q)+1)) and gcd(pow(6,(p-1)//q,p)-1,p)==1 for q in F)
    print('nested cyclic then genuine C7 affine',p,7*p)
if __name__=='__main__':run()
