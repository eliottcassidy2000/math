"""Physical readout formulas and absolute target inverse for transitive tails."""
from fractions import Fraction
from math import isqrt
from functools import lru_cache
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'moonshine_prime_blocks'),str(HERE.parent/'dihedral_height')]
from rooted_pair import rooted,TT,append,decode
from replay import dp

@lru_cache(None)
def state(T,a):
    h,e,t=rooted(T,a)[:3]
    H=dp(T); s=(1<<len(T))-1-(1<<a)
    return h,e,t,sum(H[s]) if s else 1

def move(T,a,i,j):
    V,b=append(T,a,TT(i+j+1),i+j)
    return V,b-j

def predicted(h,e,t,d,i,j):
    k=i+j
    return h+2**k*t,2**i*t,2**j*h+2**(i+1)*3**(j-1)*t-(e if i==0 else 0),h+2**(k-1)*t

def band_j(N,T):
    if T<=N:return None
    j=1
    while Fraction(T,N)>Fraction(3,2)**j:j+=1
    return None if Fraction(T,N)==Fraction(3,2)**j else j

def decode_double(p,h,t):
    for j in range(1,p.bit_length()):
        for i in range(1,p.bit_length()):
            N=h+2**(i+j)*t
            if N>=p:break
            T=2**j*h+2**(i+1)*3**(j-1)*t
            x,r=divmod(p-N,T)
            if not r and x>=2 and not x&(x-1):
                yield (i,j,x.bit_length()-1)

def run():
    inputs=(((8,1,3,6),1),((64,77,105,33,79,83,8),0))
    for inputs_pair,(h,t),p,expected in [(inputs,(5,6),421,(1,2,2)),(inputs,(5,6),613,(1,3,1))]:
        hits=list(decode_double(p,h,t))
        assert expected in hits
        i,j,r=expected
        for q,(A,a) in zip((1,7),inputs_pair):
            S=state(A,a)
            T,b=move(A,a,i,j)
            U=state(T,b)
            assert U==predicted(*S,i,j)
            assert band_j(U[0],U[2])==j
            c=Fraction(3,2)**(j-1)
            assert (U[2]-c*U[0])/(2**j-c)==q*h
            W,z=append(T,b,TT(r+1),r)
            assert state(W,z)==(q*p,U[0],q*p,U[0])
        assert all(p%d for d in range(2,isqrt(p)+1))
        print('prescribed two-stage targets',p,7*p,'choice',expected)
    count=0
    for n in range(1,4):
        for code in range(1<<(n*(n-1)//2)):
            A=decode(code,n)
            for a in range(n):
                S=state(A,a)
                for i,j in ((1,1),(1,2),(2,1),(1,3),(0,2),(0,3)):
                    B,b=move(A,a,i,j)
                    U=state(B,b)
                    assert U==predicted(*S,i,j)
                    assert band_j(U[0],U[2])==j
                    count+=1
    print(count,'independent literal four-state checks')
if __name__=='__main__':run()
