"""Two-port-synchronizing reset with autonomous lower and upper sources."""
from itertools import product
from pathlib import Path
from math import isqrt
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'moonshine_prime_blocks'),str(HERE.parent/'dihedral_height')]
from rooted_pair import rooted,append,TT,decode
from replay import dp

def src(T,a):
    assert T[a].bit_count()==1
    h,e,t,*_=rooted(T,a)
    return h,e,t-h

def reset(T,a,n):
    B=TT(n)
    U,last=append(T,a,B,n-1)
    return U,last-1

def solve(p,h,k):
    z,r=divmod(p-h,h+k)
    if r or z<4 or z&(z-1):return None
    return z.bit_length()

def run():
    # Raw lower and upper inputs have matching (h,eta), not matching e.
    for A,B,lower,upper,p in (
        ((2,4,9,3,15,31,63,127),(2,4,33,23,39,75,31,127),(5,2,1),(35,12,7),29),
        ((16,5,9,19,38,15,63,127),(16,33,11,19,102,77,143,63),(29,8,1),(203,64,7),149)):
        assert src(A,0)==lower and src(B,0)==upper
        assert all(p%d for d in range(2,isqrt(p)+1))
        for target,U,s in ((p,A,1),(7*p,B,7)):
            n=solve(target,lower[0]*s,lower[2]*s)
            assert n==3
            V,b=reset(U,0,n)
            assert src(V,b)==(target,(target-lower[0]*s)//2,lower[0]*s)
        print('synchronized exact targets',p,7*p)
    checks=0
    for n0 in range(2,5):
        for code in range(1<<(n0*(n0-1)//2)):
            B=decode(code,n0)
            for a in range(n0):
                if B[a].bit_count()!=1:continue
                h,e,d=src(B,a)
                for n in range(3,6):
                    V,b=reset(B,a,n)
                    assert src(V,b)==(h+(1<<(n-1))*(h+d),(1<<(n-2))*(h+d),h)
                    checks+=1
    print('literal reset checks',checks)
if __name__=='__main__':run()
