"""Verify literal two-copy norm forms, including the positive Pell branch."""
from pathlib import Path
from random import Random
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'moonshine_prime_blocks'))
from clone_replay import vectors,inflate,TT
from gap_replay import rand_t
from replay import profile

def run():
    data=(
      ((2,4,1),(3,1,1,1),lambda h,u:u*u+2*h*h),
      ((6,4,8,3),(5,2,2,1),lambda h,u:(u+h)**2+h*h),
      ((30,8,58,48,34,3),(23,11,15,7),lambda h,u:(2*h+3*u)**2-2*u*u),
      ((2,4,9,19,7),(9,5,3,3),lambda h,u:4*h*h+2*h*u+3*u*u),
      ((54,16,10,19,36,14),(29,15,15,7),lambda h,u:6*h*h+16*h*u+7*u*u),
    )
    r=Random(111079)
    for A,c,norm in data:
        assert dict(vectors(A,(0,1),1,len(A)-1))==dict(zip(((1,1),(1,2),(2,1),(2,2)),c))
        for k in range(1,4):
            for T in (TT(k),rand_t(k,r)):
                h,b=profile(T)
                assert profile(inflate(A,{0:T,1:T}))[0]==norm(h,b-h)
    A=data[2][0]
    print('literal Pell outputs',[profile(inflate(A,{0:TT(k),1:TT(k)}))[0] for k in (1,2,3)])
    assert 23**2-2*7**2==431
    for ind,ks,primes in ((3,(2,4),(37,709)),(4,(1,3,5),(29,461,7229))):
        A=data[ind][0]
        assert tuple(profile(inflate(A,{0:TT(k),1:TT(k)}))[0] for k in ks)==primes
        assert all(all(p%d for d in range(2,1+int(p**.5))) for p in primes)
if __name__=='__main__':run()
