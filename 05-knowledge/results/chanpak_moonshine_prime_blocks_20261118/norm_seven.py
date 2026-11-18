"""Literal input profiles and direct prime / sevenfold norm witnesses."""
from pathlib import Path
from math import isqrt
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'moonshine_prime_blocks'))
from clone_replay import inflate
from replay import profile,dp
Q=(30,8,58,48,34,3)
R=(54,16,10,19,36,14)
W=(
    (Q,(30,12,24,16,2),5,64,4001),
    (Q,(126,92,24,112,96,70,4),29,640,404177),
    (R,(126,76,56,112,98,66,4),37,752,572867)
)
def run():
    assert pow(2,21,49)==1
    assert profile((460,25,426,400,997,331,14,98,194,495))[0]==4001
    for T,U,h,b,p in W:
        assert profile(U)==(h,b)
        assert all(p%d for d in range(2,isqrt(p)+1))
        assert sum(dp(inflate(T,{0:U,1:U}))[-1])==7*p
        f=lambda h,u:(2*h+3*u)**2-2*u*u if T==Q else 6*h*h+16*h*u+7*u*u
        assert f(h,b-h)==7*p
        assert f(h,b*2**21-h)%49==7*(p%7)
        print(p,7*p,'profile',h,b)
    assert 187**2-2*59**2==28007
if __name__=='__main__':run()
