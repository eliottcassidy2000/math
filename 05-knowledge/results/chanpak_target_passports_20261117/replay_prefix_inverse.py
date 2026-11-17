from math import isqrt
from functools import cache
from replay_ordinal_digits import profile, external, ordinal, strong, end_profiles

R=(50,104,19,21,98,76,13)
MASKS=(31,47)
def need(x):
    if not x: raise AssertionError('prefix inverse replay failed')
def tt(n):
    return tuple((1<<n)-(1<<(i+1)) for i in range(n))

@cache
def F(v):
    if not v:return 1
    q,b=divmod(v,2)
    return F(q)+(q+(q&1) if b==0 else 0)

def build(v,k,digit):
    n=v.bit_length()
    tags=sum(((v>>(n-1-i))&1)<<i for i in range(n))
    return external(ordinal(tt(n),R),tags|(MASKS[digit]<<n),k)

def count(v,k):
    t=1<<k
    return 185*t-112+(284*t-210)*v+25*F(v)

def low(r,ell):
    if ell==0:return 0,0
    x=0
    for i in range(ell):x=(x<<1)|(1-(r>>i&1))
    return x,1-(r>>(ell-1)&1)

def invert(p):
    if p<283:return []
    found=[]
    k=1
    while 185*(1<<k)-87<=p:
        A=284*(1<<k)-210
        C=185*(1<<k)-112
        ell=0
        while (A<<ell)**2 < 4*p:ell+=1
        L=1<<ell
        for r in range(0,L,2):
            X,delta=low(r,ell)
            B=A*L+25*X
            Y=p-C-A*r-25*(F(r)-1)
            lo=max(0,-(-(Y-50)//(B+25)));hi=(Y-25)//B
            res=Y*pow(B,-1,25)%25
            q=lo+(res-lo)%25
            if q<=hi and (L*q+r)%2==0:
                need(q+25>hi)
                if B*q+25*(F(q)+delta*(q&1))==Y:
                    found.append((L*q+r,k))
        k+=1
    return found

def main():
    p=(175,1988,1295,1295,125,125)
    need(profile(R,31)==(*p,-609) and profile(R,47)==(*p,-595))
    for v,k in [(0,1),(2,1),(6,1),(8,2),(34,2)]:
        for s in (0,1):
            a=build(v,k,s);need(strong(a) and end_profiles(a)[0][-1]==7*(count(v,k)+2*s))
    for ell in range(8):
        for r in range(1<<ell):
            X,d=low(r,ell)
            for q in (0,1,2,3,13,128):
                need(F((q<<ell)+r)==F(q)+X*q+d*(q&1)+F(r)-1)
    expected={}
    for k in range(1,10):
        for v in range(0,500):
            if v%2==0 and count(v,k)<=12000:expected.setdefault(count(v,k),[]).append((v,k))
    for p in range(12001):need(sorted(invert(p))==sorted(expected.get(p,[])))
    for v,k in [(1024,1),(99998,2),(10**5,1),(2**11,14),(70398,9)]:
        need((v,k) in invert(count(v,k)))
    for v,k,p in [(2,1,1049),(32,1,12539),(128,3,268529)]:
        need(count(v,k)==p and all(q%d for q in (p,p+2) for d in range(2,isqrt(q)+1)))
    print('prime twins',count(2,1),count(2,1)+2)
    print('Exact inverse, literal full-boundary switch, and high-height tests: PASS')
if __name__=='__main__':main()
