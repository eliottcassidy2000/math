from math import gcd
from random import Random

# The action e=1 or 2 is T_k with 2^(k+1)=e modulo m.
def word(m,v,h,V,H):
    d=(m-1).bit_length()
    ans=[]
    def apply(e):
        nonlocal v,h
        h=(h+v+1)%m
        v=(e*v-1)%m
        ans.append(e)
    def steer(W):
        i=((v-1)*(1<<d)-(W-1))%m
        for j in range(d-1,-1,-1):
            apply(2)
            if i>>j&1:apply(1)
    steer(1)
    x,y=v,h
    sz=len(ans)
    steer(4)
    for _ in range(d):
        apply(1);apply(2)
    steer(V)
    h0=h
    del ans[sz:]
    v,h=x,y
    c=gcd(3,m)
    r=next(r for r in range(c) if (H-h0-8*d-2*r)%c==0)
    q=m//c
    z=(H-h0-8*d-2*r)//c * (pow(3//c,-1,q) if q>1 else 0) %q
    for _ in range(r+4*(d-z.bit_count())):apply(2)
    steer(4)
    for j in range(d):
        for e in ((2,1,1) if z>>j&1 else (1,2)):apply(e)
    steer(V)
    if (v,h)!=(V%m,H%m):raise AssertionError((m,v,h,V,H))
    return ans

def lift(m,a,v,h,V,H):
    E={}
    for e in (1,2):
        k=max(1,a)
        while pow(2,k+1,m)!=e%m:k+=1
        E[e]=k
    return [E[e] for e in word(m,v,h,V,H)]

if __name__=='__main__':
    rng=Random(8447)
    for m in list(range(3,402,2))+[511,4095,65535]:
        for _ in range(15):
            v,h,V,H=[rng.randrange(m) for _ in range(4)]
            u,g=v,h
            for k in lift(m,7,v,h,V,H):
                g=(g+u+1)%m;u=(pow(2,k+1,m)*u-1)%m
            if (u,g)!=(V,H):raise AssertionError
    print('logarithmic phase words checked')
