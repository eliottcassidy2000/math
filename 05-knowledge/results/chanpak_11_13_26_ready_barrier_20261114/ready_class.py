from itertools import combinations
from bridge import bridge,source,rooted_profile
from bilateral import path_count,splice,value,ready_pair
from alternate_target import strong
from math import isqrt

def ends_at(T,u):
    n=len(T);S=((1<<n)-1)^(1<<u);v=[]
    # Count only full paths ending at u, by one marked DP.
    from functools import cache
    @cache
    def end(m,i):
        t=m^(1<<i)
        if not t:return 1
        return sum(end(t,j) for j in range(n) if t>>j&1 and T[j]>>i&1)
    return end((1<<n)-1,u)

if __name__=='__main__':
    for n in range(1,6):
        pairs=tuple(combinations(range(n),2))
        for mask in range(1<<len(pairs)):
            T=[0]*n
            for k,(i,j) in enumerate(pairs):
                if mask>>k&1:T[i]|=1<<j
                else:T[j]|=1<<i
            for u in range(n):
                h,p,m=rooted_profile(T,u)
                assert m==-ends_at(T,u)
                if h==p:assert T[u].bit_count()<=1
    P=bridge(source((8,1,3,6)),0,0)
    assert strong(P) and rooted_profile(P,len(P)-1)==(37,37,-5)
    assert ready_pair(P,len(P)-1)==(5,32)
    H=value(7,(5,32),(1,4))
    assert H==7*2027 and all(2027%d for d in range(2,isqrt(2027)+1))
    N=180131;out=[]
    r=1
    while 193+230*(1<<r)<=N:
        y=1<<r;u,v=57+68*y,68+81*y
        for a in range(1,N//(u+2*v)+1,2):
            b,t=divmod(N-u*a,v)
            if t==0 and b%2==0 and b>=2*a:out.append((y,a,b))
        r+=1
    assert out==[(4,59,410)]
    print('all rooted signed endpoints through order five; ready pair (5,32); unique terminal-seven arithmetic point')
