"""Literal proof certificate for the 11-13 bridge and the full-ready line obstruction."""
from functools import cache
from itertools import combinations, permutations
from math import isqrt

def h_table(A):
    n=len(A)
    Z=1<<n
    E=[[0]*n for _ in range(Z)]
    for i in range(n): E[1<<i][i]=1
    for s in range(1,Z):
        for i,x in enumerate(E[s]):
            if not x: continue
            choices=A[i]&~s
            while choices:
                b=choices&-choices
                choices-=b
                E[s|b][b.bit_length()-1]+=x
    return [1]+[sum(row) for row in E[1:]]

def rooted(A,r):
    h=h_table(A)
    full=(1<<len(A))-1
    @cache
    def starts(S,v):
        S^=1<<v
        if not S:return 1
        m=A[v]&S
        ret=0
        while m:
            br=m&-m
            m-=br
            ret+=starts(S,br.bit_length()-1)
        return ret
    b,c=0,0
    for s in range(1,full+1):
        if not s>>r&1:continue
        z=starts(s,r)*h[full^s]
        b+=z
        c+=(-1)**s.bit_count()*z
    return h[-1],b,c

def bridge(A,x,r):
    n=len(A)
    m=n+r+1
    later=((1<<m)-1)^((1<<n)-1)
    out=[s|later for s in A]
    out[x]^=1<<(m-1)
    for i in range(n,m-1):
        out.append(((1<<m)-1)^((1<<(i+1))-1))
    return tuple(out+[1<<x])

def extends(A,d):
    return tuple(r|(((d>>i&1)^1)<<len(A)) for i,r in enumerate(A))+(d,)

def canon(A):
    groups={}
    for i,r in enumerate(A):
        groups.setdefault(r.bit_count(),[]).append(i)
    def codes(remaining, run=()):
        if not remaining:
            yield tuple(sum(((A[x]>>y)&1)<<j for j,y in enumerate(run)) for x in run)
        else:
            for p in permutations(remaining[0]):
                yield from codes(remaining[1:],run+p)
    return min(codes(tuple(groups[d] for d in sorted(groups))))

def strong(A):
    mask=(1<<len(A))-1
    for v in range(len(A)):
        reach=1<<v
        while 1:
            r=reach
            for i,a in enumerate(A):
                if reach>>i&1:r|=a
            if r==reach:break
            reach=r
        if reach!=mask:return False
    return True

def prime(n):
    return n>1 and all(n%u for u in range(2,isqrt(n)+1))

def outer(A,k):
    for _ in range(k):
        A=bridge(A,len(A)-1,0)
        A=bridge(A,len(A)-1,1)
    for _ in range(7):
        A=bridge(A,len(A)-1,0)
    return bridge(A,len(A)-1,2)

def report59():
    """Return every H=59 strong core and its complete per-root signed profile."""
    classes={(0,)}
    expected=(1,1,2,4,12,56,456,6880)
    hits={}
    for n in range(1,9):
        if n>1:
            classes={canon(extends(A,d)) for A in classes for d in range(1<<(n-1))}
        assert len(classes)==expected[n-1]
        for A in classes:
            if strong(A) and h_table(A)[-1]==59:
                hits[A]=tuple(rooted(A,i)[1:] for i in range(n))
    return hits

if __name__=='__main__':
    p=25733
    assert prime(p) and 159**2<=p<907**2
    points=[(a,(p-47*a)//56) for a in range(1,p//159+1,2)
            if (p-47*a)%112==0]
    assert points==[(59,410)]
    N=7*p
    terminal=[]
    for r in range(1,10):
        y=1<<r
        u,v=57+68*y,68+81*y
        for a in range(1,N//(u+2*v)+1,2):
            b,rem=divmod(N-u*a,v)
            if rem==0 and not b%2:
                assert b>=2*a
                terminal.append((y,a,b))
    assert terminal==[(4,59,410)]
    hit=report59()
    assert hit and min(map(len,hit))>=6
    assert not any(R in (205,410) for profiles in hit.values() for R,_ in profiles)
    print('H=59 cores:',len(hit),'by order:',{n:sum(len(A)==n for A in hit) for n in (6,7,8)},flush=True)
    print('25733 has no ready source, unrestricted over rooted tournaments',flush=True)

    T=(2,20,9,3,13)
    assert rooted(T,1)[:2]==(11,13)
    A=bridge(T,1,1)
    assert A==(98,52,105,99,109,64,2)
    assert rooted(A,6)==(37,37,-11)
    B=outer(A,0)
    assert prime(1973) and len(B)==17 and h_table(B)[-1]==7*1973
    A=bridge(T,1,7)
    assert rooted(A,len(A)-1)==(1675,1675,-11)
    assert prime(93701) and 159**2<=93701<907**2
    assert 47*11+56*1664==93701
    assert 11==93701*pow(47,-1,112)%112
    print('11,13,26; ordinary 7·1973 by literal path DP, and maximal-depth cofactor 93701')

    target=7*25733
    all_widths=[]
    r=1
    while 193+230*(1<<r)<=target:
        y=1<<r
        ca,cb=57+68*y,68+81*y
        for a in range(1,target//(ca+2*cb)+1,2):
            b,residue=divmod(target-ca*a,cb)
            if residue==0 and b%2==0 and b>=2*a:
                all_widths.append((y,a,b))
        r+=1
    assert all_widths==[(4,59,410)]
    print('All positive terminal-seven widths: sole arithmetic input (4,59,410).')
