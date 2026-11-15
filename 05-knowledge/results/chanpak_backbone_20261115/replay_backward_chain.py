"""Independent literal path DP + exact density certificate (Python 3)."""
from itertools import product
from fractions import Fraction

def triples(n):
    b = [(1,1,1),(1,2,2),(3,3,4)]
    for r in range(3,n+1):
        b.append(tuple(sum(b[r-j][i] for j in (1,2,3)) for i in range(3)))
    return b[:n+1]

def literal(gaps):
    anchors = [0]
    for m in gaps:
        anchors.append(anchors[-1]+1+m)
    n=anchors[-1]+1
    row = [((1<<n)-1)^((1<<(u+1))-1) for u in range(n)]
    for u,v in zip(anchors,anchors[1:]):
        row[u]^=1<<v
        row[v]^=1<<u
    a=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        a[1<<v][v]=1
    for s in range(1,1<<n):
        for u,cu in enumerate(a[s]):
            if not cu:
                continue
            nxt=row[u]&~s
            while nxt:
                bit=nxt&-nxt
                nxt^=bit
                a[s|bit][bit.bit_length()-1]+=cu
    return sum(a[-1])

def formula(gaps):
    b=triples(len(gaps))
    t=0
    xy=None
    for m in gaps:
        if m==0:
            t+=1
            continue
        a,v,w=b[t]
        if xy is None:
            xy=a,(1<<m)*v
        else:
            x,y=xy
            xy=a*x+v*y,(1<<m)*(v*x+w*y)
        t=0
    a,v,w=b[t]
    if xy is None:
        return a
    return a*xy[0]+v*xy[1]

def segmented(gaps):
    result=1
    start=0
    for i,j in enumerate(gaps[:-1]):
        if j and gaps[i+1]:
            result*=formula(gaps[start:i+1])
            start=i+1
    return result*formula(gaps[start:])

def upper(num,den,scale=1000000):
    """Return integer ceil(scale*(num/den)**(49/50)); no floats."""
    target=scale**50 * num**49
    d=den**49
    lo,hi=0,2*scale
    while lo+1<hi:
        m=(lo+hi)//2
        if m**50*d>=target:
            hi=m
        else:
            lo=m
    if not ((hi-1)**50*d<target<=hi**50*d):
        raise RuntimeError("rational power calculation")
    return hi

def main():
    cnt=0
    for length in range(7):
        for gs in product(range(4), repeat=length):
            if sum(gs)+length+1>13:
                continue
            h=literal(gs)
            if h!=formula(gs) or h!=segmented(gs):
                raise RuntimeError((gs,h,formula(gs),segmented(gs)))
            cnt+=1
    s=1000000
    up = Fraction(sum(upper(3,a+2*b+(1<<m)*(b+2*c),s)
                      for a,b,c in triples(23)[1:] for m in range(1,20)),s)
    z=Fraction(3,5)
    if 2**49*3**50<=5**50:
        raise RuntimeError("factor bound")
    tail = 2*z**12/(1-z)*z/(1-z)+2/(1-z)*z**20/(1-z)
    if not up+tail < Fraction(24,25):
        raise RuntimeError((up,tail))
    print(f"literal path DP: {cnt} words")
    print(f"finite rational upper bound: {up}; tail: {tail}")
    print(f"growth contraction certified < 24/25")

if __name__=='__main__':
    main()
