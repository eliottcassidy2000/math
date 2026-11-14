"""Construct and recognize physical source words in the near-boundary seven-axis family."""
from functools import cache
from math import isqrt

@cache
def C(q):
    if q<3:
        return ((1,1,1),(1,2,2),(3,3,4))[q]
    return tuple(map(sum,zip(C(q-3),C(q-2),C(q-1))))

def step(q,r,state):
    e,f,g=C(q)
    x,y=state
    return e*x+f*y,(f*x+g*y)<<r

def transfer(word):
    s=1,1<<word[0]
    q=0
    for r in word[1:]:
        if r:
            s=step(q,r,s)
            q=0
        else:
            q+=1
    assert q==0
    return s

def cone(a,b):
    return a>0 and a&1 and b%2==0 and b>=2*a

def dyadic(v):
    return v>=2 and not v&(v-1)

@cache
def decode(a,b):
    """All-seam existence, with one literal witness on success."""
    if not cone(a,b):
        return None
    if a==1:
        return (b.bit_length()-1,) if dyadic(b) else None
    d=(b-2*a)//2
    if d>0 and a>6*d:
        # All-seam geometric restriction, not a heuristic.
        inner=decode(d,(a-d)//2)
        return (inner+(0,1)) if inner is not None else None
    r=(b//a).bit_length()-1
    y=1<<r
    if b%y or 3*b>=4*y*a:
        return None
    t=b//y
    q=0
    while True:
        e,f,g=C(q)
        if e+2*f>a:
            return None
        det=e*g-f*f
        if det==0:
            if e*t==f*a:
                for x in range(1,a//(e+2*f)+1,2):
                    yy,res=divmod(a-e*x,f)
                    if not res:
                        stem=decode(x,yy)
                        if stem is not None:
                            return stem+(0,)*q+(r,)
        else:
            x,yy=g*a-f*t,e*t-f*a
            if x%det==yy%det==0:
                stem=decode(x//det,yy//det)
                if stem is not None:
                    return stem+(0,)*q+(r,)
        q+=1

def row(k):
    u,v=47,56
    for _ in range(k):
        u,v=u+4*v,2*u+4*v
    return u,v

def order_of_depth(p):
    k=0
    u,v=row(k)
    while True:
        U,V=u+4*v,2*u+4*v
        if (U+2*V)**2>p:
            return k,u,v
        k+=1
        u,v=U,V

def source_candidates(p,u,v):
    # Each dyadic band has a separate CRT stride, not the wider cone's.
    y=2
    while u+v*y<=p:
        lo=3*p//(3*u+4*v*y)+1
        hi=p//(u+v*y)
        m=v*y
        offset=p*pow(u,-1,m)%m
        for a in range(lo+(offset-lo)%m, hi+1,m):
            yield a,(p-u*a)//v
        y*=2

def target_sources(p,k):
    u,v=row(k)
    for a,b in source_candidates(p,u,v):
        w=decode(a,b)
        if w is not None:
            yield (a,b),w

def only_boundary_source(p,k,d):
    """Test exactly the declared small-defect shell d; return a physical full word."""
    u,v=row(k)
    A=u+2*v
    B,rem=divmod(p-(u+4*v)*d,2*A)
    if rem or B<2*d:
        return None
    stem=decode(d,B)
    if stem is None:
        return None
    return stem+(0,1)*(k+1)+(0,)*7+(2,)

def max_depth_family(k,stem=(1,)):
    """An actual source at maximal relaxed k, as long as d is not too large."""
    d,c=transfer(stem)
    u,v=row(k)
    A=u+2*v
    y=1
    # Multiply second coordinate by 2 via enlarging the last positive gap.
    r=0
    while (u+4*v)*d+2*A*c*y<A*A:
        r+=1
        y*=2
    p=(u+4*v)*d+2*A*c*y
    w=stem[:-1]+(stem[-1]+r,)+(0,1)
    return p,w

def adjacency(gaps):
    a=[0]
    for g in gaps:
        a.append(a[-1]+g+1)
    n=a[-1]+1
    back=set(zip(a,a[1:]))
    return tuple(sum(1<<j for j in range(n) if i!=j and ((i<j) != ((min(i,j),max(i,j)) in back)))
                 for i in range(n))

def all_paths(rows):
    @cache
    def end(S,i):
        if S==1<<i:
            return 1
        S^=1<<i
        return sum(end(S,j) for j in range(len(rows)) if S>>j&1 and rows[j]>>i&1)
    all=(1<<len(rows))-1
    return sum(end(all,j) for j in range(len(rows)))

def interval(gaps):
    P=[1]
    sums=[0]
    for g in gaps:
        sums.append(sums[-1]+g)
    for j,s in enumerate(sums):
        P.append(P[-1]+sum((1<<(s-sums[i]))*P[i] if s>sums[i]
                          else 2*P[i]*int((j-i)%2==0) for i in range(j)))
    return P[-1]

def prime(p):
    return p>=2 and all(p%q for q in range(2,isqrt(p)+1))

if __name__=='__main__':
    p=1137067
    assert prime(p)
    k,u,v=order_of_depth(p)
    assert (k,u,v)==(1,271,318)
    a=p*pow(u,-1,2*v)%(2*v)
    b=(p-u*a)//v
    assert (a,b)==(469,3176) and decode(a,b) is None
    assert (a,b) not in tuple(source_candidates(p,u,v))
    stem=(3,)+(0,)*8+(1,)
    assert transfer(stem)==(1105,2634)==(a+2*v,b-2*u)
    assert ((1105,2634),stem) in tuple(target_sources(p,k))
    literal=adjacency(stem)
    assert all_paths(literal)==3739
    full=stem+(0,1)+(0,)*7+(2,)
    assert sum(transfer(full))==interval(full)==7*p==7959469
    # Different real d stems: 1 from a single gap, 5 from a single q=1 seam,
    # 9 from initial width 2 followed by q=1.
    for stem in ((1,),(1,0,1),(2,0,1)):
        d,_=transfer(stem)
        for k in range(1,18):
            p,w=max_depth_family(k,stem)
            a,b=transfer(w)
            use,u,v=order_of_depth(p)
            assert use==k and u*a+v*b==p and a>6*d
            assert decode(a,b) is not None
            full=w+(0,1)*k+(0,)*7+(2,)
            assert only_boundary_source(p,k,d) is not None
            assert sum(transfer(full))==interval(full)==7*p
    # Replay exact threshold and converse in small ordinary words, all seams included.
    from itertools import product
    for n in range(1,5):
        for W in product(range(3),repeat=n):
            if W[0]<=0 or W[-1]<=0:
                continue
            a,b=transfer(W)
            d=(b-2*a)//2
            if a>6*d>0:
                assert W[-2:]==(0,1) and len(W)>=3
                assert transfer(W[:-2])==(d,(a-d)//2)
            assert decode(a,b) is not None
    print('literal 15-vertex prefix for prime cofactor 1137067; 51 changed-source maximal-depth constructions')
