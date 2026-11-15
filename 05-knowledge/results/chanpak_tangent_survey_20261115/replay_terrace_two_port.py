"""Independent literal subset checks of exact two-terminal terrace inversion."""
from itertools import product
from math import isqrt

I=((1,0),(0,1))
def mul(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in (0,1)) for j in (0,1)) for i in (0,1))

def corridors(n):
    vals=[(1,1,1),(1,2,2),(3,3,4)]
    while len(vals)<=n:
        vals.append(tuple(sum(x[i] for x in vals[-3:]) for i in range(3)))
    return [((a,b),(b,d)) for a,b,d in vals[:n+1]]

def profile(q, c):
    """q zero-corridor lengths; c exponents of the fixed interior positives."""
    assert len(q)==len(c)+1 and all(z>=0 for z in q)
    mats=corridors(max(q))
    p=mats[q[0]]
    for z,r in zip(c,q[1:]):
        p=mul(mul(p,((1,0),(0,2**z))),mats[r])
    # rows = selection of left bit, columns = selection of right bit.
    return p[0][0], p[1][0], p[0][1], p[1][1]

def evaluate(m,a,b):
    A,B,C,D=m
    return A + (B<<a) + (C<<b) + (D<<(a+b))

def inverse(m,n):
    A,B,C,D=m
    if n<evaluate(m,1,1): return set()
    if n==evaluate(m,1,1): return {(1,1)}
    k=(n//D).bit_length()-1
    if k<3: return set()
    E=n-A-(D<<k)
    delta=E*E-4*B*C*(1<<k)
    if delta<0: return set()
    r=isqrt(delta)
    if r*r!=delta:return set()
    roots=set()
    for s in set((E+r,E-r)):
        x,d=divmod(s,2*B)
        if d or x<2 or x&(x-1):continue
        a=x.bit_length()-1
        b=k-a
        if b>=1 and evaluate(m,a,b)==n: roots.add((a,b))
    return roots

def literal(gaps):
    anchors=[0]
    for g in gaps: anchors.append(anchors[-1]+g+1)
    n=anchors[-1]+1
    adj=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    for i,j in zip(anchors,anchors[1:]):
        adj[i]^=1<<j
        adj[j]^=1<<i
    vals=[[0]*n for _ in range(1<<n)]
    for i in range(n):vals[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,k in enumerate(vals[mask]):
            if not k:continue
            step=adj[v]&~mask
            while step:
                bit=step&-step
                vals[mask|bit][bit.bit_length()-1]+=k
                step-=bit
    return sum(vals[-1])

def check():
    literal_cases=0
    for r in (0,1,2):
        for q in product(range(3),repeat=r+1):
            for c in product((1,2),repeat=r):
                m=profile(q,c)
                A,B,C,D=m
                assert 0<A<=B<=D and 0<A<=C<=D
                for a,b in product(range(1,7),repeat=2):
                    n=evaluate(m,a,b)
                    assert inverse(m,n)=={(u,v) for u,v in product(range(1,12),repeat=2) if evaluate(m,u,v)==n}
                    gaps=[a]
                    for z,t in zip(q[:-1],c):gaps.extend([0]*z+[t])
                    gaps.extend([0]*q[-1]+[b])
                    if sum(gaps)+len(gaps)+1<=12:
                        assert literal(gaps)==n
                        literal_cases+=1
    for c in range(1,33):
        w=1<<c
        assert profile((1,2),(c,))==(3+6*w,6+6*w,3+8*w,6+8*w)
    print("literal subset controls",literal_cases,"and exact inverse controls passed")
if __name__=="__main__":check()
