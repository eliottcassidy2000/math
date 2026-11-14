"""Independent controls for the generic Syracuse ordinal compiler."""
from itertools import permutations, product
from random import Random

def ensure(c,*ctx):
    if not c: raise AssertionError(ctx)

def word_perm(a):
    q=[0]
    for v in a:q.append(q[-1]+v)
    A=q[-1]
    return tuple(sorted(range(len(q)),key=lambda i:3**i*2**(A-q[i])))

def compile_perm(p):
    m=len(p)
    b=[(3**i).bit_length()-1 for i in range(m)]
    num=[3**i*2**(b[-1]-b[i]) for i in range(m)]
    h=[0]*m
    for i,j in zip(p,p[1:]):
        h[j]=h[i]+(num[i]>num[j])
    a=tuple(b[i]-b[i-1]-h[i]+h[i-1] for i in range(1,m))
    return (a,h) if min(a,default=1)>=1 else None

def seed(a,sig):
    A=sum(a); x=1
    n=1
    for v in reversed(a):
        n+=v; q=2**n
        x=(2**v*x-sig)*pow(3,-1,q)%q
    B=2**A*3**(len(a)+1)
    return B+1+(x-B-1)%(2**(A+1))

def read(x,k,sig):
    v=[]; a=[]
    for i in range(k+1):
        v.append(x)
        if i<k:
            y=3*x+sig; div=y&-y
            a.append(div.bit_length()-1); x=y//div
    ensure(len(set(v))==len(v),v)
    return tuple(a),tuple(sorted(range(len(v)),key=v.__getitem__))

def strong(p):
    r={v:i for i,v in enumerate(p)}
    m=len(p)
    # A cut must be a chronological cut because the time edges form a path.
    return all(any(j-i>1 and r[i]<r[j]
        for i in range(k) for j in range(k,m)) for k in range(1,m))

def check(p,c):
    a,h=c
    k=len(a); b=(3**k).bit_length()-1
    ensure(word_perm(a)==p,p,a)
    ensure(sum(a)<=b+k,a)
    if strong(p):
        ensure(max(h)-min(h)<=b-k and sum(a)<=2*b-k,a,h)
    for sig in [-1,1]:
        ensure(read(seed(a,sig),k,sig)==(a,p),p,a,sig)

def main():
    counts=[]
    for m in range(1,8):
        # A cap of floor(log2(3^(m-1)))+1 preserves all comparisons.
        cap=(3**(m-1)).bit_length()
        allowed={word_perm(a) for a in product(range(1,cap+1),repeat=m-1)} if m<=6 else None
        mine=set()
        for p in permutations(range(m)):
            c=compile_perm(p)
            if c: check(p,c); mine.add(p)
        if allowed is not None: ensure(mine==allowed,m,len(mine),len(allowed))
        counts.append(len(mine))
    rng=Random(8252)
    for m in [8,16,32,64,100,250]:
        for _ in range(20):
            a=[rng.randrange(1,5*m) if rng.randrange(8)==0 else rng.randrange(1,4) for i in range(m-1)]
            p=word_perm(a)
            check(p,compile_perm(p))
    print("generic ordinal patterns orders 1..7:",counts)
    print("plus/minus literal seed, strong-budget controls and random long clocks: PASS")

if __name__=="__main__":main()
