from functools import cache
from itertools import product

def check(x):
    if not x:
        raise AssertionError('right-fiber control')

def g(bits):
    n=P=S=0
    for u in bits:
        if not u: n+=P+S
        P,S=2*P+u,u
    return n+1

def offsets(W,i):
    x=0
    for v in W[:i]: x=2*x+v
    q=int(i+2<len(W) and not W[i+2])
    r=2*W[i]+W[i+1]
    ans=(0,2+3*x-x%2-2*q,-4-3*x+x%2+2*q,0)[r]
    bias=r if i+2==len(W) else 3*W[-1]
    return ans,bias

def independent(W,i):
    h,N,P,S=1,1,0,0
    j=0
    while j<len(W):
        if j==i:
            r=2*W[j]+W[j+1]
            K=3 if r in (0,3) else 5
            h,N,P,S=3*h,3*N+(K-3)*h+3*(3-r)*P+(3-r)*S,12*P+3*r*h,r*h
            j+=2
        else:
            v=W[j]
            if not v: N+=P+S
            P,S=2*P+h*v,h*v
            j+=1
    return h,N,P,S

def build(W,i):
    layers=[(1,v) for v in W[:i]]+[(3,2*W[i]+W[i+1])]+[(1,v) for v in W[i+2:]]
    n=sum(z[0] for z in layers); a=[]; u=0
    for size,r in layers:
        b=len(a)
        for j in range(size):
            a.append(((1<<n)-(1<<(b+size))) | ((1<<(b+(j+1)%3)) if size==3 else 0))
        u|=((1<<r)-1)<<b
    return tuple(z|((1<<n) if not (u>>j)&1 else 0) for j,z in enumerate(a))+(u,)

def H(T):
    n=len(T); full=(1<<n)-1
    @cache
    def starts(s,i):
        x=s^(1<<i)
        if not x: return 1
        t=x&T[i]; z=0
        while t:
            j=(t&-t).bit_length()-1
            t^=1<<j; z+=starts(x,j)
        return z
    return sum(starts(full,i) for i in range(n))

def find(W,a,bias):
    n=len(W)
    lead=next((i for i,v in enumerate(W) if v),n)
    if a==0:
        return [i for i in range(n-1) if offsets(W,i)==(0,bias)]
    target=(0,1) if a>0 else (1,0)
    hits=set()
    for p,q in product((0,1),repeat=2):
        if a>0:
            v=a-2+p+2*q
        else:
            v=-a-4+p+2*q
        if v<0 or v%3: continue
        x=v//3
        if x%2!=p: continue
        if x:
            i=lead+x.bit_length()
        else:
            i=lead-int(target==(0,1))
        if not 0<=i<n-1 or W[i:i+2]!=target: continue
        if offsets(W,i)==(a,bias): hits.add(i)
    return sorted(hits)

literal=0
for w in range(2,12):
    for W in product((0,1),repeat=w):
        G=g(W)*3
        fibers={}
        for i in range(w-1):
            a,b=offsets(W,i)
            h,N,P,S=independent(W,i)
            check((h,N,S)==(3,G+a,b))
            fibers.setdefault(b,{}).setdefault(a,[]).append(i)
            if w<=7:
                check(H(build(W,i))==N)
                literal+=1
        for b,vals in fibers.items():
            for a,ii in vals.items():
                check(find(W,a,b)==ii)
                if a+2 in vals: check(a==0 or a==-2)
                check(not (a+2 in vals and a+4 in vals))
for w in range(2,50):
    W=tuple(1-i%2 for i in range(w-1)); W+=(W[-1],)
    V=[offsets(W,i) for i in range(w-1)]
    check(len({v[1] for v in V})==1 and len(set(V))==w-1)
    W=(0,)*(w-3)+(1,0,0) if w>=3 else (0,0)
    check(len({offsets(W,i)[0] for i in range(w-1)})==(2 if w>=3 else 1))
print(literal,'literal graphs, widths up to eleven of offset inverses, and extremal widths through 49.')
for w in range(2,10):
    for W in product((0,1),repeat=w):
        lead=next((i for i,v in enumerate(W) if v),w)
        N=(0,0)+W[lead:]
        A={offsets(W,i) for i in range(len(W)-1)}
        B={offsets(N,i) for i in range(len(N)-1)}
        check(A <= B and g(W)==g(N))
        if lead>=2: check(A==B)
        for a,b in A:
            if (a+2,b) in A: check(b%3==0)
print('Width-free right-fiber normalization and mod-three controls pass.')
W=(1,0,1,0,1,0,0)
base=tuple(sorted(H(build(W,i)) for i in range(len(W)-1)))
check(base==(187,209,215,219,221,233))
check(tuple(sorted(H(build(W+(0,1,0),i)) for i in range(len(W)-1)))==
      tuple(q+1266 for q in base))
print('Width-seven shared right continuation is exact.')
