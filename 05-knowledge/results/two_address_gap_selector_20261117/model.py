from itertools import product,permutations
from fractions import Fraction

def step(v,i,o):
    # state index = 2 L status + R occupancy; i/out from predecessor port encoded by bit.
    w=[0]*6
    for j,a in enumerate(v):
        l,r=divmod(j,2)
        w[2*(1 if i else 2)+r]+=a
        if r or o:w[2*l+1]+=a
    return tuple(w)
def read(v): return v[0]+v[1]+v[4]+v[5]
def shared(v,t):return step(v,t,t)
def transform(v,t):return tuple(read(shared(tuple(int(i==j) for i in range(6)),t)) for j in range(6))
def dual(r,t):return tuple(sum(a*b for a,b in zip(r,shared(tuple(int(i==j) for i in range(6)),t))) for j in range(6))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def law():
    # Difference functional read(w v)-read(v) should vanish on all words.
    R=(1,1,0,0,1,1)
    basis={}
    todo=[tuple(a-b for a,b in zip(dual(R,t),R)) for t in (0,1)]
    while todo:
        f=todo.pop()
        w=list(map(Fraction,f))
        for k in sorted(basis):
            z=w[k]
            w=[x-z*y for x,y in zip(w,basis[k])]
        if any(w):
            k=next(i for i,x in enumerate(w) if x)
            w=tuple(x/w[k] for x in w)
            basis[k]=w
            todo.extend(dual(w,t) for t in (0,1))
    return tuple(basis.values())
def advance(v,d,plus):
    return step(v,(d>>(plus))&1,(d>>(1-plus))&1)
def triangle(v,plus):
    # ports' outgoing flags to b,c,y
    flag=(3,2,1)
    w=[0]*6
    for p in permutations(range(3)):
        for k in range(4):
            if any((p[j]-p[j-1])%3!=1 for j in (1,2) if j!=k):continue
            for j,a in enumerate(v):
                if not a: continue
                l,r=divmod(j,2)
                if k:l=1 if (flag[p[k-1]]>>plus)&1 else 2
                if k<3:
                    if not r and not (flag[p[k]]>>(1-plus))&1:continue
                    r=1
                w[2*l+r]+=a
    return tuple(w)
def difference(prefix,suffix):
    res=[]
    for p in (0,1):
        v=(1,0,0,0,0,0)
        for c in prefix:v=advance(v,c,p)
        v=triangle(v,p)
        for c in suffix:v=advance(v,c,p)
        res.append(v)
    return tuple(a-b for a,b in zip(res[1],res[0]))
if __name__=='__main__':
    b=law();print(b)
    for pref,suf in (((2,1),()),((1,),(1,3,2))):
        d=difference(pref,suf);print(pref,suf,d,read(d),[(t,read(shared(d,t)))for t in (0,1)])
    for n in range(2,8):
        hist=[0,0]
        for w in product(range(4),repeat=n):
            if not any(x in (1,2) for x in w):continue
            for k in range(n+1):
                d=difference(w[:k],w[k:])
                if abs(read(d))==2 and all(not dot(r,d) for r in b):
                    hist[0]+=1
                    if hist[0]==1:print('witness',w[:k],w[k:],read(d),d)
                if abs(read(d))==2:hist[1]+=1
        print(n,hist,flush=True)
