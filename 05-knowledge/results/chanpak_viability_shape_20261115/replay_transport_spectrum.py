from random import Random
from collections import Counter
from itertools import product
from replay_unit_suffix import need,p,rev,T,S,Q

P="110222011212"
R="11022"
rng=Random(529)
def delta(u,v):
    return S*v-T*u+45*(p(v)-p(u))
def encode(k,bits):
    n=12*k-7
    u,v=87<<n,89<<n
    free_iter=iter(bits)
    for j,e in enumerate(map(int,P*(k-1)+R)):
        b=(next(free_iter) if e==1 else e//2)
        v|=b<<j
        u|=(e-b)<<(n-1-j)
    return u,v

rows=0
for k in range(1,31):
    size=5*k-3
    d=12*k; A=1<<d
    f0=16*(k-2)
    weights=[-32,50,50,-32,-68]*(k-1)+[-32,50]
    for b in range(size):
        bits=[0]*size
        bits[b]=1
        need(delta(*encode(k,bits))==f0+weights[b])
    configs=list(product((0,1),repeat=size)) if size<=12 else [tuple(rng.getrandbits(1) for _ in range(size)) for _ in range(200)]
    for bits in configs:
        u,v=encode(k,bits)
        need(128*v+rev(u,d)==(A-1)*(T+45)//45)
        need(128*u+rev(v,d)==(A-1)*(S+45)//45)
        need(delta(u,v)==f0+sum(a*b for a,b in zip(weights,bits)))
        rows+=1
base_k=9
U={16*a+25*b+34*c for a in range(2*base_k) for b in range(2*base_k) for c in range(base_k)}
need(set(range(248,116*base_k-75-248+1))<=U)
apery=[min(t for b in range(16) for c in range(16) if (t:=25*b+34*c)%16==r) for r in range(16)]
# b>=16 or c>=16 can always be reduced by 16 in its own residue.
need(max(apery)-16==247)
need(sum((t-i)//16 for i,t in enumerate(apery))==124)
print("all-width displacement controls",rows)
print("mod-16 Apery minima",apery)
print("width k=9 central support certified",248,116*9-75-248)
count=0
for k in range(2,8):
    n=12*k-7
    for _ in range(300):
        bits=tuple(rng.getrandbits(1) for _ in range(5*k-3))
        u,v=encode(k,bits)
        for j in range(1,k):
            low=12*j
            a1,a2=1<<low,1<<(12*k-low)
            u1,u2=divmod(u,a2)
            v2,v1=divmod(v,a1)
            def legal(w,z,d):
                A=1<<d
                return (128*z+rev(w,d)==(A-1)*(T+45)//45
                        and 128*w+rev(z,d)==(A-1)*(S+45)//45)
            marker=(bits[5*(j-1)+2:5*j]==(0,1,0))
            need(legal(u1,v1,low)==marker)
            need(legal(u2,v2,12*k-low)==marker)
            need(((v1>>(low-7))==89)==marker)
            need(((u2>>(12*k-low-7))==87)==marker)
            count+=1
print("primitive factor and both suffix-marker controls",count)
