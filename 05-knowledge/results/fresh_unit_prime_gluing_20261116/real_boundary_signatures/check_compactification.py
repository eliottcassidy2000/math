from itertools import combinations

def induced(T,order):
    return tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(order)) for i in order)

def reachable(a,s,mask):
    r=1<<s
    while True:
        z=r
        for j in range(len(a)):
            if r>>j&1: z|=a[j]&mask
        if z==r: return r
        r=z

def scc(a,verts):
    verts=list(verts)
    mask=sum(1<<i for i in verts)
    grouped={}
    for i in verts:grouped.setdefault(reachable(a,i,mask),[]).append(i)
    return [grouped[z] for z in sorted(grouped,key=int.bit_count,reverse=True)]

def Htable(a):
    n=len(a);F=(1<<n)-1
    path=[[0]*n for _ in range(1<<n)]
    H=[1]+[0]*F
    for s in range(1,1<<n):
        for j in range(n):
            if s>>j&1:
                u=s^(1<<j)
                path[s][j]=sum(path[u][i] for i in range(n) if a[j]>>i&1) if u else 1
        H[s]=sum(path[s])
    return H,path

def profile(a,r):
    h,path=Htable(a);F=(1<<len(a))-1
    E=sum(path[s][r]*h[F^s] for s in range(1<<len(a)) if s.bit_count()%2==0)
    O=sum(path[s][r]*h[F^s] for s in range(1<<len(a)) if s.bit_count()%2==1)
    return h[F],E,O,h[F^(1<<r)]

def pad(a,r,m):
    for _ in range(m):
        a=tuple(v+((i!=r)<<len(a)) for i,v in enumerate(a))+(1<<r,)
    return a

def reduce(a,r):
    blocks=scc(a,(v for v in range(len(a)) if v!=r))
    contacts=[i for i,B in enumerate(blocks) if any(a[r]>>j&1 for j in B)]
    if not contacts:return None
    l,u=contacts[0],contacts[-1]
    left=[B for B in blocks[:l] if len(B)>1]
    active=blocks[l:u+1]
    right=[B for B in blocks[u+1:] if len(B)>1]
    tail=sum(len(B)==1 for B in blocks[u+1:])
    C=induced(a,[v for B in left+active+right for v in B]+[r])
    return C,len(C)-1,tail,len(list(B for B in active if len(B)==1))

rng_count=0
for n in range(2,6):
    E0=list(combinations(range(n),2))
    for m in range(1<<len(E0)):
        a=[0]*n
        for k,(u,v) in enumerate(E0):
            if m>>k&1:a[u]|=1<<v
            else:a[v]|=1<<u
        for r in range(n):
            P=reduce(a,r)
            if P is None:continue
            C,r0,tail,ns=P
            assert profile(a,r)==profile(pad(C,r0,tail),r0)
            h,e,o,g=profile(a,r)
            K=(2*o-g-h)//2
            assert K>=0 and ns<=(3+K.bit_length()-1 if K else 2)
            rng_count+=1
print(rng_count,"literal reductions with equal full rooted profiles")

def ordinal(a,b):
    return tuple(v+(((1<<len(b))-1)<<len(a)) for v in a)+tuple(v<<len(a) for v in b)

def nearblock(a,r,b):
    n=len(a);s=((1<<len(b))-1)<<n
    return tuple(v+(s if i!=r else 0) for i,v in enumerate(a))+tuple((v<<n)|(1<<r) for v in b)

a=(2,4,1);c=(2,4,1)
B1=ordinal(ordinal((0,),c),(0,))
B2=ordinal(c,(2,0))
assert profile(nearblock(a,0,B1),0)==profile(nearblock(a,0,B2),0)
C,r,m,_=reduce(nearblock(a,0,B1),0)
assert m==3 and profile(pad(C,r,m),r)==profile(nearblock(a,0,B1),0)
print("Nontrivial suffix SCC and moved singleton factors verified.")
