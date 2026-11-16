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


def strong(a):
    F=(1<<len(a))-1
    return all(reachable(a,i,F)==F for i in range(len(a)))
data={
(167,15,-75):[
((0,1,11,19,71,31,175,63),228,668,45),
((0,1,11,35,15,87,159,63),41,167,44),
((0,1,3,135,47,79,31,119),50,334,43),
((0,1,67,23,39,15,59,127),216,668,42),
((0,1,11,19,71,31,175,63),21,167,41)],
(601,33,-315):[
((0,1,35,23,135,91,31,111),77,601,174),
((0,1,35,23,135,91,31,111),101,601,173),
((0,13,25,49,35,71,31,127),234,1202,172),
((0,17,11,19,69,31,175,63),113,601,171),
((0,1,35,39,79,211,143,31),53,601,170)]
}
for (q,g,I),rows in data.items():
  t0=max((Z//q).bit_length() for _,_,Z,_ in rows)
  for j,(D,Q,Z,K) in enumerate(rows):
    root=len(D)
    A=tuple(row+((not (Q>>i&1))<<root) for i,row in enumerate(D))+(Q,)
    assert all((A[i]>>k&1)^(A[k]>>i&1) for i in range(len(A)) for k in range(i))
    H,E,O,g0=profile(A,root)
    assert (g0,E+O-g0,H-2*O+2*g0)==(g,Z,I+2*j)
    assert H-2*O+2*g0==g-2*K
    print(q,j,A,(H,E,O,g0),end=' ')
    for t in range(t0,t0+2):
        m=t-(Z//q).bit_length()+1
        W=pad(A,root,m)
        assert strong(W)
        assert profile(W,root)==(q*(1<<t)+I+2*j,q*(1<<(t-1)),q*(1<<(t-1))+g,g)
        print('t',t,'OK',end=' ')
    print()
