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

def strong(a):
    F=(1<<len(a))-1
    return all(reachable(a,i,F)==F for i in range(len(a)))

pairs=[(217,[((0,1,67,23,71,31,171,63),30,23),((0,1,3,23,39,79,31,127),242,5)])]
for E,lr in pairs:
 for p,(a,m,g) in zip([41,43],lr):
  n=len(a)
  b=tuple(x+((not (m>>i&1))<<n) for i,x in enumerate(a))+(m,)
  assert all((b[i]>>j&1)^(b[j]>>i&1) for i in range(n+1) for j in range(i))
  assert all(not (b[i]>>i&1) for i in range(n+1))
  assert strong(b)
  h,o=7*p,E+g
  assert profile(b,n)==(h,E,o,g)
  assert profile(pad(b,n,1),n)==(h+2*E,2*E,2*E+g,g)
  assert profile(pad(b,n,2),n)==(h+6*E,4*E,4*E+g,g)
  table,_=Htable(b)
  assert table[367 if p==41 else 319]==p
  print(p,b,(h,E,o,g),'induced parent',[(s,tuple(i for i in range(n+1) if s>>i&1)) for s in range(1<<len(b)) if table[s]==p])

from math import isqrt
def pr(p):return p>1 and all(p%d for d in range(2,isqrt(p)+1))

def near_block(A,root,B):
    n=len(A);m=len(B)
    mask=((1<<m)-1)<<n
    return tuple((a|(mask if j!=root else 0)) for j,a in enumerate(A))+tuple((b<<n)|(1<<root) for b in B)
def uv(a):
    H,_=Htable(a);F=len(H)-1
    return H[-1],sum(H[s]*H[F^s] for s in range(len(H)))

rails=((256,1,67,23,71,287,427,319,30),(256,1,259,279,39,79,31,127,242))
for h,p,s in [(41,0,[0,1,2,3,5,6,8]),(43,1,[0,1,2,3,4,5,8])]:
    root=induced(rails[p],s)
    assert profile(root,6)[:2]==(h,27)
    assert profile(root,6)[2]-profile(root,6)[1]==profile(root,6)[3]
for B,u,v,targets in [
   ((12,9,2,4),5,32,(1879,1889)),
   ((26,60,25,32,8,21),23,272,(16381,16427))]:
    assert uv(B)==(u,v)
    assert targets==(62*v-21*u,62*v-19*u)
    for a,h in zip(rails,targets):
        assert pr(h) and Htable(near_block(a,8,B))[0][-1]==7*h
    print('Near-block u,v and new prime factors',u,v,targets)

B=(16,17,3,71,12,31,55)
assert uv(B)==(35,676)
p,q=432*676-31*5,432*676-29*5
assert (p,q)==(291877,291887) and pr(p) and pr(q)
assert pr(227) and pr(229)
print('Same-profile old 756 reader: sharp prime gap',p,q)

rails30=((256,273,35,279,357,395,303,351,4),(256,269,281,97,299,327,407,319,8))
for A,p,g in zip(rails30,(41,43),(77,91)):
    assert strong(A)
    assert profile(A,8)==(7*p,105,105+g,g)
    H,_=Htable(A)
    parents=[s for s,h in enumerate(H) if h==p and s.bit_count()==7]
    assert parents
    assert profile(pad(A,8,1),8)==(7*(p+30),210,210+g,g)
    print('30 rail',p,'parent masks',parents)
assert all(pr(p) for p in (71,73,1931,1933,52859,52957))
B=(2,20,17,7,105,15,47)
assert uv(B)==(49,872)
for A,I,p in zip(rails30,(11,13),(52859,52957)):
    assert 60*872+I*49==p
    # smaller literal replay; add the extra neutral singleton by the transfer law.
    assert Htable(near_block(A,8,B))[0][-1]==7*(30*872+I*49)

for A,r,E,g in [((16,5,17,7,74,31,175,63),4,38,1),
                ((2,16,3,23,165,79,31,111),1,26,11)]:
    assert profile(A,r)==(63,E,E+g,g)
    assert profile(pad(A,r,3),r)[0]==63+14*E
    cyclic=[(9+2*E*((8**j-1)//7))%49 for j in range(49)]
    assert len(set(cyclic))==49 and {x%3 for x in [9+2*E*((8**j-1)//7) for j in range(6)]}=={0,1}
assert pr(61)
