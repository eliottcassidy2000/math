from itertools import permutations, combinations
def dp(r,root=None):
 n=len(r);N=1<<n
 D=[[0]*n for _ in range(N)]
 for i in range(n) if root is None else (root,):D[1<<i][i]=1
 for s,row in enumerate(D):
  for i,h in enumerate(row):
   if not h:continue
   e=r[i]&~s
   while e:
    b=e&-e;e^=b
    D[s|b][b.bit_length()-1]+=h
 return [1]+[sum(row) for row in D[1:]]
def prof(r,root):
 h=dp(r);k=dp(r,root);N=len(h)-1
 q=[sum(k[s]*h[N^s] for s in range(1,N+1) if s.bit_count()%2==i) for i in (0,1)]
 return h[N],h[N^(1<<root)],*q
def adjoin(r,x):
 n=len(r)
 return tuple(rr|((not (x>>i&1))<<n) for i,rr in enumerate(r))+(x,)
def append(r,root):return adjoin(r,1<<root)
def join(a,b):
 n,m=len(a),len(b)
 return tuple(r|(((1<<m)-1)<<n) for r in a)+tuple(r<<n for r in b)
def reorder(r,vs):return tuple(sum(((r[a]>>b)&1)<<i for i,b in enumerate(vs)) for a in vs)
def flips(a,b):
 return tuple((i,j) for i,j in combinations(range(len(a)),2) if (a[i]^b[i])>>j&1)
def centered(n,k):
 L=1<<k
 h=2*(n//(2*L))+1
 d=n-L*h
 assert abs(d)<L
 return h,d

def main():
 base=(34,36,41,51,39,0)
 assert [prof(adjoin(base,x),6)[0] for x in (2,3,6,7)]==[27,41,43,53]
 a,b=[adjoin(base,x) for x in (3,6)]
 assert flips(a,b)==((0,6),(2,6))
 for k in range(6):
  assert prof(a,6)==(1+40*2**k,9,20*2**k,9+20*2**k)
  assert prof(b,6)==(3+40*2**k,9,20*2**k,9+20*2**k)
  assert flips(a,b)==((0,6),(2,6))
  va=[u for u in range(len(a)) if u!=3]
  vb=[u for u in range(len(b)) if u!=1]
  p,q=reorder(a,va),reorder(b,vb)
  assert prof(p,5)==(1+12*2**k,3,6*2**k,3+6*2**k)
  assert prof(q,5)==(3+12*2**k,3,6*2**k,3+6*2**k)
  assert flips(p,q)==((0,1),(0,5),(2,3))
  a,b=append(a,6),append(b,6)
 t=(12,21,24,18,1)
 a=join((2,4,1),t)
 b=(216,253,185,32,136,209,156,8)
 edits=((0,1),(0,5),(2,6),(3,6),(4,5))
 for k in range(5):
  assert prof(a,3)==(3+36*2**k,9,18*2**k,9+18*2**k)
  assert prof(b,3)==(5+36*2**k,5,18*2**k,5+18*2**k)
  assert flips(a,b)==edits
  a,b=append(a,3),append(b,3)
 for k in range(1,10):
  L=1<<k
  fs=[0,0,0]
  for d in range(1-L,L,2):
   h,d1=centered(3*(L*7+d),k)
   c=(h-21)//2
   assert c in (-1,0,1) and d1==3*d-2*L*c
   fs[c+1]+=1
  u=(L-(-1)**k)//3;v=(L+2*(-1)**k)//3
  assert fs==[u,v,u]
  for c in (1,3,5,7,11):
   for P in range(L*(7*c-1)+1,L*(7*c+1),2):
    H=2*(P//(2*c))+1
    assert centered(H,k)[0]==7
    assert abs(P-c*H)<=c-1 and (P-c*H)%2==0
 # Direct 3h+2 carry when z has a unique Hamiltonian start.
 C0=reorder((10,4,9,18,7),range(4))
 z=3
 rev=tuple(((1<<4)-1)^r^(1<<i) for i,r in enumerate(C0))
 assert dp(C0,z)[-1]==1
 b=join((0,),(10,4,9,18,7))
 new=adjoin(b,((1<<6)-1)^(1<<5))
 assert dp(b)[-1]==13 and dp(new)[-1]==41
 print('Verified rooted rows 13/15 → 41/43 and fixed old-arc repairs through five tail appends.')
 print('Verified scaled 39/41 and exact nearest-source carries through depth 9.')
if __name__=='__main__':main()
