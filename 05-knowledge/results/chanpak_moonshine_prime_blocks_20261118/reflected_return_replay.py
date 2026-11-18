from itertools import combinations
D=(54,32,2,119,6,20,55);mask=26

def tables(d):
 n=len(d);E=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
 for v in range(n): E[1<<v][v]=1
 for m in range(1,1<<n):
  for v in range(n):
   if m>>v&1:
    if m!=(1<<v):E[m][v]=sum(E[m^(1<<v)][u] for u in range(n) if d[u]>>v&1)
  h[m]=sum(E[m])
 return h,E

def profile(d,mask):
 h,e=tables(d);_,s=tables(tuple(sum((r>>i&1)<<j for j,r in enumerate(d)) for i in range(len(d))))
 f=len(h)-1;L=[0,0];R=[0,0]
 for v in range(f+1):
  j=v.bit_count()%2
  R[j]+=h[f^v]*sum(c for i,c in enumerate(e[v]) if not mask>>i&1)
  L[j]+=h[f^v]*sum(c for i,c in enumerate(s[v]) if mask>>i&1)
 a=tuple(row|(((mask>>j&1)^1)<<len(d)) for j,row in enumerate(d))+(mask,)
 return(h[-1],sum(h[i]*h[f^i] for i in range(f+1)),sum(L),sum(R),L[0]-L[1],R[0]-R[1],tables(a)[0][-1])

base=profile(D,mask)
if base!=(5,256,142,142,-2,0,91):raise RuntimeError(base)
records=[]
edges=list(combinations(range(8),2))
for k in (1,2,3):
 hits=[]
 for es in combinations(edges,k):
  d=list(D);m=mask
  for u,v in es:
   if v==7:m^=1<<u
   else:d[u]^=1<<v;d[v]^=1<<u
  p=profile(d,m)
  if p[2]==p[3] and p[5]-p[4]==2:
   if p!=base:raise RuntimeError((es,p))
   hits.append(es)
 records.append(hits)
 print('preserving edits',k,hits)
if list(map(len,records))!=[1,1,4]:raise RuntimeError(records)
Ts=[[252,189,504,128,8,24,58,112,251],[244,253,504,129,8,24,56,112,251]]

def rooted(T,v):
 vs=[i for i in range(len(T)) if i!=v]
 d=[sum((T[i]>>j&1)<<k for k,j in enumerate(vs)) for i in vs]
 m=sum((T[v]>>j&1)<<k for k,j in enumerate(vs))
 return d,m

def F(x):
 a=0;z=1;p=0
 for b in bin(x)[2:]:
  b=int(b)
  if b==0:z+=a+p
  a,p=2*a+b,b
 return z

def read(p,u,w):
 h,P,L,R,l,r,c=p
 return c+h*(F(u)+F(w)-2)+R*u+L*w+P*u*w-r*(u%2)-l*(w%2)

def make(T,v,u,w):
 d,m=rooted(T,v);s=bin(u)[2:];t=''.join(str(1-int(b)) for b in reversed(bin(w)[2:]));n=len(s)+len(d)+len(t)
 out=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if len(s)<=i<j<len(s)+len(d):winner=i if d[i-len(s)]>>(j-len(s))&1 else j
   else:winner=i
   loser=j if winner==i else i
   out[winner]|=1<<loser
 tags=sum(int(b)<<i for i,b in enumerate(s))|(m<<len(s))|sum(int(b)<<(len(s)+len(d)+i) for i,b in enumerate(t))
 return [r|(((tags>>i&1)^1)<<n) for i,r in enumerate(out)]+[tags]
ps=[]
for T in Ts:
 profiles=[profile(*rooted(T,v)) for v in range(9)]
 if any(p[2]==p[3] for p in profiles):raise RuntimeError(profiles)
 ps.append(profiles[8])
 print('all rerooted profiles',profiles)
for u,w in [(1,2),(2,1),(3,2),(2,3)]:
 hs=[tables(make(T,8,u,w))[0][-1] for T in Ts]
 if hs!=[read(p,u,w) for p in ps]:raise RuntimeError((u,w,hs))
 print('addresses',u,w,'counts',hs,'source-change',hs[1]-hs[0])
print('PASS: 3682 core edits, 18 reroots, 8 literal flank constructions')
