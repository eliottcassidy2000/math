import itertools
O=((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1))
Z=((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1))
J=((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
def mul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)) for i in range(4))
D={'':I}; mats={O:'O',Z:'Z',J:'J'};seen={I:''}; groups={}
for k in range(11):
 hits=[]
 for w,m in D.items():
  a=mul(Q,m);key=tuple(a[i][j] for i in range(4) for j in range(4) if (i,j)!=(3,3))
  if key in groups:
   w0,a0=groups[key]
   if a0!=a[3][3]:hits.append((w0,w,a0,a[3][3]))
  else:groups[key]=(w,a[3][3])
 print(k,len(D),'delta_constants',hits[:10],flush=True)
 nxt={}
 for w,m in D.items():
  for t,letter in mats.items():
   n=mul(t,m)
   if n not in seen:
    seen[n]=w+letter;nxt[w+letter]=n
 D=nxt
