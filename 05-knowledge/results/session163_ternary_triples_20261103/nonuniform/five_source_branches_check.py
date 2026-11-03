#!/usr/bin/env python3
"""Exhaust 1-vertex induced extensions of five marked sources and all ports."""
from itertools import combinations

def need(p,msg):
 if not p:raise ValueError(msg)
def H(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for m in range(1<<n):
  for i,c in enumerate(dp[m]):
   z=rows[i]&~m
   while c and z:
    b=z&-z;dp[m|b][b.bit_length()-1]+=c;z-=b
 return sum(dp[-1])
def add(rows,m):
 n=len(rows)
 return tuple(x|((1<<n) if not(m>>i&1) else 0) for i,x in enumerate(rows))+(m,)
def state(rows,p):
 order=tuple(i for i in range(len(rows)) if i!=p)+(p,)
 r=tuple(sum(((rows[a]>>b)&1)<<j for j,b in enumerate(order)) for a in order)
 h=H(r);r1=add(r,1<<(len(r)-1));h1=H(r1);r2=add(r1,1<<(len(r1)-1));h2=H(r2)
 need((h1-h)%2==0 and (h2-h1)%2==0,'observer parity')
 return h,(h1-h)//2,(h2-h1)//2
G={'U':((112,97,115,103,74,16,32),(31,14,25)),
'U2':((252,249,194,228,204,148,32,64),(95,41,75)),
'U4':((496,497,499,487,456,400,288,64,128),(97,40,75)),
'R2':((1008,1017,1011,997,776,912,816,592,128,256),(285,123,225)),
'R4':((996,993,994,999,975,912,800,576,128,256),(291,120,225))}
expected={'U':[(1,97),(33,95),(117,95),(125,97)],'U2':[(4,289)],'U4':[],'R2':[],'R4':[]}
for name,(r,source) in G.items():
 n=len(r);need(state(r,n-1)==source,f'source state {name} {state(r,n-1)} expected {source}')
 targets={(3*source[0]+2*k,3*source[1]-k,3*source[2]) for k in (1,2)}
 target_h={v[0] for v in targets};f0=H(add(r,0));f=[H(add(r,1<<i)) for i in range(n)]
 a=[x-f0 for x in f];b={(i,j):H(add(r,(1<<i)|(1<<j)))-f[i]-f[j]+f0 for i,j in combinations(range(n),2)};hits=[]
 for mask in range(1<<n):
  D=[i for i in range(n) if mask>>i&1];h=f0+sum(a[i] for i in D)+sum(b[i,j] for i,j in combinations(D,2))
  need(h==H(add(r,mask)),f'quadratic vs literal {name} {mask}')
  if h in target_h:
   hits.append((mask,h));new=add(r,mask)
   for p in range(n+1):need(state(new,p) not in targets,f'found full branch {name} {mask} {p}')
 need(hits==expected[name],f'unexpected hits {name}: {hits}')
 print(name,'target-count hits:',hits,'no full-state hits for any port')
print('PASS: all orientations and all target-count-hit ports checked')
