#!/usr/bin/env python3
import itertools,json
from pathlib import Path
D=Path(__file__).resolve().parent
L=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def check(x,s):
 if not x:raise RuntimeError(s)
def arrays(n):
 sq={a*a%n for a in range(n)}
 for c,a,b in itertools.product(sq,repeat=3):
  u=(a-c)%n;v=(b-c)%n
  g=tuple((c+s*u+t*v)%n for s,t in [(1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0)])
  if len(set(g))==9 and set(g)<=sq:
   check(len({sum(g[i] for i in line)%n for line in L})==1,'direct magic');yield g
rows=[]
for p,K in [(2,9),(3,5),(5,3),(7,3),(29,1)]:
 counts={-1:0,0:0}
 for k in range(1,K+1):
  n=p**k;N=P=U=0;desc=set()
  for g in arrays(n):
   N+=1;primitive=any(x%p for x in g);P+=primitive;U+=all(x%p for x in g)
   if not primitive:
    check(k>=2 and all(x%(p*p)==0 for x in g),'nonunit square divisibility');h=tuple(x//(p*p) for x in g);check(h not in desc,'descent injectivity');desc.add(h)
  check(N==P+counts[k-2],'content recurrence')
  if k>=2:check(desc==set(arrays(p**(k-2))),'descent onto lower universe')
  if p in (2,3):check(P==U,'special all-unit branch')
  counts[k]=N;rows.append({'p':p,'k':k,'N':N,'P':P,'U':U,'descent':len(desc)})
perms=[]
for reflected in [False,True]:
 for rotations in range(4):
  perm=[]
  for r,c in itertools.product(range(3),repeat=2):
   if reflected:c=2-c
   for _ in range(rotations):r,c=c,2-r
   perm.append(3*r+c)
  perms.append(perm)
check(len({tuple(p) for p in perms})==8,'D4 permutations')
orbits={}
for n in [35,55]:
 S=set(arrays(n));O=set()
 for g in S:
  orbit={tuple(g[i] for i in perm) for perm in perms};check(len(orbit)==8 and orbit<=S,'actual D4 free orbit');O.add(min(orbit))
 check(len(O)=={35:3,55:5}[n],'explicit orbit count');orbits[n]=len(O)
result=json.dumps({'status':'PASS','rings':rows,'actual_D4_orbits':orbits},sort_keys=True,indent=2)+'\n';path=D/'ring-results.json'
if path.exists():check(path.read_text()==result,'normal optimized equality')
else:path.write_text(result)
print('PASS: direct ring recurrence, bijective descent, primitive/all-unit and actual D4 orbits')
