from itertools import product
import json
from pathlib import Path

def lift(m,edges,b):
 T=[0]*(2*m)
 for i in range(m):T[2*i]|=1<<(2*i+1)
 for c,w in edges:
  for si in (1,-1):
   for sj in (1,-1):
    u=2*c+(si==-1);v=2*w+(sj==-1)
    if si==-b[c,w]:T[u]|=1<<v
    else:T[v]|=1<<u
 return T

def fixed(T):
 m=len(T)//2;D=[[0]*(2*m) for _ in range(1<<m)]
 for v in range(2*m):D[1<<(v//2)][v]=1
 for S in range(1,1<<m):
  available=[w for w in range(2*m) if not S>>(w//2)&1]
  for v,a in enumerate(D[S]):
   if a:
    for w in available:
     if T[v]>>w&1:D[S|1<<(w//2)][w]+=a
 return sum(D[-1][::2])

rows=[]
for name,diffs in [('R7',[1,2,3]),('P7',[1,2,4])]:
 m=8;edges=[(i,(i+d)%7) for i in range(7) for d in diffs]+[(i,7) for i in range(7)];best=None
 for signs in product((-1,1),repeat=4):
  b={(c,w):(signs[3] if w==7 else signs[diffs.index((w-c)%7)]) for c,w in edges}
  T=lift(m,edges,b);f=fixed(T);Bs=[]
  for c in range(m):
   bb={e:(-x if e[0]==c else x) for e,x in b.items()};Bs.append((f-fixed(lift(m,edges,bb)))//2)
  A=f-sum(Bs);opt=A+sum(map(abs,Bs));optb={e:(-x if Bs[e[0]]<0 else x) for e,x in b.items()};TT=lift(m,edges,optb)
  if fixed(TT)!=opt:raise RuntimeError('phase certificate')
  row={'controller':name+'+sink','signs':signs,'baseline':f,'A':A,'B':Bs,'optimized':opt,'outmasks':TT}
  rows.append(row)
  if best is None or opt>best['optimized']:best=row
 print('BEST',best)
open(Path(__file__).with_name('cyclic-controller-results.json'),'w').write(json.dumps(rows,indent=2)+'\n')
