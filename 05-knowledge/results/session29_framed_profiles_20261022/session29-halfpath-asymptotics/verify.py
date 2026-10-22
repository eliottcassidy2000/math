import itertools,json,math
from pathlib import Path
P=Path(__file__).parent
def dp(m,q,neg):
 E=list(itertools.combinations(range(m),2));T=[0]*(2*m)
 for i in range(m):T[2*i]|=1<<(2*i+1)
 for e,(i,j) in enumerate(E):
  c=j if q>>e&1 else i;b=-1 if neg>>e&1 else 1
  for si in (1,-1):
   for sj in (1,-1):
    a=2*i+(si==-1);z=2*j+(sj==-1)
    need=-b if c==i else b
    yes=(si if c==i else sj)==need
    T[a if yes else z]|=1<<(z if yes else a)
 D=[[0]*(2*m) for _ in range(1<<m)]
 for a in range(2*m):D[1<<(a//2)][a]=1
 for S in range(1,1<<m):
  for a in range(2*m):
   if not D[S][a]:continue
   for j in range(m):
    if S>>j&1:continue
    for z in (2*j,2*j+1):
     if T[a]>>z&1:D[S|1<<j][z]+=D[S][a]
 return sum(D[-1][::2]),T
out=[]
for m in (6,7):
 rows=[list(map(int,l.split())) for l in (P/f'm{m}_classes.txt').read_text().splitlines()]
 best=max(r[2] for r in rows)
 for r in rows:
  if r[2]!=best:continue
  mm,q,f,A,neg,d,*B=r
  got,T=dp(m,q,neg);assert got==f,(r,got)
  out.append(dict(m=m,q=q,negative_label_mask=neg,fixed_hp=f,A=A,B=B,outmasks=T,sigma=[i^1 for i in range(2*m)]))
(P/'witnesses.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
