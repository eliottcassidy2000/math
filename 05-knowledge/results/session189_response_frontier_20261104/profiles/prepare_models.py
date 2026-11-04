import json,itertools
from pathlib import Path
D=Path(__file__).parent
cases=json.load(open(D/'quotients.json'));qq=[r['rows'] for r in cases if not any(r['sigma'])];models=[]
for Q in qq:
 ss=[[1]*7]
 for i in range(7):
  for s in (2,3):v=[1]*7;v[i]=s;ss.append(v)
 for i,j in itertools.combinations(range(7),2):v=[1]*7;v[i]=v[j]=2;ss.append(v)
 for sizes in ss:
  blocks=[];n=0
  for s in sizes:blocks.append(list(range(n,n+s)));n+=s
  for side in range(2):
   r=[0]*n
   for i in range(7):
    for u in blocks[i]:
     for v in blocks[i]:
      if u<v:r[u]|=1<<v
    for j in range(7):
     if i!=j and bool(Q[i]>>j&1)!=bool(side):
      for u in blocks[i]:
       for v in blocks[j]:r[u]|=1<<v
   models.append(dict(Q=Q,sizes=sizes,side=side,rows=r,blocks=blocks))
json.dump(models,open(D/'models.json','w'),indent=2)
with open(D/'models.in','w') as f:
 print(len(models),file=f)
 for m in models:print(len(m['rows']),*m['rows'],file=f)
