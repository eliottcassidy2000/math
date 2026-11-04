import json
from pathlib import Path
D=Path(__file__).parent
A=json.load(open(D/'anchor_exploration.json'));B=json.load(open(D/'models.json'));qs={tuple(m['Q']) for m in B};out=[]
for a in A:
 Q=a['rows']
 if not a['blind'] or tuple(Q) in qs:continue
 sizes=[[1]*7]
 for j in range(7):v=[1]*7;v[j]=2;sizes.append(v)
 for ss in sizes:
  blocks=[];n=0
  for s in ss:blocks.append(list(range(n,n+s)));n+=s
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
   out.append(dict(Q=Q,sizes=ss,side=side,rows=r,blocks=blocks))
json.dump(out,open(D/'extra_models.json','w'))
with open(D/'extra_models.in','w') as f:
 print(len(out),file=f)
 for m in out:print(len(m['rows']),*m['rows'],file=f)
