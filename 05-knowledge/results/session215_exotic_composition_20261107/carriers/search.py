from inherited import *
import json
banks=[]
for pair,js,root,a,b in PAIRS:
 for side,rows,j in zip('AB',(a,b),js):
  init=profile(remark(rows,root)); seen={init:''}; level={init:''}; bybase={}
  for depth in range(10):
   nxt={}
   for p,w in level.items():
    x=base(p,j)
    if not w.endswith('O'):
     bybase.setdefault(x,(pair,side,w,p,j))
    if depth<9:
     for c in 'OZJ':
      q=move(p,c)
      if q not in seen: seen[q]=w+c;nxt[q]=w+c
   level=nxt
  banks.append(bybase)
  print(pair,side,len(seen),len(bybase),flush=True)
res=[]
for i in range(0,8,2):
 for k in range(1,8,2):
  for x in banks[i].keys() & banks[k].keys():
   a,b=banks[i][x],banks[k][x]
   if 'J' in a[2]+b[2]:res.append((x,a,b))
res.sort(key=lambda z:(max(len(z[1][2]),len(z[2][2])),z[0]))
json.dump(res,open('collisions.json','w'),indent=2)
print('COLLISIONS',len(res));print(res[:10])
