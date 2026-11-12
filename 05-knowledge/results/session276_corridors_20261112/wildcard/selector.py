"""Exact inverse constructor from the minimized all-order H=N language."""
from capped import build,minimize,count
from fibres import graph,ham
from collections import deque
import json

def select(N):
 vs,ts,ws,hs=build(N)
 tr,out,part=minimize(ts,[int(h==N) for h in hs],0)
 initial=part[0]
 q=deque([(initial,())]);seen=set();word=None
 while q:
  i,w=q.popleft()
  if i in seen:continue
  seen.add(i)
  if out[i]:word=w;break
  for l,j in enumerate(tr[i]):q.append((j,w+(l,)))
 witness=None
 if word is not None:
  g=graph(word)
  if count(word)!=N or ham(g)!=N:raise RuntimeError('bad witness')
  i=initial;walk=[i]
  for l in word:i=tr[i][l];walk.append(i)
  witness={'word':''.join(map(str,word)),'rows':g,'count':N,'walk':walk}
 return {'N':N,'initial':initial,'transitions':tr,'accepting':[i for i,a in enumerate(out) if a],'witness':witness,'status':'NONEMPTY' if witness else 'EMPTY'}
if __name__=='__main__':
 for N in (7,21,63):
  x=select(N)
  with open('05-knowledge/results/session276_corridors_20261112/wildcard/SELECTOR'+str(N)+'.json','w') as f:json.dump(x,f,indent=2)
  print(N,x['witness'])
