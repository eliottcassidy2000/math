from itertools import product
from math import gcd
import json
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def components(m,steps):
 unseen=set(range(m));count=0
 while unseen:
  count+=1;todo=[unseen.pop()]
  while todo:
   v=todo.pop()
   for s in steps:
    for w in ((v+s)%m,(v-s)%m):
     if w in unseen:unseen.remove(w);todo.append(w)
 return count
patterns=0
for m in (3,5,7,9):
 for levels in product(range(3),repeat=(m-1)//2):
  for d in (0,1):
   steps=[i+1 for i,l in enumerate(levels) if l>=d+1]
   c=gcd(m,*steps)
   need(components(m,steps)==c,'threshold gcd components')
  q=9
  for i in range(m):
   incoming=outgoing=0
   for s,l in enumerate(levels,1):
    block=3**(2-l);outgoing+=block;incoming+=q-block
   need((q-1)//2+incoming+outgoing==(m*q-1)//2,'stepwise balance')
  patterns+=1
# Concrete81-vertex regular tournament, prefix depths1 except step3 depth2.
m=9;k=2;levels={1:1,2:1,3:2,4:1};words=list(product(range(3),repeat=k));V=list(product(range(m),words))
def inner(u,v):
 for a,b in zip(u,v):
  if a!=b:return (b-a)%3==1
 return False
def arc(x,y):
 i,u=x;j,v=y
 if i==j:return inner(u,v)
 delta=(j-i)%m
 if delta<=m//2:return u[:levels[delta]]==v[:levels[delta]]
 l=levels[m-delta];return u[:l]!=v[:l]
need(all(sum(arc(x,y) for y in V)==40 for x in V),'actual81vertex regularity')
c0=components(m,[1,2,3,4]);c1=components(m,[3])
need((c0,c1)==(1,3),'example threshold components')
print(json.dumps({'status':'PASS','checks':checks,'depth_patterns':patterns,'example':{'vertices':81,'degree':40,'step_depths':levels,'threshold_components':[c0,c1],'declared_chart_order':3**(c0+3*c1),'owner_base':3*c1},'scope':'exact graph components and regularity; group and base formulas proved structurally'},indent=2))
