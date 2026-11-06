#!/usr/bin/env python3
from functools import cache
import json

def need(ok,why):
 if not ok:raise RuntimeError(why)
@cache
def g(u):
 if u==1:return 1
 if u%2==0:return g(u//2)
 v=(u+1)//2;return g(v)+v-v%2

def graph(u):
 bits=bin(u-1)[2:] if u>1 else ''
 n=len(bits)+2;a=[0]*n;a[0]=2
 for j,b in enumerate(bits,2):
  a[j]|=1
  if b=='1':a[1]|=1<<j
  else:a[j]|=2
  for k in range(j+1,n):a[j]|=1<<k
 return tuple(a)
def valid(a):
 n=len(a)
 return all(not(a[i]>>i&1) for i in range(n)) and all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i in range(n) for j in range(i+1,n))
def literal(a):
 n=len(a);N=1<<n
 d=[[0]*n for _ in range(N)];f=[[0]*n for _ in range(N)];cy=[[0]*n for _ in range(N)]
 for v in range(n):d[1<<v][v]=f[1<<v][v]=1
 cy[1][0]=1
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   rest=S^(1<<v)
   if rest:
    for w in range(n):
     if rest>>w&1 and a[w]>>v&1:
      d[S][v]+=d[rest][w]
      if not(v==0 and rest&2):f[S][v]+=f[rest][w]
      cy[S][v]+=cy[rest][w]
 h=sum(d[-1]);forward=sum(f[-1]);cycles=sum(cy[-1][v] for v in range(n) if a[v]&1)
 return h,forward,h-forward,cycles

def extend(a,kind):
 out=list(a)+[0];z=len(a)
 for v in range(z):
  if v==0 or(kind=='B' and v==1):out[z]|=1<<v
  else:out[v]|=1<<z
 return tuple(out)
rows=[]
for u in range(1,130):
 a=graph(u);need(valid(a),'invalid graph');h,f,r,C=literal(a)
 need((h,f,r)==(u+g(u)-u%2,u,g(u)-u%2),'sector count')
 need(C==(0 if u==1 else 1),'cycle count')
 need(len(a)==2+(u-1).bit_length(),'order')
 if u>=2:
  for kind,v in [('A',2*u),('B',2*u-1)]:
   b=extend(a,kind);need(b==graph(v),'literal graph extension')
   hb,fb,rb,cb=literal(b)
   expected=(2*u,r+u%2) if kind=='A' else (2*u-1,r+u-1)
   need((fb,rb)==expected and hb==sum(expected) and cb==1,'two-overlap arithmetic')
 rows.append({'u':u,'H':h,'forward':f,'reverse':r,'cycles':C})
need([x['H'] for x in rows[:13]]==[1,3,5,5,9,9,11,9,17,15,19,15,23],'13 hostile table')
need(all(x['H']!=13 for x in rows[:13]),'13 finite inverse')
h13=(18,8,3,5,14);need(valid(h13) and literal(h13)[0]==13,'unrestricted13 witness')
need(literal(graph(39))[0]==63,'63 repair')
for k in range(1,20):
 u=2**k+1;need(g(u)==u,'sharp approach family')
print(json.dumps({'status':'PASS','scope':'literal graph and sector controls support separately proved all-order formulas; no scalar-extraction or eventual-coverage claim','indices':{'minimum':1,'maximum':129,'count':129},'extension_graphs_count':256,'arc_forward_sector_surjective_on_tested_indices':True,'E39_rows':list(graph(39)),'E39_H':63,'H13_rows':list(h13),'H13_H':13,'hostile_table':[x['H'] for x in rows[:13]],'sharp_approach_controls_k_max':19,'largest_literal_graph_order':len(graph(258))},indent=2,sort_keys=True))
