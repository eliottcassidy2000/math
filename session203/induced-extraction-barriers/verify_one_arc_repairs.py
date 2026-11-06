#!/usr/bin/env python3
"""Literal graph/DP controls for one-old-arc surgery. No infinite coverage claim."""
import json
from pathlib import Path
BASE=Path(__file__).resolve().parent

def need(ok,why):
 if not ok:raise RuntimeError(why)
def source(u):
 bits=bin(u-1)[2:];n=len(bits)+2;a=[0]*n;a[0]=2
 for j,b in enumerate(bits,2):
  a[j]|=1
  if b=='1':a[1]|=1<<j
  else:a[j]|=2
  for k in range(j+1,n):a[j]|=1<<k
 return a
def add(a,mask):
 n=len(a);return [r|(0 if mask>>j&1 else 1<<n) for j,r in enumerate(a)]+[mask]
def valid(a):
 n=len(a)
 return all(not(a[j]>>j&1) for j in range(n)) and all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i in range(n) for j in range(i+1,n))
def H(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for j in range(n):d[1<<j][j]=1
 for S in range(1,1<<n):
  for j in range(n):
   if not d[S][j]:continue
   z=a[j]&~S
   while z:
    w=z&-z;z-=w;d[S|w][w.bit_length()-1]+=d[S][j]
 return sum(d[-1])
def strong(a):
 n=len(a)
 for v in range(n):
  seen=1<<v;old=0
  while old!=seen:
   old=seen
   for j in range(n):
    if seen>>j&1:seen|=a[j]
  if seen!=(1<<n)-1:return False
 return True
repairs=json.loads((BASE/'ONE_ARC_REPAIRS.json').read_text())
for r in repairs:
 a=source(r['u']);i,j=r['reversed_old_arc'];need(a[i]>>j&1,'specified old arc')
 a[i]^=1<<j;a[j]^=1<<i
 a=add(add(a,r['first_new_mask']),r['second_new_mask'])
 need(a==r['rows'] and valid(a),'repair graph')
 need(H(a)==r['target']==r['H'],'repair H')
 need(len(a)==r['order'] and strong(a)==r['strong'],'repair metadata')
rows391=[514,1536,2043,1011,995,963,1923,1795,1539,1024,57]
a=source(65);a[1]^=4;a[2]^=2;a=add(add(a,0),57)
need(a==rows391 and valid(a) and H(a)==391,'391 repair')
b=source(129);b[1]^=4;b[2]^=2;b=add(add(b,0),57)
need(H(b)==783 and H(b)!=6*129+1,'nearby generalization hostile')
print(json.dumps({'status':'PASS','one_arc_repair_targets':[r['target'] for r in repairs],
 'all_five_recorded_graphs_nonstrong':all(not r['strong'] for r in repairs),
 'extra_repair':{'u':65,'target':391,'reversed_old_arc':[1,2],'first_new_mask':0,'second_new_mask':57,'rows':rows391},
 'nearby_same_rule_hostile':{'u':129,'actual':783,'desired':775},
 'scope':'Exact finite repair graphs; no all-u extraction or eventual threshold proved.'},indent=2,sort_keys=True))
