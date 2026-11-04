#!/usr/bin/env python3
"""Exact endpoint-interface experiments; no external dependencies."""
import itertools,json,hashlib
from collections import deque
from pathlib import Path

def req(x,m):
 if not x: raise RuntimeError(m)
def tables(a):
 n=len(a);N=1<<n;e=[[0]*n for _ in range(N)]
 for i in range(n):e[1<<i][i]=1
 for s in range(1,N):
  for i in range(n):
   if not(s>>i&1):continue
   r=s^(1<<i)
   if r:e[s][i]=sum(e[r][j] for j in range(n) if (r>>j&1) and (a[j]>>i&1))
 return e
def reverse(a):
 n=len(a);return tuple(((1<<n)-1)^x^(1<<i) for i,x in enumerate(a))
def data(a):
 e=tables(a);b=tables(reverse(a));return sum(e[-1]),e,b
def count(a):return sum(tables(a)[-1])
def flip(a,u,v):
 a=list(a);a[u]^=1<<v;a[v]^=1<<u;return tuple(a)
def literal(a):
 return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(a)-1)) for p in itertools.permutations(range(len(a))))
def edge_uses(a,e,b,u,v):
 n=len(a);others=((1<<n)-1)^(1<<u)^(1<<v);s=others;ans=0
 while True:
  ans+=e[s|(1<<u)][u]*b[(others^s)|(1<<v)][v]
  if not s:break
  s=(s-1)&others
 return ans
def insertion(a,bits):
 n=len(a);return tuple(a[i]|((1<<n) if bits>>i&1 else 0) for i in range(n))+(((1<<n)-1)^bits,)
def orphan_formula(a,bits,e,b):
 full=(1<<len(a))-1;total=0
 for s in range(1,full):
  r=full^s
  for u in range(len(a)):
   if not(s>>u&1 and bits>>u&1):continue
   for w in range(len(a)):
    if r>>w&1 and not(bits>>w&1) and a[w]>>u&1:
     total+=e[s][u]*b[r][w]
 return total
def insertion_literal(a,bits):
 n=len(a);old=[];positions=0
 for p in itertools.permutations(range(n)):
  if all(a[p[i]]>>p[i+1]&1 for i in range(n-1)):
   old.append(p)
   positions+=sum((k==0 or bits>>p[k-1]&1) and (k==n or not(bits>>p[k]&1)) for k in range(n+1))
 new=insertion(a,bits);orphan=0
 for p in itertools.permutations(range(n+1)):
  if all(new[p[i]]>>p[i+1]&1 for i in range(n)):
   q=tuple(x for x in p if x!=n)
   orphan+=not all(a[q[i]]>>q[i+1]&1 for i in range(n-1))
 return len(old),positions,orphan

def tour(n,k):
 a=[0]*n;z=0
 for i in range(n):
  for j in range(i+1,n):
   if k>>z&1:a[i]|=1<<j
   else:a[j]|=1<<i
   z+=1
 return tuple(a)
def strong(a):
 n=len(a)
 for source in range(n):
  seen=1<<source
  while True:
   nxt=seen
   for i in range(n):
    if seen>>i&1:nxt|=a[i]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True
def box(A,B,k):
 a=list(A)+[x<<3 for x in B]+[7]
 for j in range(3):a[3+j]|=1<<6
 for i in range(3):
  for j in range(3):
   if k>>(3*i+j)&1:a[i]|=1<<(3+j)
   else:a[3+j]|=1<<i
 return tuple(a)
def prime(x):return x>=2 and all(x%d for d in range(2,int(x**.5)+1))
def main():
 checks=0;orphan_examples=[]
 for n in range(1,6):
  for k in range(1<<(n*(n-1)//2)):
   a=tour(n,k);h,e,b=data(a);req(h==literal(a),'permutation H')
   for u in range(n):
    for v in range(u+1,n):
     x,y=(u,v) if a[u]>>v&1 else (v,u)
     d=edge_uses(a,e,b,y,x)-edge_uses(a,e,b,x,y)
     req(count(flip(a,x,y))==h+d,'edge identity');checks+=1
   if n<=4:
    for bits in range(1<<n):
     h0,slots,o=insertion_literal(a,bits)
     req(o==orphan_formula(a,bits,e,b),'orphan identity')
     req(count(insertion(a,bits))==slots+o,'insertion identity')
     if o and len(orphan_examples)<3:orphan_examples.append(dict(base=list(a),bits=bits,H=h0,legal_slots=slots,orphans=o,new_H=slots+o))
 families=[];witnesses={}
 for name,A,B in [('TT',(6,4,0),(6,4,0)),('CC',(2,4,1),(2,4,1)),('CT',(2,4,1),(6,4,0))]:
  graphs=[box(A,B,k) for k in range(512)];ds=[data(a) for a in graphs];hs=[d[0] for d in ds];st=[strong(a) for a in graphs]
  for k,a in enumerate(graphs):
   h,e,b=ds[k]
   for bit in range(9):
    u,v=bit//3,3+bit%3
    if not(a[u]>>v&1):u,v=v,u
    req(hs[k^(1<<bit)]-h==edge_uses(a,e,b,v,u)-edge_uses(a,e,b,u,v),'interface edge identity');checks+=1
   if st[k] and (prime(h) or (h%7==0 and prime(h//7) and h//7!=3)):
    witnesses.setdefault(str(h),dict(family=name,interface=k,rows=list(a),strong=True,lane='prime' if prime(h) else '7p'))
  components=[];seen=set()
  for seed in range(512):
   if seed in seen:continue
   seen.add(seed);todo=[seed];comp=[]
   while todo:
    k=todo.pop();comp.append(k)
    for bit in range(9):
     j=k^(1<<bit)
     if j not in seen and abs(hs[k]-hs[j])<=2:seen.add(j);todo.append(j)
   vals=sorted({hs[k] for k in comp});req(vals==list(range(vals[0],vals[-1]+1,2)),'component interval')
   components.append(dict(size=len(comp),min=vals[0],max=vals[-1],values=vals,states=sorted(comp)))
  components.sort(key=lambda c:(-c['max']+c['min'],-c['size']))
  strong_components=[];seen_strong=set()
  for seed in range(512):
   if not st[seed] or seed in seen_strong:continue
   seen_strong.add(seed);todo=[seed];comp=[]
   while todo:
    k=todo.pop();comp.append(k)
    for bit in range(9):
     j=k^(1<<bit)
     if st[j] and j not in seen_strong and abs(hs[k]-hs[j])<=2:seen_strong.add(j);todo.append(j)
   vals=sorted({hs[k] for k in comp});req(vals==list(range(vals[0],vals[-1]+1,2)),'strong component interval')
   strong_components.append(dict(size=len(comp),min=vals[0],max=vals[-1],values=vals,states=sorted(comp)))
  strong_components.sort(key=lambda c:(-c['max']+c['min'],-c['size']))
  families.append(dict(strong_components=strong_components,name=name,states=512,strong_states=sum(st),distinct_H=sorted(set(hs)),components=components,H=hs))
 out=dict(schema=1,scope='All labelled tournaments orders1..5; all insertions into orders1..4; three 512-state labelled 3+3+1 interface cubes, not all order7 tournaments',edge_checks=checks,orphan_examples=orphan_examples,families=families,witnesses=witnesses)
 print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()
