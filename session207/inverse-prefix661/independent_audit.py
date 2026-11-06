import json,pathlib
from functools import lru_cache
P=pathlib.Path(__file__).resolve().parent

def decode(code,n):
 a=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   x,y=(i,j) if code>>k&1 else (j,i)
   a[x]|=1<<y;k+=1
 return a

def count(a):
 n=len(a)
 @lru_cache(None)
 def end(s,v):
  if s==(1<<v):return 1
  t=s^(1<<v)
  return sum(end(t,u) for u in range(n) if t>>u&1 and a[u]>>v&1)
 return sum(end((1<<n)-1,v) for v in range(n))

m={}
rows=json.loads((P/'census.json').read_text())['orders']
for r in rows:
 for h in r['spectrum']:m.setdefault(h,r['n'])
found=set()
for line in (P/'fill.out').read_text().splitlines():
 h,c,base,pat=map(int,line.split());a=decode(c,9)
 assert count(a)==h and h not in m and h not in found
 old=decode(base,8)
 assert [x&255 for x in a[:8]]==old and a[8]==pat
 for v in range(9):
  seen=1<<v
  for _ in range(9):
   for u in range(9):
    if seen>>u&1:seen|=a[u]
  assert seen==511
 m[h]=9;found.add(h)
assert found=={611,615,617,619,623,625,635,647,655}
assert set(m)==set(range(1,662,2))-{7,21}
assert {str(h):m[h] for h in sorted(m)}==json.loads((P/'minimum_orders.json').read_text())['minimum_orders']
print('PASS all nine strong order-nine witnesses, encoded augmentations, and exact 329-entry inverse table')
