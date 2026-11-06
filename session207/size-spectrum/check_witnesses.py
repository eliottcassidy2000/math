#!/usr/bin/env python3
"""Independent Python endpoint recurrence checks every published witness."""
import json,pathlib
from functools import lru_cache
P=pathlib.Path(__file__).resolve().parent
j=json.loads((P/'census.json').read_text())
def decode(code,n):
 a=[0]*n;k=0
 for i in range(n):
  for z in range(i+1,n):
   if code>>k&1:a[i]|=1<<z
   else:a[z]|=1<<i
   k+=1
 return a
def hp(a):
 n=len(a)
 @lru_cache(None)
 def end(s,v):
  if s==1<<v:return 1
  t=s^(1<<v)
  return sum(end(t,u) for u in range(n) if t>>u&1 and a[u]>>v&1)
 return sum(end((1<<n)-1,v) for v in range(n))
def ext(a):
 n=len(a);out=set()
 for p in range(1<<n):
  b=[a[i]|(0 if p>>i&1 else 1<<n) for i in range(n)]+[p]
  out.add(hp(b))
 return out
for r in j['orders']:
 for h,c in r['witness_codes'].items(): assert hp(decode(c,r['n']))==int(h)
print('Independent endpoint recurrence: all 430 witnesses verified.')
by={};collision=None
for c in range(1024):
 a=decode(c,5);h=hp(a);s=ext(a)
 if h in by and by[h][1]!=s:
  collision={'n':5,'h':h,'code_a':by[h][0],'code_b':c,'extensions_a':sorted(by[h][1]),'extensions_b':sorted(s)};break
 by[h]=(c,s)
assert collision
(P/'scalar_extension_counterexample.json').write_text(json.dumps(collision,indent=2)+'\n')
print(json.dumps(collision))
