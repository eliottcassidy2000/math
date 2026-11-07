import sys,json,itertools,random,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
MULTI=ROOT/'multi-exit'
if not MULTI.exists(): MULTI=ROOT/'session215-multi-exit'
def mod(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
multi=mod('multi',MULTI/'verify.py')
def require(b,msg):
 if not b:raise RuntimeError(msg)
def enum_h(rows):
 return sum(all(rows[a]>>b&1 for a,b in zip(p,p[1:])) for p in itertools.permutations(range(len(rows))))
def subset_h(rows):
 # A memoized path-from-first recursion, independent from the iterative end DP.
 from functools import lru_cache
 n=len(rows)
 @lru_cache(None)
 def paths(mask,v):
  rest=mask^(1<<v)
  if rest==0:return 1
  return sum(paths(rest,w) for w in range(n) if rest>>w&1 and rows[v]>>w&1)
 @lru_cache(None)
 def h(mask):return sum(paths(mask,v) for v in range(n) if mask>>v&1) if mask else 1
 return paths,h
rng=random.Random(215);tests=0
for k in range(5):
 for trial in range(15):
  A=[0]*k
  for i,j in itertools.combinations(range(k),2):
   a,b=(i,j) if rng.randrange(2) else (j,i);A[a]|=1<<b
  w=[rng.randrange(1<<k) for _ in range(8-k)]
  rows=multi.rows(A,w);actual=enum_h(rows)
  require(actual==multi.transfer(A,w),('multi',A,w,actual));tests+=1
print(json.dumps({"status":"PASS","independent_permutation_cases":tests,"vertex_order":8,"apex_range":[0,4]},indent=2,sort_keys=True))
