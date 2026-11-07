import sys,json,itertools,random,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
CARRIERS=ROOT/'carriers'
if not CARRIERS.exists(): CARRIERS=ROOT/'session215-carriers'
def mod(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
car=mod('car',CARRIERS/'inherited.py')
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
