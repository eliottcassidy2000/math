#!/usr/bin/env python3
"""Threshold42 extension of the frozen Session124 recursive square-prefix core."""
import importlib.util,json,hashlib
from pathlib import Path
from math import isqrt
HERE=Path(__file__).resolve().parent
DEFAULT_CORE=HERE.parent/'session124_recursive_square_prefix_20261031'
DEFAULT_BASE=HERE.parent/'square_sum_B_base_interval_20261019'
def threshold(l):return (85*49**l-1)//2

def load_core(core_dir):
 spec=importlib.util.spec_from_file_location('frozen_square_core',Path(core_dir)/'core.py')
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def create_core(core_dir=DEFAULT_CORE,base_dir=DEFAULT_BASE):
 global vertices,prefix
 original=load_core(core_dir);vertices=original.vertices;prefix=original.prefix
 return ExtendedCore(original,core_dir,base_dir)

class ExtendedCore:
 def __init__(self,original,core_dir,base_dir):
  self.original=original;self.old=original.Core(base_dir);self.core_dir=Path(core_dir)
  self.small=json.loads((HERE/'SMALL_BASES.json').read_text());self.used_pairs=[]
  self.used_base_chunks=self.old.used_base_chunks
 @property
 def checks(self):return self.old.checks
 def check(self,ok,msg):self.old.check(ok,msg)
 def base(self,n):
  self.check(42<=n<=2081,'extension base domain')
  if n>=53:return self.old.base(n)
  pair=self.small[str(n)];self.old.check_pair(pair);self.check(len(pair[0])==n,'small base index');return pair
 def verify_bases(self):
  self.check(set(self.small)=={str(n) for n in range(42,53)},'small base exact catalog')
  for n in range(42,2082):self.base(n)
  return {'base_pairs':2040,'small_pairs':11,'checks':self.checks}
 def build_B(self,n):
  self.check(n>=42,'extension B index')
  if n<=2081:return self.base(n)
  m=(n-24)//49;d=n-49*m;pair=self.build_B(m)
  a=json.loads((self.core_dir/'paired'/f'B_p{m%2}_d{d}.json').read_text())
  result=[list(vertices(pair,w)) for w in a['words']]
  self.check(len(result[0])==n,'B output index');self.old.check_pair(result)
  self.used_pairs.append([m%2,d]);return result
 def cycle(self,N):
  self.check(N>=2082,'direct cycle domainN>=2082');m=(N-24)//49;d=N-49*m;pair=self.build_B(m);a=json.loads((self.core_dir/'cycles'/f'cycle_p{m%2}_d{d}.json').read_text())
  depth=1
  while threshold(depth+1)<=N:depth+=1
  P=prefix(depth);seen=bytearray(N+1);initial=[];first=last=None;count=0;digest=hashlib.sha256()
  for v in vertices(pair,a['word']):
   self.check(1<=v<=N and not seen[v],'cycle vertex ownership');seen[v]=1;count+=1;digest.update((str(v)+'\n').encode())
   if len(initial)<len(P):initial.append(v)
   if first is None:first=v
   if last is not None:self.check(isqrt(last+v)**2==last+v,'cycle square')
   last=v
  self.check(count==N and sum(seen)==N,'cycle exact coverage');self.check(initial==P,'maximum certified unscaled prefix')
  self.check(last==(8 if N%2==0 else 15),'closing endpoint');self.check(isqrt(first+last)**2==first+last,'cycle closing square')
  return {'status':'passed','N':N,'depth':depth,'threshold':threshold(depth),'prefix':P,'closing_vertex':last,'cycle_decimal_lines_sha256':digest.hexdigest(),'checks':self.checks,'used_base_chunks':sorted(self.used_base_chunks),'used_new_pair_states':self.used_pairs,'direct_state':[m%2,d],'cycle_stored':False}
