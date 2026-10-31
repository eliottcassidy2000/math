#!/usr/bin/env python3
"""Portable minimal constructor: old base pairs + new paired/terminal tables only."""
import gzip,hashlib,json,struct
from functools import lru_cache
from math import isqrt
from pathlib import Path
HERE=Path(__file__).resolve().parent
DEFAULT_BASE=HERE.parent/'square_sum_B_base_interval_20261019'

def threshold(l):return (107*49**l-1)//2

def prefix(l):
 p=[1,3]
 for _ in range(l):p=[1,3,6,19]+[49*v+(-1)**(j+1)*4 for j,v in enumerate(p)]
 return p

def vertices(pair,word):
 for a in word:
  if a[0]=='single':yield a[1]
  elif a[0]=='head':yield 49-a[1]
  else:
   kind,c,i,o=a;start=int(kind=='tail');end=len(pair[i])
   indices=range(start,end) if o==1 else range(end-1,start-1,-1)
   for j in indices:yield 49*pair[i][j]+(-1)**(j+1)*c

class Core:
 def __init__(self,base_dir=DEFAULT_BASE):
  self.base_dir=Path(base_dir);self.checks=0;self.used_base_chunks=set();self.used_pairs=[]
  self.deps=json.loads((HERE/'BASE_DEPENDENCIES.json').read_text())
  raw=(self.base_dir/self.deps['manifest_file']).read_bytes()
  self.check(hashlib.sha256(raw).hexdigest()==self.deps['manifest_sha256'],'base manifest hash')
  self.check(json.loads(raw)==self.deps['chunks'],'base manifest contents')
  self.chunks=self.deps['chunks'];next_index=53
  for c in self.chunks:self.check(c['first']==next_index and c['last']>=c['first'],'contiguous base chunks');next_index=c['last']+1
  self.check(next_index==2621,'complete base coverage')
 def check(self,ok,msg):
  self.checks+=1
  if not ok:raise RuntimeError(msg)
 @lru_cache(maxsize=2)
 def chunk(self,name):
  c=next(x for x in self.chunks if x['file']==name);z=(self.base_dir/name).read_bytes()
  self.check(len(z)==c['bytes'] and hashlib.sha256(z).hexdigest()==c['sha256'],'compressed base hash')
  raw=gzip.decompress(z);self.check(len(raw)==c['decoded_bytes'] and hashlib.sha256(raw).hexdigest()==c['decoded_sha256'],'decoded base hash')
  self.check(raw[:8]==b'SSBAS1\0\0','base header');lo,hi=struct.unpack('<II',raw[8:16]);self.check((lo,hi)==(c['first'],c['last']),'base header bounds')
  self.check(len(raw)==16+(hi-lo+1)*(lo+hi-1),'base exact byte count');self.used_base_chunks.add(name);return raw
 def check_pair(self,pair,P=(1,3)):
  n=len(pair[0]);self.check(len(pair[1])==n+1,'pair lengths');colors=bytearray(n+2)
  for i,path in enumerate(pair):
   self.check(path[:len(P)]==list(P) and path[-1]==(8 if len(path)%2==0 else 10),'B prefix/terminal');seen=bytearray(len(path)+1)
   for j,v in enumerate(path):
    self.check(1<=v<=len(path) and not seen[v],'pair interval coverage');seen[v]=1
    if i==0:colors[v]=(j+1)%2+1
    elif v<=n:self.check(colors[v]==(j+1)%2+1,'pair common parity')
    if j:self.check(isqrt(path[j-1]+v)**2==path[j-1]+v,'pair square')
 def base(self,n):
  self.check(53<=n<=2620,'base index');c=next(x for x in self.chunks if x['first']<=n<=x['last']);raw=self.chunk(c['file']);lo=c['first'];offset=16+(n-lo)*(n+lo-2);pair=[]
  for length in(n,n+1):
   roots=raw[offset:offset+length-1];self.check(len(roots)==length-1,'base root count');offset+=length-1;path=[1]
   for r in roots:path.append(r*r-path[-1])
   pair.append(path)
  self.check_pair(pair);return pair
 def verify_bases(self):
  for n in range(53,2621):self.base(n)
  return {'base_pairs':2568,'base_chunks':len(self.used_base_chunks),'checks':self.checks}
 def build_B(self,n):
  self.check(n>=53,'B index')
  if n<=2620:return self.base(n)
  m=(n-24)//49;d=n-49*m;pair=self.build_B(m);a=json.loads((HERE/'paired'/f'B_p{m%2}_d{d}.json').read_text())
  result=[list(vertices(pair,w)) for w in a['words']];self.check(len(result[0])==n,'B output index');self.check_pair(result);self.used_pairs.append([m%2,d]);return result
 def cycle(self,N):
  self.check(N>=2621,'direct cycle domainN>=2621');m=(N-24)//49;d=N-49*m;pair=self.build_B(m);a=json.loads((HERE/'cycles'/f'cycle_p{m%2}_d{d}.json').read_text())
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
