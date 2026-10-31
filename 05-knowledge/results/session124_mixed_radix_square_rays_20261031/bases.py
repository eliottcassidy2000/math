import gzip,hashlib,json,struct
from functools import lru_cache
from math import isqrt
from pathlib import Path
HERE=Path(__file__).resolve().parent

class Bases:
 def __init__(self,base_dir,small_base_file):
  smallraw=Path(small_base_file).read_bytes();dep=json.loads((HERE/'DEPENDENCIES.json').read_text());
  if hashlib.sha256(smallraw).hexdigest()!=dep['small_base_sha256']:raise RuntimeError('small base hash')
  self.small=json.loads(smallraw)
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
  self.check(42<=n<=2081,'base index')
  if n<53:
   pair=self.small[str(n)];self.check_pair(pair);self.check(len(pair[0])==n,'small base index');return pair
  c=next(x for x in self.chunks if x['first']<=n<=x['last']);raw=self.chunk(c['file']);lo=c['first'];offset=16+(n-lo)*(n+lo-2);pair=[]
  for length in(n,n+1):
   roots=raw[offset:offset+length-1];self.check(len(roots)==length-1,'base root count');offset+=length-1;path=[1]
   for r in roots:path.append(r*r-path[-1])
   pair.append(path)
  self.check_pair(pair);return pair
