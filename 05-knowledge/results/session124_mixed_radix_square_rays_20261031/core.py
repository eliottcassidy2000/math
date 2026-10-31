"""Portable A/B/C schedule constructor; standard library, no solver imports."""
import hashlib,json
from math import isqrt,prod
from pathlib import Path
from bases import Bases
HERE=Path(__file__).resolve().parent
DEFAULT_BASE=HERE.parent/'square_sum_B_base_interval_20261019'
DEFAULT_SMALL=HERE.parent/'session124_square_prefix_threshold42_20261031'/'SMALL_BASES.json'
DEFAULT_Q49=HERE.parent/'session124_recursive_square_prefix_20261031'
OPS={'A':(25,12,(1,3)),'B':(25,8,(1,3,6,19)),'C':(49,4,(1,3,6,19))}
def valid_word(word):
 if not isinstance(word,str) or any(x not in OPS for x in word):raise ValueError('word must use only A, B, C')
 return word
def threshold(word):return(85*prod(OPS[x][0] for x in valid_word(word))-1)//2
def prefix(word):
 p=[1,3]
 for symbol in reversed(valid_word(word)):
  q,c,h=OPS[symbol];p=list(h)+[q*v+(-1)**(j+1)*c for j,v in enumerate(p)]
 return p
def vertices(q,pair,word):
 for a in word:
  kind,c=a[:2]
  if kind=='single':yield c
  elif kind=='head':yield q-c
  else:
   i,o=a[2:];start=int(kind=='tail');indices=range(start,len(pair[i])) if o==1 else range(len(pair[i])-1,start-1,-1)
   for j in indices:yield q*pair[i][j]+(-1)**(j+1)*c
class Constructor:
 def __init__(self,base_dir=DEFAULT_BASE,small_base_file=DEFAULT_SMALL,q49_dir=DEFAULT_Q49):
  self.bases=Bases(base_dir,small_base_file);self.q49_dir=Path(q49_dir);self.deps=json.loads((HERE/'DEPENDENCIES.json').read_text());self.checks=0;self.used_q49=False;self.hashes=json.loads((HERE/'CERTIFICATES.json').read_text())
 def need(self,ok,msg):
  self.checks+=1
  if not ok:raise RuntimeError(msg)
 def certificate(self,symbol,p,d,terminal=False):
  if symbol=='C':
   name=f'{"cycles/cycle" if terminal else "paired/B"}_p{p}_d{d}.json';raw=(self.q49_dir/name).read_bytes();self.need(hashlib.sha256(raw).hexdigest()==self.deps['Q49_sha256'][name],'Q49 certificate hash');self.used_q49=True
  else:
   name=f'{symbol}/{"cycles/C" if terminal else "paired/B"}_p{p}_d{d}.json';raw=(HERE/name).read_bytes();self.need(hashlib.sha256(raw).hexdigest()==self.hashes[name],'Q25 certificate hash')
  return json.loads(raw)
 def pair(self,n,word=''):
  valid_word(word);self.need(isinstance(n,int) and n>=42,'pair index')
  # The empty suffix uses an A-only weak-B supply, avoiding Q49 dependencies.
  symbol=word[0] if word else 'A';suffix=word[1:] if word else '';q,_,_=OPS[symbol];half=(q-1)//2
  if n<42*q+half:pair=self.bases.base(n)
  else:
   m=(n-half)//q;d=n-q*m;source=self.pair(m,suffix);a=self.certificate(symbol,m%2,d);pair=[list(vertices(q,source,w)) for w in a['words']]
  self.bases.check_pair(pair,prefix(word) if n>=threshold(word) else[1,3]);self.need(len(pair[0])==n,'pair exact index');return pair
 def cycle(self,N,word):
  valid_word(word);self.need(bool(word) and isinstance(N,int) and N>=threshold(word),'requested cycle domain')
  symbol=word[0];q,_,_=OPS[symbol];half=(q-1)//2;m=(N-half)//q;d=N-q*m;pair=self.pair(m,word[1:]);a=self.certificate(symbol,m%2,d,True);atoms=a['word'] if symbol=='C' else a['words'][0]
  seen=bytearray(N+1);P=prefix(word);initial=[];first=last=None;count=0;digest=hashlib.sha256()
  for v in vertices(q,pair,atoms):
   self.need(1<=v<=N and not seen[v],'cycle interval uniqueness');seen[v]=1;count+=1;digest.update((str(v)+'\n').encode())
   if len(initial)<len(P):initial.append(v)
   if first is None:first=v
   if last is not None:self.need(isqrt(last+v)**2==last+v,'cycle square edge')
   last=v
  self.need(count==N and sum(seen)==N,'cycle exact coverage');self.need(initial==P,'requested unscaled prefix');self.need(first==1 and last==(8 if symbol!='C' or N%2==0 else 15),'closing endpoint');self.need(isqrt(first+last)**2==first+last,'closing square')
  return {'status':'passed','N':N,'word':word,'threshold':threshold(word),'prefix':P,'cycle_decimal_lines_sha256':digest.hexdigest(),'checks':self.checks,'base_checks':self.bases.checks,'Q49_used':self.used_q49,'cycle_stored':False}
