#!/usr/bin/env python3
"""Portable checks of diagram bytes, labels and numeric invariants; not a proof replay."""
from pathlib import Path
from hashlib import sha256
from math import comb
import json, subprocess, sys, tempfile, shutil
import xml.etree.ElementTree as ET
P=Path(__file__).resolve().parent

def require(ok,msg):
 if not ok: raise ValueError(msg)
def integrity(root):
 manifest=json.loads((root/'manifest.json').read_text())
 actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
 require(actual==set(manifest),'manifest membership mismatch')
 for rel,digest in manifest.items(): require(sha256((root/rel).read_bytes()).hexdigest()==digest,'digest: '+rel)
def mul(a,b):
 out=0
 while b:
  if b&1: out^=a
  b>>=1; a<<=1
  if a&64: a^=67 # X^6 + X + 1
 return out
def power(a,n):
 out=1
 while n:
  if n&1: out=mul(out,a)
  a=mul(a,a); n>>=1
 return out
def trace(a): return a^power(a,4)^power(a,16)
def checks():
 integrity(P)
 c=json.loads((P/'counts.json').read_text())
 for b in range(7):
  v=[2]+([2]+[0]*(b-1) if b else [])
  for r in range(1,11):
   require(v==[2*3**(r-1)*comb(r,j) if j<=r else 0 for j in range(b+1)],'transfer/binomial mismatch')
   if r==10 and b in (0,1): require(sum(v)==c['r10']['b'+str(b)],'r10 count')
   v=[3*v[j]+(3*v[j-1] if j else 0) for j in range(b+1)]
 # Independent exact field arithmetic, not copied incidence tables.
 g=[power(2,e) for e in range(63)]
 require(len(set(g))==63,'primitive element')
 fibers=[{g[e+21*j] for j in range(3)} for e in range(21)]
 require(len(set.union(*fibers))==63 and all(len(f)==3 for f in fibers),'scalar fibers')
 flags={(e,l) for e in range(21) for l in range(21) if trace(g[(e-l)%63])==0}
 require(len(flags)==105 and all(sum(e==p for e,l in flags)==5 for p in range(21)), 'point flags/degrees')
 require(all(sum(l==p for e,l in flags)==5 for p in range(21)),'line degrees')
 require(all(((-l)%21,(-e)%21) in flags for e,l in flags),'typed polarity')
 require(c['native']=={'blockers':63,'boundary_images':21,'fiber_size':3},'native counts')
 require(c['typed']=={'objects':42,'flags':105,'degree':5,'lift_vertices':126,'lift_edges':945,'lift_degree':15},'typed counts')
 require(2*21*3==126 and len(flags)*3*3==945 and 5*3==15,'blow-up arithmetic')
 required={
  'digit-budget.svg':['39,366','433,026','r ≥ 1','N e_b = 0','dimension 1/2','finiteness or infinitude is NOT settled','not a coefficient congruence'],
  'boundary-incidence.svg':['63 native','21 points','42 objects','105 flags','945 edges','degree 5 × 3 = 15','ONE point, not a line','NEW relation','Not a locally bijective graph covering','NOT a collineation']}
 for name,labels in required.items():
  xml=(P/name).read_text(); root=ET.fromstring(xml)
  txt=' '.join(root.itertext())
  for label in labels: require(label in txt,'missing label: '+label)
  require(len(xml.encode())<20000,'SVG not lightweight')
  require('href=' not in xml and '<script' not in xml,'external or script content')
 with tempfile.TemporaryDirectory() as t:
  q=Path(t); shutil.copy(P/'generate.py',q/'generate.py')
  subprocess.run([sys.executable,'-B',str(q/'generate.py')],check=True)
  for name in [*required,'counts.json']: require((q/name).read_bytes()==(P/name).read_bytes(),'non-reproducible '+name)
 return {'status':'PASS','figures':2,'r10_counts':[39366,433026],'base_flags':105,'lift_edges':945,'scope':'diagram integrity and numeric/label sanity, not full source theorem replay'}
def selftest():
 with tempfile.TemporaryDirectory() as t:
  q=Path(t)/'packet'; shutil.copytree(P,q)
  integrity(q)
  for kind in ('mutation','extra','missing'):
   target=q/'digit-budget.svg'; original=target.read_bytes()
   if kind=='mutation': target.write_bytes(original+b' ')
   elif kind=='extra': (q/'extra.txt').write_text('hostile')
   else: target.unlink()
   try: integrity(q)
   except ValueError: pass
   else: raise ValueError('hostile accepted: '+kind)
   target.write_bytes(original)
   if kind=='extra': (q/'extra.txt').unlink()
 normal=subprocess.check_output([sys.executable,'-B',str(P/'replay.py')])
 optimized=subprocess.check_output([sys.executable,'-O','-B',str(P/'replay.py')])
 require(normal==optimized,'normal/optimized mismatch')
if __name__=='__main__':
 result=checks()
 if '--self-test' in sys.argv: selftest(); result['self_test']='PASS: mutation, missing, extra, normal/-O'
 print(json.dumps(result,sort_keys=True))
