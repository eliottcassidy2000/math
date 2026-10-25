#!/usr/bin/env python3
"""Portable SVG reproduction, exact arithmetic checks, and packet integrity."""
from pathlib import Path
from hashlib import sha256
import itertools,json,subprocess,sys,tempfile,shutil
import xml.etree.ElementTree as ET
P=Path(__file__).resolve().parent
LINES=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
def need(ok,msg):
 if not ok:raise ValueError(msg)
def integrity(root):
 d=json.loads((root/'manifest.json').read_text())
 files={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='manifest.json'}
 need(files==set(d),'manifest membership')
 for rel,h in d.items():need(sha256((root/rel).read_bytes()).hexdigest()==h,'digest '+rel)
def ring(n):
 Q={r*r%n for r in range(n)}
 for c,a,b in itertools.product(Q,repeat=3):
  z=[a,(3*c-a-b)%n,b,(c-a+b)%n,c,(c+a-b)%n,(2*c-b)%n,(a+b-c)%n,(2*c-a)%n]
  if len(set(z))==9 and all(x in Q for x in z):
   need(len({sum(z[k] for k in L)%n for L in LINES})==1,'ring line sum');return False
 return True

def field(row):
 p=row['p'];d=row['n'];q=row['q'];f=row['modulus_low_first']
 def digits(x):return [(x//p**i)%p for i in range(d)]
 def enc(a):return sum((x%p)*p**i for i,x in enumerate(a))
 def mul(a,b):
  aa=digits(a);bb=digits(b);z=[0]*(2*d-1)
  for i in range(d):
   for j in range(d):z[i+j]=(z[i+j]+aa[i]*bb[j])%p
  for i in range(2*d-2,d-1,-1):
   for j in range(d):z[i-d+j]=(z[i-d+j]-z[i]*f[j])%p
  return enc(z[:d])
 # Independent monic divisor test for the certificate's field modulus.
 for k in range(1,d//2+1):
  for coeff in itertools.product(range(p),repeat=k):
   g=list(coeff)+[1];rem=f[:]
   for i in range(d,k-1,-1):
    t=rem[i]
    for j in range(k+1):rem[i-k+j]=(rem[i-k+j]-t*g[j])%p
   need(any(rem[:k]),'reducible modulus')
 digit=[digits(x) for x in range(q)]
 add=[[enc([a+b for a,b in zip(digit[x],digit[y])]) for y in range(q)] for x in range(q)]
 neg=[enc([-a for a in digit[x]]) for x in range(q)]
 Q={mul(x,x) for x in range(q)}
 counts=[]
 for c in [0,1]:
  count=0
  for u in ([1] if c==0 else range(q)):
   for v in range(q):
    offsets=[u,neg[add[u][v]],v,add[neg[u]][v],0,add[u][neg[v]],neg[v],add[u][v],neg[u]]
    z=[add[c][x] for x in offsets]
    if len(set(z))==9 and all(x in Q for x in z):
     totals=[]
     for L in LINES:
      t=0
      for k in L:t=add[t][z[k]]
      totals.append(t)
     need(len(set(totals))==1,'field line sums');count+=1
  counts.append(count)
 need(counts==[row['center_zero_u1_count'],row['center_one_count']],'normalized count mismatch')
 return counts

def checks():
 integrity(P);D=json.loads((P/'countdata.json').read_text());c=D['crt'];w=c['entries']
 need(len(set(w))==9 and all(r*r%35==x for r,x in zip(c['square_roots'],w)),'CRT roots/distinctness')
 need([sum(w[k] for k in L) for L in LINES[:3]]==[10,45,80],'ordinary row sums')
 need(all(sum(w[k] for k in L)%35==10 for L in LINES),'eight modular sums')
 pairs=list(itertools.combinations(range(9),2));masks=[]
 for n in [5,7]:
  a=[x%n for x in w];need(a==c['mod'+str(n)],'local reduction')
  masks.append(sum(1<<i for i,(u,v) in enumerate(pairs) if a[u]==a[v]))
 coordinates=[[(i//3-i%3)%3 for i in range(9)],[(i//3+i%3+1)%3 for i in range(9)]]
 for n,coord in zip([5,7],coordinates):
  need(all((w[i]%n==w[j]%n)==(coord[i]==coord[j]) for i,j in pairs),'affine coordinate fibres')
 need(len(set(zip(*coordinates)))==9,'affine coordinate orthogonality')
 need([hex(m) for m in masks]==c['masks'] and masks[0]&masks[1]==0,'collision masks')
 need(all(sum(1 for i in range(9) if w[i]%5==a and w[i]%7==b)==1 for a in set(c['mod5']) for b in set(c['mod7'])),'orthogonal partitions')
 R={n:ring(n) for n in list(range(1,101))+[128,200,243]}
 need(R[5] and R[7] and not R[35],'CRT Parker values')
 need([n for n in range(1,51) if not R[n] and R[2*n]]==[50],'smallest doubling reentry')
 need(not R[200] and all(not R[n] for n in [27,81,243,128]),'ring examples')
 need(all(R[2**k] for k in range(1,7)),'dyadic threshold controls')
 fieldcounts={r['q']:field(r) for r in json.loads((P/'sources/field_rows.json').read_text())}
 need(fieldcounts=={27:[0,0],81:[6,48],243:[0,0]},'field counts')
 requirements={'crt-shared-cells.svg':['10, 45, 80','C₃₅ = C₅ ∩ C₇ = ∅','not integer magic','cell 8'], 'fields-rings-reentry.svg':['k ≥ 7','ℤ/50ℤ','ℤ/100ℤ','ℤ/200ℤ','not integer-square existence']}
 for name,labels in requirements.items():
  xml=(P/name).read_text();e=ET.fromstring(xml);text=' '.join(e.itertext())
  need('href=' not in xml and '<script' not in xml,'nonlocal SVG content')
  for label in labels:need(label in text,'missing label '+label)
 with tempfile.TemporaryDirectory() as t:
  out=Path(t);shutil.copy(P/'generate.py',out/'generate.py');subprocess.run([sys.executable,'-B',str(out/'generate.py')],check=True)
  for name in [*requirements,'countdata.json']:need((out/name).read_bytes()==(P/name).read_bytes(),'determinism '+name)
 return {'status':'PASS','figures':2,'ring_exhaustive_checks':len(R),'field_normalized_counts':fieldcounts,'crt_collision_masks':[hex(m) for m in masks],'ordinary_row_sums':[10,45,80],'scope':'Exact illustrated finite cases and deterministic bytes; all-exponent statements depend on pinned proofs.'}
def selftest():
 with tempfile.TemporaryDirectory() as t:
  q=Path(t)/'packet';shutil.copytree(P,q);integrity(q)
  target=q/'countdata.json';original=target.read_bytes()
  for case in ['mutation','missing','extra']:
   if case=='mutation':target.write_bytes(original+b' ')
   elif case=='missing':target.unlink()
   else:(q/'rogue.txt').write_text('x')
   try:integrity(q)
   except ValueError:pass
   else:raise ValueError('accepted '+case)
   target.write_bytes(original)
   if case=='extra':(q/'rogue.txt').unlink()
 a=subprocess.check_output([sys.executable,'-B',str(P/'replay.py')]);b=subprocess.check_output([sys.executable,'-O','-B',str(P/'replay.py')]);need(a==b,'normal/-O mismatch')
if __name__=='__main__':
 result=checks()
 if '--self-test' in sys.argv:selftest();result['self_test']='PASS: mutation, missing, extra, normal/-O'
 print(json.dumps(result,sort_keys=True))
