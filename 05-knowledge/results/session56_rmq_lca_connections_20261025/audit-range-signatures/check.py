#!/usr/bin/env python3
"""Independent explicit-exception audit; no assert-based gates or producer writes."""
import sys,json,hashlib,importlib.util,itertools,math
from pathlib import Path
from collections import Counter
ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent/'range-signatures'; OUT=Path(__file__).parent
spec=importlib.util.spec_from_file_location('packet',ROOT/'verify.py'); p=importlib.util.module_from_spec(spec);sys.modules['packet']=p;spec.loader.exec_module(p)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def mask(a):return sum(1<<k for k,(i,j) in enumerate(itertools.combinations(range(9),2)) if a[i]==a[j])
lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def valid(a,n):return len({sum(a[i] for i in ln)%n for ln in lines})==1 and all(x in {t*t%n for t in range(n)} for x in a)
def census(n):
 # Solve from arbitrary top row and center; no producer parametrization calls.
 Q={t*t%n for t in range(n)}; centers={}
 for e in Q:centers.setdefault(3*e%n,[]).append(e)
 h=Counter()
 for a,b,c in itertools.product(Q,repeat=3):
  s=(a+b+c)%n
  for e in centers.get(s,[]):
   g=(s-c-e)%n;i=(s-a-e)%n;d=(s-a-g)%n;f=(s-c-i)%n;h0=(s-b-e)%n
   A=(a,b,c,d,e,f,g,h0,i)
   if all(x in Q for x in A) and len({sum(A[j] for j in ln)%n for ln in lines})==1:h[mask(A)]+=1
 return h

def minset(F):return frozenset(x for x in F if all(not(y!=x and y&x==y) for y in F))
def conv(A,B):
 H=Counter()
 for a,c in A.items():
  for b,d in B.items():H[a&b]+=c*d
 return H
need(hashlib.sha256((ROOT/'sources/session55_verify.py').read_bytes()).hexdigest()=='eaf0b74d374c56971fe67d90c79f8c9d5b8deaf9e26c45f3ed0466f4aa33b5f3','source pin')
cert=json.loads((ROOT/'certificate.json').read_text()); cache={}; globals_=[]; tested=0
for row in cert['rows']:
 mods=row['mods'][row['interval'][0]:row['interval'][1]+1]; n=math.prod(mods)
 if n not in cache:cache[n]=census(n)
 expected=cache[n]; claimed=Counter({int(k,16):v for k,v in row['histogram'].items()})
 need(expected==claimed,f'independent full weighted enumeration {mods}')
 need(sum(expected.values())==row['total_arrays'],f'total {mods}')
 need(sorted(minset(expected))==sorted(map(lambda s:int(s,16),row['minimal_masks'])),f'minima {mods}')
 leaves=[p.local(q,('component',i,q)) for i,q in enumerate(row['mods'])]
 dst=p.DisjointSparseTable(leaves,p.merge);seg=p.SegmentTree(leaves,p.merge,p.IDENTITY);l,r=row['interval']
 a,b=dst.query(l,r),seg.query(l,r)
 need(a.hist==b.hist==expected,'interval implementation')
 need(a.ids==b.ids==tuple(x.ids[0] for x in leaves[l:r+1]),'component order')
 need(set(a.witnesses)==set(expected),'complete witness supports')
 for bits,local_arrays in a.witnesses.items():
  need(len(local_arrays)==len(mods),'provenance size')
  # Incremental CRT, independent from producer closed-form CRT.
  output=[]
  for column in zip(*local_arrays):
   x=0;M=1
   for residue,q in zip(column,mods):x+=M*((residue-x)*pow(M,-1,q)%q);M*=q
   output.append(x)
  need(valid(output,n) and mask(output)==bits,'global admissibility/exact mask')
  need(all(tuple(x%q for x in output)==tuple(a) and valid(a,q) for q,a in zip(mods,local_arrays)),'local witness match')
  globals_.append({'components':mods,'modulus':n,'mask':hex(bits),'entries':output,'square_roots':[next(t for t in range(n) if t*t%n==x) for x in output]})
 tested+=1
need(tested==31,'interval count')
F5=census(5);dup=conv(conv(census(2),F5),conv(F5,census(3))); true=census(30)
need(sum(true.values())==112 and sum(dup.values())==784 and 0 not in true and 0 in dup,'hostile')
need(minset(true)==frozenset((0x1000880,0x80040200)),'hostile minima')
need(0 not in census(25),'mod25 control')
for n in range(1,81):
 leaves=[(j,) for j in range(n)];dst=p.DisjointSparseTable(leaves,lambda a,b:a+b);seg=p.SegmentTree(leaves,lambda a,b:a+b,())
 for l in range(n):
  for r in range(l,n):need(dst.query(l,r)==seg.query(l,r)==tuple(range(l,r+1)),'noncommutative order')
for leaf in [p.local(2,'x'),p.local(5,'x')]:
 try:p.merge(leaf,leaf)
 except ValueError:pass
 else:raise RuntimeError('duplicate identity accepted')
try:p.OverlapSparseTable([p.local(q,i) for i,q in enumerate((2,5,3))],p.merge).query(0,2)
except ValueError:pass
else:raise RuntimeError('overlap accepted')
families=[frozenset(i for i in range(8) if f>>i&1) for f in range(256)]
for F in families:
 need((p.meet(F,F)==F)==all(a&b in F for a in F for b in F),'meet closure iff')
 A=minset(F);need((p.anti(A,A)==A)==(len(A)<=1),'antichain iff')
 for G in families:need(minset(p.meet(F,G))==p.anti(minset(F),minset(G)),'compression')
for w in itertools.product(range(4),repeat=4):
 h=Counter({i:v for i,v in enumerate(w) if v});need((conv(h,h)==h)==(sum(w)<=1),'weighted iff')
result={'status':'PASS','explicit_exception_gates':True,'intervals':tested,'distinct_moduli':sorted(cache),'full_global_witness_count':len(globals_),'noncommutative_queries':sum(n*(n+1)//2 for n in range(1,81)),'source_sha256':hashlib.sha256((ROOT/'verify.py').read_bytes()).hexdigest(),'certificate_sha256':hashlib.sha256((ROOT/'certificate.json').read_bytes()).hexdigest(),'hostile_true_count':112,'hostile_overlap_count':784}
(OUT/'global-witnesses.json').write_text(json.dumps(globals_,indent=2)+'\n');(OUT/'audit-certificate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
