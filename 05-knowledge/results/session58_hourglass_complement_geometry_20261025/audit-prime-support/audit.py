#!/usr/bin/env python3
"""Independent exact audit; only local audit output is written."""
import sys,json,hashlib,subprocess
from math import isqrt,gcd,prod
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve().parent.parent/'hourglass-prime-support'
sys.dont_write_bytecode=True
import importlib.util
spec=importlib.util.spec_from_file_location('candidate',SOURCE/'check.py');c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
def need(v,msg):
 if not v: raise RuntimeError(msg)
B=3000
spf=list(range(B+1))
for p in range(2,isqrt(B)+1):
 if spf[p]==p:
  for n in range(p*p,B+1,p):
   if spf[n]==n:spf[n]=p
def fac(n):
 d={}
 while n>1:
  p=spf[n];d[p]=d.get(p,0)+1;n//=p
 return d
def divisors(n):
 ds=[1]
 for p,a in fac(n).items():ds=[x*p**b for x in ds for b in range(a+1)]
 return ds
def valuation(x,p):
 n=0
 while x%p==0:x//=p;n+=1
 return n
def endpoints(e):
 squares={b*b:b for b in range(1,isqrt(2*e*e)+1)}
 return {e*e-a*a:(a,squares[2*e*e-a*a]) for a in range(1,e) if 2*e*e-a*a in squares}
def pyth(e):
 return {2*u*isqrt(e*e-u*u) for u in range(1,isqrt(e*e//2)+1) if (v:=isqrt(e*e-u*u))>u and u*u+v*v==e*e}
D={};all_divisor_checks=triple_checks=0;counts={};gauss={}
for e in range(1,B+1):
 D[e]=endpoints(e);need(set(D[e])==pyth(e),f'independent square enumerations {e}')
 f=fac(e);need(f==c.factor(e),'independent factorization')
 split={p:a for p,a in f.items() if p%4==1};m=prod(p**a for p,a in f.items() if p%4!=1)
 need(len(D[e])==(prod(2*a+1 for a in split.values())-1)//2,'count')
 G,ps,cm=c.gaussian(e);gauss[e]=G;need(set(G)==set(D[e]) and cm==m,'Gaussian bijection')
 need(set(D[e])=={m*m*d for d in D[e//m]},'stripping')
 for t in divisors(e):
  need({d for d in D[e] if d%(t*t)==0}=={t*t*d for d in D[e//t]},'all divisor heredity');all_divisor_checks+=1
 for d,j in G.items():
  for (p,a),k in zip(ps,j):need(valuation(d,p)==2*(a-abs(k)) if k else valuation(d,p)>=2*a,'valuation')
 passed=relations=0
 for ds in combinations(sorted(D[e]),3):
  js=[G[d] for d in ds];M=[max(abs(j[i]) for j in js) for i in range(len(ps))]
  t=prod(p**(a-mi) for (p,a),mi in zip(ps,M))
  rootgcd=e
  for d in ds:
   for r in D[e][d]:rootgcd=gcd(rootgcd,r)
  need(rootgcd==m*t,'exact common root divisor');triple_checks+=1
  ok=all(sum(abs(j[i])==mi for j in js)>=2 for i,mi in enumerate(M))
  passed+=ok
  if ok:
   reduced=e//rootgcd;rd=[d//(rootgcd*rootgcd) for d in ds]
   need(all(gcd(reduced,gcd(x,y))==1 for x,y in combinations(rd,2)),'primitive pair gcd filter')
  if ds[0]+ds[1]==ds[2]:relations+=1;need(ok,'relation must pass')
 need(relations==0,'bounded nonexistence')
 if e in (65,1105):counts[e]={'triples':len(list(combinations(D[e],3))),'filter_pass':passed}
# Independent projective lines as kernels of nonzero covectors over F_3.
def canonical(v):
 w=tuple(x%3 for x in v);lead=next(x for x in w if x)
 return tuple(x*lead%3 for x in w)
pts={canonical(v) for v in product(range(3),repeat=3) if any(v)}
lines={frozenset(v for v in pts if sum(a*b for a,b in zip(v,h))%3==0) for h in pts}
need(len(pts)==13 and len(lines)==13 and all(len(L)==4 for L in lines),'dual projective census')
need(all(sum(u in L and v in L for L in lines)==1 for u,v in combinations(pts,2)),'dual projective pair incidence')
G=gauss[1105];bypt={canonical(j):d for d,j in G.items()}
need(set(bypt)==pts,'projective bijection independently realized')
for S in product((False,True),repeat=3):
 t=prod(p for p,z in zip((5,13,17),S) if z)
 coordinate={d for j,d in bypt.items() if all(not z or j[i]==0 for i,z in enumerate(S))}
 need(coordinate=={t*t*d for d in D[1105//t]},'all coordinate subsets arithmetic deletion')
need({bypt[j] for j in pts if j[2]==0}=={582624,867000,1068144,1172184},'coordinate line hostile')
# Explicit negative controls and assertion-independent failures.
need(valuation(3000,5)==3 and gauss[65][3000][0]==0,'balanced strict inequality')
need(len(D[125])==3 and not any(x+y==z for x,y,z in combinations(sorted(D[125]),3)),'cardinality hostile')
need(counts=={65:{'triples':4,'filter_pass':4},1105:{'triples':286,'filter_pass':150}},'filter hostile counts')
for n in (0,-1):
 try:c.val(n,5)
 except RuntimeError:pass
 else:raise RuntimeError('valuation nonpositive guard failed')
try:c.require(False,'injected negative control')
except RuntimeError:pass
else:raise RuntimeError('require guard disabled')
need(sum(len(v) for v in D.values())==3172,'gap total')
# Replay original program without modifying original certificate or cache.
copy=ROOT/'candidate-replay.py';copy.write_bytes((SOURCE/'check.py').read_bytes())
replays=[]
for flag in ([],['-O']):
 r=subprocess.run([sys.executable,*flag,str(copy)],capture_output=True,text=True)
 need(r.returncode==0,'candidate replay failed');j=json.loads(r.stdout)
 need(j==json.loads((SOURCE/'certificate.json').read_text()),'certificate mismatch')
 replays.append({'mode':'optimized' if flag else 'normal','exit':r.returncode,'certificate_matches':True})
# Historical primary-source facts are retained separately from this arithmetic replay.
meta=json.loads((ROOT/'SOURCE_STATUS.json').read_text())
need(meta['metadata']['resource_type']['subtype']=='preprint','preprint status')
need(meta['metadata']['publication_date']=='2026-08-13','source date')
need(meta['doi']=='10.5281/zenodo.21910218','DOI')
source_pdf_rechecked=False
if '--source-dir' in sys.argv:
 source_dir=Path(sys.argv[sys.argv.index('--source-dir')+1])
 need(hashlib.md5((source_dir/'MagicHourGlassPrimes.pdf').read_bytes()).hexdigest()==meta['recorded_pdf_md5'],'PDF checksum')
 source_pdf_rechecked=True
out={'status':'PASS','source_pdf_rechecked':source_pdf_rechecked,'bound':B,'positive_gaps':sum(len(v) for v in D.values()),'all_divisor_checks':all_divisor_checks,'triple_root_gcd_checks':triple_checks,'tropical_counts':counts,'replays':replays,'negative_controls':'PASS','source_type':meta['metadata']['resource_type'],'source_date':meta['metadata']['publication_date'],'sha256_candidate':hashlib.sha256((SOURCE/'check.py').read_bytes()).hexdigest()}
print(json.dumps(out,indent=2))
