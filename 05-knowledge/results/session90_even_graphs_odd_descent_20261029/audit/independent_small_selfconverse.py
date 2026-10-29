#!/usr/bin/env python3
"""Audit minimal odd nonflat order from ordinary tournament banks, not pair lifts.
Uses hash-pinned previously audited S87/S88 generator; literal permutations
independently find all automorphisms and anti-involutions for each ordinary type.
"""
import argparse,sys,types,json,itertools,hashlib
from pathlib import Path
from fractions import Fraction
from collections import Counter
from independent_holonomy import mul,inv,power
FILES={'generation':('05-knowledge/results/session87_four_vertices_enumeration_20261028/check_generation.py','f5d3da961f34dc1a32eaa0c1b69fad9bda983154578ea0187c5fd6b2e2d10fae'),'blockfill':('05-knowledge/results/session88_symmetry_certificates_20261028/blockfill/check_blockfill.py','cb246f280786f6a7f8b1f829d9d9b9e0993cb2298eb399783809bedb27780bdb')}
CHECKS=0
def ck(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
def load(root):
 mods={}
 for name,(path,digest) in FILES.items():
  raw=(root/path).read_bytes();ck(hashlib.sha256(raw).hexdigest()==digest,'pinned generator hash')
  m=types.ModuleType(name);m.__file__=str(root/path)
  if name=='blockfill':
   b=types.ModuleType('baseline_sources');b.generation=mods['generation'];sys.modules['baseline_sources']=b
  exec(compile(raw,str(root/path),'exec'),m.__dict__);mods[name]=m
 return mods['blockfill']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--baseline-root',type=Path);args=ap.parse_args()
 root=args.baseline_root
 if root is None:
  candidates=list(Path(__file__).resolve().parents)+[Path(__file__).resolve().parents[1]/'repo']
  root=next((x for x in candidates if (x/FILES['generation'][0]).is_file()),None)
 ck(root is not None,'supply --baseline-root repository')
 block=load(root);bank=[()];out=[]
 for n in range(1,8):
  bank=list(block.build(bank,n-1,'blockfill')['bank'])
  if n%2==0:continue
  ps=list(itertools.permutations(range(n)));one=tuple(range(n));counts=Counter();SC=0;triangles=0;mass=Fraction();marked=Fraction()
  for W in bank:
   ds=[row.bit_count() for row in W];A=[];I=[]
   for p in ps:
    if all(ds[i]==ds[p[i]] for i in range(n)) and all(((W[i]>>j)&1)==((W[p[i]]>>p[j])&1) for i in range(n) for j in range(i)):
     A.append(p)
    if all(p[p[i]]==i and ds[i]+ds[p[i]]==n-1 for i in range(n)) and all(((W[i]>>j)&1)==((W[p[j]]>>p[i])&1) for i in range(n) for j in range(i)):
     I.append(p)
   if not I:continue
   SC+=1;tau=I[0];C=[a for a in A if mul(a,tau)==mul(tau,a)];K={mul(s,tau) for s in I}
   ck(len(A)==len(I)*len(C),'literal all-anti centralizer count')
   ck(all(mul(x,y)==mul(y,x) for x in K for y in K),'all small K commute')
   g={(t,s):power(mul(t,s),(len(A)+1)//2) for t in I for s in I}
   for u,t,s in itertools.product(I,repeat=3):
    ck(mul(g[u,t],g[t,s])==g[u,s],'literal flat triangles');triangles+=1
   mass+=Fraction(1,len(A));marked+=Fraction(len(I),len(A));counts[len(A),len(I),len(C)]+=1
  ck(SC=={1:1,3:2,5:8,7:88}[n],'self-converse count')
  out.append(dict(order=n,ordinary_classes=len(bank),selfconverse_classes=SC,triangles=triangles,ordinary_mass=str(mass),marked_mass=str(marked),symmetry_profiles=[[*k,v] for k,v in sorted(counts.items())]))
 print(json.dumps(dict(status='PASS',checks=CHECKS,source='ordinary S87/S88 banks; no pair lift used',rows=out),indent=2))
if __name__=='__main__':main()
