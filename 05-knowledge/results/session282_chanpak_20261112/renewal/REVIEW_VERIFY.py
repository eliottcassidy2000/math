#!/usr/bin/env python3
"""Independent backward-DP and matrix audit of receiver inversion."""
import importlib.util,json
from pathlib import Path
from itertools import permutations
P=Path(__file__).resolve().parent

def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
f=load('fiber',P/'receiver_fiber.py');w=load('wildcard',P.parent/'wildcard/verify.py')
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
C3=(2,4,1);D=(2,4,9,3);controls=0
for core in (C3,D):
 for a in range(3):
  for b in range(3-a):
   for p in range(1<<a):
    for q in range(1<<b):
     for cut in range(1<<len(core)):
      rec=dict(a=a,b=b,p=p,q=q,cut=cut);T=f.witness(core,rec);got=w.packet(T);c=w.packet(f.rooted(core,cut));h=c['h'];Q=c['Q']
      check(got['V']==h*(1+p)+(1<<a)*(c['V']-h+Q*q),'literal receiver law')
      check(got['H'] in {v['profile']['H'] for v in f.fiber(core,got['V'])},'inverse includes literal input');controls+=1
for r in range(10001):
 vals={v['profile']['H'] for v in f.fiber(C3,3*(r+1))};expected={3*f.unit(r)['H']+x for x in f.triangle_deltas(r)}
 check(vals==expected,'triangle exact defect support');check(not any({h,h+2,h+4}<=vals for h in vals),'triangle no-three control')
records=f.fiber(D,1906);check({r['profile']['H'] for r in records}=={195,741,819,821,823,825},'globalh5 fiber packet')
selected=[]
for target in (819,821,823,825):
 rec=next(r for r in records if r['profile']['H']==target);T=f.witness(D,rec);got=w.packet(T)
 check((got['H'],got['h'],got['V'])==(target,5,1906),'literal fourbank');selected.append(dict(H=target,order=len(T),packet=got))
# Exhaustive h5 strong core classification through order5.
def strong(T):
 n=len(T)
 for v in range(n):
  S=1<<v
  while True:
   S2=S
   for i in range(n):
    if S>>i&1:S2|=T[i]
   if S2==S:break
   S=S2
  if S!=(1<<n)-1:return False
 return True
counts={};representative=None
for n in (3,4,5):
 edges=[(i,j) for i in range(n) for j in range(i+1,n)];num=0
 for mask in range(1<<len(edges)):
  T=[0]*n
  for j,(a,b) in enumerate(edges):
   if mask>>j&1:T[a]|=1<<b
   else:T[b]|=1<<a
  T=tuple(T)
  if strong(T) and w.hp_subsets(T)[0][-1]==5:
   num+=1
   check(n==4,'noother orderstrongh5')
   check(any(all(((T[i]>>j)&1)==((D[pp[i]]>>pp[j])&1) for i in range(n) for j in range(n)) for pp in permutations(range(n))),'uniqueh5 isoclass')
 counts[n]=num
# Eight-bank: independent exact core packets + independently multiplied matrices.
def ordinal(A,B):
 n=len(A);return tuple(x|(((1<<len(B))-1)<<n) for x in A)+tuple(x<<n for x in B)
def unit(p,n):
 M=[[int(i==j) for j in range(4)] for i in range(4)]
 for i in range(n):M=w.mul(M,w.mat(w.packet((0,1) if p>>i&1 else (2,0))))
 return M
A=(6,12,24,17,3);E=(44,1,2,6,15,30);BC=ordinal(C3,D);cores={'A':A,'E':E,'BC':BC}
rows=[(54365,'A',8,57,7,4),(54367,'A',8,57,11,4),(54369,'BC',10,313,96,0),(54371,'A',2,1,3,291),(54373,'A',2,1,5,291),(54375,'E',4,9,57,40),(54377,'A',0,0,7,1165),(54379,'A',0,0,11,1165)]
outs=[]
for H,name,a,p,cut,q in rows:
 core=w.packet(f.rooted(cores[name],cut));out=w.decode(w.mul(w.mul(unit(p,a),w.mat(core)),unit(q,q.bit_length())))
 check((out['H'],out['h'],out['V'])==(H,15,139878),'eightbank core packet assembly');outs.append(out)
print(json.dumps(dict(status='PASS',literal_padding_inputs=controls,triangle_receiver_codes=10001,V1906_witnesses=selected,strong_count5_labelled_census=counts,eightbank_packets=outs,eightbank_audit_scope='independent source-core DP and matrix assembly, not full18vertex DP'),indent=2,sort_keys=True))
