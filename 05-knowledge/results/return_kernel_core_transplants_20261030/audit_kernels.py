#!/usr/bin/env python3
"""Independent Bellman audit, full-integer witnesses and 25-state paid kernels."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd
import json,subprocess,hashlib,tempfile,sys
P=Path(__file__).resolve().parent
BIN='/tmp/math-audit-excursion'
def need(ok,msg):
 if not ok:raise RuntimeError(str(msg))
def height(P0,Q,B):
 g=gcd(abs(P0-Q),B);h=0
 while g%5==0:g//=5;h+=1
 return h
def witness(row,E,scale):
 p,q,b=1,1,0;reward=0;word=row['word'];need(any(a>1 for a in word),'non-one flag retained')
 for j,a in enumerate(word):
  p,q,b=3*p,2**a*q,3*b+q;h=height(p,q,b);reward+=min(h,E)
  need((h>=E)==(j==len(word)-1),'literal primitive return')
 mod=5**(E+1);inv=pow(q,-1,mod);m=p*inv%mod;c=b*inv%mod
 need(((m-1)//5**E,c//5**E)==(row['x'],row['y']),'literal tangent endpoint')
 N=5*mod-5
 need(N*sum(word)-scale*reward==row['cost'],'literal word cost')
def kernel_costs(data):
 E=data['E'];D=data['scale'];c=[F(r['cost'],D) for r in data['rows']]
 c[24]=min(c[24],F(1,5)) # the separate all-ones primitive excursion
 # All costs are positive; exact shortest paths on the25-state additive group.
 dis=[None]*25;dis[0]=F(0);done=set()
 for _ in range(25):
  u=min((u for u in range(25) if u not in done and dis[u] is not None),key=lambda u:dis[u]);done.add(u)
  for j,cost in enumerate(c):
   v=5*((u//5+j//5)%5)+(u%5+j%5)%5;val=dis[u]+cost
   if dis[v] is None or val<dis[v]:dis[v]=val
 loop=min(dis[5*((-(j//5))%5)+(-(j%5))%5]+cost for j,cost in enumerate(c))
 need(loop==1,'exact paid zero-sum loop')
 for u in range(25):
  for j,cost in enumerate(c):
   v=5*((u//5+j//5)%5)+(u%5+j%5)%5
   need(F(v==0)+dis[v]-dis[u]<=cost,'paid tangent potential edge')
 return loop,dis

def main():
 subprocess.run(['c++','-O2','-std=c++17',str(P/'audit_excursion.cpp'),'-o',BIN],check=True)
 receipts=[];factorial_checks=0
 for E in (1,2,3,4):
  path=P/f'kernel{E}.json';data=json.loads(path.read_text());D=data['scale']
  need(len(data['rows'])==25 and all((r['x'],r['y'])==divmod(j,5) for j,r in enumerate(data['rows'])),'all ordered kernel endpoints')
  costs=','.join(str(r['cost']) for r in data['rows'])
  cmd=[BIN,str(E),str(path)+'.dist',costs,str(data['cutoff'])]
  run=subprocess.run(cmd,capture_output=True,text=True);need(run.returncode==0,run.stdout+run.stderr)
  receipt=json.loads(run.stdout)
  for row in data['rows']:witness(row,E,D)
  minrow=min(data['rows'],key=lambda r:r['cost']);minimum=F(minrow['cost'],D);need(minimum>=1,'non-one fee')
  # Full integer charge for repeated native excursions versus the factorial identity.
  for row in data['rows']:
   w=row['word'];p0,q0,b0=1,1,0;basecharge=0
   for a in w:
    p0,q0,b0=3*p0,2**a*q0,3*b0+q0;basecharge+=height(p0,q0,b0)
   endpoint=height(p0,q0,b0)
   for repeat in (2,5):
    p0,q0,b0=1,1,0;total=0
    for a in w*repeat:
     p0,q0,b0=3*p0,2**a*q0,3*b0+q0;total+=height(p0,q0,b0)
    vfact=0;k=5
    while k<=repeat:vfact+=repeat//k;k*=5
    need(total==repeat*basecharge+vfact,'exact first-return factorial charge');factorial_checks+=1
  loop,potential=kernel_costs(data)
  # Independent integer all-ones clock and next-layer tangent.
  T=2*5**(E-1);p,q=3**T,2**T;b=p-q;inv=pow(q,-1,5**(E+1))
  need(((p*inv%5**(E+1)-1)//5**E,(b*inv%5**(E+1))//5**E)==(4,4),'all-ones tangent')
  need(F(data['numerator']*T-D*((5**E-1)//4),D)==F(1,5),'all-ones fee')
  # Hash the generated lower-bound certificate without loading it all into RAM.
  h=hashlib.sha256()
  with Path(str(path)+'.dist').open('rb') as f:
   for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
  receipt.update(min_nonone=str(minimum),min_word=minrow['word'],min_clock=sum(minrow['word']),
                 min_tangent=[minrow['x'],minrow['y']],shortest_zero_sum=str(loop),potential=list(map(str,potential)),certificate_sha256=h.hexdigest())
  receipts.append(receipt);print(E,minimum,'edges',receipt['bellman_edges'],'tail',receipt['tail_margin'],flush=True)
 # Corruption control: a reachable non-root potential cannot be arbitrarily increased.
 src=P/'kernel1.json.dist';lines=src.read_text().splitlines();idx=next(i for i,s in enumerate(lines) if s.split()[:3]!=['0','1','0'] and int(s.split()[-1])<10**15)
 parts=lines[idx].split();parts[-1]=str(int(parts[-1])+1000000);lines[idx]=' '.join(parts)
 with tempfile.TemporaryDirectory(prefix='math-kernel-audit-') as tmp:
  corrupt=Path(tmp)/'corrupt.dist';corrupt.write_text('\n'.join(lines)+'\n');d=json.loads((P/'kernel1.json').read_text())
  bad=subprocess.run([BIN,'1',str(corrupt),','.join(str(r['cost']) for r in d['rows']),str(d['cutoff'])],capture_output=True,text=True)
  need(bad.returncode!=0,'corrupted potential rejected')
 # Costs need not increase with refinement: fourth layer is a required hostile.
 need(F(receipts[3]['min_nonone'])<F(receipts[2]['min_nonone']),'nonmonotone fee control')
 out=dict(status='independently checked finite certificates; all-word minima at E=1,2,3,4',rows=receipts,
          native_witnesses=100,tangent_potential_edges=2500,factorial_checks=factorial_checks,corruption_rejected=True)
 (P/'kernel_audit.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
