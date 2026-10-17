#!/usr/bin/env python3
"""Independent scalar integer-clock replay of every source at terminal depth64."""
import argparse,json,hashlib,struct
from pathlib import Path
from collections import defaultdict


def need(x,m):
 if not x:raise RuntimeError(m)
def seed(x):
 while x%8==5:x=(x-1)//4
 return x

def run(path):
 data=json.loads(Path(path).read_text());out=[]
 for row in data['runs']:
  N=row['N'];depth=row['rows'][-1]['depth'];endhash=hashlib.sha256();clockhash=hashlib.sha256();mass=defaultdict(int);count=0;actual_steps=0
  for source in range(1,4*N+1,2):
   if source%8==5:continue
   h=0;m=source
   while m<=4*N:
    h+=(-3*(0 if m%3==0 else 1 if m%3==1 else -1)) if m<=2*N else (2 if m%3==2 else -1)
    m=4*m+1
   x=source;clock=0
   for d in range(depth+1):
    s=seed(x)
    if s==1:
     while x!=1:x=x//2 if x%2==0 else (3*x+1)//2;clock+=1
    if d==depth:break
    if x!=1:
     x=(3*x+1)//2;clock+=1
     while x%2==0:x//=2;clock+=1
   endhash.update(struct.pack('<q',x));clockhash.update(struct.pack('<q',clock));actual_steps+=clock;count+=1
   mass[seed(x)]+=h*(-1 if clock%2 else 1)
  last=row['rows'][-1];root=mass.pop(1,0);res=sum(abs(v) for v in mass.values())
  need(endhash.hexdigest()==last['source_endpoint_sha256'] and clockhash.hexdigest()==last['source_clock_sha256'],'every original source endpoint and actual integer clock')
  need(root==last['root_current'] and res==last['residual_l1'],'independent terminal signed mass')
  out.append(dict(N=N,source_count=count,shortcut_steps_replayed=actual_steps,root_current=root,residual_l1=res,source_endpoint_sha256=endhash.hexdigest(),source_clock_sha256=clockhash.hexdigest()))
 return dict(status='FINITE-EXACT independent scalar replay, all sources at depth64',runs=out)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('input');a=p.parse_args();print(json.dumps(run(a.input),sort_keys=True))
