#!/usr/bin/env python3
"""Exact arbitrary-gadget, arbitrary-order image of the rooted start-one class."""
from collections import deque
from pathlib import Path
import json,sys

def require(x,msg):
 if not x:raise RuntimeError(msg)

def gadgets(bound):
 # Index L is exactly the root-ending two-cover count of the physical gadget.
 out=[None]*(bound+1)
 for L in range(1,bound+1):
  if L==1:g,Q,B=1,1,1
  elif L%2==0:g,Q,oldB,old_delta=out[L//2];B=2*oldB
  else:
   ell=(L+1)//2;oldg,oldQ,oldB,old_delta=out[ell]
   g,Q,B=oldg+ell-ell%2,oldQ+oldB,2*oldB
  out[L]=(g,Q,B,int(L==1))
 return out

def word(L):
 if L==1:return ()
 t=L-1
 return tuple(1-((t>>i)&1) for i in range(t.bit_length()))

def transfer(x,L,g):
 h,e,R=x;H,Q,B,delta=g
 return (H*h-(L%2)*e+L*R,h-delta*e,Q*h+B*R)

def bounded(bound):
 gs=gadgets(bound);seed=(1,1,1);seen={seed};depth={seed:0};q=deque([seed]);values={1};edges=0;maxdepth=0
 while q:
  x=q.popleft()
  for L in range(1,bound//x[2]+1):
   y=transfer(x,L,gs[L])
   if y[0]>bound:continue
   edges+=1;values.add(y[0])
   if y[2]<=bound and y not in seen:
    seen.add(y);depth[y]=depth[x]+1;maxdepth=max(maxdepth,depth[y]);q.append(y)
 return {'bound':bound,'gadgets':bound,'states':len(seen),'edges':edges,'max_depth':maxdepth,'missing_odd':sorted(set(range(3,bound+1,2))-values),'values':sorted(values)}

if __name__=='__main__':
 B=int(sys.argv[1]) if len(sys.argv)>1 else 5001
 out=bounded(B)
 known={501:(2727,3900,9),5001:(81585,117980,13)}
 if B in known:require((out['states'],out['edges'],out['max_depth'])==known[B],'closed image changed')
 dest=Path(__file__).resolve().parent/f'unique_start_{B}.json'
 dest.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='values'},sort_keys=True,indent=2))
