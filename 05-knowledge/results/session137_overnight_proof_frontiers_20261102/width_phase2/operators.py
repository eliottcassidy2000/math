#!/usr/bin/env python3
"""Positive unary width-three modular operator probe; standard-library only."""
from math import factorial,gcd
from collections import deque
import json,hashlib

def require(ok,msg):
 if not ok:raise RuntimeError(msg)

def precisions(m):
 out=[];c=1
 while True:
  d=factorial(c)*factorial(c-1)
  if d%m==0:return tuple(out)
  out.append(m//gcd(m,d));c+=1

def step(state,mods,op):
 p=(0,)+state+(0,)
 def f(k):return p[k] if 0<=k<len(p) else 0
 out=[]
 for j,mod in enumerate(mods,1):
  if op=='P':v=j*f(j)+f(j-1)
  elif op=='C':v=f(j-2)+(2*j-1)*f(j-1)+j*(j+2)*f(j)+j*(j+1)*f(j+1)
  elif op=='D':v=f(j-3)+3*(j-1)*f(j-2)+(3*j*j+j-3)*f(j-1)+j*(j*j+6*j-2)*f(j)+j*(j+1)*(2*j+1)*f(j+1)
  else:raise ValueError(op)
  out.append(v%mod)
 return tuple(out)

def evaluate_word(word,m):
 mods=precisions(m);state=(1,)+(0,)*(len(mods)-1)
 for op in word:state=step(state,mods,op)
 return state

def search(m,maxstates=20000):
 mods=precisions(m);one=(1,)+(0,)*(len(mods)-1)
 states=[one];parents=[(-1,'')];seen={one:0};hits={one[0]:0}
 target=set(range(m)) if m%2 else set(range(1,m,2))
 q=0;terminal=''
 while q<len(states) and len(hits)<len(target):
  state=states[q]
  for op in 'PCD':
   nxt=step(state,mods,op)
   if nxt not in seen:
    if len(states)>=maxstates:terminal='budget';break
    i=len(states);seen[nxt]=i;states.append(nxt);parents.append((q,op));hits.setdefault(nxt[0],i)
  if terminal:break
  q+=1
 if not terminal:terminal='covered' if len(hits)==len(target) else 'complete orbit'
 def word(i):
  out=[]
  while i:
   i,op=parents[i];out.append(op)
  return ''.join(reversed(out))
 witnesses={h:word(i) for h,i in sorted(hits.items())}
 for h,w in witnesses.items():require(evaluate_word(w,m)[0]==h,'reconstructed word')
 return {'modulus':m,'precisions':mods,'status':terminal,'states_seen':len(states),'residues':sorted(hits),
         'missing':sorted(target-set(hits)),'max_word_length':max(map(len,witnesses.values())),
         'witnesses':witnesses}

def main():
 moduli=list(range(2,129))+[243,256,343,625,729,1024,2187,3125,4096]
 rows=[]
 for m in moduli:
  r=search(m)
  if m>32:
   raw=json.dumps(r.pop('witnesses'),sort_keys=True).encode()
   r['witnesses_sha256']=hashlib.sha256(raw).hexdigest()
  rows.append(r)
 print(json.dumps({'status':'PASS','scope':'Positive exact unary witnesses only; no all-modulus inference','rows':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()
