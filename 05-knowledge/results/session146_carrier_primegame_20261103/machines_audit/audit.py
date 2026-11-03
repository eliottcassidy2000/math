#!/usr/bin/env python3
"""Independent interpreter and index-pair oracle for the finite FRACTRAN compiler."""
import importlib.util,itertools,json,math,sys
from fractions import Fraction
from pathlib import Path

def require(v,msg):
 if not v:raise RuntimeError(msg)
p=Path(__file__).resolve().parent.parent/'machines'/'check.py'
s=importlib.util.spec_from_file_location('producer',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)

def trial(P,vector):
 z=P['initial']
 for prime,b in zip(P['bits'],vector):
  if b:z*=prime
 trace=[]
 while True:
  for n,(a,b) in enumerate(P['rules']):
   value=Fraction(z*a,b)
   if value.denominator==1:
    trace.append(n);z=value.numerator;break
  else:break
  require(len(trace)<=4*len(P['pairs']),'nonterminating or bound failure')
 data=1
 for prime,b in zip(P['bits'],vector):
  if b:data*=prime
 require(z % data==0,'missing bit factor')
 control=z//data
 if control==P['reject']:got=None
 else:
  require(control in P['accepts'],'wrong halt control or extra data')
  got=P['pairs'][P['accepts'].index(control)]
 return got,len(trace)

def oracle(N,T,v):
 candidates=[]
 for a in range(N):
  for b in range(a,N):
   if 2*a+2*b+6==T and v[a] and v[b]:candidates.append((a,b))
 return min(candidates) if candidates else None

def main():
 count=0
 for N in range(1,9):
  for T in range(4*N+5):
   p=m.compile_program(N,T)
   for v in itertools.product((False,True),repeat=N):
    got,steps=trial(p,v);require(got==oracle(N,T,v),'truth table mismatch');count+=1
 p=m.compile_program(1,6);p['rules'][0],p['rules'][1]=p['rules'][1],p['rules'][0]
 require(trial(p,(True,))[0] is None,'priority hostile not detected')
 print(json.dumps({'status':'PASS','scope':'Independent fraction interpreter and independently enumerated first-pair oracle applied to producer-generated programs; not an independent program generator.','cases':count,'primitive_step_bound':'at most 4m','priority_hostile_detected':True,'assert_gates_used':False},indent=2,sort_keys=True))
if __name__=='__main__':main()
