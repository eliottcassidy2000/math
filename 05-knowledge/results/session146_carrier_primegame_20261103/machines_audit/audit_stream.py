#!/usr/bin/env python3
"""Literal integer FRACTRAN interpreter, independent of valuation-vector runner."""
import importlib.util,json,math
from pathlib import Path

def need(v,msg):
 if not v:raise RuntimeError(msg)
H=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('producer',H.parent/'machines'/'check.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
head=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11))
lists={'owner':head+((15,2),(1,7),(55,1)),'guy':head+((15,14),(15,2),(55,1))}
def direct(N,which):
 z=2;bits=[False]*N;cap=2*N+2
 for step in range(2000000):
  if z>0 and z&(z-1)==0:
   e=z.bit_length()-1
   if e>=3 and e%2 and (e-3)//2<N:bits[(e-3)//2]=True
  n=z;e=0
  while n%2==0:e+=1;n//=2
  if e>=cap and all(z%q for q in (3,5,11,13,17,19,23,29)):
   return bits,step
  for a,b in lists[which]:
   if z%b==0:z=z//b*a;break
  else:raise RuntimeError('unexpected halt')
 raise RuntimeError('budget exhausted; inconclusive')
def prime(n):
 return n>=2 and not any(n%d==0 for d in range(2,math.isqrt(n)+1))
for variant in ('owner','guy'):
 need(tuple(m.PG_LISTS[variant])==lists[variant],'list differs from inherited source')
 for N in range(1,6):
  bits,steps=direct(N,variant);actual,steps2=m.primegame_bits(N,variant)
  need(bits==actual==[prime(2*i+3) for i in range(N)],'stream discrepancy')
  need(steps==steps2,'stop step differs')
print(json.dumps({'status':'PASS','integer_interpreter_cases':10,'lists':'both 14-rule published lists','prime_oracle':'trial division','stop':'first qualifying state-free configuration at or past 2N+2','assert_gates_used':False},indent=2,sort_keys=True))
