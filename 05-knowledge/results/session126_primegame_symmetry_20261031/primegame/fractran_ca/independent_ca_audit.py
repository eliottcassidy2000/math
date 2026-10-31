#!/usr/bin/env python3
"""Independent auditor: own Fraction semantics, codec, observer, dense sweep.
No use of the production macro driver, observer, head_delta, or scheduler in
macro replays. The implemented local_rule is the sole transition under test.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import prod
from pathlib import Path
from random import Random
import json
import primegame_ca as ca

ROOT=Path(__file__).resolve().parent
PS=(2,3,5,7,11,13,17,19,23,29)
RACK=tuple(Fraction(x) for x in ('17/91','78/85','19/51','23/38','29/33','77/29','95/23','77/19','1/17','11/13','13/11','15/2','1/7','55/1'))
T=4096;MARK=1024;PASSIVE=2048;READY=1
checks=0

def check(x,why):
 global checks
 checks+=1
 if not x:raise RuntimeError(why)

def factor(n):
 exp=[]
 for p in PS:
  a=0
  while n%p==0:n//=p;a+=1
  exp.append(a)
 return exp,n

def integer_step(n):
 for i,f in enumerate(RACK):
  y=f*n
  if y.denominator==1:return y.numerator,i
 raise RuntimeError('missing final instruction')

# Construct fixed finite head-state numbering independently from fractions.
OPS=[];HEADS=[('READY',-1,-1)]
for i,f in enumerate(RACK):
 den,r=factor(f.denominator);check(r==1,'den support')
 num,r=factor(f.numerator);check(r==1,'num support')
 check(all(x<=1 for x in den+num),'not squarefree')
 ops=[(j,-1) for j,e in enumerate(den) if e]+[(j,1) for j,e in enumerate(num) if e]
 OPS.append(ops)
 for j,(_,sign) in enumerate(ops):
  HEADS.append(('SEEK',i,j))
  if sign<0:HEADS.append(('DECLEFT',i,j))
  HEADS.append(('RETURN',i,j))
CODES={s:i+1 for i,s in enumerate(HEADS)}
ERROR=T*(1+len(HEADS))
check(len(HEADS)==93 and ERROR+1==385025,'finite alphabet cardinality')
check(tuple(HEADS)==ca.HEAD_STATES,'head-state enumeration disagreement')
check(tuple(tuple(x) for x in OPS)==ca.OPS,'operation enumeration disagreement')

def head_ref(data,code):
 phase,i,j=HEADS[code-1]
 origin=bool(data&MARK)
 if phase=='READY':
  if not origin:return None
  for i,f in enumerate(RACK):
   den,_=factor(f.denominator)
   if all(not e or data&(1<<k) for k,e in enumerate(den)):break
  return data,CODES['SEEK',i,0],0
 k,sign=OPS[i][j];bit=1<<k
 nextcode=READY if j+1==len(OPS[i]) else CODES['SEEK',i,j+1]
 if phase=='RETURN':
  return (data,nextcode,0) if origin else (data,code,-1)
 if phase=='DECLEFT':
  if sign!=-1 or not data&bit:return None
  return (data^bit,nextcode,0) if origin else (data^bit,CODES['RETURN',i,j],-1)
 if data&bit:return data,code,1
 if sign==-1:
  return None if origin else (data,CODES['DECLEFT',i,j],-1)
 return (data|bit,nextcode,0) if origin else (data|bit,CODES['RETURN',i,j],-1)

def local_ref(l,c,r):
 cells=[l,c,r]
 if ERROR in cells:return ERROR
 decoded=[divmod(x,T) for x in cells]
 heads=[i for i,(h,d) in enumerate(decoded) if h]
 if not heads:return c
 if len(heads)>1:return ERROR
 loc=heads[0];h,d=decoded[loc]
 v=head_ref(d,h)
 if v is None:return ERROR
 d2,h2,direction=v
 if loc==1:return d2+(T*h2 if direction==0 else 0)
 if loc+direction==1:return decoded[1][1]+T*h2
 return c

# Exhaust the 380928 head/tape pairs, checking both own head semantics and all
# three possible relative head positions in a fixed tape neighborhood.
head_pairs=0
for h in range(1,94):
 for data in range(T):
  ref=head_ref(data,h)
  check(ca.head_delta(data,h)==ref,'finite head semantics')
  for loc in range(3):
   triple=[37,73,19];triple[loc]=data+T*h
   check(ca.local_rule(*triple)==local_ref(*triple),'single-head local lift')
  head_pairs+=1

# Own dense full-radius sweep, ignoring production scheduler/candidate sets.
def full_sweep(cfg):
 if not cfg:return {}
 lo,hi=min(cfg)-1,max(cfg)+1
 answer={}
 for x in range(lo,hi+1):
  v=ca.local_rule(cfg.get(x-1,0),cfg.get(x,0),cfg.get(x+1,0))
  if v:answer[x]=v
 return answer

def encode(n,origin=0):
 exp,r=factor(n);cfg={origin:MARK+T*READY}
 for k,a in enumerate(exp):
  for x in range(origin,origin+a):cfg[x]=cfg.get(x,0)|(1<<k)
 for j in range(r.bit_length()):
  if r&(1<<j):cfg[origin-1-j]=PASSIVE
 return cfg

def decode_ready(cfg,origin):
 # All READY boundaries must match this independent canonical shape exactly.
 heads=[(x,v//T) for x,v in cfg.items() if v//T]
 if heads!=[(origin,READY)]:return None
 exp=[]
 for k in range(10):
  locs=sorted(x-origin for x,v in cfg.items() if v%T&(1<<k))
  check(locs==list(range(len(locs))),'nonunary READY track')
  exp.append(len(locs))
 spectator=0
 for x,v in cfg.items():
  if v%T&PASSIVE:
   check(x<origin,'spectator crossed marker');spectator|=1<<(origin-1-x)
 n=spectator*prod(p**a for p,a in zip(PS,exp))
 check(cfg==encode(n,origin),'noncanonical READY state')
 return n

micro_replay=0
def macro(cfg,n,origin=0):
 global micro_replay
 y,i=integer_step(n);exp,_=factor(n)
 touched=[k for k,sgn in OPS[i]]
 expected_time=1+len(touched)+2*sum(exp[k] for k in touched)
 check(expected_time<=10*max(exp)+6,'duration bound')
 for t in range(1,expected_time+1):
  cfg=full_sweep(cfg);micro_replay+=1
  check(ERROR not in cfg.values(),'valid encoding developed error')
  z=decode_ready(cfg,origin)
  if z is not None:
   check(t==expected_time,'premature READY')
   check(z==y,'Fraction/CA mismatch')
   return cfg,y,t,i
 raise RuntimeError('missing expected READY boundary')

# Full literal sweep to first five prime-output events, counting every pure
# power rather than silently discarding unexpected nonprime exponents.
cfg=encode(2);n=2;time=0;events=[]
for tick in range(1,2376):
 cfg,n,duration,i=macro(cfg,n);time+=duration
 if n>1 and not n&(n-1):events.append([n.bit_length()-1,tick,time])
check(events==[[2,19,130],[3,69,537],[5,281,2677],[7,710,7942],[11,2375,34603]],'milestones')

# Every ordered instruction, special cycle, arbitrary positive integer inputs.
for i,f in enumerate(RACK):
 cfg,y,t,j=macro(encode(f.denominator),f.denominator)
 check(j==i,'instruction not covered')
cycle=encode(55);cfg=cycle;n=55;durations=[]
for _ in range(2):cfg,n,t,i=macro(cfg,n);durations.append(t)
check(cfg==cycle and durations==[5,5],'55/65 CA cycle length')
seed_cases=[]
for seed in (1,7,15,16,17,23,25,55,343,60,31,37,31*37,2**2*11**3*17**2*31):
 cfg=encode(seed,-4);n=seed
 for _ in range(4):cfg,n,t,i=macro(cfg,n,-4)
 seed_cases.append(seed)
rng=Random(1262026)
for _ in range(80):
 exp=[rng.randrange(6) for _ in PS]
 spectator=rng.choice((1,31,37,31*37,41**2,43*47))
 n=spectator*prod(p**e for p,e in zip(PS,exp));origin=rng.randrange(-15,16)
 cfg=encode(n,origin)
 for _ in range(2):cfg,n,t,i=macro(cfg,n,origin)

# Malformed yet in-alphabet configurations: total local semantics, scheduler
# equivalence for canonical finite-support maps, and translation equivariance.
for _ in range(10000):
 triple=[rng.randrange(ERROR+1) for _ in range(3)]
 check(ca.local_rule(*triple)==local_ref(*triple),'arbitrary triple semantics')
for _ in range(500):
 cfg={x:rng.randrange(ERROR+1) for x in range(-8,9) if rng.randrange(2)}
 cfg={x:v for x,v in cfg.items() if v}
 out=full_sweep(cfg)
 check(out==ca.ca_step_full(cfg),'production full sweep mismatch')
 check(out==ca.ca_step(cfg),'sparse scheduler mismatch')
 shift=rng.randrange(-10,11)
 shifted={x+shift:v for x,v in cfg.items()}
 check(full_sweep(shifted)=={x+shift:v for x,v in out.items()},'translation equivariance')
for locs in ((0,1),(0,2),(-1,1)):
 cfg={x:MARK+T*READY for x in locs}
 check(ERROR in full_sweep(cfg).values(),'head collision error')
check(full_sweep({0:ERROR})=={-1:ERROR,0:ERROR,1:ERROR},'error propagation')
# Zero-valued keys are redundant storage, not extra CA states. Semantically
# sparse/full agree; dictionary equality needs canonical map convention.
zero_caveat=ca.ca_step({0:0})!=ca.ca_step_full({0:0})
result={'status':'PASS','method':'Own Fraction arithmetic, canonical encoder/READY observer, dense min-to-max radius-one sweep; separate head/local oracle; no production macro-driver use','alphabet_size':ERROR+1,'head_states':len(HEADS),'head_tape_pairs_exhausted':head_pairs,'single_head_neighborhoods_checked':3*head_pairs,'seed2_events_prime_macro_micro':events,'cycle55_macro_durations':durations,'literal_CA_microsteps_executed':micro_replay,'special_seeds_checked':seed_cases,'random_general_integer_seeds':80,'random_alphabet_triples':10000,'arbitrary_finite_support_configurations':500,'explicit_checks':checks,'zero_key_representation_caveat_present':zero_caveat,'scope':'All-positive-integer theorem follows from guard selection and unary operation proof; finite checks cover each rule, multi-control seeds, spectators and shifts. No metagraph conjugacy is asserted.'}
(ROOT/'independent_ca_audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
