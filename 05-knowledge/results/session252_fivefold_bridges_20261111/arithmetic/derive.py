#!/usr/bin/env python3
"""Small controls for new Session252 arithmetic constructions; no repo writes."""
import json, sys, importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
spec=importlib.util.spec_from_file_location('cpcompiler', ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)

def step(v,b,m=None):
 g,n,e=v
 w=(g+n+e,2*n,0) if b==0 else (g,2*n+1,1)
 return tuple(x%m for x in w) if m else w

def word_state(word,m=None):
 v=(1,0,0)
 for b in word:v=step(v,int(b),m)
 return v

def order2(m):
 n=1;x=2%m
 while x!=1:x=2*x%m;n+=1
 return n

# Four-coordinate transfer and exact block formulas, not the inherited graph census.
for n in range(256):
 if not (word_state(bin(n)[2:])==(c.g(n),n,n%2)):raise AssertionError("arithmetic control failed")
 for L in range(1,17):
  z=word_state(bin(n)[2:]+'0'*L)
  o=word_state(bin(n)[2:]+'1'*L)
  if not (z==(c.g(n)+(2**L-1)*n+n%2,2**L*n,0)):raise AssertionError("arithmetic control failed")
  if not (o==(c.g(n),2**L*n+2**L-1,1)):raise AssertionError("arithmetic control failed")

# CRT/ray controls using modular state iteration avoid giant target integers.
fermat=[2**(2**i)+1 for i in range(5)]
M=2**32-1
ray=[]
for s in range(1,8):
 m=2**(2**s)-1;L=2**s;delta=4*L*pow(7,-1,m)%m
 for k in [0,1,2,3,7,19]:
  v=word_state('100'*(L*k),m)
  if not (v==( (1+delta*k)%m,0,0)):raise AssertionError("arithmetic control failed")
 ray.append({'s':s,'modulus':str(m),'block_repetitions':L,'drift':str(delta)})
for t in range(35):
 v=word_state('100'*t)
 if not (7*v[1]==4*(8**t-1)):raise AssertionError("arithmetic control failed")
 if not (49*v[0]==25+24*8**t+28*t):raise AssertionError("arithmetic control failed")

# New direct bridge to historical Lucas/Gilbreath addresses.
def W(t):
 out=1
 for i in range(t.bit_length()):
  if t>>i&1:out*=2**(2**i)+1
 return out
lucas=[]
for t in range(32):
 n=W(t)
 if not (n==sum(2**j for j in range(t+1) if j&t==j)):raise AssertionError("arithmetic control failed")
 # zero insertion, bit duplication (leading zero is harmless)
 if not (int(''.join('0'+b for b in bin(n)[2:]),2)==W(2*t)):raise AssertionError("arithmetic control failed")
 if not (int(''.join(b+b for b in bin(n)[2:]),2)==W(2*t+1)):raise AssertionError("arithmetic control failed")
 lucas.append({'t':t,'address':n,'g':c.g(n)})

# Compact real tournament certificates at selected Lucas addresses.
cores=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]
witnesses=[]
for j,core in enumerate(cores):
 O=1;s=c.signature(core,O)
 for t in [0,1,2,3]:
  A=0;B=W(t);H=c.count(s,A,B);rows=c.construct(core,O,A,B)
  if not (c.hp_subsets(rows)[-1]==H):raise AssertionError("arithmetic control failed")
  witnesses.append({'core':j,'O':O,'A':A,'B':B,'H':H,'rows':rows})

ogg=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]
report={'status':'FINITE-EXACT controls of newly derived algebra, not all-height coverage',
 'fermat_periods':dict(zip(fermat,[order2(p) for p in fermat])),
 'ogg_odd_periods':{p:order2(p) for p in ogg if p!=2},
 'genus_zero_pentad':[2,3,5,7,13],
 'genus_zero_odd_periods':{p:order2(p) for p in [3,5,7,13]},
 'crt_modulus':M,'crt_idempotents':{p:(M//p)*pow(M//p,-1,p)%M for p in fermat},
 'modular_surjectivity_rays':ray,'lucas_addresses':lucas,'graph_witnesses':witnesses}
Path(__file__).with_name('results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['graph_witnesses','lucas_addresses']},indent=2))
