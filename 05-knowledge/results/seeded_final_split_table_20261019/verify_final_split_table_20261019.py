#!/usr/bin/env python3
"""Replay every final split-lift state, using only the Python standard library."""
import hashlib,json
from math import isqrt
from pathlib import Path
ROOT=Path(__file__).resolve().parent
checks=0

def check(ok,message):
 global checks
 checks+=1
 if not ok:raise RuntimeError(message)

path=ROOT/'final_split_glue_certificates_20261019.json'
raw=path.read_bytes();table=json.loads(raw)
expected={(p,d,a,b)for p in(0,1)for d in range(24,73)for a in(-1,1)for b in(-1,1)}
seen_states=set()
for data in table:
 p,d,o0,o1=(data[k]for k in('parity','residue','order0','order1'))
 state=(p,d,o0,o1)
 check(state in expected and state not in seen_states,'state range or duplication')
 seen_states.add(state)
 check(data['cutparity']==1,'2 must have odd input position')
 choices={c:int(49+(-1)**(p+1)*c<=d)for c in range(-24,25)}
 j=int(d>=57);c10=-8*(-1)**(p+j)
 check(choices[c10]==j,'the 482 port must select the input ending in 10')
 endpoints={};labels=[]
 for c in range(-24,25):
  i=choices[c];length_parity=(p+i)%2;end=10 if i==j else 8
  first,last=(98-c,343+c)if(o0 if i==0 else o1)==1 else(343+c,98-c)
  endpoints[(c,0)]=(3381-c,first)
  endpoints[(c,1)]=(last,49*end+(-1)**length_parity*c)
  labels.extend(endpoints[(c,0)]+endpoints[(c,1)])
 check(len(set(labels+list(range(1,25))))==220,'all ports and singleton vertices must be distinct')
 seen_segments=set();seen_singles=set();glues=[];previous=None
 for atom in data['word']:
  check(atom[0]in('single','segment'),'unknown atom type')
  if atom[0]=='single':
   check(len(atom)==2,'singleton schema')
   v=atom[1];check(1<=v<=24 and v not in seen_singles,'singleton coverage')
   seen_singles.add(v);a=b=v
  else:
   check(len(atom)==5,'segment schema')
   _,c,i,s,orientation=atom
   check((c,s)in endpoints and(c,s)not in seen_segments,'segment coverage')
   seen_segments.add((c,s))
   check(i==choices[c] and orientation in(-1,1),'forced input choice or orientation')
   a,b=endpoints[(c,s)]
   if orientation==-1:a,b=b,a
  if previous is not None:
   check(isqrt(previous+a)**2==previous+a,'nonsquare external connector')
   glues.append((previous,a))
  previous=b
 check(seen_segments==set(endpoints),'missing segment')
 check(seen_singles==set(range(1,25)),'missing singleton')
 check(data['word'][0]==['single',1]and data['word'][-1]==['single',8],'closing endpoints')
 check(isqrt(1+8)**2==1+8,'closing square')
 check((482,3362)in glues or(3362,482)in glues,'missing seed connector')
 check(endpoints[(19,0)][0]==3362,'wrong protected-prefix segment')
 check(49*127+19==6242 and 3362+6242==98**2,'protected internal seed edge')
 check({69,127}.isdisjoint({2,7}),'split cannot cut protected prefix')
check(seen_states==expected,'incomplete state table')
print(json.dumps({'status':'passed','states':len(table),'checks':checks,
 'certificate_bytes':len(raw),'certificate_sha256':hashlib.sha256(raw).hexdigest(),
 'scope':'All 392 conditional final-lift states, independent of input length and internal order.',
 'conclusion':'Every admissible input pair produces a square-sum Hamiltonian cycle containing 482-3362-6242.'},sort_keys=True,indent=2))
