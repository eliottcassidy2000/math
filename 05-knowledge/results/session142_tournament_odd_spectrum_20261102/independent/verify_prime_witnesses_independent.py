#!/usr/bin/env python3
"""Independent adjacency reconstruction and native-coordinate certificate.
No producer imports. Full order29 subset DP is deliberately not claimed.
"""
from pathlib import Path
from itertools import permutations,combinations
from math import isqrt
import json
BASE=Path(__file__).resolve().parent.parent/'prime80557'

def need(c,m):
 if not c:raise RuntimeError(m)
def words(a,vs):return [p for p in permutations(vs) if all(a[u][v] for u,v in zip(p,p[1:]))]
def native(a,cut):
 n=len(a);vs=set(range(n));ps=words(a,range(n));h=len(ps)
 g=sum(p[-1] in cut for p in ps)
 q=0
 for mask in range(1,1<<n):
  qs={v for v in range(n) if mask>>v&1}
  q+=len(words(a,vs-qs))*sum(p[0] in cut for p in words(a,qs))
 full=[row+[int(i not in cut)] for i,row in enumerate(a)]+[[int(i in cut) for i in range(n)]+[0]]
 return (h,len(words(full,range(n+1))),q,g)
def step(state,bit):
 h,f,q,g=state
 return (h,f,2*q+h,h) if bit else (h,f+q+g,2*q,0)
def local(tag):
 if tag=='S':return [[0]]
 need(tag=='D','unknown block')
 rs=[8,1,3,6]
 return [[(r>>j)&1 for j in range(4)] for r in rs]
def reconstruct(recipe):
 blocks=[];owner=[];cuts=[]
 for tag,mask in recipe:
  a=local(tag);need(type(mask)==int and 0<=mask<1<<len(a),'local cut')
  block=list(range(len(owner),len(owner)+len(a)));blocks.append(block)
  owner.extend((len(blocks)-1,j) for j in range(len(a)))
  cuts.extend((mask>>j)&1 for j in range(len(a)))
 n=len(owner);out=[[0]*(n+1) for _ in range(n+1)]
 for i,(bi,li) in enumerate(owner):
  for j,(bj,lj) in enumerate(owner):
   out[i][j]=local(recipe[bi][0])[li][lj] if bi==bj else int(bi<bj)
  out[i][n]=1-cuts[i];out[n][i]=cuts[i]
 return out,blocks

def valid(a):return all(not a[i][i] and all(a[i][j] in (0,1) and a[i][j]+a[j][i]==1 for j in range(i)) for i in range(len(a)))
def strong(a):
 for v in range(len(a)):
  seen={v};front={v}
  while front:
   front={j for i in front for j in range(len(a)) if a[i][j]}-seen;seen|=front
  if len(seen)<len(a):return False
 return True

def one_apex_partition_count(a):
 m=len(a)-1
 need(m>=1 and all(a[i][j]==int(i<j) for i in range(m) for j in range(m)), 'parent not transitively ordered')
 bits=a[m][:m]
 empty_before=bits[0];empty_after=1-bits[-1]
 forward=sum((1-bits[i])*bits[i+1] for i in range(m-1))
 backward=sum((1-bits[i])*bits[j]*2**(i-j-1) for i in range(m) for j in range(i))
 return empty_before+empty_after+forward+backward

def controls():
 tests=0
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   a=[[0]*n for _ in range(n)]
   for k,(i,j) in enumerate(combinations(range(n),2)):a[i][j]=(code>>k)&1;a[j][i]=1-a[i][j]
   for mask in range(1<<n):
    cut={i for i in range(n) if mask>>i&1};state=native(a,cut)
    for bit in (0,1):
     b=[row+[1] for row in a]+[[0]*(n+1)];newcut=cut|({n} if bit else set())
     need(native(b,newcut)==step(state,bit),'literal S0/S1 native transition');tests+=1
 return tests

def main():
 out=[]
 for target in (80557,100000007):
  cert=json.loads((BASE/f'witness_{target}.json').read_text());recipe=cert['word'];a,blocks=reconstruct(recipe)
  rows=[sum(b<<j for j,b in enumerate(row)) for row in a]
  need(valid(a) and strong(a),'actual tournament legality or strongness')
  need(rows==cert['rows'] and blocks==cert['blocks'] and len(a)==cert['order'],'certificate graph reconstruction')
  tag,mask=recipe[0];seed=local(tag);state=native(seed,{i for i in range(len(seed)) if mask>>i&1});ledger=[list(state)]
  for tag,bit in recipe[1:]:
   need(tag=='S' and bit in (0,1),'certificate outside audited singleton-tail grammar')
   state=step(state,bit);ledger.append(list(state))
  need(state[1]==target and cert['target']==target,'compositional target mismatch')
  need(all(target%d for d in range(2,isqrt(target)+1)),'target not prime')
  text=(BASE/f'adjacency_{target}.txt').read_text().splitlines()
  need(text[0]==f'{len(a)} {target}' and text[1:]==[''.join(map(str,row)) for row in a],'adjacency text differs')
  direct=None
  if all(tag=='S' for tag,_ in recipe):
   direct=one_apex_partition_count(a)
   need(direct==target,'direct max-before/min-after graph count')
  # Reverse one actual apex edge, keeping a valid tournament, and recount.
  altered=[list(x) for x in recipe];altered[0][1]^=1
  other,_=reconstruct(altered);need(valid(other),'hostile is not a tournament')
  tag,mask=altered[0];seed=local(tag)
  other_state=native(seed,{i for i in range(len(seed)) if mask>>i&1})
  for _,bit in altered[1:]:other_state=step(other_state,bit)
  need(other_state[1]!=target,'changed apex edge escaped target-count check')
  out.append({'target':target,'order':len(a),'strong':True,'prime':True,'literal_seed_native_state':ledger[0],
    'native_state_order':['h','F','q','g'],'native_state_ledger':ledger,'verified_H':state[1],'direct_one_apex_graph_count':direct,'hostile_reversed_apex_edge_H':other_state[1],
    'large_graph_count_method':'exact native-coordinate composition from a literal seed; not full subset DP'})
 print(json.dumps({'status':'PASS','literal_small_transition_controls':controls(),'witnesses':out},indent=2,sort_keys=True))
if __name__=='__main__':main()
