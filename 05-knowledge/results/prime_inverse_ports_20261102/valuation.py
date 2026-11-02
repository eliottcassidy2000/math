#!/usr/bin/env python3
"""Literal block profiles, weighted valuation gates and two-event residues."""
from functools import lru_cache
from itertools import permutations,product
from pathlib import Path
import json,random
from research import all_h,circ,need,strong
HERE=Path(__file__).resolve().parent
@lru_cache(None)
def data(rows):
 n=len(rows);H=0;starts=[0]*n;ends=[0]*n;W=[[0]*n for _ in rows]
 def path(w):return all(rows[u]>>v&1 for u,v in zip(w,w[1:]))
 for w in permutations(range(n)):
  if path(w):H+=1;starts[w[0]]+=1;ends[w[-1]]+=1
  for k in range(1,n):
   if path(w[:k]) and path(w[k:]):W[w[k-1]][w[k]]+=1
 return H,starts,ends,W
@lru_cache(None)
def profile(rows,mask):
 H,starts,ends,W=data(rows);n=len(rows)
 E=2*H+sum(map(sum,W));s=sum(starts[v] for v in range(n) if mask>>v&1)
 t=sum(ends[v] for v in range(n) if not(mask>>v&1))
 Q=sum(starts[v]+sum(W[u][v] for u in range(n)) for v in range(n) if mask>>v&1)
 P=sum(ends[u]+sum(W[u][v] for v in range(n)) for u in range(n) if not(mask>>u&1))
 graph=list(rows)+[mask]
 for v in range(n):
  if not(mask>>v&1):graph[v]|=1<<n
 f=all_h(graph)[-1]
 return H,E,P,Q,s,t,f

def step(z,b):
 F,q,g,h=z;H,E,P,Q,s,t,f=b
 return H*F+P*q+(H-s)*g+(f-H)*h,E*q+Q*h,(H-t)*h,H*h

def append(rows,B,mask):
 # mark remains vertex0; old parent precedes the new block
 n=len(rows);m=len(B);out=list(rows)+[x<<n for x in B]
 for u in range(1,n):out[u]|=((1<<m)-1)<<n
 for v in range(m):
  if mask>>v&1:out[0]|=1<<(n+v)
  else:out[n+v]|=1
 return out

def vp(n,p):
 if n==0:return 10**6
 k=0
 while n%p==0:k+=1;n//=p
 return k

def weighted_bound(bs,p):
 weights=sorted(min(vp(b[0],p),vp(b[1],p)) for b in bs)
 return sum(weights[:-2])

def audit():
 triangle=(2,4,1);diamond=(8,1,3,6)
 c5=tuple(circ(5,(0,0)));c7=tuple(circ(7,(0,0,0)));p7=tuple(circ(7,(0,0,1)))
 # Select the other order-seven orbit profile, independently of a set-name assumption.
 for bits in product([0,1],repeat=3):
  candidate=tuple(circ(7,bits))
  if data(candidate)[0]==189:p7=candidate;break
 blocks=[(0,),triangle,diamond,c5,c7,p7]
 bases=[(B,mask,profile(B,mask)) for B in blocks for mask in range(1<<len(B))]
 primes=[2,3,5,7,11];rng=random.Random(20261102);checks=0;direct=0;records=[]
 for _ in range(2500):
  word=[rng.choice(bases) for _ in range(rng.randrange(1,25))];z=(1,0,0,1)
  for _,_,b in word:z=step(z,b)
  for p in primes:
   need(vp(z[0],p)>=weighted_bound([b for _,_,b in word],p),'weighted valuation');checks+=1
 for _ in range(120):
  word=[];size=1
  while size<11:
   B,mask,b=rng.choice(bases)
   if size+len(B)>12:break
   word.append((B,mask,b));size+=len(B)
  z=(1,0,0,1);graph=[0]
  for B,mask,b in word:z=step(z,b);graph=append(graph,B,mask)
  need(all_h(graph)[-1]==z[0],'literal adjacency block action');direct+=1
 # Exact two-event formula; all cuts for p=3,5; all cardinality choices and
 # 64 actual cut pairs for each ordered choice of the two p=7 cyclic types.
 pair_controls=0;residue_examples=[]
 for p,family in [(3,[triangle]),(5,[c5]),(7,[c7,p7])]:
  for A in family:
   for B in family:
    HA,EA,*_=profile(A,0);HB,EB,*_=profile(B,0)
    u1,d1=HA//p,(EA-HA)//p;u2,d2=HB//p,(EB-HB)//p
    pairs=list(product(range(1<<p),repeat=2)) if p<=5 else [(rng.randrange(1<<p),rng.randrange(1<<p)) for _ in range(64)]
    for a,b in pairs:
     Aprof=profile(A,a);Bprof=profile(B,b)
     for middle in [[],[(0,),0],[(0,),1],[diamond,1]]:
      z=step((1,0,0,1),Aprof);Ew=1
      if middle:z=step(z,profile(tuple(middle[0]),middle[1]));Ew=profile(tuple(middle[0]),middle[1])[1]
      z=step(z,Bprof)
      predicted=-a.bit_count()*b.bit_count()*(d1*d2*Ew if middle else d1*d2+u1*u2)
      need((z[0]-predicted)%p==0,'two-event residue');pair_controls+=1
    residue_examples.append({'p':p,'left_H':HA,'right_H':HB,'u_pair':[u1,u2],'d_pair':[d1,d2],
       'adjacent_gram_mod_p':(d1*d2+u1*u2)%p,'separated_q_product_mod_p':d1*d2%p})
 # Sharpness: C1 C1 C0^r, r>=1.
 z=step(step((1,0,0,1),profile(triangle,1)),profile(triangle,1));need(z==(41,45,3,9),'sharp seed')
 graph=append(append([0],triangle,1),triangle,1)
 for r in range(1,51):
  z=step(z,profile(triangle,0));expected=3**r*(45*4**r-1)
  need(z[0]==expected and vp(z[0],3)==r,'sharp family')
  if r<=2:
   graph=append(graph,triangle,0);need(all_h(graph)[-1]==expected and strong(graph),'sharp actual graph');direct+=1
  records.append({'r':r,'F':z[0],'v3':r,'vertices':7+3*r})
 # One-zero terminal lower bound from insertion runs, for actual small sources.
 terminal_checks=0
 for B in blocks[:4]:
  H,_,_,_=data(B)
  for mask in range(1,1<<len(B)):
   z=step((1,0,0,1),profile(B,mask));terminal=step(z,profile((0,),0))[0]
   need(terminal>=(mask.bit_count()+2)*H,'terminal h budget');terminal_checks+=1
 # A genuinely heavier strong SCC: the circulant on9 with steps8,2,6,5.
 R=circ(9,(1,0,1,1));hs=all_h(R);H=hs[-1];E=sum(hs[S]*hs[511^S] for S in range(512))
 need((H,E)==(3159,48924) and strong(R),'heavy singular block')
 heavy=[];z=(1,0,0,1)
 for cut in [1,3,7]:
  a=cut.bit_count();u=H//9;d=(E-H)//9
  child=list(R)+[cut]
  for v in range(9):
   if not(cut>>v&1):child[v]|=1<<9
  b=(H,E,(9-a)*d,a*d,a*u,(9-a)*u,all_h(child)[-1]);heavy.append(b);z=step(z,b)
 need(weighted_bound(heavy,3)==4 and vp(z[0],3)>=4,'weighted strict improvement')
 heavy_record={'rows':R,'H':H,'E':E,'weight_at3':4,'cuts':[1,3,7],
               'output':z[0],'output_v3':vp(z[0],3),'weighted_bound':4,'unweighted_bound':1,'vertices':28}
 out={'status':'PASS','heavy_block_control':heavy_record,'literal_profile_bank':len(bases),'weighted_word_prime_checks':checks,'full_adjacency_checks':direct,
      'two_event_residue_checks':pair_controls,'terminal_budget_controls':terminal_checks,
      'two_event_types':residue_examples,'sharp_family':records,
      'scope':'Weighted valuation and two-event formulas are proved for arbitrary words; finite controls do not prove scalar entry.'}
 (HERE/'valuation.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('sharp_family','two_event_types')},indent=2))
if __name__=='__main__':audit()
