#!/usr/bin/env python3
"""Independent physical controls and exact finite certificates. No assert statements."""
from itertools import combinations, permutations
from collections import deque
from pathlib import Path
import argparse,json
from interface import shapes,transitions,step,H
from compile_smooth import compile_graph,compile_word,factor_parts
ROOT=Path(__file__).resolve().parent

def require(v,msg):
 if not v:raise RuntimeError(msg)

def tournament(rows):
 n=len(rows)
 require(all(0<=r<1<<n for r in rows),'row range')
 for i in range(n):
  require(not(rows[i]>>i&1),'loop')
  for j in range(i):require(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1,'tournament edge')

def birth(rows,ports,S,drop):
 # Independent graph construction, deliberately does not call interface.physical.
 n=len(rows);out=rows+[0];targets={ports[i] for i in range(2) if S>>i&1}
 for i in range(n):
  if i in targets:out[n]+=1<<i
  else:out[i]+=1<<n
 q=ports+[n];del q[drop]
 return out,q

def realize(word):
 rows=[2,0];ports=[0,1]
 for S,d in word:rows,ports=birth(rows,ports,S,d)
 return rows,ports

def count_h(rows):
 n=len(rows);N=1<<n;f=[[0]*n for _ in range(N)]
 for v in range(n):f[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   if not(mask>>v&1):continue
   rem=mask^(1<<v)
   if rem:f[mask][v]=sum(f[rem][u] for u in range(n) if rem>>u&1 and rows[u]>>v&1)
 return sum(f[-1])

def perm_h(rows):
 return sum(all(rows[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(rows))))

def cover_profile(rows,ports):
 # Literal vertex permutations and all cuts. Start-label order gives each
 # unordered path cover exactly one representative.
 n=len(rows);ix={v:i for i,v in enumerate(ports)};catalog=shapes(2);pos={s:i for i,s in enumerate(catalog)};out=[0]*27
 for p in permutations(range(n)):
  bad=sum((not(rows[p[i]]>>p[i+1]&1))<<i for i in range(n-1))
  for cuts in range(1<<(n-1)):
   if bad&~cuts:continue
   comp=[];last=0
   for i in range(n):
    if i==n-1 or cuts>>i&1:
     comp.append((ix.get(p[last],2),ix.get(p[i],2)));last=i+1
   starts=[s for s,e in comp]
   if starts.count(2)>1 or starts!=sorted(starts):continue
   sig=tuple(comp);require(sig in pos,'physical signature missing')
   out[pos[sig]]+=1
 return tuple(out)

def det_bareiss(a):
 a=[list(map(int,r)) for r in a];n=len(a);prev=1;sign=1
 for k in range(n-1):
  if a[k][k]==0:
   p=next((i for i in range(k+1,n) if a[i][k]),None)
   if p is None:return 0
   a[k],a[p]=a[p],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    z=a[i][j]*pivot-a[i][k]*a[k][j]
    require(z%prev==0,'Bareiss exact division')
    a[i][j]=z//prev
  for i in range(k+1,n):a[i][k]=0
  prev=pivot
 return sign*a[-1][-1]

def strong(rows):
 n=len(rows)
 for r in range(n):
  seen=1<<r;old=-1
  while old!=seen:
   old=seen
   for v in range(n):
    if seen>>v&1:seen|=rows[v]
  if seen!=(1<<n)-1:return False
 return True

def saturation(bound):
 cap=bound+1;seed=cover_profile([2,0],[0,1]);seen={seed:0};q=deque([seed]);support={1}
 while q:
  st=q.popleft();dep=seen[st]
  for S in range(4):
   for d in range(3):
    ns=tuple(min(cap,v) for v in step(st,2,S,d));h=H(ns,2)
    if h>bound:continue
    support.add(h)
    if ns not in seen:seen[ns]=dep+1;q.append(ns)
 return {'bound':bound,'states':len(seen),'max_shortest_word':max(seen.values()),'support':sorted(support),'missing_odds':sorted(set(range(1,bound+1,2))-support)}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--bound',type=int,default=21);parser.add_argument('--quick',action='store_true');args=parser.parse_args()
 require(args.bound>=1,'positive bound')
 require(len(shapes(2))==27,'shape count')
 data=json.loads((ROOT/'automaton.json').read_text());require(data['shapes']==[[list(x) for x in s] for s in shapes(2)],'catalog mismatch')
 for g in data['generators']:
  require(g['transitions']==[[list(x) for x in r] for r in transitions(2,g['S'],g['drop'])],'matrix mismatch')
 tested=0;profiles=0
 if not args.quick:
  # Every labelled tournament of orders 2,3,4, every ordered pair of ports,
  # and every one of the twelve physical operations. Outputs have order <=5.
  for n in range(2,5):
   edges=list(combinations(range(n),2))
   for bits in range(1<<len(edges)):
    rows=[0]*n
    for j,(a,b) in enumerate(edges):
     if bits>>j&1:rows[a]|=1<<b
     else:rows[b]|=1<<a
    for ports in permutations(range(n),2):
     ports=list(ports);st=cover_profile(rows,ports);profiles+=1
     require(H(st,2)==count_h(rows),'input path count')
     for S in range(4):
      for d in range(3):
       out,p=birth(rows,ports,S,d);pred=step(st,2,S,d)
       require(pred==cover_profile(out,p),'physical transition')
       require(H(pred,2)==count_h(out),'output path count');tested+=1
 wits=json.loads((ROOT/'port_witnesses.json').read_text())
 for key,w in wits.items():
  rows,ports=realize(w['word']);tournament(rows)
  require(rows==w['rows'] and ports==w['ports'],'witness graph')
  require(count_h(rows)==int(key)==perm_h(rows),'witness independent counts')
  st=tuple(data['seed'])
  for S,d in w['word']:st=step(st,2,S,d)
  require(st==cover_profile(rows,ports),'witness complete interface')
 interaction=json.loads((ROOT/'interaction_witness.json').read_text())
 for key,target in [('first',41),('second',43)]:
  rows,h,ports,out=interaction[key];st=cover_profile(rows,ports)
  require(count_h(out)==target,'interaction scalar')
  h0=H(st,2)
  # Extract each individual old four-coordinate profile from cover states.
  for root in (0,1):
   e=sum(c for c,s in zip(st,shapes(2)) if len(s)==1 and s[0][1]==root)
   R=sum(c for c,s in zip(st,shapes(2)) if (len(s)==1 and s[0][0]==root) or (len(s)==2 and any(a==root for a,b in s)))
   d=sum(c for c,s in zip(st,shapes(2)) if len(s)==2 and (root,root) in s)
   require((h0,e,R,d)==(15,3,21,5),'individual profile hostile')
  cross=sum(c for c,s in zip(st,shapes(2)) if len(s)==2 for a,b in [(0,1),(1,0)] if any(e==a for v,e in s) and any(v==b for v,e in s) and not any(v==b and e==a for v,e in s))
  require(cross==51-target,'cross interaction')
 # Exact 17-dimensional upper bound and independent physical lower bound.
 cert=json.loads((ROOT/'rank_certificate.json').read_text());O=cert['observable_basis'];require(len(O)==17,'basis size')
 hout=[int(len(s)==1) for s in shapes(2)]
 for row,word in zip(O,cert['suffixes']):
  obs=hout[:]
  for S,d in reversed(word):obs=[sum(obs[j]*m for j,m in t) for t in transitions(2,S,d)]
  require(obs==row,'observable basis words')
 for gi,(S,d) in enumerate((S,d) for S in range(4) for d in range(3)):
  A=[[int(z) for z in r] for r in cert['reduced_transitions'][gi]]
  for i in range(17):
   left=[sum(O[i][j]*m for j,m in t) for t in transitions(2,S,d)]
   right=[sum(A[i][j]*O[j][c] for j in range(17)) for c in range(27)]
   require(left==right,'observable span closure')
 C=[]
 for suf in cert['suffixes']:
  C.append([count_h(realize(pre+suf)[0]) for pre in cert['prefixes']])
 require(C==cert['continuation_matrix'],'physical continuation matrix')
 determinant=det_bareiss(C);require(determinant==-2097152==cert['determinant'],'rank determinant')
 # Two explicit all-order 1/2/3-exponential families, controlled at m<=5.
 rays={67:(68,26,-29),173:(228,-36,-19)};ray_checks=0
 for h,(a,b,c) in rays.items():
  w=wits[str(h)];rows=w['rows'];ports=w['ports'];require(strong(rows),'strong seed')
  st=cover_profile(rows,ports)
  for m in range(1,6):
   rows,ports=birth(rows,ports,3,2);st=step(st,2,3,2);expected=a*3**m+b*2**m+c
   require(count_h(rows)==H(st,2)==expected,'mixed clock formula control')
   require(strong(rows),'strong family');ray_checks+=1
 # Product closure, all finite seed graphs, and the unbounded target compiler.
 bank=json.loads((ROOT/'smooth_seeds.json').read_text())
 for key,w in bank.items():
  rows,ports=realize(w['word']);require(rows==w['rows'] and ports==w['ports'],'smooth seed graph')
  require(count_h(rows)==int(key),'smooth seed Hamilton count')
 for a,b in [(3,13),(49,5),(13,343),(343,161)]:
  wa=bank[str(a)]['word'];wb=bank[str(b)]['word'];A,_=realize(wa);B,_=realize(wb)
  out,_=realize(wa+[[0,0],[0,0]]+wb);na=len(A)
  require(len(out)==na+len(B),'product order')
  for i in range(na):require(out[i]==A[i]+(((1<<len(B))-1)<<na),'product old block')
  for i in range(len(B)):require(out[na+i]==B[i]<<na,'product new block')
 for N in (1,39,245,343,441):
  graph=compile_graph(N);require(count_h(graph['rows'])==N,'smooth compiled count')
 for N in (7**5*3**4*13*23,49**5*19**3,23**20,3**50):
  parts=factor_parts(N);prod=1
  for p in parts:prod*=p
  require(prod==N and all(str(p) in bank for p in parts),'large smooth factorization')
 for N in (7,21,2,29,0):
  rejected=False
  try:factor_parts(N)
  except ValueError:rejected=True
  require(rejected,'compiler domain guard')
 rec=saturation(args.bound)
 known={21:(4697,12),67:(58098,16),173:(485591,19)}
 if args.bound in known:
  require((rec['states'],rec['max_shortest_word'])==known[args.bound],'census count')
  require(rec['missing_odds']==[7,21],'bounded support')
 print(json.dumps({'physical_input_profiles':profiles,'physical_birth_checks':tested,'selected_physical_witnesses':len(wits),'rank_determinant':determinant,'observable_dimension':17,'strong_mixed_clock_checks':ray_checks,'smooth_seed_checks':len(bank),'ordered_product_controls':4,'census':rec},sort_keys=True,indent=2))
if __name__=='__main__':main()
