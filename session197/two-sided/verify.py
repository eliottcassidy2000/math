#!/usr/bin/env python3
from itertools import combinations,permutations
from collections import deque
from math import factorial
from pathlib import Path
import argparse,json
from core import require,wrap,operator,count_h,cover_polynomial,det,transitive,strong,feedback_number
from rightward import shapes,step,H,transitions
HERE=Path(__file__).resolve().parent

def independent_profile(rows,ports):
 # Canonicalize by actual first vertices, unlike rightward.brute's port labels.
 n=len(rows);label={v:i for i,v in enumerate(ports)};pos={s:i for i,s in enumerate(shapes(2))};counts=[0]*27
 for word in permutations(range(n)):
  bad=sum((not(rows[word[i]]>>word[i+1]&1))<<i for i in range(n-1))
  for cut in range(1<<(n-1)):
   if bad&~cut:continue
   paths=[];start=0
   for i in range(n):
    if i==n-1 or cut>>i&1:paths.append((word[start],word[i]));start=i+1
   actual_starts=[a for a,b in paths]
   if actual_starts!=sorted(actual_starts):continue
   sig=tuple(sorted((label.get(a,2),label.get(b,2)) for a,b in paths))
   if sum(a==2 for a,b in sig)>1:continue
   counts[pos[sig]]+=1
 return tuple(counts)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--cpp');args=parser.parse_args()
 matrices=json.loads((HERE/'automaton.json').read_text())
 for g in matrices['generators']:
  require(g['transitions']==[[list(a) for a in r] for r in transitions(2,g['S'],g['drop'])],'rightward matrix')
 physical_cases=0
 for n in range(5):
  edges=list(combinations(range(n),2))
  for bits in range(1<<len(edges)):
   T=[0]*n
   for k,(a,b) in enumerate(edges):
    if bits>>k&1:T[a]|=1<<b
    else:T[b]|=1<<a
   require(operator(cover_polynomial(T))==cover_polynomial(wrap(T)),'full polynomial wrapper')
   if n:
    require(strong(wrap(T)),'strong wrapper')
    require(feedback_number(wrap(T))==feedback_number(T)+1,'feedback number increment')
   physical_cases+=1
 # Physical continuation matrix on Tr_2,...,Tr_6 and 0,...,4 wrappers.
 N=5;C=[]
 for n in range(2,N+2):
  T=transitive(n);row=[]
  for m in range(N):row.append(count_h(T));T=wrap(T)
  C.append(row)
 expected=factorial(N)
 for j in range(1,N+1):expected*=factorial(j)*factorial(j-1)
 require(det(C)==expected==1194393600,'infinite-rank finite control')
 # Matching graph on the central root is exactly m disjoint edges.
 T=[0]
 for m in range(1,8):
  T=wrap(T);O=[v for v in range(1,len(T)) if T[0]>>v&1];I=[v for v in range(1,len(T)) if T[v]>>0&1]
  cross={(a,b) for a in O for b in I if T[a]>>b&1}
  require(cross=={(2*j,2*j-1) for j in range(1,m+1)},'central exact matching')
 # New exact physical scalar families, not numerical interpolation.
 family_checks=0
 for n in range(1,8):
  T=transitive(n)
  T=wrap(wrap(T));require(count_h(T)==2*3**n+5*2**n-1,'double wrapper family');family_checks+=1
  T=wrap(T);require(count_h(T)==6*4**n+28*3**n+13*2**n-11,'triple wrapper family');family_checks+=1
 T=[122,112,67,102,108,4,32]
 left=[r|((1<<7) if i in (1,4) else 0) for i,r in enumerate(T)]+[109]
 require(left==[122,240,67,102,236,4,32,109] and count_h(left)==151,'leftward H151')
 witness149=json.loads((HERE/'escape149.json').read_text());rows=T[:];ports=witness149['initial_ports'][:]
 for side,S,d in witness149['word']:
  n=len(rows);chosen={ports[i] for i in range(2) if S>>i&1};out=rows+[0]
  for i in range(n):
   if (i in chosen) if side=='R' else (i not in chosen):out[n]|=1<<i
   else:out[i]|=1<<n
  ports.append(n);del ports[d];rows=out
 require(rows==[378,112,67,102,108,4,32,127,254]==witness149['rows'] and ports==witness149['ports'],'leftward H149 physical word')
 require(count_h(rows)==149,'leftward H149 count')
 for graph,target in [(rows,149),(left,151)]:
  direct=sum(all(graph[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(graph))))
  require(direct==target,'independent permutation child count')
  require(strong(graph),'strong child witness')
 extensions={149:[],151:[]}
 for mask in range(128):
  out=[r|((1<<7) if not(mask>>i&1) else 0) for i,r in enumerate(T)]+[mask];h=count_h(out)
  if h in extensions:extensions[h].append(mask)
 require(extensions=={149:[],151:[109]},'complete one-vertex orientation census')
 # Exact arbitrary-depth rightward exclusions for every initial port pair.
 certificate=json.loads((HERE/'rightward_certificates.json').read_text());cases=certificate['cases']
 require(certificate['rows']==T and certificate['bound']==151,'rightward certificate metadata')
 expected_pairs=list(combinations(range(7),2));require([tuple(c['ports']) for c in cases]==expected_pairs,'all 21 pair choices')
 cpp=json.loads(Path(args.cpp).read_text()) if args.cpp else None
 summary=[];union=set();total=0
 for case in cases:
  ports=case['ports'];actual=independent_profile(T,ports)
  require(list(actual)==case['uncapped_seed'],'independent physical profile')
  initial=tuple(min(152,x) for x in actual);require(list(initial)==case['capped_seed'],'capped profile')
  states={tuple(z) for z in case['states']};require(len(states)==len(case['states']) and initial in states,'unique state inventory')
  reached={initial};q=deque([initial]);support=set()
  while q:
   z=q.popleft();support.add(H(z,2));require(H(z,2)<=151,'retained bound')
   for S in range(4):
    for d in range(3):
     nz=tuple(min(152,x) for x in step(z,2,S,d))
     if H(nz,2)>151:continue
     require(nz in states,'exact closure')
     if nz not in reached:reached.add(nz);q.append(nz)
  require(reached==states and sorted(support)==case['support'],'exact reachability/support')
  require(149 not in support and 151 not in support,'rightward target exclusion')
  if cpp:
   cc=next(c for c in cpp if c['ports']==ports)
   require(cc['states']==case['states'] and cc['support']==case['support'],'full independent C++ state match')
  union|=support;total+=len(states);summary.append({'ports':ports,'states':len(states),'support':sorted(support)})
 require(union=={49,85,91,115,137,147},'union support')
 print(json.dumps({'full_polynomial_physical_controls':physical_cases,'rank_matrix':C,'rank_determinant':expected,'all_order_family_controls':family_checks,'H151_leftward_witness':left,'H149_leftward_witness':rows,'feedback_increment_controls':physical_cases-1,'one_vertex_extensions':extensions,'rightward_pair_cases':summary,'total_rightward_states':total,'independent_cpp_full_state_match':bool(cpp)},sort_keys=True,indent=2))
if __name__=='__main__':main()
