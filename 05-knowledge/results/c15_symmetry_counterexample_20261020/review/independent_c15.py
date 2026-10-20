#!/usr/bin/env python3
"""Independent C15 n=10 verification; imports no producer code.
All 512 fixed graphs/tournaments, exact twin quotient, unique strong modules.
"""
import itertools,json,math
from pathlib import Path
N=10
PAIRS=list(itertools.combinations(range(N),2)); INDEX={e:k for k,e in enumerate(PAIRS)}
BLOCKS=[tuple(range(3)),tuple(range(3,8)),(8,),(9,)]
QPAIRS=list(itertools.combinations(range(4),2))
IDENTITY=tuple(range(N)); HGEN=(1,2,0,4,5,6,7,3,8,9)
def need(x,msg):
 if not x: raise AssertionError(msg)
def graph_from(q,ia=0,ib=0):
 E=set()
 if ia:E.update(itertools.combinations(BLOCKS[0],2))
 for i,j in itertools.combinations(range(5),2):
  d=(j-i)%5; orb=0 if d in (1,4) else 1
  if ib>>orb&1:E.add((i+3,j+3))
 for bit,(a,b) in enumerate(QPAIRS):
  if q>>bit&1:E.update(tuple(sorted((x,y))) for x in BLOCKS[a] for y in BLOCKS[b])
 return frozenset(E)
def tournament(q,ia=0,ib=0):
 arcs=set()
 for i,j in itertools.combinations(range(3),2):
  forward=(j-i)%3==1
  if ia:forward=not forward
  arcs.add((i,j) if forward else (j,i))
 for i,j in itertools.combinations(range(5),2):
  d=(j-i)%5; orb=0 if d in (1,4) else 1; forward=d in (1,2)
  if ib>>orb&1:forward=not forward
  arcs.add((i+3,j+3) if forward else (j+3,i+3))
 for bit,(a,b) in enumerate(QPAIRS):
  for x in BLOCKS[a]:
   for y in BLOCKS[b]:arcs.add((y,x) if q>>bit&1 else (x,y))
 return frozenset(arcs)
def gaction(E,p):return frozenset(tuple(sorted((p[x],p[y]))) for x,y in E)
def taction(A,p):return frozenset((p[x],p[y]) for x,y in A)
def sign(E,p):return (-1)**sum(p[x]>p[y] for x,y in E)
def mask(E):return sum(1<<INDEX[e] for e in E)
def quotient(E):
 neigh=[frozenset(y for y in range(N) if tuple(sorted((x,y))) in E) for x in range(N)]
 cells={}
 for x,v in enumerate(neigh):cells.setdefault(v,[]).append(x)
 cells=list(cells.values());k=len(cells);sizes=[len(c) for c in cells]
 adj=[[int(tuple(sorted((cells[i][0],cells[j][0]))) in E) if i!=j else 0 for j in range(k)] for i in range(k)]
 return cells,sizes,adj
def graph_type_and_sign(E):
 cells,sizes,A=quotient(E);k=len(cells)
 can=min(tuple(sizes[p[i]] for i in range(k))+tuple(A[p[i]][p[j]] for i,j in itertools.combinations(range(k),2)) for p in itertools.permutations(range(k)))
 auts=[];odd_witness=None
 for p in itertools.permutations(range(k)):
  if any(sizes[i]!=sizes[p[i]] for i in range(k)):continue
  if any(A[i][j]!=A[p[i]][p[j]] for i in range(k) for j in range(k)):continue
  lift=list(range(N))
  for i in range(k):
   for x,y in zip(cells[i],cells[p[i]]):lift[x]=y
  need(gaction(E,lift)==E,'lift is graph automorphism')
  sg=sign(E,lift);auts.append(sg)
  if sg<0:odd_witness=lift
 return can,all(s>0 for s in auts),odd_witness,math.prod(math.factorial(s) for s in sizes)*len(auts),sizes
def module(A,S):
 S=set(S)
 return all(len({(v,s) in A for s in S})==1 for v in range(N) if v not in S)
def strong(A,S):
 S=set(S)
 for v in S:
  seen={v};todo=[v]
  while todo:
   x=todo.pop()
   for y in S-seen:
    if (x,y) in A:seen.add(y);todo.append(y)
  if seen!=S:return False
 return True
def prime(A,S):
 S=set(S)
 for k in range(2,len(S)):
  for W in itertools.combinations(sorted(S),k):
   if all(len({(v,w) in A for w in W})==1 for v in S-set(W)):return False
 return True
def swap_q_t(q):
 p=(0,1,3,2);arcs={(b,a) if q>>k&1 else (a,b) for k,(a,b) in enumerate(QPAIRS)}
 A={(p[x],p[y]) for x,y in arcs}
 return sum(1<<k for k,(a,b) in enumerate(QPAIRS) if (b,a) in A)
def swap_q_g(q):
 p=(0,1,3,2);E={tuple(sorted((p[a],p[b]))) for k,(a,b) in enumerate(QPAIRS) if q>>k&1}
 return sum(1<<k for k,e in enumerate(QPAIRS) if e in E)
def main():
 allg={};allt=set();excluded=0
 for ia,ib,q in itertools.product(range(2),range(4),range(64)):
  E=graph_from(q,ia,ib);A=tournament(q,ia,ib)
  need(gaction(E,HGEN)==E,'H-fixed graph');need(taction(A,HGEN)==A,'H-fixed tournament')
  allg[mask(E)]=(ia,ib,q);allt.add(A)
  if ia or ib:
   w=list(IDENTITY)
   if ia:w[0],w[1]=w[1],w[0]
   elif ib==3:w[3],w[4]=w[4],w[3]
   else:
    for i in range(5):w[3+i]=3+((-i)%5)
   need(gaction(E,w)==E and sign(E,w)==-1,'internal odd witness');excluded+=1
 need(len(allg)==512 and len(allt)==512,'full fixed universes')
 graph_types={};normalizer_even=set();extra_odd=[];records=[]
 single_swap=list(IDENTITY);single_swap[8],single_swap[9]=9,8
 for q in range(64):
  E=graph_from(q);canon,even,w,aut,sizes=graph_type_and_sign(E)
  qcan=min(q,swap_q_g(q));ne=not(gaction(E,single_swap)==E and sign(E,single_swap)<0)
  if ne:normalizer_even.add(qcan)
  if ne and not even:extra_odd.append(q)
  if canon in graph_types:need(graph_types[canon][0]==qcan,'no hidden graph merging')
  graph_types[canon]=(qcan,even)
  records.append({'q':q,'canonical_q_under_single_swap':qcan,'graph_mask':mask(E),'false_twin_sizes':sizes,'aut_order':aut,'normalizer_even':ne,'globally_even':even,'odd_witness':w})
 tournaments=set();mod_checks=0
 for q in range(64):
  A=tournament(q)
  C3=[S for S in itertools.combinations(range(N),3) if module(A,S) and strong(A,S)]
  R5=[S for S in itertools.combinations(range(N),5) if module(A,S) and strong(A,S) and prime(A,S)]
  need(C3==[BLOCKS[0]],'unique C3 strong module');need(R5==[BLOCKS[1]],'unique prime strong5 module')
  need(swap_q_t(q)!=q,'singleton swap free on tournaments')
  tournaments.add(min(q,swap_q_t(q)));mod_checks+=1
 even_types=sum(v[1] for v in graph_types.values())
 need((len(tournaments),len(graph_types),len(normalizer_even),even_types)==(32,40,32,31),'final counts')
 need(len(extra_odd)==1,'one extra odd presentation')
 q=extra_odd[0];E=graph_from(q);need(all(sum(x in e for e in E)==5 for x in range(10)),'K5,5 degrees')
 cells,sizes,A=quotient(E);need(sizes==[5,5] and A==[[0,1],[1,0]],'exact K5,5 quotient')
 out={'H_generator':HGEN,'H_order':15,'H_orbit_sizes':[3,5,1,1],
  'labelled_H_fixed_graphs':len(allg),'labelled_H_fixed_tournaments':len(allt),
  'internal_odd_graphs_excluded':excluded,'normalized_tournament_presentations':64,
  'unique_module_presentations_checked':mod_checks,'tournament_types':len(tournaments),
  'remaining_graph_presentations':64,'remaining_graph_types':len(graph_types),
  'normalizer_even_graph_types':len(normalizer_even),'globally_even_graph_types':even_types,
  'extra_odd_graph_q':q,'extra_odd_graph_name':'K5,5','extra_odd_graph_edges':[list(e) for e in sorted(E)],
  'graph_records':records}
 Path('independent_c15.json').write_text(json.dumps(out,indent=2)+'\n')
 print('C15 <= S10 of orbit shape (3,5,1,1):')
 print('  labelled H-fixed tournaments =512; normalized presentations =64; actual types =32')
 print('  labelled H-fixed graphs =512; internal odd exclusions =448; surviving presentations =64')
 print('  surviving graph types =40; normalizer-even types =32; globally sign-even types =31')
 print('  unique additional odd graph =K5,5; quotient-mask =',q)
 print('  all 64 tournament presentations have unique C3 and prime R5 modules')
 print('ALL INDEPENDENT CHECKS PASSED')
if __name__=='__main__':main()
