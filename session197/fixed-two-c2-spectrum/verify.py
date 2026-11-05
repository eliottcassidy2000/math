"""Reproduce the all-length language certificate and literal graph controls."""
from pathlib import Path
from itertools import product,combinations,permutations
import json
from core import inventories,cycles,literal_h,birth,require,adjacency_map,pair_response
from families import g,scalar,word,prefix_graph
from spectrum import Q,q_word,from_word,inverse
from control_core import brute,cycle_readout,capped_step
import certificate

BASE=Path(__file__).parent

def hc_literal(T):
 n=len(T)
 if n<3:return 0
 return sum(all(T[a]>>b&1 for a,b in zip((0,)+p,p+(0,)))for p in permutations(range(1,n)))
def hp(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for M in range(1,1<<n):
  for v in range(n):
   if D[M][v]:
    for w in range(n):
     if not M>>w&1 and T[v]>>w&1:D[M|1<<w][w]+=D[M][v]
 return sum(D[-1])
def transitive_after(T,deleted):
 V=[v for v in range(len(T))if v not in deleted]
 return all(not((T[a]>>b&1 and T[b]>>c&1 and T[c]>>a&1)or(T[b]>>a&1 and T[c]>>b&1 and T[a]>>c&1))for a,b,c in combinations(V,3))
def profile4(T,r):
 I=inventories(T);return I['H'],I['e'][r],I['R'][r],I['allH'][((1<<len(T))-1)^(1<<r)],I['s'][r]
def main():
 saved=json.loads((BASE/'CONTROL.json').read_text())
 regenerated=json.loads(json.dumps(certificate.result))
 require(saved==regenerated,'saved control certificate differs')
 nwords=literal=fullvectors=0
 for length in range(7):
  for W in product(range(4),repeat=length):
   T=[2,0];v=brute(T,[0,1]);q=certificate.cls[0];accept=0
   for c in W:
    C=cycle_readout(v,c);T,_=birth(T,[0,1],c,2)
    require(min(C,3)==min(cycles(T),3),'physical cycle readout')
    v=capped_step(v,c,2,3);q,accept=certificate.Q[q][c]
   actual=cycles(T);require(bool(accept)==(actual==2),'controller physical acceptance')
   if length<=3:
    require(v==tuple(min(x,3)for x in brute(T,[0,1])),'full physical capped vector');fullvectors+=1
   if length<=5:require(actual==hc_literal(T),'literal cycle permutations');literal+=1
   phi=(3,1,2,0);U=from_word([phi[c]for c in W[::-1]]);f=[1,0]+list(range(len(T)-1,1,-1))
   require(all((T[a]>>b&1)==(U[f[b]]>>f[a]&1)for a in range(len(T))for b in range(len(T))if a!=b),'physical converse word map')
   nwords+=1
 qseeds=qgraphs=qliteral=0
 for u in range(2,34):
  D,p=g(u),u%2;T=from_word(q_word(u,1));m=len(bin(u-1)[2:])
  require(profile4(T,1)==(3*D+6*u-2*p-2,D+2*u-p-1,3*D+7*u-2*p-3,D+3*u-p-2,2),'Q seed profile')
  A=prefix_graph(u);A,_=birth(A,[0,1],2,2)
  # J adds b beating only the first clone, which is old port0.
  A,_=birth(A,[0,1],1,2)
  f=[m+2,0]+list(range(2,m+2))+[m+3,1]
  require(adjacency_map(A,T,f)and f[-1]==1,'Q exact marked graph map')
  require(cycles(T)==2,'Q seed cycles');qseeds+=1
 for u in range(2,10):
  for v in range(1,9):
   T=from_word(q_word(u,v));h=Q(u,v)
   require(hp(T)==h and cycles(T)==2,'Q scalar/cycles')
   require(transitive_after(T,{0,1})and h>=2*u+v*(4*u+1),'Q feedback/product bound')
   if len(T)<=8:require(literal_h(T)==h,'Q literal paths');qliteral+=1
   qgraphs+=1
 terminals=0
 for u in range(2,18):
  D,p=g(u),u%2;A=[1 if c=='1'else 3 for c in bin(u-1)[2:]];B=[2 if c=='1'else 3 for c in bin(u-1)[2:]]
  require(pair_response(from_word(A),0,1)==(D+u-p,D-p+1,D+3*u-1,u+p),'A pair response')
  require(pair_response(from_word(B+[1]),0,1)==(D+2*u,u+1,2*D+6*u-p-1,u+p),'B1 pair response')
  for k in range(1,5):
   T=from_word(A+[0]*k+[3]);U=from_word(B+[1]+[0]*k+[3])
   require(hp(T)==scalar('N',u,2**k)and cycles(T)==2,'A terminal reduction')
   require(hp(U)==scalar('M',u,2**k)and cycles(U)==2,'B terminal reduction');terminals+=2
  require(hp(from_word(A+[0,3]))==scalar('N',u,2),'M empty suffix reduction')
 require(hp(from_word([0,1,3]))==15 and cycles(from_word([0,1,3]))==2,'regular five')
 raw={}
 for u in range(2,103):
  for v in range(1,103):
   h=Q(u,v)
   if h<=101:raw.setdefault(h,[]).append((u,v))
 for h in range(1,102):
  got=[(x['u'],x['v'])for x in inverse(h)if x['family']=='Q']
  require(got==raw.get(h,[]),'Q independent rectangular inverse check')
 exclusions={str(h):inverse(h)for h in (533,821,1861)}
 require(all(not x for x in exclusions.values()),'restricted spectrum exclusion')
 boundary=[2,4,9,19,39,335,159,319,607,1279,511]
 require(hp(boundary)==533 and cycles(boundary)==2,'H533 boundary graph')
 require(not any(transitive_after(boundary,set(S))for S in combinations(range(11),2)),'H533 no feedback pair')
 require(transitive_after(boundary,{2,5,8}),'H533 feedback triple')
 examples=[]
 for h in (13,15,29,55,231,247):
  witnesses=inverse(h);require(witnesses,'positive inverse control')
  for x in witnesses:
   T=from_word(x['word']);require(hp(T)==h and cycles(T)==2,'inverse witness count')
  examples.append(dict(H=h,witnesses=witnesses))
 out=dict(status='PASS',closed_raw_vectors=len(certificate.V),controller_classes=len(certificate.Q),regex_nfa_states=len(certificate.nfa),complete_language_product_states=len(certificate.queue),physical_words_through_length6=nwords,literal_cycle_words=literal,full_capped_vector_checks=fullvectors,Q_seed_marked_isomorphisms=qseeds,Q_physical_graphs=qgraphs,Q_literal_path_checks=qliteral,terminal_reduction_graphs=terminals,inverse_rectangle_bound=101,exclusions=exclusions,boundary=dict(rows=boundary,H=533,C=2,feedback_triple=[2,5,8],feedback_pairs=0),examples=examples)
 text=json.dumps(out,sort_keys=True,indent=2)+'\n';(BASE/'VALIDATION.json').write_text(text);print(text,end='')
if __name__=='__main__':main()
