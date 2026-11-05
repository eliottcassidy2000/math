from pathlib import Path
from itertools import combinations
import json,hashlib
from core import *
from bilinear import *
P=Path(__file__).resolve().parent

def hc(T):
 n=len(T);N=1<<n;A=[[0]*n for _ in range(N)];C=[[0]*n for _ in range(N)]
 for v in range(n):A[1<<v][v]=1
 C[1][0]=1
 for S in range(1,N):
  for u in range(n):
   if not S>>u&1:continue
   outs=T[u]&~S
   while outs:
    bit=outs&-outs;outs-=bit;v=bit.bit_length()-1
    A[S|bit][v]+=A[S][u];C[S|bit][v]+=C[S][u]
 return sum(A[-1]),sum(C[-1][u]for u in range(n)if T[u]&1)

def graphs(n):
 E=list(combinations(range(n),2))
 for bits in range(1<<len(E)):
  T=[0]*n
  for i,(a,b)in enumerate(E):T[a if bits>>i&1 else b]|=1<<(b if bits>>i&1 else a)
  yield T

def profile(T,r):
 p=inventories(T)
 return (p['H'],p['e'][r],p['R'][r],p['allH'][((1<<len(T))-1)^(1<<r)]),p['s'][r]

def main():
 arbitrary=literal=lefts=0;small_right=[]
 for n in range(1,4):
  for G in graphs(n):
   p=inventories(G)
   for r in range(n):small_right.append((G,r,(p['H'],p['s'][r],p['L'][r],p['allH'][((1<<n)-1)^(1<<r)])))
 for n in range(2,5):
  for T in graphs(n):
   p=inventories(T)
   for a,b in combinations(range(n),2):
    h,E,R,K=pair_response(T,a,b)
    require(E<=h and K<=R and R-K>=h,'general pair inequalities')
    for G,r,(g,s,L,d)in small_right:
     A,_=graft(T,{a,b},G,r);H,C=hc(A)
     require(H==g*h-s*E+L*R-d*K,'arbitrary pair graft')
     require(H>=h*L,'general hL bound')
     if len(G)>=2:require(C==(p['s'][a]+p['s'][b])*inventories(G)['e'][r],'cycle product')
     if len(A)<=6:require(literal_h(A)==H,'literal arbitrary graft');literal+=1
     arbitrary+=1
 for n in range(1,5):
  for V in graphs(n):
   p=inventories(V)
   for c in range(n):
    if p['s'][c]!=1:continue
    s,_=profile(V,c)
    T,b=glue(V,{c},[0])
    if n>=2:
     require(pair_response(T,b,c)==left_response('N',s),'adjacent response');require(cycles(T)==1,'adjacent C1');lefts+=1
    T,u=glue(T,{b},[0]);require(pair_response(T,u,c)==left_response('M',s),'two-step response');require(cycles(T)==1,'two-step C1');lefts+=1
 basis=((1,0,1,0),(0,0,1,0),(1,1,1,0),(1,0,1,1))
 for kind in('N','M'):
  for i,v in enumerate(basis):
   for j,w in enumerate(basis):require(value(kind,v,w)==MATRICES[kind][i][j],'exact matrix expansion')
 B=1861;S,bank,pairs=support(B)
 require(len(bank)==2048 and len(S)==884,'bounded source/image inventory')
 recipe_controls=0
 for q,recipe in source_bank(40).items():
  T,r=build_source(recipe);got,start=profile(T,r);require(got==q and start==1,'R-bounded physical oracle');recipe_controls+=1
 factors=[];maxorder=0
 for H,record in sorted(S.items()):
  T,r=realize(record);require(hc(T)==(H,2),'all support graphs')
  f=factor_graph(T)
  if f['kind']!='R5':
   a,sa=profile(f['left_rows'],f['left_root']);b,sb=profile(f['right_rows'],f['right_root'])
   require(sa==sb==1 and value(f['kind'],a,b)==H,'physical reverse factor profiles')
  if H in(13,15,41,533,821,1861):factors.append(dict(H=H,rows=T,record=record,factor=f))
  maxorder=max(maxorder,len(T))
 # Literal finite completeness control independent of the constructed examples:
 # every labelled C2 graph through six is factored back to its exact graph.
 labelled=0
 for n in range(3,7):
  for T in graphs(n):
   if cycles(T)!=2:continue
   f=factor_graph(T);labelled+=1
   if f['kind']!='R5':
    a,sa=profile(f['left_rows'],f['left_root']);b,sb=profile(f['right_rows'],f['right_root'])
    require(sa==sb==1 and value(f['kind'],a,b)==inventories(T)['H'],'complete small graph factorization')
 # The retained graft entry is not generally a start-two return.
 q=source_bank(20);V,c=build_source(q[1,0,2,1]);W,a=build_source(q[5,2,6,3]);T,r=construct('N',V,c,W,a);z=inventories(T)
 require((z['H'],z['s'][r],cycles(T))==(41,6,2),'nonreturning marked-root hostile')
 result=dict(status='PASS',bound=B,arbitrary_pair_grafts=arbitrary,literal_graft_counts=literal,left_response_controls=lefts,physical_oracle_controls=recipe_controls,source_quadruples=len(bank),bilinear_pairs_evaluated=pairs,complete_C2_values=len(S),all_support_graph_checks=len(S),max_witness_order=maxorder,complete_labelled_C2_factorizations=labelled,selected_factors=factors,marked_root_hostile=dict(rows=T,root=r,H=41,start=6,C=2))
 (P/'SUPPORT.json').write_text(json.dumps(dict(bound=B,support=S),indent=2,sort_keys=True)+'\n')
 encoded=json.dumps(result,indent=2,sort_keys=True)+'\n';(P/'VALIDATION.json').write_text(encoded);print(encoded,end='')
if __name__=='__main__':main()
