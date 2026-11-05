from pathlib import Path
from itertools import permutations,combinations
from collections import Counter
import json,hashlib,time
from core import require,tournament,inventories,cycles,literal_h
from c3 import *
ROOT=Path(__file__).resolve().parent

def labelled(n):
 edges=list(combinations(range(n),2))
 for mask in range(1<<len(edges)):
  T=[0]*n
  for q,(a,b)in enumerate(edges):T[a if mask>>q&1 else b]|=1<<(b if mask>>q&1 else a)
  yield T

def words(T):
 n=len(T);out=[]
 def visit(p,S):
  if len(p)==n:out.append(tuple(p));return
  for v in range(n):
   if not S>>v&1 and(not p or T[p[-1]]>>v&1):visit(p+[v],S|1<<v)
 visit([],0);return out

def cycle_words(T):
 return [p for p in words(T)if p[0]==0 and T[p[-1]]&1]

def canon(T):
 cs=cycle_words(T);require(bool(cs),'canonicalized weak graph')
 return min(tuple(induced(T,p[q:]+p[:q]))for p in cs for q in range(len(T)))

def hc(T):
 # Direct Hamilton-word DFS, separate from subset endpoint recurrence.
 n=len(T);full=(1<<n)-1;H=0;C=0
 def go(v,S,first):
  nonlocal H,C
  if S==full:
   H+=1
   if first==0 and T[v]&1:C+=1
   return
  B=T[v]&~S
  while B:
   z=B&-B;B-=z;go(z.bit_length()-1,S|z,first)
 for v in range(n):go(v,1<<v,v)
 return H,C

def main():
 records={};table=0;arbitrary=0;literal=0;cone=0;source_classes={}
 for n in range(1,6):
  for V in labelled(n):
   p=inventories(V)
   for c in range(n):
    s=(p['H'],p['e'][c],p['R'][c],p['allH'][((1<<n)-1)^(1<<c)])
    if n<=4:
     for kind in TYPE_I:
      T,S=left_graph(kind,V,c)
      require(actual_response(T,S)==response(kind,s),'universal left response')
      table+=1
    if p['s'][c]==1:
     path=next(q for q in words(V)if q[0]==c)
     canonV=tuple(induced(V,path));source_classes[(n,canonV)]=(list(canonV),0,s)
     if n>1:
      D,E,U,d=s;require(D-E>=d and U>=D+d-E and U>=2*d,'physical source cone');cone+=1
 # Fully arbitrary right grafts, not only ending-one sources.
 for n in range(1,4):
  for T in labelled(n):
   p=inventories(T)
   for mask in range(1,1<<n):
    S={a for a in range(n)if mask>>a&1};h,E,R,K=actual_response(T,S)
    for m in range(1,4):
     for G in labelled(m):
      q=inventories(G)
      for a in range(m):
       A,r=graft(T,S,G,a);actual=inventories(A)['H'];d=q['allH'][((1<<m)-1)^(1<<a)]
       require(actual==h*q['H']-E*q['s'][a]+R*q['L'][a]-K*d,'arbitrary graft scalar')
       if m>1:require(cycles(A)==sum(p['s'][b]for b in S)*q['e'][a],'general graft cycle')
       arbitrary+=1
 # Complete source graph classes through order five; preserve graph classes, not only responses.
 require(Counter(n for n,T in source_classes)==Counter({1:1,2:1,3:2,4:5,5:13}),'source-class inventory')
 classes={n:set()for n in range(5,9)};templates=Counter();physical=0
 src=list(source_classes.values())
 for kind in TYPE_I+tuple('J'+str(k)for k in range(4)):
  for V,c,v in src:
   for W,a,w in src:
    if not admissible(kind,v,w):continue
    order=len(V)+len(W)+(3 if kind=='C'else 2)
    if order>8:continue
    T,r=construct(kind,V,c,W,a);H,C=hc(T)
    require(tournament(T)and H==value(kind,v,w)and C==3,'literal C3 construction')
    classes[len(T)].add(canon(T));classes[len(T)].add(canon(converse(T)))
    templates[kind]+=1;physical+=1
    if len(T)<=6:require(H==inventories(T)['H']and C==cycles(T),'subset/literal agreement');literal+=1
 # Exhaustive labelled universe through six, independent of source enumeration.
 allclasses={5:set(),6:set()};labelled_c3=Counter()
 for n in (5,6):
  for T in labelled(n):
   if cycles(T)==3:
    allclasses[n].add(canon(T));labelled_c3[n]+=1
  require(allclasses[n]==classes[n],'complete labelled C3 source union')
 require({n:len(x)for n,x in classes.items()}=={5:1,6:11,7:38,8:124},'C3 grammar coefficient controls')
 # Exact bounded inverse with physical recipes.
 support_,bank,attempts=support(501);maxorder=0
 for h,record in sorted(support_.items()):
  T,r=realize(record);a,c=hc(T)
  require(a==h and c==3,'physical inverse witness')
  maxorder=max(maxorder,len(T))
 for kind in TYPE_I+tuple('J'+str(k)for k in range(4)):
  for v in bank:
   for w in bank:
    if not admissible(kind,v,w):continue
    H=value(kind,v,w)
    if kind in TYPE_I:
     h,E,R,K=response(kind,v)
     require(H>=h*w[2]and H>=R and h>=2*v[2],'Type I source bounds')
    elif v!=ONE and w!=ONE:
     x=(v[0]-v[1]-v[3],v[2]-v[0],v[1],v[3]);y=(w[0]-w[1]-w[3],w[2]-w[0],w[1],w[3]);k=int(kind[1:])
     require(H==sum(x[i]*POSITIVE[k][i][j]*y[j]for i in range(4)for j in range(4))and H>=3*v[2]*w[2],'Type II positive bound')
    require(v[2]<=H//3 and w[2]<=H//3,'uniform source cutoff')
    if H<=501:require(H in support_,'pruned inverse missed target')
 out=dict(status='PASS',universal_left_response_controls=table,arbitrary_grafts=arbitrary,source_cone_controls=cone,source_graph_class_counts={str(n):sum(k[0]==n for k in source_classes)for n in range(1,6)},literal_template_graphs=physical,template_counts=dict(templates),subset_literal_controls=literal,labelled_C3_counts=dict(labelled_c3),C3_class_counts={str(n):len(v)for n,v in classes.items()},inverse_bound=501,source_states=len(bank),pair_attempts=attempts,support_size=len(support_),max_witness_order=maxorder,support=sorted(support_),source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted(ROOT.glob('*.py'))})
 (ROOT/'SUPPORT.json').write_text(json.dumps(support_,sort_keys=True,indent=2)+'\n')
 (ROOT/'VALIDATION.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
