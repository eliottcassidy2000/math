"""Independent literal vertex-permutation cut polynomials with physical labels."""
from pathlib import Path
import json,hashlib,math
D=Path(__file__).resolve().parent

def need(c,s):
 if not c:raise RuntimeError(s)
def op(R):return [((1<<len(R))-1)^(1<<i)^r for i,r in enumerate(R)]
def transitive(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def substitute(Q,children):
 starts=[];p=0
 for C in children:starts.append(p);p+=len(C)
 R=[0]*p;sets=[]
 for i,C in enumerate(children):
  sets.append(list(range(starts[i],starts[i]+len(C))))
  for j,row in enumerate(C):
   R[starts[i]+j]=row<<starts[i]
   for k,K in enumerate(children):
    if Q[i]>>k&1:R[starts[i]+j]|=((1<<len(K))-1)<<starts[k]
 return R,sets

def profile(R):
 n=len(R);B=max(32,(math.factorial(n)*(1<<(n-1))).bit_length()+1);mask=(1<<B)-1;dp=[[0]*n for _ in range(1<<n)];total=[0]*(1<<n)
 for S in range(1,1<<n):
  rest=S
  while rest:
   b=rest&-rest;rest-=b;v=b.bit_length()-1;M=S^b
   if not M:dp[S][v]=1;continue
   out=total[M]<<B;left=M
   while left:
    x=left&-left;left-=x;u=x.bit_length()-1
    if R[u]>>v&1:out+=dp[M][u]
   dp[S][v]=out
  total[S]=sum(dp[S])
 return [(total[-1]>>(B*k))&mask for k in range(n)]

def augment(R,a,c):
 n=len(R);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(R)]+[(((1<<n)-1)^a)|(1<<(n+1)),((1<<n)-1)^c]

def main():
 C4=[4,1,2,7];C6=substitute([2,4,1],[transitive(2)]*3)[0]
 cases=[('ordered_mixed4',[2,0],[transitive(4),C4]),('ordered_unequal46',[2,0],[C4,C6]),('cycle_nontransitive444',[2,4,1],[C4,transitive(4),C4])]
 records=[];fullcoords=0;nonzero=0;ambients=[]
 for label,Q,children in cases:
  T,sets=substitute(Q,children);U,other=substitute(op(Q),[op(C) for C in children]);need(sets==other,'fixed physical child sets');childprof=[profile(C) for C in children]
  for C,q in zip(children,childprof):need(q==profile(op(C)),'exact corresponding child profile equality')
  queries=[(a,c) for a in range(1<<len(Q)) for c in range(1<<len(Q))] if len(Q)==2 else [(0,0),(1,0),(0,1),(1,2),(2,1),(3,4),(7,0),(0,7)]
  for a,c in queries:
   A=sum(1<<v for i,S in enumerate(sets) if a>>i&1 for v in S);C=sum(1<<v for i,S in enumerate(sets) if c>>i&1 for v in S);TA=augment(T,A,C);UA=augment(U,A,C);p=profile(TA);q=profile(UA);delta=[x-y for x,y in zip(p,q)];need((delta[1]-2*delta[0])%16==0 and all(x%16==0 for x in delta[2:]),'literal complete transfer congruence')
   records.append(dict(label=label,quotient_rows=Q,source_children=children,target_children=[op(X) for X in children],child_profiles=childprof,child_sets=sets,source_rows=T,target_rows=U,core_query=[a,c],physical_query=[A,C],port_names=[len(T),len(T)+1],source_augmented_rows=TA,target_augmented_rows=UA,source_complete_q=p,target_complete_q=q,complete_delta=delta));fullcoords+=2*len(p);nonzero+=delta[0]%16!=0
  # One actually nonzero query in a cyclic exterior, if one occurs.
  candidates=[r for r in records if r['label']==label and r['complete_delta'][0]%16]
  if candidates:
   r=candidates[0];R=[2,4,1];AT,_=substitute(R,[T,[0],[0]]);AU,_=substitute(R,[U,[0],[0]]);n=len(T);A,C=r['physical_query'];A|=1<<(n+1);C|=1<<(n+1);pa=profile(augment(AT,A,C));pb=profile(augment(AU,A,C));d=pa[0]-pb[0];mult=profile(substitute(R,[transitive(2),[0],[0]])[0])[0];need(mult%2==1 and (d-mult*r['complete_delta'][0])%16==0 and d%16,'literal ambient odd multiplier and survival')
   ambients.append(dict(label=label,ambient_rows=R,physical_old_source_rows=AT,physical_old_target_rows=AU,physical_query=[A,C],exterior_names=[n,n+1],port_names=[n+2,n+3],source_augmented_rows=augment(AT,A,C),target_augmented_rows=augment(AU,A,C),source_complete_q=pa,target_complete_q=pb,ambient_H_TT2=mult,source_core_query=r['core_query'],source_local_delta=r['complete_delta'][0],actual_ambient_delta=d));fullcoords+=2*len(pa)
 outputs={'physical_controls.json':records,'ambient_controls.json':ambients}
 for name,obj in outputs.items():
  p=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if p.exists():need(p.read_text()==txt,'retained literal physical bank stable')
  else:p.write_text(txt)
 return dict(status='PASS',method='Independent literal endpoint subset recurrence with exact packed cut-polynomial coefficients; no producer imports.',local_queries=len(records),ambient_queries=len(ambients),nonzero_local_mod16_queries=nonzero,complete_exact_cover_coordinates=fullcoords,maximum_augmented_order=max(len(r['source_augmented_rows']) for r in records+ambients),source_sha256=hashlib.sha256((D/'physical_audit.py').read_bytes()).hexdigest(),output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
