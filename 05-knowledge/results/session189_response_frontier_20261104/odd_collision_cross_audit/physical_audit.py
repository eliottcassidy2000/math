"""Literal complete-cover finite differences and exact conversed-child restoration."""
from pathlib import Path
from functools import lru_cache
import json,hashlib,math
from algebra_audit import paths,kernel,weights
D=Path(__file__).resolve().parent

def need(c,s):
 if not c:raise RuntimeError(s)
def op(R):return [((1<<len(R))-1)^(1<<i)^r for i,r in enumerate(R)]
def tt(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def sub(Q,children):
 R=[];sets=[]
 for C in children:sets.append(list(range(len(R),len(R)+len(C))));R.extend([0]*len(C))
 for i,C in enumerate(children):
  shift=sets[i][0]
  for v,row in enumerate(C):R[shift+v]=(row<<shift)|sum(1<<u for j,S in enumerate(sets) if Q[i]>>j&1 for u in S)
 return R,sets

def profile(R):
 n=len(R);B=(math.factorial(n)*(1<<(n-1))).bit_length()+1;mask=(1<<B)-1;at=[[0]*n for _ in range(1<<n)];total=[0]*(1<<n)
 for S in range(1,1<<n):
  rem=S
  while rem:
   b=rem&-rem;rem-=b;v=b.bit_length()-1;M=S^b
   if not M:at[S][v]=1;continue
   z=total[M]<<B;left=M
   while left:
    x=left&-left;left-=x;u=x.bit_length()-1
    if R[u]>>v&1:z+=at[M][u]
   at[S][v]=z
  total[S]=sum(at[S])
 return [(total[-1]>>(B*j))&mask for j in range(n)]
def augment(R,a,c):
 n=len(R);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(R)]+[(((1<<n)-1)^a)|(1<<(n+1)),((1<<n)-1)^c]
CACHE={};QUERIES=[]
def actual(R,a,c):
 key=(tuple(R),a,c)
 if key not in CACHE:
  rows=augment(R,a,c);p=profile(rows);CACHE[key]=len(QUERIES);QUERIES.append(dict(old_rows=R,physical_query=[a,c],port_names=[len(R),len(R)+1],augmented_rows=rows,complete_q=p))
 return CACHE[key],QUERIES[CACHE[key]]['complete_q']
def quartic(R,p,q,r,group,base=0):
 out=[0]*(len(R)+2);refs=[]
 for bits in range(16):
  a=base|((1<<p) if bits&1 else 0)|((1<<q) if bits&2 else 0);c=base|((1<<r) if bits&4 else 0)|(sum(1<<i for i in group) if bits&8 else 0)
  ix,value=actual(R,a,c);sign=(-1)**(4-bits.bit_count());refs.append(ix)
  for j,x in enumerate(value):out[j]+=sign*x
 return out,refs

def module(R,S):
 mask=sum(1<<v for v in S);full=(1<<len(R))-1
 return all(((R[u]^R[v])&(full^mask))==0 for u in S for v in S)
def induced(R,S):return [sum(1<<j for j,v in enumerate(S) if R[u]>>v&1) for u in S]
def mul(a,b):
 n=len(a);out=[0]*n
 for S in range(n):
  A=S
  while True:
   out[S]+=a[A]*b[S^A]
   if A==0:break
   A=(A-1)&S
 return out

def capacity(R):
 n=len(R);G,C,A,H=paths(R);mask=(1<<64)-1;N=(1<<n)-1
 return [[sum((C[u][S]&mask)*(A[v][N^S]&mask) for S in range(N+1)) for v in range(n)] for u in range(n)]
def exterior(R,into):
 n=len(R);N=1<<n;G,C,A,H=paths(R);mask=(1<<64)-1
 C=[[v&mask for v in row] for row in C];A=[[v&mask for v in row] for row in A];G=[[v&mask for v in row] for row in G]
 alpha=[1]+[0]*(N-1);beta=alpha.copy();tau=[0]*N;ell=[[0]*N for _ in R];rho=[[0]*N for _ in R]
 for S in range(1,N):
  for v in range(n):
   if into>>v&1:alpha[S]+=C[v][S]
   else:beta[S]+=A[v][S]
  for u in range(n):
   for v in range(n):
    z=G[S][u*n+v]
    if into>>v&1:ell[u][S]+=z
    if not(into>>u&1):rho[v][S]+=z
    if not(into>>u&1) and into>>v&1:tau[S]+=z
 return alpha,beta,tau,A,C,ell,rho

def general_controls():
 records=[];outside=[2,4,1]
 for inside,into in [([2,4,1],1),([10,4,9,2],3),([4,1,2,7],5),([10,4,9,2],2)]:
  m=len(inside);T=inside.copy()+[row<<m for row in outside]
  for u in range(m):
   for j in range(3):
    if into>>j&1:T[m+j]|=1<<u
    else:T[u]|=1<<(m+j)
  U=T.copy()
  for u,row in enumerate(op(inside)):U[u]=(U[u]&~((1<<m)-1))|row
  X,Y=kernel(inside,0,1,2,paths(inside));alpha,beta,tau,A,C,ell,rho=exterior(outside,into);unit=[1]+[0]*7;powers=[unit]
  for j in range(m):powers.append(mul(powers[-1],tau))
  for i in range(3):
   K=[a-b for a,b in zip(mul(alpha,A[i]),mul(beta,C[i]))];gamma=mul(mul(alpha,beta),[a-b for a,b in zip(ell[i],rho[i])]);pred=sum(x*mul(K,powers[j])[-1] for j,x in enumerate(X))+sum(y*mul(gamma,powers[j])[-1] for j,y in enumerate(Y));s,si=quartic(T,0,1,2,[m+i]);t,ti=quartic(U,0,1,2,[m+i]);need(s[0]-t[0]==pred,'general mixed two-channel formula with nontransitive child')
   records.append(dict(inside_rows=inside,exterior_rows=outside,into_mask=into,physical_marks=[0,1,2,m+i],source_rows=T,target_rows=U,X=X,derivative=Y,first_channel=K,second_channel=gamma,tau=tau,source_complete_quartic=s,target_complete_quartic=t,expected_scalar_contrast=pred,source_query_ids=si,target_query_ids=ti))
 return records

def grouped_controls():
 records=[]
 for m,ambient in [(3,True),(5,False),(5,True),(7,True)]:
  Q=[2,4,1];children=[tt(m),tt(2),tt(2)];T,sets=sub(Q,children);V,_=sub(Q,[op(children[0]),children[1],children[2]]);U,_=sub(op(Q),children);core_n=len(T);root=list(range(core_n));base=0;outside=[]
  if ambient:
   T,_=sub(Q,[T,[0],[0]]);V,_=sub(Q,[V,[0],[0]]);U,_=sub(Q,[U,[0],[0]]);base=1<<(core_n+1);outside=[core_n,core_n+1]
  for S in [root]+sets:need(module(T,S),'all declared original source modules')
  # Root plus each original child gives quotient-only converse, on literal names.
  toggled=T.copy()
  for i in range(len(T)):
   for j in range(i+1,len(T)):
    if sum(i in S and j in S for S in [root]+sets)%2:toggled[i]^=1<<j;toggled[j]^=1<<i
  need(toggled==U,'literal simultaneous original-module XOR target')
  contracted=[sets[0][0]]+list(range(m,len(T)));C=capacity(induced(T,contracted));K={v:C[0][j]-C[j][0] for j,v in enumerate(contracted) if j};grouped=[]
  for group in sets[1:]:
   combo=[0]*(len(T)+2);individual=[]
   for r,w in weights(m).items():
    p,pi=quartic(T,0,1,r,group,base);q,qi=quartic(U,0,1,r,group,base);v,vi=quartic(V,0,1,r,group,base);delta=[a-b for a,b in zip(p,q)];need(delta==[a-b for a,b in zip(p,v)],'all ordered-cover saturated-quartic quotient/local restoration')
    for j,x in enumerate(delta):combo[j]+=w*x
    individual.append(dict(r=r,weight=w,source_quartic=p,quotient_target_quartic=q,local_target_quartic=v,source_query_ids=pi,quotient_query_ids=qi,local_query_ids=vi))
   need(combo[0]==2*sum(K[v] for v in group),'odd binomial full ambient contracted capacity identity');grouped.append(dict(group=group,individual=individual,complete_combined_contrast=combo,contracted_group_skew=sum(K[v] for v in group)))
  need(grouped[0]['complete_combined_contrast'][0]-grouped[1]['complete_combined_contrast'][0]>0,'ordered witness positive grouped contrast')
  qx=profile(induced(T,sets[1]));qy=profile(induced(T,sets[2]));need(qx==qy and all(qx[j]**2>=qx[j-1]*qx[j+1] for j in range(1,len(qx)-1)),'actual identical log-concave complete witness profiles')
  records.append(dict(m=m,ambient=ambient,source_rows=T,local_target_rows=V,quotient_target_rows=U,child_sets=sets,root_names=root,exterior_names=outside,base_physical_mask=base,contracted_physical_names=contracted,complete_capacity_matrix=C,grouped=grouped,witness_profile=qx))
 return records

def anchor_controls():
 P=D.parent/'collision_next'
 if not P.exists():P=D.parent/'session189-collision-next'
 data=json.loads((P/'controls.json').read_text());records=[]
 for r in data['cases']:
  T=r['core_rows'];U=r['core_target_rows'];blocks=r['child_blocks'];A=blocks[0];need(induced(U,A)==op(induced(T,A)),'distinguished child literal converse')
  for S in blocks:need(module(T,S) and module(U,S) and profile(induced(T,S))==profile(induced(U,S)),'exact corresponding complete child profile')
  source=sub(r['ambient_quotient_rows'],[T,[0],[0]])[0];target=sub(r['ambient_quotient_rows'],[U,[0],[0]])[0];restored=sub(r['ambient_quotient_rows'],[op(T),[0],[0]])[0]
  ids=[]
  for q in r['queries']:
   z=[]
   for R in [source,target,restored]:idx,p=actual(R,q['a_mask'],q['c_mask']);z.append((idx,p))
   need(z[1][1]==z[2][1],'complete augmented cover profile restoration at the identical physical query');need([p[0] for idx,p in z]==[q['source'],q['target'],q['whole_core_converse']],'all retained producer scalar controls')
   ids.append(dict(physical_query=[q['a_mask'],q['c_mask']],source_query_id=z[0][0],target_query_id=z[1][0],restored_query_id=z[2][0],complete_source_minus_target=[a-b for a,b in zip(z[0][1],z[1][1])]))
  need(any(z['complete_source_minus_target'][0]%4==2 for z in ids),'conversed child ambient anchor survives')
  records.append(dict(name=r['name'],core_source_rows=T,core_target_rows=U,child_sets=blocks,anchor=r['anchor'],query_records=ids))
 return records

def main():
 general=general_controls();grouped=grouped_controls();anchors=anchor_controls();outputs={'general_mixed_controls.json':general,'all_cover_grouped_controls.json':grouped,'all_cover_anchor_controls.json':anchors,'literal_queries.json':QUERIES}
 for name,obj in outputs.items():
  p=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if p.exists():need(p.read_text()==txt,'retained physical cross-controls stable')
  else:p.write_text(txt)
 return dict(status='PASS',general_nontransitive_mixed_contrasts=len(general),grouped_odd_source_contexts=len(grouped),conversed_child_contexts=len(anchors),conversed_child_actual_queries=sum(len(x['query_records']) for x in anchors),literal_augmented_query_instances=len(QUERIES),exact_complete_cover_coordinates=sum(len(q['complete_q']) for q in QUERIES),maximum_augmented_order=max(len(q['augmented_rows']) for q in QUERIES),source_pins={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['physical_audit.py','algebra_audit.py']},output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
