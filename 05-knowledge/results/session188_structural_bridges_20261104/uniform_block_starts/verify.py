"""Exact uniform-block aggregate-start tomography and arbitrary-third controls."""
from pathlib import Path
from itertools import combinations,product
import json,hashlib,random,math
ROOT=Path(__file__).resolve().parent
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows(n,code):
 a=[0]*n
 for bit,(u,v) in enumerate(combinations(range(n),2)):
  x,y=(v,u) if code>>bit&1 else(u,v);a[x]|=1<<y
 return a
def op(T):return [r^((1<<len(T))-1)^(1<<u) for u,r in enumerate(T)]
def paths(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)]
 for S in range(1,1<<n):
  rem=S
  while rem:
   bit=rem&-rem;rem-=bit;u=bit.bit_length()-1;rest=S^bit
   if not rest:D[S][u]=1;continue
   nxt=T[u]&rest;z=0
   while nxt:
    b=nxt&-nxt;nxt-=b;z+=D[rest][b.bit_length()-1]
   D[S][u]=z
 return D
def hp(T):return sum(paths(T)[-1])
def extend(T,a,c):
 n=len(T);out=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:out[u]|=1<<n
  else:out[n]|=1<<u
  if c>>u&1:out[u]|=1<<(n+1)
  else:out[n+1]|=1<<u
 return out
def profile(T):
 n=len(T);P=paths(T);H=[sum(r) for r in P];H[0]=1;d=[[0]*(n+1) for _ in H];d[0][0]=1
 for S in range(1,1<<n):
  A=S
  while A:
   for k in range(1,n+1):d[S][k]+=H[A]*d[S^A][k-1]
   A=(A-1)&S
 return d[-1]
def tt(n):return rows(n,0)
def ttq(n):
 q=[0,1]
 for m in range(2,n+1):q=[0]+[k*((q[k] if k<len(q) else 0)+q[k-1]) for k in range(1,m+1)]
 return q
def blow(Q,blocks):
 offs=[sum(len(b) for b in blocks[:i]) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,B in enumerate(blocks):
  for u,r in enumerate(B):
   T[offs[i]+u]=r<<offs[i]
   for j,C in enumerate(blocks):
    if Q[i]>>j&1:T[offs[i]+u]|=((1<<len(C))-1)<<offs[j]
 return T,offs
def record(T,U,I):
 s=sum(z for i,z in enumerate(paths(T)[-1]) if I>>i&1);t=sum(z for i,z in enumerate(paths(U)[-1]) if I>>i&1)
 v=[hp(extend(A,a,c)) for A in(T,U) for a,c in((I,0),(0,I))];d=[v[0]-v[2],v[1]-v[3]]
 need(d[0]-2*d[1]==2*(s-t),'physical uniform block identity')
 return {'physical_mask':I,'start_counts':[s,t],'T_queries':v[:2],'U_queries':v[2:],'contrasts':d}
def add(*ps):
 z={}
 for p in ps:
  for e,c in p.items():z[e]=z.get(e,0)+c
 return {e:c for e,c in z.items() if c}
def mul(*ps):
 z={(0,0,0):1}
 for p in ps:
  w={}
  for e,c in z.items():
   for f,d in p.items():
    g=tuple(a+b for a,b in zip(e,f));w[g]=w.get(g,0)+c*d
  z={e:c for e,c in w.items() if c}
 return z
def neg(p,k=-1):return {e:k*c for e,c in p.items()}
def symbolic():
 one={(0,0,0):1};h={(1,0,0):1};r={(0,1,0):1};c={(0,0,1):1}
 def f(r,c):return add(mul(add(one,c),add(one,h)),mul(c,add(one,h,neg(r))))
 def g(r,c):return mul(add(one,c),add(one,h,neg(r)))
 df=add(f(r,c),neg(f(c,r)));dg=add(g(r,c),neg(g(c,r)))
 need(df==mul(neg(add(c,neg(r)),2),add(one,h)),'first query polynomial')
 need(dg==mul(add(c,neg(r)),add(neg(one,2),h)),'second query polynomial')
 need(add(df,neg(dg,-2))==neg(add(r,neg(c)),2),'aggregate identity polynomial')
 return {'first_terms':[[*e,v] for e,v in sorted(df.items())],'second_terms':[[*e,v] for e,v in sorted(dg.items())]}
def main():
 poly=symbolic();singleton=[]
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   T=rows(n,code);U=op(T)
   for I in range(1<<n):singleton.append({'n':n,'code':code,**record(T,U,I)})
 rng=random.Random(188149);module=[]
 for Q,blocks in [([2,4,1],[tt(2),[2,4,1],tt(1)]),([20,1,2,7,14],[tt(2),[0],[2,4,1],[0],[0]]),([2,4,1],[tt(2),tt(3),[8,1,3,6]])]:
  T,offs=blow(Q,blocks);B=[op(A) if rng.randrange(2) else A for A in blocks];U,_=blow(op(Q),B)
  need([profile(a) for a in blocks]==[profile(a) for a in B],'profile preservation')
  for mask in range(1<<len(blocks)):
   I=sum(((1<<len(A))-1)<<o for j,(A,o) in enumerate(zip(blocks,offs)) if mask>>j&1)
   module.append({'quotient':Q,'child_rows':blocks,'new_child_rows':B,'block_union':mask,**record(T,U,I)})
 ratios=0
 for a in range(1,16):
  for b in range(a+1,17):
   A=ttq(a);B=ttq(b)
   at=lambda q,k:q[k] if k<len(q) else 0
   for k in range(1,18):need(at(B,k)*at(A,k-1)-at(B,k-1)*at(A,k)>=0,'Stirling ordered minor');ratios+=1
   need(B[2]-(A[2] if len(A)>2 else 0)==2**b-2**a,'strict second minor')
 arbitrary=[]
 for a,b in((1,2),(2,3),(3,2),(2,4)):
  for C in(tt(2),tt(3),[2,4,1],[8,1,3,6]):
   if a+b+len(C)>9:continue
   blocks=[tt(a),tt(b),C];T,offs=blow([2,4,1],blocks);U,_=blow([4,1,2],blocks);I=((1<<len(C))-1)<<offs[2];z=record(T,U,I)
   A,B,D=ttq(a),ttq(b),profile(C);at=lambda q,k:q[k] if k<len(q) else 0
   expected=sum(D[k]*(at(A,k)*at(B,k-1)-at(A,k-1)*at(B,k)) for k in range(1,len(D)))
   need(z['start_counts'][0]-z['start_counts'][1]==expected,'arbitrary-third run formula')
   need(abs(expected)>=D[2]*abs(2**a-2**b)>0,'strict quantitative separation');need(any(z['contrasts']),'actual two-query separation')
   arbitrary.append({'transitive_sizes':[a,b],'third_rows':C,'third_profile':D,'signed_start_skew':expected,**z})
 # Boundary: equal transitive sizes give zero start skew; the conclusion is not inferred there.
 blocks=[tt(2),tt(2),tt(2)];T,o=blow([2,4,1],blocks);U,_=blow([4,1,2],blocks);boundary=record(T,U,48);need(boundary['start_counts'][0]==boundary['start_counts'][1],'equal-size boundary')

 # Unequal-profile obstruction to using only aggregate starts.
 blocks=[tt(2),tt(2),[2,4,1]];T,offs=blow([2,4,1],blocks);U,_=blow([4,1,2],blocks)
 blind=[record(T,U,((1<<len(A))-1)<<o) for A,o in zip(blocks,offs)]
 need(all(z['start_counts'][0]==z['start_counts'][1] for z in blind),'all-three-start-skew blindness')
 other=[hp(extend(T,3,12)),hp(extend(U,3,12))];need(other==[1959,2043],'nonzero physical query despite start blindness')
 # Original root-containing families, including two crossing TT3 intervals.
 blocks=[tt(2),tt(3),[2,4,1]];T,offs=blow([2,4,1],blocks);N=len(T);full=(1<<N)-1;parts=[3,12,24,224];root_records=[]
 for M in parts:
  need(all(T[u]&M in(0,M) for u in range(N) if not(M>>u&1)),'selected set is original module')
 for choices in range(16):
  selected=[full]+[M for j,M in enumerate(parts) if choices>>j&1];U=list(T)
  for M in selected:
   for u in range(N):
    if M>>u&1:U[u]^=M^(1<<u)
  h,h2=hp(T),hp(U);need(h2>=h,'root-family scalar monotonicity control')
  if h==h2:
   z=record(T,U,224);need(any(z['contrasts']),'equal-scalar root-family uniform separation')
  else:z={'T_queries':[hp(extend(T,224,0)),hp(extend(T,0,224))],'U_queries':[hp(extend(U,224,0)),hp(extend(U,0,224))]}
  root_records.append({'selected_original_modules':selected,'H':[h,h2],**z})
 return {'status':'PASS','symbolic_polynomial_identities':poly,'all_labelled_converse_union_controls':singleton,'general_quotient_profile_preserving_controls':module,'stirling_minor_checks':ratios,'arbitrary_third_controls':arbitrary,'equal_size_boundary':boundary,'all_start_skews_zero_control':{'child_rows':[tt(2),tt(2),[2,4,1]],'uniform_start_records':blind,'other_query':[3,12],'other_values':other},'root_family_controls':root_records,'scope':'Universal two-uniform-query start-skew identity; unequal transitive pair plus arbitrary nontrivial third child separates. No arbitrary ambient lift or full reconstruction theorem.','source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py']}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
