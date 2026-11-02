"""Independent ordered-set-partition audit of the ordinal profile dichotomy.

No producer imports. H on subsets is computed by backward last-endpoint DP.
Unordered path partitions are counted by forcing the least vertex into the
first chosen part; multiplying the r-part coefficient by r! orders paths.
This is not the producer's permutation/separator polynomial recurrence.
"""
from pathlib import Path
from itertools import product,permutations
from functools import lru_cache
from math import factorial,comb
from collections import Counter
import json,hashlib
SOURCE=Path(__file__).resolve().parents[1]/'zero_rigidity'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def valid(a):
 n=len(a);return all(not(a[i]>>i&1) for i in range(n)) and all(((a[i]>>j)&1)+((a[j]>>i)&1)==1 for i in range(n) for j in range(i))
def subset_h(a):
 n=len(a);N=1<<n;end=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 incoming=[sum(1<<i for i in range(n) if a[i]>>j&1) for j in range(n)]
 for mask in range(1,N):
  bits=mask
  while bits:
   bit=bits&-bits;bits-=bit;j=bit.bit_length()-1;rest=mask^bit
   if not rest:end[mask][j]=1
   else:
    pred=rest&incoming[j];s=0
    while pred:
     b=pred&-pred;pred-=b;s+=end[rest][b.bit_length()-1]
    end[mask][j]=s
  h[mask]=sum(end[mask])
 return h
@lru_cache(None)
def profile(a):
 need(valid(a),'not tournament');n=len(a);h=subset_h(a);g=[None]*(1<<n);g[0]=[1]
 for mask in range(1,1<<n):
  least=mask&-mask;other=mask^least;p=[0]*(mask.bit_count()+1);sub=other
  while True:
   part=sub|least;rest=other^sub;weight=h[part]
   for r,value in enumerate(g[rest]):p[r+1]+=weight*value
   if sub==0:break
   sub=(sub-1)&other
  g[mask]=p
 return tuple(factorial(r)*g[-1][r] for r in range(1,n+1))
def blow(q):
 k=len(q);a=[0]*(3*k)
 for i in range(k):
  for t in range(3):
   a[3*i+t]|=1<<(3*i+(t+1)%3)
   for j in range(k):
    if q[i]>>j&1:a[3*i+t]|=7<<(3*j)
 return tuple(a)
def ordinal(k):return blow(tuple(sum(1<<j for j in range(i+1,k)) for i in range(k)))
def delete(a,r):
 vs=[i for i in range(len(a)) if i!=r]
 return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs),vs
def add(a,S):
 n=len(a);return tuple(row|((1<<n) if not(S>>i&1) else 0) for i,row in enumerate(a))+(S,)
def canonical(k,r,t):return t[r]==1 and all(t[i]==(0 if i<r else 3) for i in range(k) if i!=r)
def monotone(t):return not any(t[i]>0 and t[j]<3 for i in range(len(t)) for j in range(i+1,len(t)))
def low(k,r,t):
 m=t[r];external=t[:r]+t[r+1:];canonext=all(t[i]==(0 if i<r else 3) for i in range(k) if i!=r)
 return (m==2 and canonext) or (m==0 and all(t[i]==0 for i in range(r)) and monotone(external)) or (m==3 and all(t[i]==3 for i in range(r+1,k)) and monotone(external))
def reverse_state(k,r,t):return k-1-r,tuple(((3-x) if i!=r else (3-x)//2+2*((3-x)%2)) for i,x in reversed(list(enumerate(t))))
def core_type(k,r,t):
 if canonical(k,r,t):return 'canonical'
 if low(k,r,t):return 'low'
 if t[r]==1:return 'cyclic-positive'
 if t[r]==3:
  r,t=reverse_state(k,r,t)
 need(t[r] in (0,2),'pair case reduction')
 if t[r]==2:
  later=next((i for i in range(r+1,k) if t[i]<3),None)
  if later is None:r,t=reverse_state(k,r,t);later=next((i for i in range(r+1,k) if t[i]<3),None)
  need(later is not None,'missing middle-pair defect');return (2,0,(2,t[later]))
 earlier=next((i for i in range(r) if t[i]>0),None)
 if earlier is not None:return (2,1,(t[earlier],0))
 pair=next(((i,j) for i in range(r+1,k) for j in range(i+1,k) if t[i]>0 and t[j]<3),None)
 need(pair is not None,'missing nonmonotone pair');i,j=pair;return (3,0,(0,t[i],t[j]))
def literal_state(a,vs,r,S):
 k=len(a)//3;t=[]
 for i in range(k):
  if i==r:t.append(((S>>vs.index(3*i+1))&1)+2*((S>>vs.index(3*i+2))&1))
  else:t.append(sum((S>>vs.index(3*i+j))&1 for j in range(3)))
 return tuple(t)
def make(k,r,t):
 T=ordinal(k);B,vs=delete(T,3*r);chosen=[]
 for i,x in enumerate(t):
  chosen.extend(([3*i+j+1 for j in range(2) if x>>j&1]) if i==r else [3*i+j for j in range(x)])
 S=sum(1<<vs.index(v) for v in chosen);return T,add(B,S)
def signed_difference(U,T):return tuple(x-y for x,y in zip(profile(U),profile(T)))
def literal_cover(a):
 n=len(a);f=[0]*n
 for word in permutations(range(n)):
  mandatory=sum(not(a[word[i]]>>word[i+1]&1) for i in range(n-1))
  for r in range(mandatory+1,n+1):f[r-1]+=comb(n-1-mandatory,r-1-mandatory)
 return tuple(f)
def c3_context(T,A,B):
 blocks=[T,A,B];off=[0,len(T),len(T)+len(A)];n=sum(map(len,blocks));rows=[0]*n
 for i,block in enumerate(blocks):
  following=(i+1)%3
  for v,row in enumerate(block):rows[off[i]+v]=(row<<off[i])|(((1<<len(blocks[following]))-1)<<off[following])
 return tuple(rows)
def strong(a):
 n=len(a)
 for start in range(n):
  reached=1<<start
  while True:
   new=reached
   for i in range(n):
    if reached>>i&1:new|=a[i]
   if new==reached:break
   reached=new
  if reached!=(1<<n)-1:return False
 return True

def main():
 cert_rows=[json.loads(line) for line in (SOURCE/'profile_core_certificate.out').read_text().splitlines() if line.startswith('{')]
 expected={(x['k'],x['root'],tuple(x['cut_profile'])):x for x in cert_rows};need(len(expected)==15,'certificate count')
 cores=[(2,0,(2,t)) for t in range(3)]+[(2,1,(t,0)) for t in range(1,4)]+[(3,0,(0,s,t)) for s in range(1,4) for t in range(3)]
 core_results=[]
 for key in cores:
  T,U=make(*key);diff=signed_difference(U,T);r=expected[key]
  need(tuple(r['core_masks'])==U,('core construction',key));need(diff==tuple(r['difference']),('core vector',key));need(min(diff)>=0 and any(diff),('core positivity',key))
  if key[0]==2:need(profile(U)==literal_cover(U),'six-vertex literal profile')
  core_results.append({'k':key[0],'root':key[1],'cut_profile':key[2],'difference':diff})
 lower=tuple(min(r['difference'][i] for r in core_results if r['k']==3) for i in range(9));need(lower==(4,648,8088,35040,67920,60480,20160,0,0),'nine-vertex lower envelope')
 leaf_count=0
 for line in (SOURCE/'profile_leafs.txt').read_text().splitlines():
  z=list(map(int,line.split()));n=z[0];need(len(z)==1+2*n,'leaf format');need(profile(tuple(z[1:1+n]))==tuple(z[1+n:]),'leaf entire profile');leaf_count+=1
 need(leaf_count==17,'all leaf certificates')
 counts=Counter();zero_counts=[];literal_cuts=0;context_checks=0
 for k in range(1,4):
  T=ordinal(k);f=profile(T);need(f[0]==3**k,'ordinal H');need(f[-1]==factorial(3*k),'last cover coordinate')
  for r in range(k):
   zeros=0
   for t in product(range(4),repeat=k):
    _,U=make(k,r,t);d=signed_difference(U,T);kind=core_type(k,r,t);counts[str(kind if isinstance(kind,str) else 'positive-core')]+=1
    if kind=='canonical':need(all(x==0 for x in d),'canonical profile')
    elif kind=='low':need(max(d)<=0 and d[0]<0,('low domination',k,r,t,d))
    else:need(min(d)>=0 and (d[0]>0 or (d[0]==0 and d[1]>0)),('high domination',k,r,t,d))
    if k<=2:
     for A,B in [((0,),(0,)),((2,0),(0,)),((2,4,1),(2,4,1))]:
      before=subset_h(c3_context(T,A,B))[-1];after=subset_h(c3_context(U,A,B))[-1];need((after==before)==canonical(k,r,t),'restricted strong context zero');need((after-before)*(d[0] if d[0] else d[1])>0 if any(d) else after==before,'strict response direction');context_checks+=1
   B,vs=delete(T,3*r)
   for S in range(1<<len(B)):
    literal_cuts+=1;t=literal_state(T,vs,r,S);_,rep=make(k,r,t);h=subset_h(add(B,S))[-1];need(h==profile(rep)[0],'triangle-rotation quotient H')
    predicted=canonical(k,r,t)
    if t[r]==2:
     for near,w in [(r-1,1),(r+1,2)]:
      if 0<=near<k:predicted |= all(t[i]==(w if i==near else (0 if i<r else 3)) for i in range(k) if i!=r)
    need((h==3**k)==predicted,'ordinal zero classification');zeros+=h==3**k
   need(zeros==(1 if k==1 else (4 if r in (0,k-1) else 7)),'labelled zero count');zero_counts.append([k,r,zeros])
 reductions=Counter();seen=set();total=0
 for k in range(1,9):
  for r in range(k):
   for t in product(range(4),repeat=k):
    c=core_type(k,r,t);total+=1
    if isinstance(c,tuple):need(c in expected,'unlisted core');seen.add(c);reductions['finite_core']+=1
    else:reductions[c]+=1
 need(seen==set(expected),'fifteen core reachability')
 # Actual strong-Q4 hostile, rebuilt from pair-bit code convention.
 q=[0]*4;code=4
 for i in range(4):
  for j in range(i+1,4):
   if code&1:q[i]|=1<<j
   else:q[j]|=1<<i
   code>>=1
 T=blow(tuple(q));B,vs=delete(T,0);need(strong(tuple(q)) and strong(T) and strong(B),'strong hostile source')
 hostile_expected={r['cut']:r for r in (json.loads(x) for x in (SOURCE/'strong_profile_hostile.out').read_text().splitlines() if x.startswith('{'))}
 hostile=[]
 for cut in [867,995,996]:
  U=add(B,cut);need(strong(U),'hostile completion strong');d=signed_difference(U,T);e=hostile_expected[cut];need(tuple(e['parent_masks'])==B and tuple(e['child_masks'])==U,'hostile masks');need(d==tuple(e['difference']),'hostile full profile');need(min(d)<0<max(d),'mixed signs');need(d[0]!=0,'not a zero counterexample')
  actual_context_gap=subset_h(c3_context(U,(0,),(0,)))[-1]-subset_h(c3_context(T,(0,),(0,)))[-1];need(actual_context_gap==3*d[0]+d[1]==e['C3_singleton_context_gap'],'literal singleton-context difference')
  hostile.append({'cut':cut,'original_H':profile(T)[0],'new_H':profile(U)[0],'difference':d,'singleton_context_gap':actual_context_gap,'strong':True})
 result={'status':'PASS','method':'backward Hamilton-subset DP plus minimum-vertex anchored unordered set-partition DP; no producer imports','finite_core_profiles':15,'complete_leaf_profiles':leaf_count,'nine_vertex_lower_envelope':lower,'all_block_profile_cases_k_le3':sum(counts.values()),'block_profile_class_counts':dict(counts),'literal_cut_cases_k_le3':literal_cuts,'labelled_zero_counts':zero_counts,'literal_cyclic_context_checks':context_checks,'bounded_core_case_reduction_states_k_le8':total,'reduction_categories':dict(reductions),'all_fifteen_core_types_reached':len(seen),'hostiles':hostile,'core_certificate_sha256':hashlib.sha256(json.dumps(core_results,sort_keys=True,separators=(',',':')).encode()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
