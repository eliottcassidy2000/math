"""Independent full96 prime-five profile certificate and label transport audit."""
from pathlib import Path
from itertools import combinations,permutations
from math import comb
import json,subprocess,tempfile,random,hashlib
root=Path(__file__).resolve().parent;source=root.parent/'prime_five_triples'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows(n,code):
 a=[0]*n
 for bit,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(j,i) if code>>bit&1 else(i,j);a[u]|=1<<v
 return tuple(a)
def op(a):return tuple(((1<<len(a))-1)^r^(1<<i) for i,r in enumerate(a))
def transport(a,p):return tuple(sum(1<<j for j,v in enumerate(p) if a[u]>>v&1) for u in p)
def prime(a):return all(any(a[u]&S not in(0,S) for u in range(5) if not S>>u&1) for S in range(1,31) if S.bit_count()>1)
def strong(a):
 R=[r|(1<<i) for i,r in enumerate(a)]
 for k in range(5):
  for i in range(5):
   if R[i]>>k&1:R[i]|=R[k]
 return all(r==31 for r in R)
def qwords(a):
 n=len(a);q=[0]*(n+1)
 for w in permutations(range(n)):
  bad=sum(not(a[u]>>v&1) for u,v in zip(w,w[1:]));good=n-1-bad
  for k in range(bad+1,n+1):q[k]+=comb(good,k-1-bad)
 return q
def blow(Q,blocks):
 out=[0]*15
 for i,B in enumerate(blocks):
  for u,r in enumerate(B):
   out[3*i+u]=r<<(3*i)
   for j in range(5):
    if Q[i]>>j&1:out[3*i+u]|=7<<(3*j)
 return tuple(out)
def extend(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:R[u]|=1<<n
  else:R[n]|=1<<u
  if c>>u&1:R[u]|=1<<(n+1)
  else:R[n+1]|=1<<u
 return R
def recover(values):
 d=list(values)
 for bit in range(5):
  for S in range(32):
   if S>>bit&1:d[S]-=d[S^(1<<bit)]
 for S in range(32):
  need(d[S]%(1<<S.bit_count())==0,'profile coefficient divisibility');d[S]//=1<<S.bit_count()
 return [sum(d[T]*(-1)**((T^S).bit_count()) for T in range(32) if T&S==S) for S in range(32)]
def main():
 data=json.loads((source/'receipt.json').read_text());recs=data['representative_certificates'];reps=[tuple(r['quotient']) for r in recs]
 types={1:[],3:[]};triple=[]
 for code in range(8):
  A=rows(3,code);q=qwords(A);need(q in([0,1,6,6],[0,3,6,6]),'triple profile classification');types[q[1]].append(A);triple.append(q)
 need([len(types[1]),len(types[3])]==[6,2],'labelled triple counts')
 labelled=[];counts=[0]*3
 for code in range(1024):
  Q=rows(5,code)
  if not prime(Q) or not strong(Q):continue
  matches=[]
  for i,R in enumerate(reps):
   pp=next((p for p in permutations(range(5)) if transport(Q,p)==R),None)
   if pp is not None:matches.append((i,pp))
  need(len(matches)==1,'quotient classification not exactly one class');i,p=matches[0];counts[i]+=1
  need(transport(op(Q),p)==op(reps[i]),'opposite quotient used different label map')
  for mask in range(32):
   moved=sum(((mask>>p[j])&1)<<j for j in range(5))
   restored=sum(((moved>>j)&1)<<p[j] for j in range(5));need(restored==mask,'profile type transport')
  labelled.append({'code':code,'representative':i,'canonical_to_physical_children':p,'physical_query':[7<<(3*p[0]),7<<(3*p[1])]})
 need(len(labelled)==264 and counts==[120,120,24],'exact prime-five census')
 jobs=[];meta=[];rng=random.Random(188131)
 for i,R in enumerate(reps):
  for mask in range(32):
   h=[3 if mask>>j&1 else 1 for j in range(5)]
   p=list(range(5));rng.shuffle(p)
   moved=[0]*5
   for u in range(5):moved[p[u]]=sum(1<<p[v] for v in range(5) if R[u]>>v&1)
   blocks=[None]*5
   for u in range(5):blocks[p[u]]=rng.choice(types[h[u]])
   need(transport(moved,p)==R,'literal label transport setup')
   before=len(jobs)
   for Q in(R,op(R)):
    jobs.append('R '+' '.join(map(str,Q))+' 1 2 '+' '.join(map(str,h)))
   canonical_blocks=[types[h[j]][0] for j in range(5)]
   for Q in(R,op(R)):
    T=extend(blow(Q,canonical_blocks),7,56);jobs.append('P 17 '+' '.join(map(str,T)))
   for Q in(tuple(moved),op(moved)):
    T=extend(blow(Q,blocks),7<<(3*p[0]),7<<(3*p[1]));jobs.append('P 17 '+' '.join(map(str,T)))
   meta.append((i,mask,p,blocks,before))
 with tempfile.TemporaryDirectory() as d:
  exe=Path(d)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(root/'count.cpp'),'-o',str(exe)],check=True)
  values=list(map(int,subprocess.check_output([str(exe)],input=(str(len(jobs))+'\n'+'\n'.join(jobs)+'\n').encode()).split()))
 need(len(values)==len(jobs),'batch output count')
 certificates=[];by_rep=[[[],[]] for _ in reps]
 for i,mask,p,blocks,start in meta:
  got=values[start:start+6];expected=recs[i]['profile_pattern_queries'][mask];want=[expected['T'],expected['U']]
  need(got==want*3,'run/physical/transported values mismatch');need(got[0]!=got[1],'profile pattern not separated')
  by_rep[i][0].append(got[0]);by_rep[i][1].append(got[1])
  certificates.append({'representative':i,'profile_mask':mask,'canonical_to_physical':p,'actual_named_child_rows':blocks,'run_counts':got[:2],'physical_counts':got[2:4],'transported_physical_counts':got[4:]})
 for i,rec in enumerate(recs):
  need(recover(by_rep[i][0])==rec['T_profile_polynomial'],'complete T profile polynomial')
  need(recover(by_rep[i][1])==rec['U_profile_polynomial'],'complete U profile polynomial')
  gap=min(abs(a-b) for a,b in zip(*by_rep[i]));need(gap==rec['minimum_absolute_gap'],'minimum gap')
 return {'status':'PASS','labelled_prime_quotients':264,'class_sizes':counts,'triple_profiles':triple,'complete_profile_patterns':96,'independent_weighted_run_counts':192,'literal_physical_counts':384,'profile_polynomial_coefficients':192,'label_transport_records':labelled,'all96_certificates':certificates,'producer_pins':{f:hashlib.sha256((source/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','literal.cpp','receipt.json','receipt_optimized.json']},'scope':'Finite exact prime-five quotient with triple children; no other class inferred.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
