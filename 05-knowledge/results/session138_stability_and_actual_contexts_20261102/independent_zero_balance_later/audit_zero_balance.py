"""Independent tridiagonal-context / likelihood-ratio continuation audit.

Imports only the already-audited independent graph/profile kernel, not the
producer's audit.py or any zero_balance_later algorithm. The context kernel
is rebuilt by literal forced cyclic run words. Side profiles use the positive
monomial/Stirling expansion, not alternating finite differences.
"""
from pathlib import Path
from itertools import product
from math import comb,factorial
from collections import Counter
import importlib.util,json,hashlib
BASE=Path(__file__).resolve().parents[1];SOURCE=BASE/'zero_balance_later'
KERNEL=BASE/'independent/audit_profile_dichotomy.py'
spec=importlib.util.spec_from_file_location('independent_profile_kernel',KERNEL);G=importlib.util.module_from_spec(spec);spec.loader.exec_module(G)
def need(c,why):
 if not c:raise RuntimeError(why)
def lift(v):
 return [0]+[r*((v[r] if r<len(v) else 0)+(v[r-1] if r-1<len(v) else 0)) for r in range(1,len(v)+1)]
def logconcave(v):return all(v[i]*v[i]>=v[i-1]*v[i+1] for i in range(1,len(v)-1))
def mlr(v,w):return all(w[j]*v[i]>=w[i]*v[j] for i in range(1,len(v)) for j in range(i+1,len(v)))
def cyclic_words(r,s,t):
 want=(r,s,t);n=r+s+t
 return sum(tuple(sum((start+k)%3==j for k in range(n)) for j in range(3))==want for start in range(3))
def code_tournament(n,c):
 a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if c&1:a[i]|=1<<j
   else:a[j]|=1<<i
   c//=2
 return tuple(a)
def scalar(A,K,B):return sum(A[i]*K[i][j]*B[j] for i in range(1,len(K)) for j in range(1,len(K)))
def single_crossing(h):
 seen_negative=False
 for v in h[1:]:
  if v<0:seen_negative=True
  if v>0 and seen_negative:return False
 return True

def main():
 data=json.loads((SOURCE/'mixed_profile_inputs.json').read_text());M=tuple(data['original_masks']);old=G.profile(M);need(old[0]==241137,'original module H');rows=[];dim=13
 # Literal module substitutions validate the universal run-word law on all
 # labelled tournaments of orders one through three, in each of three slots.
 small=[code_tournament(n,c) for n in range(1,4) for c in range(1<<(n*(n-1)//2))];contexts=0
 for X,A,B in product(small,repeat=3):
  fx,fa,fb=map(G.profile,(X,A,B));expected=sum(fx[r-1]*fa[s-1]*fb[t-1]*cyclic_words(r,s,t) for r in range(1,len(fx)+1) for s in range(1,len(fa)+1) for t in range(1,len(fb)+1));actual=G.subset_h(G.c3_context(X,A,B))[-1];need(actual==expected,'literal run-word substitution');contexts+=1
 # The positive Stirling expansion of (m^3+2m)^a, for the prefix detected by K.
 maxa=64;S=[[0]*(dim+1) for _ in range(3*maxa+1)];S[0][0]=1
 for n in range(1,len(S)):
  for r in range(1,dim+1):S[n][r]=r*S[n-1][r]+S[n-1][r-1]
 profiles={a:[0]+[factorial(r)*sum(comb(a,j)*2**(a-j)*S[a+2*j][r] for j in range(a+1)) for r in range(1,dim+1)] for a in range(1,maxa+1)}
 v=[0,9,126,684,1656,1800,720];need(profiles[2][:7]==v and not any(profiles[2][7:]),'a2 base');need(logconcave(v),'base LC');slacks=[v[r]**2-v[r-1]*v[r+1] for r in range(1,len(v)-1)]
 recursive=[1]+[0]*dim
 for a in range(1,maxa+1):
  L1=lift(recursive);L3=lift(lift(L1));recursive=[L3[r]+2*L1[r] for r in range(dim+1)];need(recursive==profiles[a],('multiplication recurrence',a))
  if a>=2:need(mlr(v,profiles[a]),('side profile MLR',a))
 # Generic lemma controls, not a replacement for the algebraic proof.
 lc_count=0;lc_lifts=0;multiplier_checks=0
 for n in range(1,6):
  for values in product(range(1,7),repeat=n):
   base=[0]+list(values)
   if not logconcave(base):continue
   lc_count+=1;curr=base;power=[base]
   for e in range(1,7):
    nxt=lift(curr);need(logconcave(nxt),'L preserves LC');need(mlr(curr,nxt),'one-step MLR');need(mlr(base,nxt),'iterated MLR');power.append(nxt);curr=nxt;lc_lifts+=1
   for coeff in [(1,0,1),(0,2,0,1),(3,0,0,0,0,0,1)]:
    w=[sum(c*(power[e][r] if r<len(power[e]) else 0) for e,c in enumerate(coeff)) for r in range(len(power[len(coeff)-1]))];need(mlr(base,w),'nonnegative polynomial multiplier MLR');multiplier_checks+=1
 hostile=[0,2,1,2];need(not logconcave(hostile) and not mlr(hostile,lift(hostile)),'log-concavity assumption hostile')
 for inputrow in data['rows']:
  cut=inputrow['cut'];U=tuple(inputrow['child_masks']);d=[0]+[x-y for x,y in zip(G.profile(U),old)]+[0,0];need(d[1:13]==inputrow['difference'],'fixed actual profile difference')
  literal=[[0]*(dim+1) for _ in range(dim+1)];closed=[[0]*(dim+1) for _ in range(dim+1)]
  for s in range(1,dim+1):
   for t in range(1,dim+1):
    literal[s][t]=sum(d[r]*cyclic_words(r,s,t) for r in range(1,13));closed[s][t]=(d[s-1]+3*d[s]+d[s+1]) if s==t else (d[min(s,t)]+d[max(s,t)] if abs(s-t)==1 else 0)
  need(literal==closed,'tridiagonal kernel from actual run words');K=literal;base_gap=scalar(profiles[2],K,profiles[2]);thin= [sum(K[i][j]*profiles[1][j] for j in range(1,dim+1)) for i in range(dim+1)]
  if cut==996:
   need(all(x>=0 for row in K for x in row) and K[1][1]==135618,'arbitrary-side positive kernel')
   # At singleton sides the bound is equality, recounted on actual order14 graphs.
   ambient=G.subset_h(G.c3_context(U,(0,),(0,)))[-1]-G.subset_h(G.c3_context(M,(0,),(0,)))[-1];need(ambient==135618,'sharp arbitrary-side lower bound')
  else:
   boundary={867:4,995:5}[cut];need(all(K[i][j]>=0 for i in range(1,boundary) for j in range(1,dim+1)),'positive prefix rows');need(all(K[i][j]<=0 for i in range(boundary+1,dim+1) for j in range(1,dim+1)),'negative suffix rows');need(K[boundary][boundary-1]>0>K[boundary][boundary] and K[boundary][boundary+1]<0,'mixed row');need(base_gap=={867:-373209154588134,995:-93517858646046}[cut],'negative base exact');need(thin[1]>0 and min(thin)>=0,'thin arbitrary-side positive functional')
   for a in range(2,maxa+1):
    h=[sum(K[i][j]*profiles[a][j] for j in range(1,dim+1)) for i in range(dim+1)];need(single_crossing(h) and all(x<=0 for x in h[7:]),'crossing support');need(sum(v[i]*h[i] for i in range(1,7))<0,'second comparison base')
  signs=Counter()
  for a in range(1,maxa+1):
   for b in range(1,maxa+1):
    value=scalar(profiles[a],K,profiles[b]);sign=1 if value>0 else -1 if value<0 else 0;want=1 if cut==996 or min(a,b)==1 else -1;need(sign==want,('finite sign control',cut,a,b,value));signs[str(sign)]+=1
  rows.append({'cut':cut,'gap_2_2':base_gap,'gap_1_1':scalar(profiles[1],K,profiles[1]),'kernel_diagonal':[K[i][i] for i in range(1,dim+1)],'kernel_upper':[K[i][i+1] for i in range(1,dim)],'thin_profile_functional':thin[1:],'finite_parameter_sign_counts_1_to64':dict(signs)})
 out={'status':'PASS','independence':'Imports prior independent profile kernel only; no producer algorithms imported','literal_three_module_contexts':contexts,'fixed_actual_module_profile_pairs':3,'literal_cyclic_run_word_kernel_entries':3*dim*dim,'positive_stirling_profiles':maxa,'logconcave_sequences_tested':lc_count,'L_logconcavity_MLR_steps':lc_lifts,'nonnegative_multiplier_checks':multiplier_checks,'base_logconcavity_slacks':slacks,'nonlogconcave_hostile':{'v':hostile[1:],'Lv':lift(hostile)[1:],'MLR':False},'rows':rows,'input_sha256':hashlib.sha256((SOURCE/'mixed_profile_inputs.json').read_bytes()).hexdigest(),'independent_kernel_sha256':hashlib.sha256(KERNEL.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()
