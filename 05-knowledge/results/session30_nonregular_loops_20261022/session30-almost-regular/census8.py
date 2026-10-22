"""Exact mass certification of full regular9 and almostregular8 universes."""
import numpy as np,itertools,functools,math,json,pathlib,random
from experiment import rot,flip,bits
ROOT=pathlib.Path(__file__).resolve().parent
@functools.lru_cache(None)
def degree_count(d):
 n=len(d)
 if not n:return 1
 if min(d)<0 or max(d)>=n or sum(d)!=n*(n-1)//2:return 0
 total=0
 # Remove distinguished last vertex, choose its outneighbors.
 for out in itertools.combinations(range(n-1),d[-1]):
  out=set(out);rem=tuple(sorted(d[i]-(i not in out) for i in range(n-1)))
  total+=degree_count(rem)
 return total

def signatures(S):
 Q=S@S;R=Q@Q;V=R@Q
 return [(int(S[i].sum()),int(R[i,i]),int(V[i,i])) for i in range(len(S))]

def iso_count(A,B,stop=False):
 n=len(A);sa=signatures(A);sb=signatures(B)
 if sorted(sa)!=sorted(sb):return 0
 candidates=[[j for j in range(n) if sb[j]==sa[i]] for i in range(n)];order=sorted(range(n),key=lambda i:len(candidates[i]));mapping={};total=0
 def rec(k,used):
  nonlocal total
  if k==n:total+=1;return stop
  i=order[k]
  for j in candidates[i]:
   if used>>j&1:continue
   if any(A[i,h]!=B[j,t] for h,t in mapping.items()):continue
   mapping[i]=j
   if rec(k+1,used|1<<j):return True
   del mapping[i]
  return False
 rec(0,0);return total

def main():
 total=degree_count((4,)*9);print('exact_labelled_regular9',total,flush=True);bank=[];mass=0
 def insert(S):
  nonlocal mass
  if any(iso_count(S,T,True) for T in bank):return
  aut=iso_count(S,S);bank.append(S.copy());mass+=math.factorial(9)//aut;print('class',len(bank),'aut',aut,'mass',mass,'of',total,flush=True)
 for signs in itertools.product([-1,1],repeat=4):insert(rot(9,list(signs)))
 insert(np.kron(rot(3),np.ones((3,3),int))+np.kron(np.eye(3,dtype=int),rot(3)))
 S=rot(9);rng=random.Random(30809);trial=0
 while mass<total:
  S=flip(S,rng,10);insert(S);trial+=1
  if trial%1000==0:print('trials',trial,flush=True)
 assert mass==total
 near=[]
 for S in bank:
  for v in range(9):
   A=np.delete(np.delete(S,v,0),v,1)
   if not any(iso_count(A,B,True) for B in near):near.append(A)
 nearrecords=[dict(bits=bits(A),aut=iso_count(A,A)) for A in near];nmass=sum(math.factorial(8)//r['aut'] for r in nearrecords);assert nmass==total
 record=dict(labelled_degree_dp=total,regular9=[dict(bits=bits(A),aut=iso_count(A,A)) for A in bank],almostregular8=nearrecords,regular9_mass=mass,almostregular8_mass=nmass,chain_trials=trial)
 json.dump(record,open(ROOT/'census8_bank.json','w'),indent=2);print('COMPLETE',len(bank),len(near),'labelled',nmass,flush=True)
if __name__=='__main__':main()
