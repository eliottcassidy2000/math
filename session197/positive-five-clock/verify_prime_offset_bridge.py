#!/usr/bin/env python3
from verify import hp_subsets,require
import itertools,math,json,functools

def tournament(n,bits):
 T=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if bits>>k&1:T[i]|=1<<j
  else:T[j]|=1<<i
 return T

def profile(T):
 h=hp_subsets(T);full=len(h)-1
 return (len(T),h[-1],sum(h[s]*h[full^s] for s in range(full+1)))

@functools.lru_cache(None)
def prime(q):return q>=2 and all(q%d for d in range(2,math.isqrt(q)+1))

def main():
 profiles=set();h5={};core=[8,1,3,6];graphs=iso=0
 for n in range(1,6):
  for bits in range(1<<(n*(n-1)//2)):
   T=tournament(n,bits);p=profile(T);profiles.add(p);graphs+=1
   require(p[2]-3*p[1]!=0,'boundary factor never zero in census')
   if p[2]-3*p[1]<0:require(p==(1,1,2),'only negative boundary')
   if p[1]==5:
    h5[n]=h5.get(n,0)+1;require(p[2]==32*2**(n-4),'H5 two-bin profile')
    if n==4:
     maps=[per for per in itertools.permutations(range(4)) if all(bool(T[i]>>j&1)==bool(core[per[i]]>>per[j]&1) for i,j in itertools.combinations(range(4),2))]
     require(len(maps)==1,'unique rigid H5 core');iso+=1
 require(h5=={4:24,5:240},'labelled H5 small inventory')
 cases=prime_cases=0;constants=set();five_amplitudes=set()
 for P,U,R,V in itertools.product(sorted(profiles),repeat=4):
  nP,hP,bP=P;nU,hU,bU=U;nR,hR,bR=R;nV,hV,bV=V
  A=2*bU*bV*(bP-hP)*(bR-hR)
  B=bV*(bR-hR)*(4*hP-bP)*hU+bU*(bP-hP)*(4*hR-bR)*hV
  C=hU*hV*(bP-3*hP)*(bR-3*hR)
  cases+=1
  if B==0 and C>=5 and prime(C):
   prime_cases+=1;constants.add(C)
   require(P==(2,1,4) and R==(2,1,4),'prime-offset P/R descent')
   require((hU==C and hV==1) or(hV==C and hU==1),'prime-offset count descent')
   if C==5:
    require(A%1152==0,'five amplitude divisibility');power=A//1152
    require(power>0 and power&(power-1)==0,'five amplitude power of two');five_amplitudes.add(A)
 print(json.dumps({'status':'PASS','physical_graphs':graphs,'physical_profile_types':len(profiles),'coefficient_cases':cases,'pure_prime_offset_cases':prime_cases,'prime_offsets_seen':sorted(constants),'H5_labelled_inventory':h5,'H5_core_isomorphism_checks':iso,'five_amplitudes_in_finite_profile_bank':sorted(five_amplitudes),'all_order_claim':'Prime q>=5 occurs as a pure positive offset in the quotient iff H=q is already realized by an input module.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
