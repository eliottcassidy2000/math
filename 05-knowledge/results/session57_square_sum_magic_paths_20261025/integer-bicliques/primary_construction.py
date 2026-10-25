#!/usr/bin/env python3
import math,json,sys
from checker import require,roots_witness
SOURCE='https://arxiv.org/html/2508.07806v1'
def family(u,v):
 A=[9*(18*u*u+50*u*v+37*v*v)**2,36*(u*u+13*u*v+16*v*v)**2,36*(9*u*u+53*u*v+64*v*v)**2,4*(3*u*u-25*u*v-48*v*v)**2]
 B=[0,20*(u-4*v)*(11*u+19*v)*(u+4*v)*(u+v),36*(3*u+11*v)*(7*u+12*v)*(3*u+4*v)*(u+v),756*(u-v)*(3*u+2*v)*(u+3*v)*(u+2*v)]
 return A,B

def normalize_optimize(A,B):
 if len(set(A))<4 or len(set(B))<4:return None
 sums=[a+b for a in A for b in B]
 require(min(sums)>=0,'negative purported square')
 roots=[math.isqrt(s) for s in sums];require(all(r*r==s for r,s in zip(roots,sums)),'formula is not square')
 h=math.gcd(*roots);require(h>0,'zero family')
 b0=min(B);AA=sorted((a+b0)//(h*h) for a in A);BB=sorted((b-b0)//(h*h) for b in B)
 require(all((a+b0)%(h*h)==0 for a in A) and all((b-b0)%(h*h)==0 for b in B),'normalization not integral')
 lo=1;hi=min(AA)-1
 if hi<lo:return None
 center=(max(AA)-max(BB))//2;anchor=max(lo,min(hi,center))
 # Objective is V-shaped, so optimum is nearest permitted integer to center or its ceiling.
 # There are at most 16 forbidden integer collision values; radius 17 suffices at boundaries too.
 candidates={max(lo,min(hi,anchor+j)) for j in range(-17,19)}
 best=None
 for t in candidates:
  L=[a-t for a in AA];R=[b+t for b in BB]
  if set(L)&set(R):continue
  item=(max(L+R),t,L,R)
  if best is None or item<best:best=item
 if best is None:return None
 N,t,L,R=best
 return {'N':N,'normalization_root_gcd':h,'normalization_b0':b0,'normalized_A':AA,'normalized_B':BB,'gauge_t':t,'left':L,'right':R,'roots':roots_witness(L,R),'left_reflection_defect':L[0]+L[3]-L[1]-L[2],'right_reflection_defect':R[0]+R[3]-R[1]-R[2]}

def run(B=200):
 old=normalize_optimize([324,54756,119716,264196],[0,79200,227205,1258560])
 best=None;count=0;valid=0
 for u in range(-B,B+1):
  for v in range(-B,B+1):
   if math.gcd(u,v)!=1:continue
   count+=1;result=normalize_optimize(*family(u,v))
   if result is None:continue
   valid+=1
   if best is None or (result['N'],u,v)<(best['N'],*best['parameters']):best=dict(result,parameters=[u,v])
 return {'status':'PASS','source':SOURCE,'attribution':'Jean Lagrange family, reproduced by Ajai Choudhry (2025), equation (1.1)','old_example':old,'parameter_universe':{'u_min':-B,'u_max':B,'v_min':-B,'v_max':B,'restriction':'gcd(u,v)=1'},'primitive_pairs_checked':count,'valid_positive_disjoint_pairs':valid,'best_in_bounded_family':best,'claim_scope':'Bounded parameter-family minimum after exact root-gcd normalization and optimal integral gauge, not global K4,4 minimum.'}
if __name__=='__main__':
 try: print(json.dumps(run(int(sys.argv[1]) if len(sys.argv)>1 else 200),indent=2))
 except Exception as e:print(json.dumps({'status':'FAIL','type':type(e).__name__,'message':str(e)}),file=sys.stderr);sys.exit(1)
