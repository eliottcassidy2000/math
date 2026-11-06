"""Independent exact actual-index, run, CDF, and rectangle controls.
No author imports. Decimal or floating-point arithmetic is not used.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from bisect import bisect_left,bisect_right
from functools import lru_cache
import hashlib,json

def need(x,m):
 if not x:raise RuntimeError(m)
def reverse(A,k):return int(format(A,'0'+str(k)+'b')[::-1],2)if k else 0
@lru_cache(None)
def g(u):
 # Independent weighted-bit formula, previously independently proved from sources.
 need(u>=1,'positive scalar index');m=u-1
 return 1+sum((1-((m>>i)&1))*((m>>(i+1))+((m>>(i+1))&1))for i in range(m.bit_length()))
def source(k,A):return (1<<(k+1))-reverse(A,k)
def profile(k,A):
 D,E,U,d=1,0,2,1
 for i in range(k):
  if(A>>i)&1:D,E,U,d=D-E+U,d,2*U-d,d
  else:D,E,U,d=D,0,2*U,d
 return D,E,U,d

def closed_runs(k,A):
 out=[];i=0
 while i<k:
  if not(A>>i)&1:i+=1;continue
  start=i
  while i<k and(A>>i)&1:i+=1
  length=i-start
  # Midpoint t at the run start, calculated from its entire preceding prefix.
  s=F(2*(A% (1<<start))+1,1<<(start+1))
  increment=F(((1<<length)-1)*(2*(A%(1<<start))+1)+(1<<(start+1)),1<<(start+length+1))
  need(0<s<=F(1,2)and 0<increment<=F(3,4),'run start/positive increment')
  need(increment==s+(1-s)/(1<<length),'closed run increment')
  out.append((start,length,s,increment))
 return out

def cdf_extrema(slopes):
 C=Counter(slopes);N=len(slopes);before=0;positive=negative=F(0)
 for x,multiplicity in sorted(C.items()):
  negative=max(negative,x-F(before,N))
  # Strict cutoff cannot move above one. Handle that endpoint explicitly.
  if x<1:positive=max(positive,F(before+multiplicity,N)-x)
  before+=multiplicity
 negative=max(negative,F(C.get(F(1),0),N))
 return positive,negative
rational=lambda x:str(x.numerator)+'/'+str(x.denominator)
words=runcount=constructors=0;details=[]
for k in range(15):
 B=1<<k;delta=F(3*((k+1)//2),4);eps=delta/(B+1);slopes=[];indices=[];maxerr=F(0)
 for A in range(B):
  u=source(k,A);D,E,U,d=profile(k,A);t=F(2*A+1,2*B);e=g(u)-u*t;runs=closed_runs(k,A)
  need((D,E,U,d)==(g(u),u%2,u,1),'actual full source profile')
  need(B<u<=2*B and u>=2,'exact original index interval')
  need(e==sum((r[3]for r in runs),F(0)),'exact whole-word run sum')
  need(0<=e<=F(3*len(runs),4)<=delta<=F(3*(k+1),8),'full run-count bound')
  slope=F(g(u),u);need(t<=slope<=t+eps and 0<=slope<=1,'one-sided slope displacement')
  indices.append(u);slopes.append(slope);maxerr=max(maxerr,e);words+=1;runcount+=len(runs)
 need(sorted(indices)==list(range(B+1,2*B+1)),'bijective dyadic index labeling')
 pos,neg=cdf_extrema(slopes)
 need(pos<=F(1,2*B)and neg<=F(1,2*B)+eps,'one-sided strict-cutoff CDF bounds')
 for j in range(3*B+1):
  theta=F(j,3*B);A=min(B-1,j//3);u=source(k,A);slope=F(g(u),u)
  need(theta-F(1,2*B)<=slope<=theta+F(1,2*B)+eps,'actual approximate-index constructor at thirds')
  constructors+=1
 need(source(k,0)==2*B and g(2*B)==1,'all-O source')
 if k:
  need(source(k,B-1)==B+1 and g(B+1)==B+1,'all-Z source')
  need(g(B+1)-(B+1)*F(2*B-1,2*B)==F(B+1,2*B),'terminal all-Z exact error')
 else:need(slopes==[F(1,2)]and maxerr==delta==0,'empty word boundary')
 details.append(dict(k=k,source_count=B,max_error=rational(maxerr),positive_CDF_error=rational(pos),negative_CDF_error=rational(neg)))
# Exact CDF sandwich at rational cutoffs, including 0,1 and shifted thresholds.
for k in range(7):
 B=1<<k;eps=F(3*((k+1)//2),4*(B+1));ss=sorted(F(g(source(k,A)),source(k,A))for A in range(B));tt=[F(2*A+1,2*B)for A in range(B)]
 for z in set(ss+tt+[F(j,3*B)for j in range(3*B+1)]):
  need(bisect_left(tt,z-eps)<=bisect_left(ss,z)<=bisect_left(tt,z),'strict cutoff sandwich at endpoints')
# Dyadic-block APs and real-cutoff rounding using thirds, not author half-grid controls.
ap_count=0
for k in range(8):
 B=1<<k
 for ell in range(k+1):
  M=1<<ell;step=B//M
  for q in range(B//M):
   vals=sorted(reverse(A,k)for A in range(q*M,(q+1)*M));offset=reverse(q,k-ell);eta=F(offset,step)
   need(vals==list(range(offset,offset+M*step,step))and 0<=eta<1,'reversed aligned block')
   for h in range(3*B+1):
    y=F(h,3*B);x=M*y-eta;ceiling=-((-x.numerator)//x.denominator);actual=sum(3*r<h for r in vals)
    need(actual==ceiling and 0<=ceiling<=M,'exact untruncated AP ceiling')
    need(-eta<=actual-M*y<1-eta and abs(actual-M*y)<1,'strict AP count error')
    ap_count+=1
# Rectangle extrema by sorted selected slices and left/right limits; no 2D prefix sum.
joint=[]
for k in range(10):
 B=1<<k;rev=[reverse(A,k)for A in range(B)];neg=pos=aligned=both=0
 for a in range(B+1):
  strict=sorted(rev[:a]);inclusive=sorted(rev[:min(a+1,B)])
  for b in range(B+1):
   left=bisect_left(strict,b)
   neg=max(neg,a*b-B*left)
   pos=max(pos,B*bisect_right(inclusive,b)-a*b)
   aligned=max(aligned,B*bisect_right(strict,b)-a*b)
   both=max(both,abs(B*left-a*b))
 real=F(max(neg,pos),B*B);first=F(max(neg,aligned),B*B);grid=F(both,B*B)
 if k:need(first<=F(k,B)and real<=F(k+1,B),'grid versus real rectangle bound')
 else:need((real,first,grid)==(1,1,0),'singleton origin rectangle boundary')
 if k==1:need(real==F(3,4)>F(k,B),'hostile to unrounded k/B real bound')
 joint.append(dict(k=k,real=rational(real),first_aligned=rational(first),both_aligned=rational(grid)))
# Fixed exact real hostile, not just a limiting informal example.
need(F(1)-F(9,16)**2>F(1,2),'cutoffs 9/16 give discrepancy 175/256')
# Independent complete small pair controls of the exact J residual.
bilinear=0
for k in range(7):
 for ell in range(7):
  Dk=F(3*((k+1)//2),4);Dl=F(3*((ell+1)//2),4)
  for A in range(1<<k):
   u=source(k,A);tu=F(2*A+1,2*(1<<k));eu=g(u)-u*tu
   for C in range(1<<ell):
    v=source(ell,C);tv=F(2*C+1,2*(1<<ell));ev=g(v)-v*tv
    n=g(u)*g(v)+u*v-(u%2)*(v%2);res=F(n+(u%2)*(v%2))-u*v*(1+tu*tv)
    need(res==u*tu*ev+v*tv*eu+eu*ev and 0<=res<=u*Dl+v*Dk+Dk*Dl,'exact J residual and bound')
    bilinear+=1
print(json.dumps(dict(status='PASS',physical_words=words,Z_runs=runcount,constructor_checks=constructors,
 dyadic_block_real_cutoff_checks=ap_count,J_pair_controls=bilinear,slope_details=details,joint_details=joint,
 fixed_real_hostile=dict(k=1,x='9/16',y='9/16',discrepancy='175/256',false_bound='1/2'),
 scope='Actual approximate slopes and exact residual bounds only; no exact target occupancy.',
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
