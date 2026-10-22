from fractions import Fraction as F
from math import comb
import json

def require(p,s):
 if not p:raise RuntimeError(s)
def v2(n):return (n&-n).bit_length()-1
def orbit(n,h):
 x=n;A=0;C=0;first=None
 for j in range(1,h+1):
  y=3*x+1;a=v2(y);C=3*C+2**A;A+=a;x=y>>a
  if x<n and first is None:first=j
 return x,A,C,first
rows=[]
for h in range(1,201):
 P=3**h;s0=P.bit_length();require(2**(s0-1)<P<2**s0,'strict threshold')
 beta=F(3,2)**h-1;threshold=beta/(1-F(P,2**s0))
 bad=sum((F(comb(s-1,h-1),2**s) for s in range(h,s0)),F(0))
 # Independent negative binomial / binomial stopping identity.
 other=sum((F(comb(s0-1,j),2**(s0-1)) for j in range(h,s0)),F(0))
 require(bad==other,'two exact tail formulas')
 require(bad**5<=F(65536,84375)**h,'Chernoff fifth-power bound')
 if h<=20 or h in (50,100,200):rows.append({'h':h,'first_contracting_clock':s0,'bad_density':str(bad),'endpoint_descent_density':str(1-bad),'uniform_source_threshold':str(threshold)})
# Literal sources, finite exceptional set explicitly permitted, no convergence test.
checked=0;exceptions=[];monotonic_hostile=None
for k in range(3,1000,2):
 n=(3**k+1)//4
 for h in range(1,31):
  x,A,C,first=orbit(n,h);P=3**h;s0=P.bit_length();beta=F(3,2)**h-1
  require(F(C,2**A)<=beta,'uniform carry bound')
  threshold=beta/(1-F(P,2**s0))
  if n>threshold:require((x<n)==(A>=s0),'endpoint iff above threshold')
  elif (x<n)!=(A>=s0):exceptions.append({'k':k,'h':h,'n':n})
  if first is not None and x>=n and monotonic_hostile is None:
   monotonic_hostile={'k':k,'n':n,'h':h,'first_descent':first,'endpoint':x}
  checked+=1
# Guarded lower tails: source mod4 fixes first valuation either1 or>=2.
guard_rows=[]
for h in range(1,25):
 s0=(3**h).bit_length();bad=F(0);bad_first_one=F(0)
 for s in range(h,s0):
  bad+=F(comb(s-1,h-1),2**s)
  count=(1 if s==1 else 0) if h==1 else (comb(s-2,h-2) if s>=h else 0)
  bad_first_one+=F(count,2**s)
 q_one=2*bad_first_one;q_other=2*(bad-bad_first_one)
 require(0<=q_one<=2*bad and 0<=q_other<=2*bad,'guarded tail domination')
 guard_rows.append({'h':h,'bad_given_k_3mod4':str(q_one),'bad_given_k_1mod4':str(q_other)})
require(F(11,16)>F(21,32),'endpoint densities not a CDF')
print(json.dumps({'status':'PASS','h_exact_tail_checks':200,'literal_endpoint_checks':checked,'density_rows':rows,'finite_exceptions_observed':exceptions,'endpoint_not_first_descent_witness':monotonic_hostile,'guarded_tail_rows':guard_rows,'scope':'Exponent density of finite descent, NOT convergence or uniform first-descent bound'},indent=2))
