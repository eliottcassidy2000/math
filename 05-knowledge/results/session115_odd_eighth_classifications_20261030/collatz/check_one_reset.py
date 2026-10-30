"""Exact controls for the one-reset classification; standard library only."""
from fractions import Fraction as F
from math import comb
import json
INF=float('inf')
checks={}
def gate(ok,label):
 checks[label]=checks.get(label,0)+1
 if not ok:raise RuntimeError(label)
def v5(x):
 if x==0:return INF
 k=0
 while x%5==0:x//=5;k+=1
 return k
def coeff(w):
 P=Q=1;B=0
 for a in w:P,Q,B=3*P,2**a*Q,3*B+Q
 return P,Q,B
def params(W):
 P,Q,B=W;m=P*pow(Q,-1,5)%5
 t=next(t for t in (1,2,4) if pow(m,t,5)==1)
 r=min(v5(P-Q),v5(B)) if t==1 else v5(P**t-Q**t)
 return t,r
def depth(W,H):
 P,Q,B=W;R,S,C=H;t,r=params(W)
 c=v5((P-Q)*C-B*(R-S))
 if t>1:
  if pow(R*pow(S,-1,5)%5,t,5)!=1:return 0
  z=v5(R**t-S**t)
  return min(c,z if z<r else INF)
 s,h=v5(P-Q),v5(B)
 a=v5(R-S) if s<=h else v5(C)
 return max(0,min(c-r,a if a<r else INF))
def mass(t,r,D):
 if D==INF:return F(4*r+1,4*t)
 if D<=r:return F(D,t)
 return F(r,t)+F(1,4*t)*(1-F(1,5**(D-r)))
def reference_rate(w,u):
 W=coeff(w);t,r=params(W)
 return sum((mass(t,r,depth(W,coeff(w[:j]+u))) for j in range(1,len(w)+1)),F(0))
def words(A,p=()):
 yield p
 for a in range(1,A+1):yield from words(A-a,p+(a,))
def formula(L):
 r=v5(9**(L+1)+16*4**L);z=F(4*r+1,16);h=2;e=1
 while h<=L:
  z+=F(L//h,4*5**max(0,e-r));h*=5;e+=1
 return z
def phase_good(L,j):
 h=2
 while h<=L:
  if j%h>L%h:return False
  h*=5
 return True
def maximizing_entry(L,u):
 w=(1,)*L+(2,)
 for j in range(L+1):
  full=(1,)*j+u
  if len(full)%len(w)==0 and full==w*(len(full)//len(w)):
   return phase_good(L,j)
 return False
def modpair(w,q):
 P,Q,B=coeff(w);inv=pow(Q,-1,q);return P*inv%q,B*inv%q
def orbit(pair,q):
 m,b=pair;x,y=1,0;out=set()
 while (x,y) not in out:out.add((x,y));x,y=m*x%q,(m*y+b)%q
 gate((x,y)==(1,0),'cyclic_closure')
 return out
def literal_layer_mean(L,u,E):
 w=(1,)*L+(2,);q=5**E;M=orbit(modpair(w,q),q);T=len(M)
 m,b=modpair(u,q);total=0
 for k in range(T):
  for a in w:
   inv=pow(2**a,-1,q);m,b=3*m*inv%q,(3*b+1)*inv%q
   total+=min(v5(m-1),v5(b),E)
 return F(total,T)
def run():
 rows=[];best=(F(0),0);a=9;b=16
 for L in range(5001):
  r=v5(a+b);z=formula(L)
  gate(8*r<=16+11*L,'elementary_stall_bound')
  gate(z<=F(L,2)+F(9,16),'uniform_mean_bound')
  gate(z/F(L+2)<F(1,2),'strict_clock_half')
  if z/F(L+2)>best[0]:best=z/F(L+2),L
  if L<=25:rows.append({'L':L,'stall_depth':r,'max_mean':str(z),'max_clock_rate':str(z/F(L+2))})
  a*=9;b*=4
 for L in range(41):
  w=(1,)*L+(2,)
  gate(reference_rate(w,())==formula(L),'general_rate_formula')
  for e in range(1,5):
   q=5**e;A=orbit(modpair((1,),q),q);M=orbit(modpair(w,q),q)
   gate(A&M=={(1,0)},'trivial_intersection')
   gate(len(A)==2*5**(e-1),'all_one_order')
 for L in range(13):
  w=(1,)*L+(2,);mx=formula(L)
  for u in words(9):
   z=reference_rate(w,u)
   gate(z<=mx,'all_small_entries_bounded')
   gate((z==mx)==maximizing_entry(L,u),'all_small_entry_equality')
 for L in range(500):
  K,d=divmod(L,2);count=0
  for j in range(L+1):
   a,e=divmod(j,2);ok=phase_good(L,j)
   gate(ok==(e<=d and comb(K,a)%5!=0),'lucas_equivalence')
   count+=ok
  product=1;k=K
  while k:product*=k%5+1;k//=5
  gate(count==(d+1)*product,'lucas_count')
 for L in range(9):
  for E in (1,2,3):
   r=v5(9**(L+1)+16*4**L)
   exact=sum((F(1+L//(2*5**(e-1)),4*5**max(0,e-r)) for e in range(1,E+1)),F(0))
   gate(literal_layer_mean(L,(),E)==exact,'literal_empty_mean')
   for j in range(L+1):
    u=(1,)*(L-j)+(2,)
    expected=sum((F(1+j//(2*5**(e-1))+(L-j)//(2*5**(e-1)),4*5**max(0,e-r)) for e in range(1,E+1)),F(0))
    gate(literal_layer_mean(L,u,E)==expected,'literal_completed_entry_mean')
 # Unique lifting of the order-stall address, not a natural-average claim.
 ell=0;mod=1;lifts=[]
 for e in range(1,9):
  q=5**e
  candidates=[x for x in range(mod) if (9*pow(9*pow(4,-1,q)%q,x,q)+16)%q==0] if e<=5 else [ell+c*(mod//5) for c in range(5)]
  if e>5:candidates=[x for x in candidates if (9*pow(9*pow(4,-1,q)%q,x,q)+16)%q==0]
  gate(len(candidates)==1,'unique_stall_lift')
  ell=candidates[0];lifts.append({'depth':e,'modulus':mod,'residue':ell});mod*=5
 return {'status':'PASS','gates':checks,'total_gates':sum(checks.values()),'L_max_rate_scan':5000,'maximum_scanned_clock_rate':str(best[0]),'maximizer_L':best[1],'examples':rows,'stall_lifts':lifts}
if __name__=='__main__':print(json.dumps(run(),indent=2))
