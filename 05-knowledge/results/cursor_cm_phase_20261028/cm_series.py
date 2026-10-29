"""Exact eta-series, reduced-form and character controls; certified exponential intervals."""
from pathlib import Path
from math import comb,gcd,isqrt
from fractions import Fraction as F
from decimal import Decimal as D, localcontext, ROUND_FLOOR, ROUND_CEILING
import json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,msg):
 global CHECKS
 CHECKS+=1
 if not b:raise RuntimeError(msg)
def mul(a,b,N):
 out=[0]*(N+1)
 for i,x in enumerate(a):
  for j,y in enumerate(b[:N+1-i]):out[i+j]+=x*y
 return out
def product_odd(N,sign):
 a=[1]+[0]*N
 for j in range(1,N+1,2):
  b=[0]*(N+1)
  for k in range(min(24,N//j)+1):b[j*k]=comb(24,k)*sign**k
  a=mul(a,b,N)
 return a
def reciprocal(a):
 b=[1]
 for n in range(1,len(a)):b.append(-sum(a[j]*b[n-j] for j in range(1,n+1)))
 return b
def forms(Dsc):
 out=[]
 for a in range(1,isqrt(abs(Dsc)//3)+2):
  for b in range(-a,a+1):
   if (b*b-Dsc)%(4*a):continue
   c=(b*b-Dsc)//(4*a)
   if a>c or gcd(gcd(a,abs(b)),c)!=1 or ((abs(b)==a or a==c) and b<0):continue
   out.append((a,b,c))
 return out
def kronecker(a,n):
 if n==0:return int(abs(a)==1)
 value=1
 while n%2==0:
  n//=2
  if a%2==0:return 0
  if a%8 in (3,5):value=-value
 a%=n
 while a:
  while a%2==0:
   a//=2
   if n%8 in (3,5):value=-value
  a,n=n,a
  if a%4==n%4==3:value=-value
  a%=n
 return value if n==1 else 0
def trace(t,k):
 a,b=2,t
 if k==0:return a
 for _ in range(1,k):a,b=b,t*b+a
 return b

def quadratic_mul(x,y,d):return x[0]*y[0]+d*x[1]*y[1],x[0]*y[1]+x[1]*y[0]

def mobius(n):
 p=2;out=1
 while p*p<=n:
  if n%p==0:
   n//=p;out=-out
   if n%p==0:return 0
  p+=1
 return -out if n>1 else out

def product_value(n,sign):
 # A numerical cross-check of cited identities, not an algebraic proof.
 lo,hi=exp_pi_sqrt(n);q=2/(lo+hi)
 value=D(2)**(-D(1)/4)*q**(-D(1)/24)
 for j in range(1,101,2):value*=1+sign*q**j
 return value
def atan_bounds(q,N=110):
 value=sum((F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N)),F())
 error=F((-1)**N,(2*N+1)*q**(2*N+1))
 return min(value,value+error),max(value,value+error)
def exp_pi_sqrt(n):
 a,b=atan_bounds(5);c,d=atan_bounds(239)
 lo,hi=16*a-4*d,16*b-4*c
 with localcontext() as ctx:
  ctx.prec=90;ctx.rounding=ROUND_FLOOR;pl=D(lo.numerator)/D(lo.denominator)
  ctx.rounding=ROUND_CEILING;pu=D(hi.numerator)/D(hi.denominator)
  s=D(n).sqrt();sl=s.next_minus();su=s.next_plus()
  ctx.rounding=ROUND_FLOOR;xlo=pl*sl
  ctx.rounding=ROUND_CEILING;xhi=pu*su
  # Decimal sqrt and exp are correctly rounded; one neighboring representable
  # number on each side encloses that rounding as well as input uncertainty.
  return xlo.exp().next_minus(),xhi.exp().next_plus()
def main():
 N=32;minus=product_odd(N,-1);plus=product_odd(N,1);inv=reciprocal(minus)
 need(plus==[(-1)**i*x for i,x in enumerate(minus)],'odd product sign substitution')
 # q*T_2A = P_minus + 24q +4096q^2/P_minus.
 moon=minus[:];moon[1]+=24
 for i in range(2,N+1):moon[i]+=4096*inv[i-2]
 need(moon[:8]==[1,0,4372,96256,1240002,10698752,74428120,431529984],'known 2A prefix')
 # Independent eta quotient: prod(1+q^j)^(-24).
 eta=[1]+[0]*N
 for j in range(1,N+1):
  f=[0]*(N+1)
  for k in range(N//j+1):f[k*j]=(-1)**k*comb(23+k,k)
  eta=mul(eta,f,N)
 need(eta==minus,'odd-product versus all-product eta quotient')
 doubled=[minus[i//2] if i%2==0 else 0 for i in range(N+1)]
 need(mul(minus,plus,N)==doubled,'sign-pair product is nome doubling')
 qderivative=[i*x for i,x in enumerate(minus)]
 logcoeff=[-x for x in mul(qderivative,inv,N)]
 for n in range(1,N+1):
  ds=[d for d in range(1,n+1) if n%d==0]
  need(logcoeff[n]==24*sum(d for d in ds if d%2),'logarithmic derivative retains odd divisor weights')
  need(sum(mobius(n//d)*logcoeff[d] for d in ds)==n*(24 if n%2 else 0),'Mobius recovers eta exponents')
 classrows=[]
 for n in (23,29,31,47,71,79,65,37):
  Ds=(-n if n%4==3 else -4*n,-8*n);counts=[len(forms(x)) for x in Ds]
  classrows.append({'n':n,'discriminants':Ds,'class_numbers':counts,'forms':[forms(x) for x in Ds]})
 need([r['class_numbers'] for r in classrows]==[[3,4],[6,2],[3,8],[5,8],[7,4],[5,8],[8,4],[2,10]],'CM parameter doubling table')
 ranking=[]
 for p,h in ((23,3),(31,3),(47,5),(71,7),(79,5)):
  need(len(forms(-p))==h,'fundamental Paley class number')
  excess=sum((p-d)*kronecker(-p,d) for d in range(1,p))
  need(excess==p*h,'Paley ranking excess equals p*h')
  ranking.append({'p':p,'h':h,'excess':excess})
 near=[]
 for n,t,power,offset,expected in ((22,2,12,24,2508952),(37,12,6,-24,199148648),(58,5,12,24,24591257752)):
  integer=64*trace(t,power)+offset
  need(integer==expected,'exact quadratic-unit trace integer')
  lo,hi=exp_pi_sqrt(n)
  with localcontext() as ctx:
   ctx.prec=90;ctx.rounding=ROUND_FLOOR;lower=D(integer)-hi
   ctx.rounding=ROUND_CEILING;upper=D(integer)-lo
   need(0<lower<=upper<1,'certified below-integer exponential')
   near.append({'n':n,'trace_parameter':t,'unit_exponent':power,'integer':integer,
                'exp_interval':[str(lo),str(hi)],'deficit_interval':[str(lower),str(upper)]})
 # Raw algebraic identities of the trace forget u versus -1/u (even power).
 need(trace(2,12)==trace(-2,12) and trace(12,6)==trace(-12,6),'even trace loses conjugate branch')
 # D_12(x) with x=6+sqrt(65), D_0=2,D_1=x,D_j=xD_(j-1)-D_(j-2).
 a,b=(2,0),(6,1)
 for j in range(1,12):
  z=quadratic_mul((6,1),b,65);a,b=b,(z[0]-a[0],z[1]-a[1])
 need(b==(28123794942082,3488327441280),'130 reciprocal trace stays quadratic irrational')
 lo130,hi130=exp_pi_sqrt(130)
 with localcontext() as ctx:
  ctx.prec=90
  frac=(lo130-D(int(lo130)),hi130-D(int(lo130)))
  need(D('.57')<frac[0]<=frac[1]<D('.58'),'130 simple-surds hostile is not near an integer')
  s5,s13,s65=D(5).sqrt(),D(13).sqrt(),D(65).sqrt()
  formulas=[(130,-1,((2+s5)*(3+s13)/2).sqrt()),
            (65,1,(((1+s5)/2)*((3+s13)/2))**(D(1)/4)*(((9+s65)/8).sqrt()+((1+s65)/8).sqrt()).sqrt()),
            (22,-1,(1+D(2).sqrt()).sqrt()),(37,1,(6+D(37).sqrt())**(D(1)/4)),
            (58,-1,((5+D(29).sqrt())/2).sqrt())]
  numeric=[]
  for n,sign,value in formulas:
   rel=abs(product_value(n,sign)/value-1)
   need(rel<D('1e-75'),'cited Ramanujan radical matches defining product')
   numeric.append({'n':n,'sign':sign,'relative_error':str(rel)})
 out={'status':'FINITE-EXACT series/forms/traces; VERIFIED certified-decimal exponential intervals and numerical radical cross-checks',
      'q_times_2A_coefficients':moon,'CM_doubling':classrows,'Paley_rankings':ranking,'near_integers':near,
      '130_hostile':{'trace_coefficients':b,'exp_interval':[str(lo130),str(hi130)],'fractional_interval':[str(x) for x in frac]},
      'radical_cross_checks':numeric,'checks':CHECKS}
 (HERE/'cm_series.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'checks':CHECKS,'CM_class_numbers':[r['class_numbers'] for r in classrows],
                   'near_integers':[{k:v for k,v in r.items() if k in ('n','integer','deficit_interval')} for r in near]},indent=2))
if __name__=='__main__':main()
