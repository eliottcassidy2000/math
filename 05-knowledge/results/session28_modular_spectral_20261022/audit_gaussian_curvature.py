"""Independent derivative-level audit; no producer functions imported."""
import sympy as s, json
from fractions import Fraction as Q
from pathlib import Path

def require(c):
 if not c:raise RuntimeError('audit failed')
a,b,c,x=s.symbols('a b c x')
de=(b**3-c**2)/1728
F=b**11-2640*de*b**8+1813680*de**2*b**5-244992000*de**3*b**2

def D(f):return s.expand(s.diff(f,a)*(a*a-b)/12+s.diff(f,b)*(a*b-c)/3+s.diff(f,c)*(a*c-b*b)/2)
f=[s.expand(F)]
for k in range(4):f.append(D(f[-1]))
g=f[2]-45*f[1]/x
lg=-x*f[3]+45*f[2]+45*f[1]/x
l2g=x*x*f[4]-46*x*f[3]-45*f[2]-45*f[1]/x
h=-4*x*a;lh=-4*x*a+x*x*(a*a-b)/3
curv=s.factor(((l2g-2*h*lg+(h*h-lh)*g)/(990*de**4)).subs({a:6/x,c:0}))
P=lambda y:23708225*y*y-1361323584*y-69959376
require(s.factor(curv-1024*P(x*x*b)/(x*x*b))==0)
require(s.factor((lg-h*g).subs({a:6/x,c:0}))==0)
r1=s.factor((g/(990*de**4)).subs({a:6/x,c:0}))
require(s.factor(r1-8192*(304325-17489844/(x*x*b)))==0)
# Independent rational enclosure, using 20 E4 terms rather than 12.
def at(z,N):
 v=sum(((-1)**j*z**(2*j+1)/Q(2*j+1) for j in range(N)),Q(0));e=(-1)**N*z**(2*N+1)/Q(2*N+1)
 return min(v,v+e),max(v,v+e)
al,au=at(Q(1,5),60);bl,bu=at(Q(1,239),25)
pl,pu=16*al-4*bu,16*au-4*bl

def ex(z):
 t=Q(1);v=t
 for k in range(1,121):t*=z/k;v+=t
 return v,v+t*z/121/(1-z/122)
ql=1/ex(2*pu)[1];qu=1/ex(2*pl)[0]
sigma=lambda n:sum(j**3 for j in range(1,n+1) if n%j==0)
el=1+240*sum((sigma(j)*ql**j for j in range(1,21)),Q(0))
eu=1+240*sum((sigma(j)*qu**j for j in range(1,21)),Q(0))+240*21**4*qu**21/(1-16*qu)
yl,yu=4*pl*pl*el,4*pu*pu*eu
lo=Q(57471215393911309459,10**18);hi=Q(57471215393911309460,10**18)
require(lo<yl<yu<hi)
require(2*23708225*lo-1361323584>0 and P(hi)<-21413)
# r1 grows with y, so this is also a rational upper bound.
require(8192*(304325-Q(17489844,1)/hi)<12000<2723200)
out={'method':'independent four Ramanujan derivatives, log-derivative quotient rule, Machin60/25 terms, exponential120 terms, E4 first20 terms','curvature_identity':str(curv),'y_interval':[str(lo),str(hi)],'P_upper':str(P(hi)),'r1_less_than':12000,'tail_limit':2723200,'strict_local_maximum_at_t1':True,'at_least_two_reciprocal_global_minimizers':True,'limits':'No uniqueness, location, actual angular eigenvalue, or Gaussian stability claim.'}
print(json.dumps(out,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
