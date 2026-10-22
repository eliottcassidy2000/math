import sympy as s,json
from pathlib import Path
E2,E4,E6,x=s.symbols('E2 E4 E6 x', real=True)
Delta=(E4**3-E6**2)/1728
F=E4**11-2640*Delta*E4**8+1813680*Delta**2*E4**5-244992000*Delta**3*E4**2
D=lambda p:s.expand(s.diff(p,E2)*(E2**2-E4)/12+s.diff(p,E4)*(E2*E4-E6)/3+s.diff(p,E6)*(E2*E6-E4**2)/2)
N=D(D(F))-45/x*D(F)
r=N/(990*Delta**4)
# x=2*pi*t; at t=1 E2=6/x, E6=0.
r1=s.factor(r.subs({E2:6/x,E6:0})); print(r1)
# Ordinary x derivative along imaginary axis = -D on modular functions.
T=lambda p:s.diff(p,x)-D(p)
# D supports rational functions too, but expands huge intermediate expressions.
# Use quotient numerator: r'= [N' + 4*E2*N]/(990 Delta^4), because DDelta=E2Delta.
Np=s.diff(N,x)-D(N)+4*E2*N
Npp=s.diff(Np,x)-D(Np)+4*E2*Np
assert s.factor(Np.subs({E2:6/x,E6:0})) == 0
curv=s.factor((x*x*Npp/(990*Delta**4)).subs({E2:6/x,E6:0}))
print('logt curvature at1:',curv)
out={'r1':str(r1),'logt_curvature_at1':str(curv),'x_at1':'2*pi','E4_at_i':'3*Gamma(1/4)^8/(64*pi^6)'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
from fractions import Fraction as Q
from decimal import Decimal,localcontext

def atan_bounds(z,N):
 p=sum(((-1)**k*z**(2*k+1)/ (2*k+1) for k in range(N)),Q(0));term=(-1)**N*z**(2*N+1)/(2*N+1)
 return min(p,p+term),max(p,p+term)
a,b=atan_bounds(Q(1,5),50);c,d=atan_bounds(Q(1,239),20)
pl,pu=16*a-4*d,16*b-4*c

def exp_bounds(z,N=100):
 term=Q(1);sm=term
 for k in range(1,N+1):term*=z/k;sm+=term
 nxt=term*z/(N+1)
 return sm,sm+nxt/(1-z/(N+2))
el,eu=exp_bounds(2*pl);fl,fu=exp_bounds(2*pu)
ql,qu=1/fu,1/el
sigma=lambda n:sum(k**3 for k in range(1,n+1) if n%k==0)
El=1+240*sum((sigma(k)*ql**k for k in range(1,13)),Q(0))
Eu=1+240*sum((sigma(k)*qu**k for k in range(1,13)),Q(0))+240*13**4*qu**13/(1-16*qu)
yl,yu=4*pl**2*El,4*pu**2*Eu
poly=lambda y:23708225*y*y-1361323584*y-69959376
assert yl>57 and yu<58
assert 2*23708225*yl-1361323584>0
assert poly(yu)<0
fmt=lambda v: format(Decimal(v.numerator)/Decimal(v.denominator),'.24f')
with localcontext() as ctx:
 ctx.prec=40
 out['rational_certificate']={'method':'Machin alternating arctangents; positive exponential Taylor with geometric remainder; E4 first12 terms with sigma3(n)<=n^4 tail','pi_bounds':[fmt(pl),fmt(pu)],'y_bounds_rounded_display':[fmt(yl),fmt(yu)],'y_rational_enclosure_denominator':10**18,'y_rational_enclosure_numerators':[str((yl*10**18).__floor__()),str((yu*10**18).__ceil__())],'curvature_polynomial_upper':fmt(poly(yu)),'negative_curvature_proved':True}
import mpmath as m
m.mp.dps=60
rf=s.lambdify((E2,E4,E6,x),r,'mpmath')
rows=[]
for ts in ['0.99','0.999','1','1.001','1.01','1.05','1.1']:
 t=m.mpf(ts);q=m.exp(-2*m.pi*t)
 sig=lambda k,p:sum(j**p for j in range(1,k+1) if k%j==0)
 e2=1-24*sum(sig(k,1)*q**k for k in range(1,101))
 e4=1+240*sum(sig(k,3)*q**k for k in range(1,101))
 e6=1-504*sum(sig(k,5)*q**k for k in range(1,101))
 rows.append({'t':ts,'r_numerical':m.nstr(rf(e2,e4,e6,2*m.pi*t),30)})
out['nearby_numerical_samples']=rows
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
