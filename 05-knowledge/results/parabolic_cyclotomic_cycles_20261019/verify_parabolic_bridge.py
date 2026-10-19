import sympy as s
from fractions import Fraction
import json
x,t,T,c,z=s.symbols('x t T c z')
def req(ok,msg):
 if not ok:raise RuntimeError(msg)
g=t**3+t**2-2*t-1;D=2*t*(t+1)
pts=[(t**3+2*t*t+t+1)/D,(t**3-t-1)/D,-(t**3+2*t*t+3*t+1)/D]
C=-(t**6+2*t**5+4*t**4+8*t**3+9*t*t+4*t+1)/(4*t*t*(t+1)**2)
lam=8*s.prod(pts);trace=sum(pts)
for i in range(3):req(s.cancel(pts[i]**2+C-pts[(i+1)%3])==0,'cycle')
req(s.cancel(C+s.Rational(7,4)+g*g/(4*t*t*(t+1)**2))==0,'cusp square')
req(s.cancel(lam-1+(t*t+t+1)**3*g/(t**3*(t+1)**3))==0,'multiplier')
req(s.cancel(trace**2+trace+C+2)==0,'trace')
req(s.cancel(lam-8*(C*(trace+1)+1))==0,'trace multiplier')
G=x**3-T*x*x-(T*T+2*T+3)*x+T**3+2*T*T+3*T+1
req(s.factor(s.discriminant(G,x))==(4*T*T+6*T+9)**2,'discriminant')
req(s.cancel(s.prod(x-p for p in pts)-G.subs(T,trace))==0,'point polynomial')
sig=t*t-2; sig2=s.rem(sig.subs(t,sig),g,t)
req(sig2==-t*t-t+1,'inverse Galois direction')
req(s.rem((-t-s.Rational(1,2))**2-s.Rational(7,4)-(-sig2-s.Rational(1,2)),g,t)==0,'parabolic conjugacy on roots')
req(s.expand(z**3*g.subs(t,z+1/z))==sum(z**i for i in range(7)),'cyclotomic identity')
req(s.discriminant(g,t)==49,'g discriminant')
control={'t':'5/4','c':str(s.cancel(C.subs(t,s.Rational(5,4)))),'points':[str(s.cancel(p.subs(t,s.Rational(5,4)))) for p in pts],'multiplier':str(s.cancel(lam.subs(t,s.Rational(5,4))))}
req(control['c']=='-226801/129600','control c');req(control['multiplier']=='5605019/5832000','control multiplier')
req(0<s.Rational(control['multiplier'])<1,'attracting')
print(json.dumps({'status':'PASS','control':control,'scope':'Exact symbolic identities; denominator theorem separately cited'},indent=2))
