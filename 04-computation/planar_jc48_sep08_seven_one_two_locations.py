#!/usr/bin/env python3
"""Exact controls for two complete (7,1) DG boundary placements.

All original finite positions are symbolic. Rational coordinate maps
retain the actual volume. The theorem is polynomial, with sharp rational
constant-L controls; the finite1/infinity7 placement is separate.
"""
from hashlib import sha256
import sympy as S
GATES=[]
def need(label,pred):
    if not bool(pred):raise RuntimeError(label)
    GATES.append(label)
def zero(label,val):need(label,S.cancel(val)==0)
def jac(f,g,x,y):return S.diff(f,x)*S.diff(g,y)-S.diff(f,y)*S.diff(g,x)
u,t,p,a,b,d,lam,mu,e,v,z,h,f,A,k=S.symbols('u t p a b d lam mu e v z h f A k')

def global_h(label,N,P,Q):
    for j,val in enumerate([Q,P-2*Q*(u+p)**2,N-P*(u+p)**2+Q*(u+p)**4]):
        need(label+str(j),S.Poly(S.expand(val),u).degree()<=4)
def global_l(label,M,R):
    for j,val in enumerate([R,M-R*(u+p)**2]):
        need(label+str(j),S.Poly(S.expand(val),u).degree()<=2)

# The complete post-first-jet source rows, with all original p retained.
N=u**7*(u-1)
P=2*u**6-(4*p+2)*u**5+a*u**4+b*u**3
Q=u**4-(4*p+1)*u**3+(a+4*p*p)*u*u+(b-2*p*a+4*p*p)*u+d
M=lam*u**3*(u-1);R=lam*(u*u-(1+2*p)*u)+e
global_h('all-finite H',N,P,Q);global_l('all-finite L',M,R)
H=N*t*t+P*t+Q;L=M*t+R
vsource=u+u**3*t-2*p
tsub=(v+2*p-u)/u**3
zero('all-finite exact H carrier',H.subs(t,tsub)-
     (u*u*v*v+((a-4*p)-v)*u*v+b*v+d+2*b*p))
zero('all-finite exact L carrier',L.subs(t,tsub)-(lam*(u-1)*v+e-2*p*lam))
zero('all-finite polynomial first jet',S.diff(H*H+L,u).subs(u,0).subs(lam,2*b*d)+4*p*d*(a+b-2*p))
zero('original v Jacobian',jac(u,vsource,u,t)-u**3)
zero('z-v volume transport',jac(u*v,v,u,v)-v)

# Exact rational source field in (z,v), then (F,H).
h0=z*z+A*z+k
HH=h0+v*(b-z);FF=HH*HH+lam*(z-v)
den=h*h-f-lam*(A+b)
q=b*(h*h-f)+lam*(k-h)
zz=q/den
ww=(h*h-f)**2-lam*A*(h*h-f)+lam*lam*(k-h)
vv=ww/(lam*den)
zero('rational inverse H',HH.subs({z:zz,v:vv},simultaneous=True)-h)
zero('rational inverse F',FF.subs({z:zz,v:vv},simultaneous=True)-f)
zero('rational inverse v identity',vv-zz-(h*h-f)/lam)
zero('actual volume in f-h field',
     vv**2/zz**3*jac(zz,vv,f,h)+ww*ww/(lam*lam*q**3))
zero('actual dominance derivative',
     S.diff(h*h+lam*((A+b)*z+k-h)/(b-z),z)-lam*(k+b*(A+b)-h)/(b-z)**2)

# Generic-residue expansion. Keep the A term in B0; omitting it is false.
Ph=h*h-lam*h/b+lam*k/b
w=S.symbols('w')
B1=2*(h-k)/b-A
B0=(h-k)**2/b**2-A*(h-k)/b-(h-k)
zero('entire residue expansion',
     ww.subs(f,Ph-w)-(w*w+lam*B1*w+lam*lam*B0))
zero('proper and polynomial differential pieces',
     (w*w+lam*B1*w+lam*lam*B0)**2/w**3-
     (w+2*lam*B1+lam*lam*(B1*B1+2*B0)/w+
      2*lam**3*B1*B0/w**2+lam**4*B0*B0/w**3))
def DP(expr,P,var):return S.diff(expr,var)/S.diff(P,var)
res=(B1*B1+2*B0)/S.diff(Ph,h)+DP(2*lam*B1*B0/S.diff(Ph,h),Ph,h)+\
    DP(DP(lam*lam*B0*B0/S.diff(Ph,h),Ph,h),Ph,h)/2
zero('unavoidable all-finite residue slope',S.limit(res/h,h,S.oo)-3/b**2)

# Finite seven, original infinity one: full source rows after the gate.
N7=u**7;P7=2*u**5+a*u**4+b*u**3
Q7=u**3+a*u*u+(b-2*p*a-4*p*p)*u+d
M7=u**3*(lam*u+mu);R7=lam*u*u+(mu-2*p*lam)*u+e
global_h('finite-seven H',N7,P7,Q7);global_l('finite-seven L',M7,R7)
H7=N7*t*t+P7*t+Q7;L7=M7*t+R7
zero('finite-seven exact H carrier',H7.subs(t,tsub)-
     (u*v*(v+a+4*p)+b*v+d+2*b*p))
zero('finite-seven exact L carrier',L7.subs(t,tsub)-
     (lam*u*v+mu*v+e+2*p*mu))
zero('finite-seven polynomial first jet',
     S.diff(H7*H7+L7,u).subs(u,0).subs(mu,-2*b*d)+2*p*(lam+2*d*(a+2*p)))

# The literal quadratic field over C(f,v), including lambda=0.
hzero=b*v+k;delta=lam/(v+A);Pbase=hzero*hzero+mu*v+e
uf=(h-hzero)/(v*(v+A))
zero('second-location actual volume',
     uf**(-3)*S.diff(uf,h)-v*v*(v+A)**2/(h-hzero)**3)
disc=delta*delta+4*(f+delta*hzero-mu*v-e)
zero('genuine quadratic discriminant f coefficient',S.diff(disc,f)-4)
# Universal full trace, reconstructed by a companion matrix.
BB,CC=S.symbols('BB CC')
Zmat=S.Matrix([[0,-CC],[1,-BB]])
fulltrace=S.trace((Zmat**3*(2*Zmat+BB*S.eye(2))).inv())
zero('actual full trace sign and factor',fulltrace-(1/CC**2-BB*BB/CC**3))
weight=v*v*(v+A)**2
L1=weight*(2*hzero+delta)**2/S.diff(Pbase,v)
zero('second-location nonzero degree mismatch',
     S.limit((S.diff(L1,v)-2*weight)/v**4,v,S.oo)-8)

# b=0 forces mu=0 by the polynomial first jet, so lambda!=0.
z0=(f-h*h)/lam;v0=(h-k)/z0-A
zero('b0 actual inverse H',z0*(v0+A)+k-h)
zero('b0 actual inverse F',(z0*(v0+A)+k)**2+lam*z0-f)
eta=(h-k-A*z0)**2/(lam*z0**6)
zero('b0 original volume transport',v0*v0/z0**3*jac(z0,v0,f,h)-eta)
rho=S.symbols('rho',nonzero=True)
proper=(lam*(h-k)-A*(rho*rho-h*h))**2/(rho*rho-h*h)**6
expected=-(80*A*A*rho**4+140*A*k*lam*rho*rho+63*k*k*lam*lam-7*lam*lam*rho*rho)/(512*rho**11)
zero('complete b0 generic residue',S.residue(proper,h,rho)-expected)
need('b0 residue first forces A zero',S.expand(-512*rho**11*expected).coeff(rho,4)==80*A*A)
need('b0 residue then contradicts lambda nonzero',S.expand((-512*rho**11*expected).subs(A,0)).coeff(rho,2)==-7*lam*lam)

# Actual rational sharpness in both exact partitions, every finite p.
Va=u+u**3*t-2*p
Ha=u*(u-1)*Va**2+d
Ga=-(8*u*u+4*u+3)/(15*u**3*Va)
global_h('all-finite sharp H',u**7*(u-1),2*u**4*(u-1)*(u-2*p),u*(u-1)*(u-2*p)**2+d)
zero('all-finite sharp rational H mate',jac(Ha,Ga,u,t)-1)
zero('all-finite sharp rational F mate',jac(Ha*Ha+e,Ga/(2*Ha),u,t)-1)
Hb=u*Va**2+d;Gb=1/(5*u**3*Va)
global_h('finite-seven sharp H',u**7,2*u**4*(u-2*p),u*(u-2*p)**2+d)
zero('finite-seven sharp rational H mate',jac(Hb,Gb,u,t)-1)
zero('finite-seven sharp rational F mate',jac(Hb*Hb+e,Gb/(2*Hb),u,t)-1)
print('Seven-one two locations: complete all-finite and finite7/infinity1 polynomial exclusions')
print('Actual fields and volumes; generic residues; all-p rational constant-L hostiles retained')
print('Finite1/infinity7 is a separate OPEN location')
print('Exact gates:',len(GATES))
print('Semantic SHA256:',sha256('\n'.join(GATES).encode()).hexdigest())
