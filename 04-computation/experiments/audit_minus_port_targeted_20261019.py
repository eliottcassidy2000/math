import sympy as s
from fractions import Fraction as F

x,h=s.symbols('x h')
checks=[]
def equal(a,b,label):
    if s.cancel(a-b) != 0: raise RuntimeError(label)
    checks.append(label)
def true(q,label):
    if not q: raise RuntimeError(label)
    checks.append(label)

A=x**4+2*x**2-8*x+1
B=3*x**6-12*x**5+17*x**4-23*x**2+12*x-5
C=x**5-2*x**4+2*x**3+4*x**2+x-2
D=x**6-4*x**5+3*x**4+16*x**3-29*x**2+20*x+1
N=s.expand(C*D); E=s.expand(A*B)
wc=[1,-8,27,-112,743,-3368,8293,-8000,-11894,44912,-39122,-48288,150510,-132112,-7414,100544,-66331,24216,-26297,21392,-6469,824,1]
W=sum(c*x**(22-i) for i,c in enumerate(wc))
f=lambda t: 4*t*(t*t-1)/(t*t+1)**2
k=(-x*x+2*x+1)/(x*x+1);K=k*k
q=((x*x-4*x+1)/(x*x+1),-8*x*(x-1)/(x*x+1)**2)
u0=-x; z0=-(x*x+1)
X0=(k*z0+K+2*u0+K*u0*u0)/(2*u0*u0)
Y0=(X0-K)*(X0*u0-1)/k
equal(X0,q[0],'negative-sheet anchor X')
equal(Y0,q[1],'negative-sheet anchor Y')
equal(q[1]**2,(q[0]-K)*(q[0]**2-1),'anchor on cubic')
equal(K,1-f(x),'anchor square')
# Independently use slope/intersection/reflection formulas on y²=x³-Kx²-x+K.
lam=s.factor((3*q[0]**2-2*K*q[0]-1)/(2*q[1]))
X2=s.factor(lam**2+K-2*q[0]);Y2=s.factor(lam*(q[0]-X2)-q[1])
equal(X2,(x*x+1)*(x*x-4*x+5)/(4*(x-1)**2),'doubling X')
lam3=s.factor((Y2-q[1])/(X2-q[0]));X3=s.factor(lam3**2+K-X2-q[0]);Y3=s.factor(lam3*(q[0]-X3)-q[1])
u3=s.factor((k*Y3+X3-K)/(X3*(X3-K)))
z3=s.factor((2*u3*u3*X3-K-2*u3-K*u3*u3)/k)
equal(u3,-N/E,'3Q inverse parameter')
equal(z3,-W/E**2,'3Q inverse z')
equal((x*x+1)**2*W**2,(x*x-2*x-1)**2*(N*N+E*E)**2+4*N*E*(N*N-E*E)*(x*x+1)**2,'integer polynomial square identity')
# Reciprocal plus section, recreated from companion note, with independent coefficient list.
Ap=x**4-8*x**3+2*x**2+1
Bp=5*x**6-12*x**5+23*x**4-17*x**2+12*x-3
Cp=2*x**5-x**4-4*x**3-2*x*x+2*x-1
Dp=x**6+20*x**5-29*x**4+16*x**3+3*x*x-4*x+1
Np=x*Ap*Bp;Ep=Cp*Dp
zc=[1,824,-6469,21392,-26297,24216,-66331,100544,-7414,-132112,150510,-48288,-39122,44912,-11894,-8000,8293,-3368,743,-112,27,-8,1]
Z=sum(c*x**(22-i) for i,c in enumerate(zc))
equal(N/E,1/(Np/Ep).subs(x,1/x),'v_minus = 1/u_plus(1/t)')
equal(W,x**22*Z.subs(x,1/x),'W = reciprocal Z')
equal(x**11*Np.subs(x,1/x),-E,'reciprocal plus numerator')
equal(x**11*Ep.subs(x,1/x),-N,'reciprocal plus denominator')
equal(f(1/x),-f(x),'spread under reciprocal')

EP=s.Poly(E.subs(x,h+1),h)
LP=s.Poly(s.cancel((N-x*E)/(x-1)).subs(x,h+1),h)
GP=s.Poly((5*E-4*N).subs(x,h+1),h)
for p,coeff,label in [(EP,[32,32,-80,-128,-16,48,40,64,50,18,3],'E shifted coefficients'),(LP,[64,256,480,384,96,0,-48,-96,-60,-16,-2],'displacement coefficients'),(GP,[32,-352,-1232,-1728,-1040,-272,-152,96,178,58,-5,-4],'endpoint coefficients')]:
    true([p.nth(i) for i in range(p.degree()+1)]==coeff,label)

def neg_tail(p):return sum(F(-int(p.nth(i)),16**i) for i in range(1,p.degree()+1) if p.nth(i)<0)
true(neg_tail(EP)==F(1409,4096)<1,'E lower bound >31')
true(neg_tail(LP)==F(1777281,549755813888)<1,'displacement coefficient lower bound >63')
true(neg_tail(GP)==F(119849021866005,4398046511104)<28,'endpoint lower bound >4')
equal(s.diff(f(x),x),-4*(x**4-6*x*x+1)/(x*x+1)**3,'spread derivative')
true(f(F(5,4))==F(720,1681)<F(1,2),'spread endpoint cap')
true(f(F(17,16))==F(35904,297025)>F(1,10),'fallback cutoff')
true(2*(100*(4+5)//7)==256,'fallback signed-even endpoint')
true(len(list(range(-256,0,2)))+len(list(range(2,257,2)))==256,'fallback candidate count')
true(s.degree(N,x)==11 and s.degree(E,x)==10,'minus polynomial degrees')
true(sum(abs(c) for c in s.Poly(N,x).all_coeffs())==446,'N coefficient norm')
true(sum(abs(c) for c in s.Poly(E,x).all_coeffs())==760,'E coefficient norm')
true(446**2+760**2==776516,'cusp height constant')
print('PASS',len(checks),'targeted exact controls')
for line in checks:print('-',line)
