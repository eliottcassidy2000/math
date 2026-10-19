import sympy as s
from fractions import Fraction

t,h=s.symbols('t h')
f=lambda x:4*x*(x*x-1)/(x*x+1)**2
k=(t*t+2*t-1)/(t*t+1); K=k*k
A=t**4-8*t**3+2*t*t+1
B=5*t**6-12*t**5+23*t**4-17*t*t+12*t-3
C=2*t**5-t**4-4*t**3-2*t*t+2*t-1
D=t**6+20*t**5-29*t**4+16*t**3+3*t*t-4*t+1
N=s.expand(t*A*B); E=s.expand(C*D)
Z=s.Poly.from_list([1,824,-6469,21392,-26297,24216,-66331,100544,-7414,-132112,150510,-48288,-39122,44912,-11894,-8000,8293,-3368,743,-112,27,-8,1],t).as_expr()
passed=[]
def ck(name,truth):
    if not truth: raise RuntimeError(name)
    passed.append(name)
def eq(name,a,b): ck(name,s.cancel(a-b)==0)
def add(p,q):
    x,y=p; X,Y=q
    if p==q: slope=s.cancel((3*x*x-2*K*x-1)/(2*y))
    else: slope=s.cancel((Y-y)/(X-x))
    xx=s.factor(slope*slope+K-x-X)
    return xx,s.factor(-y+slope*(x-xx))
P=(s.Integer(0),k)
Q=((t*t+4*t+1)/(t*t+1),8*t*(t+1)/(t*t+1)**2)
for label,(x,y) in [('P',P),('Q',Q)]:eq(label+' on curve',y*y,(x-K)*(x*x-1))
eq('k square',k*k,1+f(t))
eq('swap f',f((t+1)/(t-1)),f(t))
eq('swap k',k.subs(t,(t+1)/(t-1)),k)
eq('derivative',s.diff(f(t),t),-4*(t**4-6*t*t+1)/(t*t+1)**3)
V=add(Q,P)
eq('Vx',V[0],(t*t-4*t+1)/(t*t+1))
eq('Vy',V[1],-8*t*t*(t-1)/(t*t+1)**2)
S=add(V,V)
eq('Sx',S[0],(t*t+1)*(5*t*t-4*t+1)/(4*t*t*(t-1)**2))
sy=-(t*t+2*t-1)*(3*t*t-2*t+1)*(t**4-6*t**3+2*t-1)/(8*t**3*(t-1)**3*(t*t+1))
eq('Sy',S[1],sy)
R=add(Q,S); x,y=R
u=s.factor((k*y+x-K)/(x*(x-K)))
z=s.factor((2*u*u*x-K-2*u-K*u*u)/k)
eq('R inverse u',u,N/E)
eq('R inverse z',z,Z/E**2)
eq('polynomial square',(t*t+1)**2*Z**2,(t*t+2*t-1)**2*(N*N+E*E)**2-4*N*E*(N*N-E*E)*(t*t+1)**2)
shift=lambda p:s.Poly(p.subs(t,1+h),h)
eq('shift A',A.subs(t,1+h),h**4-4*h**3-16*h*h-16*h-4)
eq('shift B',B.subs(t,1+h),5*h**6+18*h**5+38*h**4+72*h**3+76*h*h+40*h+8)
eq('shift C',C.subs(t,1+h),2*h**5+9*h**4+12*h**3-8*h-4)
eq('shift D',D.subs(t,1+h),h**6+26*h**5+86*h**4+120*h**3+92*h*h+40*h+8)
ck('C tail',2*s.Rational(1,16)**5+9*s.Rational(1,16)**4+12*s.Rational(1,16)**3==s.Rational(1609,524288))
L=s.Poly(s.cancel(-(N-t*E)/(t-1)).subs(t,1+h),h)
ck('full displacement coefficients',[L.nth(i) for i in range(12)]==[64,448,1440,2976,4704,6048,6000,4176,1860,476,58,2])
G=shift(4*N-5*E)
ck('full endpoint coefficients',[G.nth(i) for i in range(12)]==[32,-160,-2000,-7392,-15632,-22512,-23576,-17744,-8878,-2554,-293,10])
ck('exact endpoint negative-tail bound',sum(-G.nth(i)*s.Rational(1,16)**i for i in range(1,11))==s.Rational(21856812814021,1099511627776)<20)
ck('f cutoff',f(s.Rational(17,16))==s.Rational(35904,297025)>s.Rational(1,10))
ck('f cap',f(s.Rational(5,4))==s.Rational(720,1681)<s.Rational(1,2))
ck('N degree and norm',s.degree(N,t)==11 and sum(abs(c) for c in s.Poly(N,t).all_coeffs())==760)
ck('E degree and norm',s.degree(E,t)==11 and sum(abs(c) for c in s.Poly(E,t).all_coeffs())==446)
ck('coprime N E',s.degree(s.gcd(N,E),t)==0)
ck('height coefficient',760**2+446**2==776516)
ck('budget',2*(100*(4+5)//7)==256 and len([q for q in range(-256,257,2) if q])==256)
ck('S double pole',s.limit((t-1)**2*S[0],t,1)==1)
ck('S triple pole',s.limit((t-1)**3*S[1],t,1)==1)
ck('formal parameter leading',s.limit((-S[0]/S[1])/(t-1),t,1)==-1)
ck('P Q V reductions',tuple(c.subs(t,1) for c in P)==(0,1) and tuple(c.subs(t,1) for c in Q)==(3,4) and tuple(c.subs(t,1) for c in V)==(-1,0))
zeta=lambda w:s.simplify((w-s.sqrt(2))/(w+s.sqrt(2)))
eq('zeta P',zeta(-1),-(1+s.sqrt(2))**2)
eq('zeta Q',zeta(2),(1+s.sqrt(2))**-2)
eq('zeta T',zeta(0),-1)
print('\n'.join(passed)); print('PASS:',len(passed),'independent exact controls')
