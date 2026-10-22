import sympy as s
u,z=s.symbols('u z');n=1/u
a=z*z*(n+1)/n**2;c=z*z/n**2
g0=a/(1-a);g1=c/(1-c);h0=-a/(1+a);h1=-c/(1+c);o0=(z/n)/(1+a)
gd=g0+2*(g1-g0)/n;hd=h0+2*(h1-h0)/n
gs=2*(g1-g0)/n;hs=2*(h1-h0)/n
x1=n*(2*gd*hd-s.Rational(3,2)*gd**2)+n*gd*(hd**2-s.Rational(2,3)*gd**2)+n*(n-2)/2*gs*(hs**2-o0**2-s.Rational(2,3)*gs**2)-n*(hd-gd)**2*g1
T=-n*g1;L=-2*n*(hd-gd)*g1
for name,v in [('X1',x1),('T',T),('L',L),('A1',x1+T+L)]:print(name,s.series(v,u,0,3))
a1=s.series(x1+T+L,u,0,3).removeO().expand()
X0=z**2*(1-u)+z**6*u**2/s.Integer(3)
D=lambda f:s.expand(z*s.diff(f,z))
cumulants=[];tmp=X0
for j in range(4):tmp=D(tmp);cumulants.append(tmp)
c1,c2,c3,c4=cumulants
b1=c1;b2=c2+c1*c1;b3=c3+3*c1*c2+c1**3;b4=c4+4*c1*c3+3*c2*c2+6*c1*c1*c2+c1**4
p1=(b2-b1)/2;p2=(3*b4-2*b3-3*b2+2*b1)/24
shift=(D(D(a1))+(2*c1-1)*D(a1))/2
Y1=s.series(a1+p1*u,u,0,3).removeO().expand()
Y2=s.series((7*z**4+26*z**6)*u*u+a1*a1/2+(a1*p1+shift)*u+p2*u*u,u,0,3).removeO().expand()
print('Y1(z)',s.factor(Y1))
print('Y2(z)',s.factor(Y2))
print('Y1(1)',Y1.subs(z,1));print('Y2(1)',Y2.subs(z,1))
k1=s.series(s.exp(X0.subs(z,1)-1)*(1+Y1.subs(z,1)),u,0,3)
k2=s.series(s.exp(X0.subs(z,1)-1)*(1+Y1.subs(z,1)+Y2.subs(z,1)),u,0,3)
print('kappa1/e',k1);print('kappa2/e',k2)
assert k2.removeO()==1-s.Rational(5,2)*u+s.Rational(65,8)*u*u
assert Y2.subs(z,1)==s.Rational(35,24)*u*u
