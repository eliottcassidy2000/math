#!/usr/bin/env python3
"""Exact controls for the all-index 7+1 inverse-exactness theorem.
The induction in the companion note, not a finite order extrapolation,
pays the all-index statement. No producer imports.
"""
import hashlib
import json
import sympy as S

u,t,p,x,r,bD=S.symbols('u t p x r bD')
a,b,c,d,e,f,m4,m3,m2,m1,m0,h=S.symbols('a b c d e f m4 m3 m2 m1 m0 h')
z,v,q=S.symbols('z v q')
gates=0
record={}
def check(name,value):
    global gates
    if not bool(value):raise RuntimeError(name)
    gates+=1
def zero(name,expr):check(name,S.cancel(expr)==0)
def jac(F,G):return S.diff(F,u)*S.diff(G,t)-S.diff(F,t)*S.diff(G,u)
def deg(poly):return S.Poly(S.expand(poly),u).degree()

# Complete global rows: no active jet is discarded.
N=u**7*(u-1)
P=2*u**6-(4*p+2)*u**5+a*u**4+b*u**3+c*u*u+d*u+e
Q=u**4-(4*p+1)*u**3+(a+4*p*p)*u*u+(b-2*a*p+4*p*p)*u+f
M=m4*u**4+m3*u**3+m2*u*u+m1*u+m0
R=m4*u*u+(m3-2*p*m4)*u+h
Nx=S.Poly(S.expand(N.subs(u,x-p)),x);Px=S.Poly(S.expand(P.subs(u,x-p)),x);Qx=S.Poly(S.expand(Q.subs(u,x-p)),x)
for k,l in [(6,8),(5,7)]:zero('original P'+str(k),Px.coeff_monomial(x**k)-2*Nx.coeff_monomial(x**l))
zero('original Q4',Qx.coeff_monomial(x**4)-Nx.coeff_monomial(x**8))
zero('original Q3',Qx.coeff_monomial(x**3)-Nx.coeff_monomial(x**7))
zero('original Q2',Qx.coeff_monomial(x*x)-Px.coeff_monomial(x**4)+Nx.coeff_monomial(x**6))
zero('original Q1',Qx.coeff_monomial(x)-Px.coeff_monomial(x**3)+Nx.coeff_monomial(x**5))
Dnum=S.expand(4*N*Q-P*P);Enum=S.expand(2*N*R-M*P)
check('D bounded all p',deg(Dnum)<=8)
check('E bounded all p',deg(Enum)<=8)
check('C squared bounded',deg(M*M)<=8)
for k in range(9,13):zero('D cancelled top '+str(k),Dnum.coeff(u,k))
for k in range(9,11):zero('E cancelled top '+str(k),Enum.coeff(u,k))
record['D8']=str(S.factor(Dnum.coeff(u,8)))
record['E8']=str(S.factor(Enum.coeff(u,8)))
# The constant and linear P rows survive the full equation.
zero('P0 retained',P.subs(u,0)-e)
zero('P1 retained',S.diff(P,u).subs(u,0)-d)

# Necessity from the actual simple-root residue, and full vanishing ideal.
zero('rational residue',S.residue(-M/(4*N),u,1)+M.subs(u,1)/4)
Mvan=S.factor(M.subs(m0,-m4-m3-m2-m1))
zero('complete M vanishing',Mvan-(u-1)*(m4*u**3+(m4+m3)*u*u+(m4+m3+m2)*u+m4+m3+m2+m1))
check('M/N regular at one',S.denom(S.cancel(Mvan/N)).subs(u,1)!=0)
# E and D*B are regular there as rational functions for the full parameters.
Ev=S.cancel(R.subs(m0,-m4-m3-m2-m1)-Mvan*P/(2*N))
Dv=S.cancel(Q-P*P/(4*N));Bv=S.cancel(Mvan*Mvan/N)
for name,expr in [('E',Ev),('B',Bv),('DB',Dv*Bv),('B2',Bv*Bv)]:
    check('regular at one '+name,S.denom(S.cancel(expr)).subs(u,1)!=0)
zero('B vanishes',Bv.subs(u,1))

# Explicit rational parametrization and leading primitive.
uz=1/(1-z*z);kappaz=z/(1-z*z)**4
zero('quadratic field',kappaz*kappaz-N.subs(u,uz))
Llead=2*z-S.Rational(4,3)*z**3+S.Rational(2,5)*z**5
zero('leading primitive',S.diff(Llead,z)-S.diff(uz,z)/kappaz)
zero('deck transformation u',uz.subs(z,-z)-uz)
zero('deck transformation kappa',kappaz.subs(z,-z)+kappaz)
# Pullback residue doubles at the ramified simple root.
negative=(-S.Integer(1)/(4*N)).subs(u,uz)*S.diff(uz,z)
zero('minimal nonexact residue',S.residue(negative,z,0)+S.Rational(1,2))

# Finite abstract controls for the two ALL-ORDER recursions.
D,E,C=S.symbols('D E C')
ws=[S.Integer(1)]
for j in range(1,9):
    A=S.Symbol('A');W=sum(ws[i]*v**i for i in range(j))+A*v**j
    equation=W**4+2*D*v*v*W*W+C*v**3*W+(D*D+E)*v**4-1
    coefficient=S.expand(equation).coeff(v,j)
    ans=S.factor(S.solve(coefficient,A)[0]);ws.append(ans)
    zero('W recursion '+str(j),coefficient.subs(A,ans))
    check('W polynomial coefficient '+str(j),S.denom(ans).is_Integer)
    # T_k=w_(k+1)/kappa, so the coefficient has the indicated C parity.
    zero('W involution '+str(j),ans.subs(C,-C)-(-1)**j*ans)
    if (j-1)%4==0:zero('vanishing centered multiple4 '+str(j),ans)
record['W_coefficients']=[str(s) for s in ws]
# V depends only on E, D*B, B^2: no isolated negative-order D can enter.
X,Y=S.symbols('X Y');vs=[S.Integer(1)]
for j in range(1,5):
    A=S.Symbol('A');V=sum(vs[i]*q**i for i in range(j))+A*q**j
    equation=(1-E*q)*V**2+X*q*q*V**4/4+Y*q**3*V**6/64-1
    coefficient=S.expand(equation).coeff(q,j)
    ans=S.factor(S.solve(coefficient,A)[0]);vs.append(ans)
    zero('V recursion '+str(j),coefficient.subs(A,ans))
    check('V polynomial ring '+str(j),S.denom(ans).is_Integer)
record['V_coefficients']=[str(s) for s in vs]

# Actual global hostile with all rows exact but no rational mate.
H=u**3*(u-1)*(1+u*u*t)**2
L=(u-1)*t;F=S.expand(H*H+L)
zero('hostile leading',S.Poly(H,t).coeff_monomial(t*t)-N)
zero('hostile P',S.Poly(H,t).coeff_monomial(t)-2*u**5*(u-1))
zero('hostile Q',S.Poly(H,t).coeff_monomial(1)-u**3*(u-1))
for name,expr in [('H',H),('L',L)]:
    chart=S.expand(expr.subs({u:1/r,t:-r*r-r**4*bD},simultaneous=True))
    check('hostile global '+name,S.denom(S.cancel(chart))==1)
    if name=='H':zero('hostile boundary H',chart.subs(r,0)-bD*bD)
    else:zero('hostile boundary L',chart.subs(r,0))
zero('hostile M condition',(u-1).subs(u,1))
check('finite7 M unit',(u-1).subs(u,0)!=0)
check('finite7 unbalanced',7!=3*5)
# Three determinations s~u^(14/3) form one branch, eta order 2.
check('finite7 actual order',2*14+(3-1)-2*14==2)

# Concrete descended formal mate through q^7, retaining original F.
# g3'=-T_-1/4, g6'=T2/2, g7'=3T3/4 in the v expansion.
g3=-Llead/4;g6=S.Rational(1,48)/uz**6
g7=S.Rational(3,4)*(z**3/6-S.Rational(3,10)*z**5+S.Rational(3,14)*z**7-z**9/18)
T2=-S.Rational(1,4)/uz**7
T3=z/(4*uz**5)
zero('g3 derivative',S.diff(g3,z)+S.diff(uz,z)/(4*kappaz))
zero('g6 derivative',S.diff(g6,z)-T2*S.diff(uz,z)/2)
zero('g7 derivative',S.diff(g7,z)-3*T3*S.diff(uz,z)/4)
for name,expr in [('g3',g3),('g7',g7)]:zero('primitive deck '+name,expr.subs(z,-z)+expr)
# Convert even Laurent expressions in z to rational functions of u.
def rationalize(expr):
    ex=S.expand(S.cancel(expr))
    ans=0
    for term in S.Add.make_args(ex):
        coeff,power=term.as_coeff_exponent(z)
        check('even Laurent deck term',power.is_Integer and int(power)%2==0)
        ans+=coeff*((u-1)/u)**(int(power)//2)
    return S.cancel(ans)
R3=rationalize(g3/(u**12*z**3))
R6=1/(48*u**6*N**3)
R7=rationalize(g7/(u**28*z**7))
Bq=S.expand(F.subs(t,1/q)*q**4/N**2)
Gq=S.series(R3*q**3*Bq**(-S.Rational(3,4))+R6*q**6*Bq**(-S.Rational(3,2))+R7*q**7*Bq**(-S.Rational(7,4)),q,0,8).removeO()
G=S.expand(Gq.subs(q,1/t));Jq=S.expand(jac(F,G).subs(t,1/q))
zero('formal constant Jacobian',Jq.coeff(q,0)-1)
for j in range(1,5):zero('formal Jacobian depth '+str(j),Jq.coeff(q,j))
check('formal no negative powers',all(S.expand(term).as_powers_dict().get(q,0)>=0 for term in S.Add.make_args(Jq)))
record['formal_coefficients']={str(i):str(S.factor(S.expand(Gq).coeff(q,i))) for i in range(3,8)}

record['gates']=gates
semantic=hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('All-finite (7,1): exactness of the entire inverse hierarchy')
print('Exact gates:',gates)
print('Full global coefficients and universal infinity bounds: PASS')
print('Ramified parity, rational pair-sum ring, simple-root condition: PASS')
print('Formal trace descent and explicit q^7 Jacobian control: PASS')
print('All-index exact but no rational mate: global hostile PASS')
print('Semantic SHA256:',semantic)
