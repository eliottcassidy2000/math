#!/usr/bin/env python3
"""Exact symbolic controls for the complete binary (5,2,1) entry.
No bounded mate degree or numerical parameter census is used.
"""
import hashlib
import json
import sympy as S

u,t,p,a,b,c,d,lam,mu,A=S.symbols('u t p a b c d lam mu A')
v,z,D,E,K,NN,MM=S.symbols('v z D E K NN MM')
gates=0
record={}
def check(name, value):
    global gates
    if not bool(value): raise RuntimeError(name)
    gates+=1
def zero(name,expr): check(name,S.cancel(expr)==0)
def jac(f,g): return S.diff(f,u)*S.diff(g,t)-S.diff(f,t)*S.diff(g,u)
def res(expr): return S.factor(S.residue(S.cancel(expr),u,0))
def rem(expr,poly,var=A):
    num,den=S.fraction(S.cancel(expr))
    check('remainder scalar denominator',not den.has(var))
    return S.factor(S.rem(num,poly,var)/den)

# Full shifted section constraints, before the active jets are imposed.
x,r,bD=S.symbols('x r bD');P1,P0,M4,M3,M2,M1,M0,R0=S.symbols('P1 P0 M4 M3 M2 M1 M0 R0')
N=u**5*(u-1)
P=a*u**4+b*u**3+c*u*u+P1*u+P0
Q=(a-1)*u*u+(b-2*p*a+4*p+1)*u+d
M=M4*u**4+M3*u**3+M2*u*u+M1*u+M0
R=M4*u*u+(M3-2*p*M4)*u+R0
Nx=S.Poly(S.expand(N.subs(u,x-p)),x);Px=S.Poly(S.expand(P.subs(u,x-p)),x)
Qx=S.Poly(S.expand(Q.subs(u,x-p)),x)
zero('original Q2',Qx.coeff_monomial(x*x)-Px.coeff_monomial(x**4)+Nx.coeff_monomial(x**6))
zero('original Q1',Qx.coeff_monomial(x)-Px.coeff_monomial(x**3)+Nx.coeff_monomial(x**5))
for name,f in [('H',N*t*t+P*t+Q),('L',M*t+R)]:
    fi=S.expand(f.subs({u:1/r-p,t:-r*r-r**4*bD},simultaneous=True))
    for j in range(-4,0):zero('global '+name+str(j),fi.coeff(r,j))
    check('global polynomial '+name,S.denom(S.cancel(fi))==1)

# T2 simple-root condition and full moving-root first rows.
zero('full T2 simple residue',S.residue(M/N,u,1)-M.subs(u,1))
zero('complete factored M', (M4*u**4+M3*u**3+M2*u*u).subs(M2,-M4-M3)-u*u*(u-1)*(M4*u+M4+M3))
P=u*u*(a*u*u+b*u+c)
M=u*u*(u-1)*(lam*u+mu);R=lam*u*u+(mu-lam-2*p*lam)*u
H=N*t*t+P*t+Q;L=M*t+R
zero('T2 simple condition',M.subs(u,1))
zero('T2 finite residues',res(M/N)+S.residue(M/N,u,1))
w=S.symbols('w');Fw=S.expand((H*H+L).subs(t,w/u**2))
f=(c*w+d)**2-mu*w
g=2*(c*w+d)*(-w*w+b*w+b-2*p*a+4*p+1)+(mu-lam)*(w+1)-2*p*lam
zero('full constant row',Fw.coeff(u,0)-f)
zero('full first row',Fw.coeff(u,1)-g)
zero('cubic mismatch',S.expand(g).coeff(w,3)+2*c)
zero('n2 quadratic',S.expand(g.subs(c,0)).coeff(w,2)+2*d)
zero('n2 linear',S.expand(g.subs({c:0,d:0})).coeff(w,1)-mu+lam)

# Independent direct quartic inverse to v^6. Here w=sqrt(N)*v*y.
cs=[S.Integer(1)]
for j in range(1,7):
    h=S.Symbol('h'); W=sum(cs[i]*v**i for i in range(j))+h*v**j
    poly=W**4+2*D*v*v*W*W+K*v**3*W+(D*D+E)*v**4-1
    co=S.expand(poly).coeff(v,j)
    hvalue=S.factor(S.solve(co,h)[0]);cs.append(hvalue)
    zero('direct inverse order '+str(j),co.subs(h,hvalue))
zero('T1 numerator',cs[2]+D/2)
zero('T2 numerator',cs[3]+K/4)
zero('T5 numerator',cs[6]+(2*D**3+4*D*E+K*K)/32)
record['direct_inverse']=[str(q) for q in cs]

# Ferrari opposite-pair coefficients, checked by literal substitution.
Vco=[S.Integer(1),E/2,(3*E*E-D*MM*MM/NN)/8,
     (40*E**3*NN**2-40*D*E*MM**2*NN-MM**4)/(128*NN**2),
     7*(2*D*D*MM**4-20*D*E*E*MM*MM*NN+10*E**4*NN**2-E*MM**4)/(256*NN**2)]
VV=sum(Vco[j]*z**j for j in range(5))
for j in range(1,5):
    # Only coefficients needed for this order are expanded.
    Vj=sum(Vco[i]*z**i for i in range(j+1))
    poly=(1-E*z)*Vj**2+D*MM*MM*z*z*Vj**4/(4*NN)+MM**4*z**3*Vj**6/(64*NN**2)-1
    zero('Ferrari order '+str(j),S.expand(poly).coeff(z,j))

def rows(P,Q,M,R):
    D0=S.cancel(Q-P*P/(4*N));E0=S.cancel(R-M*P/(2*N))
    rr={}
    for j in range(1,5):
        term=(-MM*Vco[j]/(4*NN)).subs({D:D0,E:E0,MM:M,NN:N},simultaneous=True)
        rr[4*j+2]=res(term)
    return D0,E0,rr

# n=2: all-parameter residues through T18 and exact ideal reduction.
P2=u**3*((A+2)*u-2-4*p-2*A)
Q2=(A+1)*u*u+(-1-2*A-2*p*A-4*p)*u
M2=lam*u*u*(u*u-1);R2=lam*(u*u-2*p*u)
D2,E2,rr2=rows(P2,Q2,M2,R2)
qA=A*A+6*p*A+12*p*p
zero('n2 T6 after tuning',rr2[6])
zero('n2 T10',rr2[10]+lam**3*qA/32)
zero('n2 T14 mod',rem(rr2[14],qA)+5*lam**4*(lam+32*p*p*(A+3*p))/512)
lam_value=-32*p*p*(A+3*p)
zero('n2 T18 final',rem(rr2[18].subs(lam,lam_value)+77*lam_value**5*p**3*(A+4*p)/64,qA))
zero('n2 contradictory quadratic',qA.subs(A,-4*p)-4*p*p)
check('n2 p0 cannot lambda',lam_value.subs(p,0)==0)
record['n2_residues']={str(k):str(v) for k,v in rr2.items()}

# Before tuning, retain independently extracted T6 in both n strata.
Pg=u**3*(a*u+b);Qg=(a-1)*u*u+(b-2*p*a+4*p+1)*u+d
Mg=u*u*(u-1)*(lam*u+mu);Rg=lam*u*u+(mu-lam-2*p*lam)*u
Eg=S.cancel(Rg-Mg*Pg/(2*N));Dg=S.cancel(Qg-Pg*Pg/(4*N))
R6=res(-Mg*Eg/(8*N));R10=res(-Mg*(3*Eg*Eg-Dg*Mg*Mg/N)/(32*N))
zero('general T6',R6-lam*(2*a*mu+b*lam+4*lam*p+2*lam-4*mu)/16)
zero('n3 T10',R10.subs(mu,0)+d*lam**3/32)

# n=3: the two algebraic infinity residues, before and after T1.
P3=u**3*((A+2)*u-4*p-2)
Q3=(A+1)*u*u-(2*p*(A+2)+1)*u
M3=lam*u**3*(u-1);R3=lam*u*(u-1-2*p)
D3=S.cancel(Q3-P3*P3/(4*N));E3=S.cancel(R3-M3*P3/(2*N))
zero('n3 centered D',D3+u*(A*A*u*u-8*A*p+16*p*p)/(4*(u-1)))
zero('n3 centered E',E3+A*lam*u*u/2)
# y^2=u(u-1); at infinity y=epsilon*u*sqrt(1-1/u).
zz,eps=S.symbols('zz eps');y=eps/zz*S.sqrt(1-zz)
T1=(A*A*u*u-8*A*p+16*p*p)/(8*S.Symbol('y')**3)
infT1=S.series(T1.subs({u:1/zz,S.Symbol('y'):y},simultaneous=True)*(-1/zz**2),zz,0,0).removeO()
zero('T1 infinity epsilon+',infT1.subs(eps,1).coeff(zz,-1)+A*A/8)
zero('T1 infinity epsilon-',infT1.subs(eps,-1).coeff(zz,-1)-A*A/8)
D30=D3.subs(A,0);E30=E3.subs(A,0)
T5=-(2*D30**3+4*D30*E30+M3*M3/N)/(32*u*u*S.Symbol('y'))
infT5=S.series(T5.subs({u:1/zz,S.Symbol('y'):y},simultaneous=True)*(-1/zz**2),zz,0,0).removeO()
zero('T5 infinity epsilon+',infT5.subs(eps,1).coeff(zz,-1)-lam*lam/32)
zero('T5 infinity epsilon-',infT5.subs(eps,-1).coeff(zz,-1)+lam*lam/32)

# Exact same-partition rational boundary with constant L, every p.
kappa,e=S.symbols('kappa e');Z=u*u*t+1;HH=u*(u-1)*Z*Z+kappa
GH=-(1+2*u)/(3*u*u*Z)
zero('sharp H mate',jac(HH,GH)-1)
zero('sharp F mate',jac(HH*HH+e,GH/(2*HH))-1)
zero('sharp leading',S.Poly(HH,t).coeff_monomial(t*t)-N)
Hinf=S.expand(HH.subs({u:1/r-p,t:-r*r-r**4*bD},simultaneous=True))
check('sharp global',S.denom(S.cancel(Hinf))==1)
zero('sharp P',S.Poly(HH,t).coeff_monomial(t)-2*u**3*(u-1))

# Minimal staged hostiles: earlier exact coefficients alone do not close n2.
zero('n2 earlier-pass control',rr2[10].subs({A:0,p:0}))
check('n2 T14 hostile',rr2[14].subs({A:0,p:0,lam:1})!=0)
zero('n3 T1-pass control',infT1.subs({A:0,eps:1}).coeff(zz,-1))
check('n3 T5 hostile',infT5.subs({eps:1,lam:1}).coeff(zz,-1)!=0)

record['gates']=gates
semantic=hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('Binary (5,2,1): complete inverse-coefficient obstruction')
print('Exact gates:',gates)
print('Full shifted global rows and moving residues: PASS')
print('n2: T6,T10,T14,T18 exact symbolic contradiction: PASS')
print('n3: T6,T10 and algebraic T1,T5 residues: PASS')
print('All-location polynomial scope; constant-L rational sharpness: PASS')
print('Semantic SHA256:',semantic)
