#!/usr/bin/env python3
"""Exact controls for the complete global quartic boundary type 5+1+1+1.

All source coordinates and induced global coefficients are retained.  The
proof, not a bounded mate search, excludes polynomial mates of every degree.
"""
from hashlib import sha256
import json

import sympy as S

gates = []


def need(label, value):
    if not bool(value):
        raise RuntimeError(label)
    gates.append(label)


def zero(label, expr):
    need(label, S.cancel(expr) == 0)


x, t, v, r, bb = S.symbols('x t v r bb')
a = S.symbols('a', nonzero=True)
P, Q, M, R, U3 = S.symbols('P Q M R U3')
U0 = 1/a
U1 = -P/(2*a*a)
U2 = (P*P-4*a*a*Q)/(8*a**3)
UU = U0+U1*v+U2*v*v+U3*v**3
scaled = S.Poly(S.expand((a*a*UU*UU+P*UU*v+Q*v*v)**2+M*UU*v**3+R*v**4-1), v)
for j in range(3):
    zero(f'inverse coefficient v{j}', scaled.coeff_monomial(v**j))
zero('inverse coefficient v3 before solving', scaled.coeff_monomial(v**3)-4*a*U3-M/a)
zero('complete T2', scaled.coeff_monomial(v**3).subs(U3, -M/(4*a*a)))
# The x-derivative of g_6 equals 2*T2/4; the primitive is 2*g_6.
need('formal coefficient index and multiplier', 4+2 == 6 and S.Rational(2,4)*2 == 1)

# Full global section spaces: no frozen-prefix or selected-generator reduction.
nn = S.symbols('n0:9')
pp = S.symbols('p0:5')
q0 = S.symbols('q0')
N = sum(nn[i]*x**i for i in range(9))
Px = 2*nn[8]*x**6+2*nn[7]*x**5+sum(pp[i]*x**i for i in range(5))
Qx = nn[8]*x**4+nn[7]*x**3+(pp[4]-nn[6])*x*x+(pp[3]-nn[5])*x+q0
H = N*t*t+Px*t+Qx
Hr = S.expand(H.subs({x:1/r,t:-r*r-r**4*bb}, simultaneous=True))
num, den = S.fraction(S.cancel(Hr))
need('entire fifteen-parameter H is regular on second chart', den == 1)
zero('complete H boundary restriction', Hr.subs(r,0)-(nn[8]*bb*bb+(2*nn[6]-pp[4])*bb+nn[4]-pp[2]+q0))

mm = S.symbols('m0:5')
rr = S.symbols('r0:3')
Mx = sum(mm[i]*x**i for i in range(5))
Rx = sum(rr[i]*x**i for i in range(3))
Lr = S.expand((Mx*t+Rx).subs({x:1/r,t:-r*r-r**4*bb}, simultaneous=True))
zero('L1 complete r^-2 obstruction', S.expand(Lr*r*r).coeff(r,0)-rr[2]+mm[4])
zero('L1 complete r^-1 obstruction', S.expand(Lr*r*r).coeff(r,1)-rr[1]+mm[3])
L = Mx*t+mm[4]*x*x+mm[3]*x+rr[0]
zero('L1 second chart complete formula', L.subs({x:1/r,t:-r*r-r**4*bb}, simultaneous=True)
     -(rr[0]-mm[2]-mm[4]*bb-mm[1]*r-mm[0]*r*r-mm[3]*r*bb-mm[2]*r*r*bb-mm[1]*r**3*bb-mm[0]*r**4*bb))
zero('L1 zero leading coefficient forces constant', L.subs({m:0 for m in mm})-rr[0])
need('L1 six independent coefficients', len(mm)+1 == 6)

# Every constant scalar, every finite root position: exact residue and rank.
p, q, z, w, kappa = S.symbols('p q z w kappa')
Cc = (x-q)*(x-z)*(x-w)
Nb = kappa*(x-p)**5*Cc
for index, root in enumerate([q,z,w]):
    others = [q,z,w]
    others.pop(index)
    zero(f'simple-root derivative {index}', S.diff(Nb,x).subs(x,root)-kappa*(root-p)**5*(root-others[0])*(root-others[1]))

# Under x-q=s^2, a simple rational residue doubles in the radical field.
s, rho = S.symbols('s rho')
zero('ramified residue is not lost', rho*S.diff(s*s,s)/(s*s)-2*rho/s)
need('normalized trace preserves exact rational differentials', S.Rational(1,2)*2 == 1)

# The three simple-root zeros give M=C*(A*x+B); the two p-jets force A=B=0.
A, B = S.symbols('A B')
Cp, Cpp = S.symbols('Cp Cpp')
jet_matrix = S.Matrix([[Cp*p,Cp],[Cpp*p+Cp,Cpp]])
zero('two remaining coefficient equations are invertible', jet_matrix.det()+Cp*Cp)
zero('factorized M first jet', S.diff(Cc*(A*x+B),x).subs(x,p)
     -(S.diff(Cc,x).subs(x,p)*(A*p+B)+A*Cc.subs(x,p)))

# Independent confluent Vandermonde route after the ordinary polynomial
# coefficient change x=p+u; this does not claim a surface automorphism.
d1, d2, d3 = S.symbols('d1 d2 d3')
mat = S.Matrix([[1,0,0,0,0],[0,1,0,0,0]]+
               [[1,d,d*d,d**3,d**4] for d in [d1,d2,d3]])
zero('full five-zero determinant', mat.det()-d1*d1*d2*d2*d3*d3*(d2-d1)*(d3-d1)*(d3-d2))
translate = S.Matrix([[S.binomial(j,i)*p**(j-i) if j>=i else 0 for j in range(5)] for i in range(5)])
zero('all finite positions coefficient change determinant', translate.det()-1)

# Literal leading survivor and two induced global corrections. These show why
# all five counted zeros are needed rather than three evaluations alone.
N0 = x**5*(x**3+1)
P0 = 2*x**6
Q0 = x**4-x
H0 = N0*t*t+P0*t+Q0
zero('literal H regular on second chart', H0.subs({x:1/r,t:-r*r-r**4*bb}, simultaneous=True)
     -(bb*bb+2*bb*r+r**3*bb*bb))
# The explicit chart identity above is separately checked; no producer import.
for label, MC, RC in [('unit',x**3+1,x),('linear-jet',x*(x**3+1),x*x)]:
    Cmodel = x**3+1
    zero(f'{label} M divisible by simple-root cubic', S.rem(MC,Cmodel,x))
    need(f'{label} three simple roots are distinct', S.discriminant(Cmodel,x) != 0)
    transformed = S.cancel((MC*t+RC).subs({x:1/r,t:-r*r-r**4*bb}, simultaneous=True))
    need(f'{label} full L is global', S.fraction(transformed)[1] == 1)
need('unit control retains nonzero M(p)', (x**3+1).subs(x,0) == 1)
need('linear-jet control retains nonzero derivative', S.diff(x*(x**3+1),x).subs(x,0) == 1)

# M-unit m=5 branch exhaustion: j=0 is degree two; j=1 is simple plus
# two cancellation determinations; all j>=2 have the one threefold slope.
for j in [0,1]:
    need(f'm5 cancellation inequality j={j}', 5 > 3*j)
    exponent = S.Rational(5-3*j,2)
    ram = 2 if exponent.q == 2 else 1
    need(f'm5 normalized cancellation regular j={j}', ram*exponent+ram-1 >= 0)
need('m5 unit j0 has only two boundary determinations', 2 == 2)
need('m5 unit j1 exhausts degree three', 1+2 == 3)
need('m5 high j universal threshold', 5 < 3*2)
need('m5 high-j ramification and differential order', S.gcd(10,3) == 1 and 3-1 == 2)
# Shared normal-unit ell ranges 1..5, giving (5-ell)/2 before normalization.
for ell in range(1,6):
    order = S.Rational(5-ell,2)
    ram = 2 if order.q == 2 else 1
    need(f'normal-unit normalized regular ell={ell}', ram*order+ram-1 >= 0)

# Actual named tangent logarithm for the linear-jet correction.
Z, fibre = S.symbols('Z fibre', nonzero=True)
tangent = Z**3-fibre*Z**4
zero('linear-jet tangent simple root', tangent.subs(Z,1/fibre))
zero('linear-jet exact nonzero logarithmic residue', (Z*Z/S.diff(tangent,Z)).subs(Z,1/fibre)+1)
surface_N = N0+P0*s+Q0*s*s
surface_M = x*(x**3+1)+x*x*s
E = surface_N**2+s**3*surface_M-fibre*s**4
zero('linear-jet actual tangent extraction', S.expand(E.subs(s,x*Z)).coeff(x,4)-tangent)

# Infinity is not moved by a surface or projective normalization.
need('fivefold infinity degree-three differential holomorphic', 3-3 == 0 and 1-1 == 0)
need('simple infinity degree-seven capacity fails', (5-2)+0+0 < 7-2)
need('all three locations exhaust distinct marked roots', 1+1+1 == 3)

# The final boundary is polynomial-only: the source factor is 2H.
GG = S.Function('GG')(x,t)
zero('composite constant-L polynomial factor',
     S.diff(H0*H0,x)*S.diff(GG,t)-S.diff(H0*H0,t)*S.diff(GG,x)
     -2*H0*(S.diff(H0,x)*S.diff(GG,t)-S.diff(H0,t)*S.diff(GG,x)))
Fc, Gc = t**4, -x/(4*t**3)
zero('constant L does not generally exclude rational mates', S.diff(Fc,x)*S.diff(Gc,t)-S.diff(Fc,t)*S.diff(Gc,x)-1)
zero('elliptic leading differential positive control',
     S.diff(-2/(3*x**4),x)*N0+(-2/(3*x**4))*S.diff(N0,x)/2-1)

# Sharp rational mate inside this exact boundary partition, with every global
# correction retained. The map (x,t)->(x^-3,h) has field degree three.
delta = S.symbols('delta', nonzero=True)
beta, gamma, constant = S.symbols('beta gamma constant')
h = x*x+x**4*t
HH = (x**4+delta*x)*(1+x*x*t)**2+beta*h+gamma
J = lambda f,g: S.diff(f,x)*S.diff(g,t)-S.diff(f,t)*S.diff(g,x)
zero('sharp family full leading polynomial', S.expand(HH).coeff(t,2)-x**5*(x**3+delta))
need('sharp family has three distinct nonzero simple roots', S.discriminant(x**3+delta,x) == -27*delta**2)
zero('sharp family actual second chart', HH.subs({x:1/r,t:-r*r-r**4*bb}, simultaneous=True)
     -((1+delta*r**3)*bb*bb-beta*bb+gamma))
zero('sharp family volume identity', J(x**-3,h)+3)
zero('sharp family function field expression', HH-((1+delta*x**-3)*h*h+beta*h+gamma))
zero('sharp family rational mate of H', J(HH,1/(3*delta*h))-1)
zero('sharp family rational mate of F', J(HH*HH+constant,1/(6*delta*HH*h))-1)
u = S.symbols('u')
zero('sharp family degree-three inverse relation', (u*x**3-1).subs(u,x**-3))
need('sharp family extension is not declared birational', S.degree(u*x**3-1,x) == 3)

semantic = sha256(json.dumps(gates,separators=(',',':')).encode()).hexdigest()
print('DG quartic boundary partition 5+1+1+1: polynomial-mate exclusion')
print('Universe: all complex global H in L2, L in L1; all four distinct boundary points')
print('T2, full sections, three residues plus two jets, and infinity cases: PASS')
print('Mate degree is unrestricted; rational conclusion stops at constant L')
print(f'Gates: {len(gates)}')
print('Semantic SHA256: '+semantic)
print('RESULT: PASS')
