#!/usr/bin/env python3
"""Independent finite algebra controls; no import from the producer checker.

No test is a proof of uniqueness of embeddings or of the universal obstruction.
Exceptions, not assert, keep every check active under python -O.
"""
import argparse
import json
from pathlib import Path
import sympy as S

x, y, W, Z, U, V, t = S.symbols('x y W Z U V t')
checks = 0
families = []

def check(name, proposition):
    global checks
    checks += 1
    if not bool(proposition):
        raise RuntimeError(name)

def equal(name, lhs, rhs=0):
    check(name, S.cancel(lhs-rhs) == 0)

def jac(p, q, a, b):
    return S.diff(p, a)*S.diff(q, b)-S.diff(p, b)*S.diff(q, a)

def integral(poly, lo, hi):
    # Direct coefficient primitive. This works for nonsymmetric endpoints too.
    primitive = sum(co*y**(mon[0]+1)/S.Integer(mon[0]+1)
                    for mon,co in S.Poly(S.expand(poly),y).terms())
    return S.expand(primitive.subs(y,hi)-primitive.subs(y,lo))

def split(poly):
    even, odd = S.Integer(0), S.Integer(0)
    for (power,), co in S.Poly(S.expand(poly),y).terms():
        if power % 2:
            odd += co*W**((power-1)//2)
        else:
            even += co*W**(power//2)
    return S.expand(even),S.expand(odd)

def in_ideal(poly, gens):
    gb=S.groebner(gens,W,x,order='lex')
    return S.expand(gb.reduce(S.expand(poly))[1]) == 0

def in_ring(poly, gens):
    return in_ideal(split(poly)[1],gens)

# Independent ring-membership and conductor controls. Test polynomial ideal
# membership over Q[x,W], not rational-function divisibility over Q(x)[W].
models = [
    ('constant_node', W-1),
    ('moving_simple', W**2+(x*x+1)*W+2),
    ('nonmonic_no_empty', x*W**3+W+2),
    ('common_zero', 2+x*W+x*x*W**2),
    ('prime_to_composite_drop', x*W**5+W**4+2),
    ('repeated_generic', (W-2)**2*(x*W+1)),
]
for name,g in models:
    gy=g.subs(W,y*y)
    w,z=y*y,y*gy
    equal(name+': relation',z*z-w*g.subs(W,w)**2)
    for idx,(ee,oo) in enumerate([(0,0),(1,0),(0,1),(g,0),(0,g),
                                  (g*(x+W),g*(W*W-x)),
                                  (g+1,g),(g,g+1)]):
        f=S.sympify(ee).subs(W,y*y)+y*S.sympify(oo).subs(W,y*y)
        rhs=in_ideal(ee,[g]) and in_ideal(oo,[g])
        lhs=in_ring(f,[g]) and in_ring(y*f,[g])
        check(name+f': conductor characterization {idx}',lhs == rhs)
        e,o=split(f)
        equal(name+f': reconstruct even odd {idx}',f,e.subs(W,y*y)+y*o.subs(W,y*y))
    # Explicit inverse witnesses prove entire R-generation, including x=0.
    # (P,Q)=(w+x z+z^2, -3z+x^2) has inverse z=(Q-x^2)/-3,
    # w=P-xz-z^2 over R, not only over C(x).
    p=w+x*z+z*z
    q=-3*z+x*x
    zi=(q-x*x)/-3
    wi=p-x*zi-zi*zi
    equal(name+': inverse z over R',zi,z)
    equal(name+': inverse w over R',wi,w)
    for xi in [0,1,-2]:
        equal(name+f': specialized inverse z {xi}',zi.subs(x,xi),z.subs(x,xi))
        equal(name+f': specialized inverse w {xi}',wi.subs(x,xi),w)
        check(name+f': fiber g never zero {xi}',g.subs(x,xi)!=0)
    families.append({'kind':'normal_form_generation','name':name,'g':str(g)})

# Principal recognition boundary: I=(x,W) gives an even/odd ring, but its
# conductor (x,y^2) is not principal. GCD=1 and 1 is absent from the ideal.
check('nonprincipal gcd',S.gcd(x,y*y)==1)
check('nonprincipal proper ideal',not in_ideal(1,[x,W]))
for f in [x,y*y,x*y,y**3,x+y**2,x*y+y**3,1,y]:
    e,o=split(f)
    check('nonprincipal conductor '+str(f),
          (in_ring(f,[x,W]) and in_ring(y*f,[x,W])) ==
          (in_ideal(e,[x,W]) and in_ideal(o,[x,W])))
# Reject the invalid odd generator shortcut using y*h: h is not in (y*h).
h=x+y*y+1
check('odd conductor impossible quotient',S.rem(h,y*h,y)!=0)
families.append({'kind':'nonprincipal_recognition_boundary','ideal':'(x,W)'})

# Same conductor support/multiplicity does not fix the equalizer algebra or
# guarantee a quadratic normalization coordinate. These are explicit abstract
# equalizer rings, not claimed to be Keller source algebras or two-generated.
roots=[-2,-1,1,2]
C=S.prod(y-r for r in roots)
standard_pairs=[(-1,1),(-2,2)]
cross_pairs=[(1,2),(-1,-2)]
for k in range(5):
    for a,b in standard_pairs+cross_pairs:
        equal(f'equalizer conductor witness {k},{a},{b}',
              (C*y**k).subs(y,a)-(C*y**k).subs(y,b))
for a,b in standard_pairs:
    equal('symmetric quadratic '+str((a,b)),a*a-b*b)
b=S.symbols('b')
cross_eq=[S.expand((a*a+b*a)-(c*c+b*c)) for a,c in cross_pairs]
check('cross pairing forbids monic quadratic',S.solve(cross_eq,[b])==[])
for a,beta in cross_pairs:
    equal('cross equalizer cubic witness '+str((a,beta)),
          (y**3-7*y).subs(y,a)-(y**3-7*y).subs(y,beta))
families.append({'kind':'same_conductor_different_pairing','conductor':str(C)})

# Full Gauss-ring membership checks. Nonunit leading coefficients are allowed.
for name,g in models:
    for lam in [-3,2]:
        L=x*W+W**2-3
        D=lam+g*L
        check(name+f': global scalar congruence {lam}',in_ideal(D-lam,[g]))
        equal(name+f': fraction field remainder {lam}',S.rem(D,g,W),lam)
        for xi in [0,1]:
            equal(name+f': exceptional scalar remainder {lam},{xi}',
                  S.rem(D.subs(x,xi),g.subs(x,xi),W),
                  lam if S.degree(g.subs(x,xi),W)>0 else 0)
    families.append({'kind':'gauss_scalar_descent_control','name':name})
# All-collision hypothesis is essential for the scalar extension argument:
# D=-W has common generic value 1/x on g=xW+1, but the x=0 fiber is empty.
equal('empty collision rational pole',S.rem(-W-1/x,x*W+1,W))
equal('empty collision cleared identity',x*(-W)-1,-(x*W+1))
equal('empty collision fiber', (x*W+1).subs(x,0),1)
# Primitivity is required for integral divisibility, not for division over K.
g=x*(W+1)
equal('nonprimitive K divisibility',S.rem(W+1,g,W))
check('nonprimitive R divisibility fails',not in_ideal(W+1,[g]))
# Units on a root cover can vary and need not be common values across nodes.
g=W*W-x*W+1
equal('cover unit inverse',S.rem(W*(x-W),g,W),1)
check('cover unit not scalar',S.degree(S.rem(W,g,W),W)==1)
families.append({'kind':'descent_failure_boundaries'})

# Ambient representative independence at every singular point, including the
# generic prime-degree -> exceptional composite-degree transition.
for name,g in models:
    rel=Z*Z-W*g*g
    p0=W+x*Z+Z*Z
    q0=-3*Z+x*x
    p=p0+(x*W+Z+1)*rel
    q=q0+(W-Z*x)*rel
    equal(name+': ambient determinant',jac(p0,q0,W,Z),-3)
    scalar=S.expand(jac(p,q,W,Z).subs(Z,0))
    check(name+': representatives scalar ideal',in_ideal(scalar+3,[g]))
    for xi in [0,2]:
        gi=g.subs(x,xi)
        equal(name+f': special representative scalar {xi}',
              S.rem((scalar+3).subs(x,xi),gi,W))
    families.append({'kind':'ambient_representatives','name':name})

# Moving-boundary identity tested with unequal moving speeds and nonzero,
# nonconstant common endpoint values. Direct total differentiation is distinct
# from the symmetric-root symbolic calculation in the producer checker.
for u,v in [(x,x+1),(x*x,2*x+3),(-x-1,x*x+2)]:
    vanish=(y-u)*(y-v)
    p=x*x+(y+x)*vanish
    q=x**3+(y*y-x)*vanish
    area=integral(p*S.diff(q,y),u,v)
    jarea=integral(jac(p,q,x,y),u,v)
    equal('moving endpoint P equality '+str((u,v)),p.subs(y,v),p.subs(y,u))
    equal('moving endpoint Q equality '+str((u,v)),q.subs(y,v),q.subs(y,u))
    equal('moving endpoint total identity '+str((u,v)),S.diff(area,x),jarea)
    boundary=(p*(S.diff(q,x)+S.diff(q,y)*S.diff(v,x))).subs(y,v)-\
             (p*(S.diff(q,x)+S.diff(q,y)*S.diff(u,x))).subs(y,u)
    equal('moving endpoint full boundary '+str((u,v)),boundary)
    families.append({'kind':'moving_endpoint_cancellation','u':str(u),'v':str(v)})
# Deliberately omit endpoint equality: the boundary correction must be nonzero.
p,q=y*y,x*y
area=integral(p*S.diff(q,y),-x,x)
bulk=integral(jac(p,q,x,y),-x,x)
check('boundary trap is nonzero',S.cancel(S.diff(area,x)-bulk)!=0)

# Coefficient contradiction normalization, independently obtained by dividing
# a symbolic H-kappa and comparing its lowest/highest coefficients.
c,k=S.symbols('c k',nonzero=True)
for d in range(1,9):
    a=S.symbols('a1:'+str(d+1))
    da=S.symbols('da1:'+str(d+1))
    g=c+sum(a[i]*W**(i+1) for i in range(d))
    H=sum(da[i]*W**(i+2)/S.Integer(2*i+5) for i in range(d))
    A,B=S.symbols('A B')
    residual=S.Poly(S.expand(H-k-(A*W+B)*g),W)
    solved=S.solve([residual.nth(0),residual.nth(1)],[A,B])
    equal(f'd{d}: quotient constant',solved[B],-k/c)
    equal(f'd{d}: quotient slope',solved[A],k*a[0]/c**2)
    equal(f'd{d}: leading obstruction',residual.nth(d+1).subs(solved),
          da[-1]/S.Integer(2*d+3)-k*a[0]*a[-1]/c**2)
    # If a1=0, vanishing consecutive coefficients recursively kills a2,... .
    for i in range(2,d+1):
        expr=residual.nth(i).subs(solved)
        prior={a[h]:0 for h in range(i-1)}
        prior.update({da[h]:0 for h in range(i-1)})
        equal(f'd{d}: zero branch step {i}',expr.subs(prior),k*a[i-1]/c)
    families.append({'kind':'coefficient_recurrence','degree':d})

# Fixed-cylinder area transport with orientation/scaling and nonlinear target
# automorphisms; generation has explicit polynomial inverses.
for multiplicity in [1,2,3,4]:
    for scale in [-2,3]:
        w=y*y
        z=y*(y*y-1)**multiplicity
        p=z+x*w*w
        q=scale*w+(x*x+1)*p*p
        zi=p-x*((q-(x*x+1)*p*p)/scale)**2
        wi=(q-(x*x+1)*p*p)/scale
        equal(f'fixed{multiplicity},{scale}: inverse w',wi,w)
        equal(f'fixed{multiplicity},{scale}: inverse z',zi,z)
        area=integral(p*S.diff(q,y),-1,1)
        base=integral(w*S.diff(z,y),-1,1)
        equal(f'fixed{multiplicity},{scale}: area orientation',area,-scale*base)
        equal(f'fixed{multiplicity},{scale}: integral Jacobian',
              integral(jac(p,q,x,y),-1,1))
        families.append({'kind':'fixed_conormal_orientation','multiplicity':multiplicity,'scale':scale})
# A chosen collision can have zero loop period: the fixed theorem does not
# need any nonvanishing hypothesis on its area invariant.
g=7*W**2-10*W+3
w=y*y; z=y*g.subs(W,w)
equal('zero-period collision plus',z.subs(y,1))
equal('zero-period collision minus',z.subs(y,-1))
equal('zero-period base area',integral(w*S.diff(z,y),-1,1))
p=w+x*z*z; q=z+x*x
equal('zero-period generated area',integral(p*S.diff(q,y),-1,1))
families.append({'kind':'zero_fixed_period','g':str(g)})

# Polynomial source shear is global and does not invert a nonunit in R.
Y=S.symbols('Y')
b=x*x-3*x; dd=2*x+1
quadratic=y*y+b*y+dd
equal('monic quadratic shear',(y+b/2)**2-quadratic,b*b/4-dd)
p=x+y*y; q=y+p*p
newp=p.subs(y,Y-b/2); newq=q.subs(y,Y-b/2)
equal('source shear preserves Jacobian',jac(newp,newq,x,Y),
      jac(p,q,x,y).subs(y,Y-b/2))
families.append({'kind':'global_monic_quadratic_shear'})

# Generic-only generation fails exactly at x=0 and gives a varying area.
w=y*y-1; z=y*w
area=integral(w*S.diff(x*z,y),-1,1)
equal('localized-only generation area',area,8*x/S.Integer(15))
check('localized-only generation odd missing',split((x*z).subs(x,0))[1]==0)

# Repeated generic root: both partials vanish and opposite collision values
# prevent any nonzero base Jacobian for every member pair, not just generators.
r=x+1
g=(W-r*r)**2*(W+2)
z=y*g.subs(W,y*y)
p=x+z; q=y*y+x*z*z
J=jac(p,q,x,y)
for sign in [-1,1]:
    equal('repeated root value '+str(sign),g.subs(W,r*r))
    equal('repeated root partial x '+str(sign),S.diff(g,x).subs(W,r*r))
    equal('repeated root partial W '+str(sign),S.diff(g,W).subs(W,r*r))
    equal('repeated root Jacobian '+str(sign),J.subs(y,sign*r),2*sign*r)
families.append({'kind':'repeated_moving_root','r':'(x+1)^2'})

parser=argparse.ArgumentParser()
parser.add_argument('--output',default='independent_results.json')
args=parser.parse_args()
report={'status':'PASS','checks':checks,'family_count':len(families),'families':families,
        'scope':'Finite exact identities and explicit hostile boundaries, not all-embedding or all-pair proof.'}
Path(args.output).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
print(f'PASS: {checks} independent exact checks; {len(families)} declared controls.')
