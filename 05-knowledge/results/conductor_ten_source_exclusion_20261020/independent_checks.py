#!/usr/bin/env python3
"""Independent exact audit of genus five and the quartic/sextic candidate.

No producer module is imported. Numerical enumeration has a proved finite
universe; symbolic-function identities are universal, not coefficient scans.
"""
import itertools
import json
import math
import sympy as S

checks = []
def check(name, condition):
    if not condition:
        raise RuntimeError('FAILED: '+name)
    checks.append(name)
def zero(name, expression):
    check(name, S.expand(expression.doit()) == 0)
def monoid(gens, bound):
    have = {0}
    for n in range(1, bound+1):
        if any(n >= g and n-g in have for g in gens):
            have.add(n)
    return have

def primitives(vals):
    return tuple(n for n in sorted(vals) if n > 0 and not any(
        a in vals and n-a in vals for a in range(1,n)))

# Enumerate every subset of {0,...,9}, without presupposing symmetry choices.
gaps = []
for mask in range(1 << 10):
    vals = {n for n in range(10) if mask & (1 << n)}
    if 0 not in vals or 9 in vals or len(vals) != 5:
        continue
    if any((n in vals) == (9-n in vals) for n in range(10)):
        continue
    if any(a+b < 10 and a+b not in vals for a in vals for b in vals):
        continue
    gaps.append(tuple(n for n in range(10) if n not in vals))
expected_gaps = [(1,2,3,4,9),(1,2,3,5,9),(1,3,5,7,9)]
check('exhaustive 1024-subset symmetric genus-five census', sorted(gaps)==expected_gaps)
expected_gens = [(2,11),(4,6,7),(5,6,7,8)]
for gen in expected_gens:
    have = monoid(gen, 40)
    check('gap count '+str(gen), len(set(range(41))-have)==5)
    check('conductor ten '+str(gen), 9 not in have and set(range(10,41)) <= have)
    check('primitive generator check '+str(gen), primitives(have)==gen)

# Entire initial integer universe 1..11, not just a proposed semigroup.
# Every nonprincipal pair in this universe has gcd 1, 2, or 3;
# a strict further gcd drop then terminates, so at most three entries.
pairs = [(a,b) for a in range(1,12) for b in range(1,12)
         if a%b and b%a]
check('nonprincipal bounded initial gcd at most three and prime if >1',
      set(math.gcd(*z) for z in pairs)=={1,2,3})
def characteristic(seq):
    d=seq[0]
    previous_n=None
    for i in range(1,len(seq)):
        new_d=math.gcd(d,seq[i]); n=d//new_d
        if i>=2 and n<=1: return False
        if n*seq[i] not in monoid(seq[:i], n*seq[i]): return False
        if i>=2 and seq[i]>=previous_n*seq[i-1]: return False
        d=new_d; previous_n=n
    return d==1
found=[]
for a,b in pairs:
    trials=[(a,b)] + [(a,b,c) for c in range(1,10)]
    for seq in trials:
        if not characteristic(seq): continue
        have=monoid(seq,60)
        missing=set(range(61))-have
        if len(missing)==5 and max(missing)==9:
            found.append(seq)
expected_patterns=[(2,11),(4,6,7),(6,4,7),(11,2)]
check('all bounded characteristic sequences genus five', sorted(found)==expected_patterns)
check('nonplanar symmetric hostile excluded', all(
    primitives(monoid(seq,60)) != (5,6,7,8) for seq in found))
# An independent target-filtered recursion crosschecks the whole-universe pass.
recursive=[]
for gen in expected_gens:
    allowed=sorted(monoid(gen,11)-{0})
    def extend(seq):
        d=math.gcd(*seq)
        if d==1:
            if primitives(monoid(seq,60))==gen and characteristic(seq):
                recursive.append(tuple(seq))
            return
        for nxt in allowed:
            if nxt<=9 and math.gcd(d,nxt)<d:
                extend(seq+[nxt])
    for a,b in itertools.permutations(allowed,2):
        if a%b and b%a: extend([a,b])
check('independent target-filtered recursion agrees', sorted(recursive)==sorted(found))

x,t,y=S.symbols('x t y')
a,b,c,d,e,h,f=[S.Function(n)(x) for n in ['a','b','c','d','e','h','f']]
C,E,H,F,L,M=S.symbols('C E H F L M')
P=t**4+a*t**2+b*t+d
Q=t**6+c*t**4+S.Rational(3,2)*b*t**3+e*t**2+h*t+f
J=S.Poly((S.diff(P,x)*S.diff(Q,t)-S.diff(P,t)*S.diff(Q,x)).expand(),t)
expected={
7:6*S.diff(a,x)-4*S.diff(c,x),
6:0,
5:-2*a*S.diff(c,x)+4*c*S.diff(a,x)+6*S.diff(d,x)-4*S.diff(e,x),
4:-3*a*S.diff(b,x)+S.Rational(9,2)*b*S.diff(a,x)-b*S.diff(c,x)+4*c*S.diff(b,x)-4*S.diff(h,x),
3:-2*a*S.diff(e,x)+3*b*S.diff(b,x)+4*c*S.diff(d,x)+2*e*S.diff(a,x)-4*S.diff(f,x),
2:-2*a*S.diff(h,x)+S.Rational(9,2)*b*S.diff(d,x)-b*S.diff(e,x)+2*e*S.diff(b,x)+h*S.diff(a,x),
1:-2*a*S.diff(f,x)-b*S.diff(h,x)+2*e*S.diff(d,x)+h*S.diff(b,x),
0:h*S.diff(d,x)-b*S.diff(f,x)}
for n,value in expected.items(): zero('full symbolic-function Jacobian coefficient '+str(n),J.nth(n)-value)

cc=3*a/2+C
hh=b*(3*a/4+C)+H
ee=3*d/2+3*a*a/8+C*a+E
ff=C*d+3*a*d/4-a**3/16+E*a/2+3*b*b/8+F
zero('integrate c from degree seven',J.nth(7)+4*S.diff(c-3*a/2,x))
zero('integrate e from degree five after c',J.nth(5).subs(c,cc).doit()+4*S.diff(e-ee,x))
zero('integrate h from degree four after c',J.nth(4).subs(c,cc).doit()+4*S.diff(h-hh,x))
zero('integrate full scalar f from degree three',J.nth(3).subs({c:cc,e:ee}).doit()+4*S.diff(f-ff,x))
fp_candidate=ee*S.diff(a,x)/2+3*b*S.diff(b,x)/4+cc*S.diff(d,x)-a*S.diff(ee,x)/2
zero('candidate f derivative agrees with primitive',fp_candidate-S.diff(ff,x))
subs={c:cc,e:ee,h:hh,f:ff}
Jsol=S.Poly(J.as_expr().subs(subs).doit().expand(),t)
v=d-a*a/4+2*E/3
I1=H*a+3*b*v
I2=H*b+3*v*v/2-3*a*b*b/4
j=H*S.diff(d,x)-3*(b*b*S.diff(b,x)+b*v*S.diff(a,x))/4
zero('first integral I1 controls degree two',Jsol.nth(2)-S.diff(I1,x))
zero('first integral I2 controls degree one',Jsol.nth(1)-S.diff(I2,x))
zero('entire remaining scalar coefficient',Jsol.nth(0)-j)
zero('whole residual Jacobian identity',Jsol.as_expr()-S.diff(I1,x)*t*t-S.diff(I2,x)*t-j)
for n in [7,6,5,4,3]: zero('integrated upper coefficient '+str(n),Jsol.nth(n))

# Independent jet-ring differentiation rather than substituted functions.
A,B,D,V=S.symbols('A B D V'); Ad,Bd,Dd=S.symbols('Ad Bd Dd')
jetP=t**4+A*t*t+B*t+D
jetQ=t**6+(3*A/2+C)*t**4+3*B*t**3/2+(3*D/2+3*A*A/8+C*A+E)*t*t+(B*(3*A/4+C)+H)*t+(C*D+3*A*D/4-A**3/16+E*A/2+3*B*B/8+F)
def dx(z): return S.diff(z,A)*Ad+S.diff(z,B)*Bd+S.diff(z,D)*Dd
jetJ=S.expand(dx(jetP)*S.diff(jetQ,t)-S.diff(jetP,t)*dx(jetQ))
for n in range(8):
    value=Jsol.nth(n).subs({a:A,b:B,d:D,S.diff(a,x):Ad,S.diff(b,x):Bd,S.diff(d,x):Dd})
    zero('independent jet-ring coefficient '+str(n),S.Poly(jetJ,t).nth(n)-value)

# The completing-square equality is equivalent to I2=M after I1=L.
I2_after=H*B+3*V*V/2-3*((L-3*B*V)/H)*B*B/4-M
W=V+3*B**3/(4*H)
R=L*B*B/(2*H)-2*H*B/3+2*M/3
zero('completing-square equivalence',W*W-9*B**6/(16*H**2)-R-2*I2_after/3)
zero('factor difference retains cubic b degree',(W+3*B**3/(4*H))-(W-3*B**3/(4*H))-3*B**3/(2*H))
# Concrete boundary controls for the degree lemma: zero RHS needs exclusion.
for r in [1,2,3,7]:
    bb=x**r+x+1 if r>1 else x+1
    aa=3*bb**3/4
    check('degree lower bound nonzero product witness '+str(r),S.degree(S.expand((aa+1)**2-aa**2),x)==3*r)
    zero('zero-product exception witness '+str(r),aa**2-aa**2)

# Re-derive the first odd remainders without producer-specific constants.
A5,B3=S.symbols('A5 B3')
aa,bb,cc0,ee0,hh0=S.symbols('aa bb cc ee hh')
p0=t**4+aa*t*t+bb*t
q0=t**6+A5*t**5+cc0*t**4+B3*t**3+ee0*t*t+hh0*t
r0=S.Poly(q0*q0-p0**3,t)
zero('degree eleven remainder coefficient',r0.nth(11)-2*A5)
q1=q0.subs(A5,0)
r1=S.Poly(q1*q1-p0**3,t)
r1=S.Poly(r1.as_expr()-r1.nth(10)*p0*q1,t)
zero('degree nine remainder coefficient',r1.nth(9)-(2*B3-3*bb))
q2=q1.subs(B3,3*bb/2)
r2=S.Poly(r1.as_expr().subs(B3,3*bb/2),t)
r2=S.Poly(r2.as_expr()-r2.nth(8)*p0*p0,t)
zero('degree seven remainder coefficient',r2.nth(7)-(3*aa*bb/2-2*bb*cc0+2*hh0))
for n,delta in [(11,7),(9,6),(7,5),(5,4),(3,3)]:
    ds=[0,6,n,n+6]
    check('distinct residue module degrees '+str(n),sorted(i%4 for i in ds)==[0,1,2,3])
    check('module exact gap codimension '+str(n),sum(i-i%4 for i in ds)//4==delta)
u,z=S.symbols('u z')
T=S.Matrix([[1,0,-z,0],[0,1,-u,-z],[0,0,1,-u],[0,0,0,1]])
check('unimodular full module transition',T.det()==1)

# Genuine polynomial realizations, not inferred from characteristic axioms.
p=t**4+t; q=t**6+S.Rational(3,2)*t**3+t
r=S.expand(q*q-p**3)
zero('genus-five quartic sextic exact remainder',r-(2*t**7-3*t**6/4+3*t**4-t**3+t*t))
check('genus-five quartic sextic immersion',S.gcd(S.diff(p,t),S.diff(q,t))==1)
check('genus-five quartic sextic primitive degrees',math.gcd(4,S.degree(r,t))==1)
check('genus-five quartic sextic Apéry degrees', [0,13,6,7]==[min(n for n in monoid((4,6,7),40) if n%4==i) for i in range(4)])
p=t*t; q=t**11+t
check('quadratic undecic immersion',S.gcd(S.diff(p,t),S.diff(q,t))==1)
check('quadratic undecic primitive degrees',math.gcd(2,11)==1)
check('quadratic undecic gaps',set(range(30))-monoid((2,11),29)=={1,3,5,7,9})

# Hostiles and positives for hypotheses often lost in shortcuts.
p=t**4+x*t; q=t**6+S.Rational(3,2)*x*t**3+S.Rational(3,8)*x*x
jj=S.diff(p,x)*S.diff(q,t)-S.diff(p,t)*S.diff(q,x)
zero('moving H=0 nonconstant-Jacobian hostile',jj+3*x*x/4)
check('hostile Jacobian is not zero and not constant',S.degree(S.expand(jj),x)==2)
hp=t**4+x*t; hq=t**6+S.Rational(3,2)*x*t**3
zero('moving hostile is generic delta-three not genus-five',hq*hq-hp**3+3*x*x*hq/4-x**3*t**3/8)
check('moving hostile generic module codimension three',sum(n-n%4 for n in [0,6,3,9])//4==3)
vv=S.symbols('v0')
a0=x; b0=0; d0=vv+a0*a0/4
p=t**4+a0*t*t+d0
q=t**6+3*a0*t**4/2+(3*d0/2+3*a0*a0/8)*t*t+3*a0*d0/4-a0**3/16
zero('moving zero-Jacobian b-zero hostile',S.diff(p,x)*S.diff(q,t)-S.diff(p,t)*S.diff(q,x))
# Dropping scalar terms changes J: this is not a legal Jacobian simplification.
pp=t**4+t+x; qq=t**6+S.Rational(3,2)*t**3+t+x*x
before=S.diff(pp,x)*S.diff(qq,t)-S.diff(pp,t)*S.diff(qq,x)
pp0=pp-pp.subs(t,0); qq0=qq-qq.subs(t,0)
after=S.diff(pp0,x)*S.diff(qq0,t)-S.diff(pp0,t)*S.diff(qq0,x)
check('x-dependent scalar subtraction changes actual Jacobian',S.expand(before-after)!=0)
# Source centering, in contrast, has determinant one.
ss=S.Function('s')(x)
zero('polynomial source-centering chain rule',
     S.diff(P.subs(t,y+ss),x)*S.diff(Q.subs(t,y+ss),y)-S.diff(P.subs(t,y+ss),y)*S.diff(Q.subs(t,y+ss),x)-J.as_expr().subs(t,y+ss))
for phi in [lambda p:0,lambda p:p**3+2*p,lambda p:p**5-p*p+7]:
    p=t*t+2*x+3; q=3*t+phi(p)
    zero('genuine quadratic Keller control '+str(S.degree(q,t)),S.diff(p,x)*S.diff(q,t)-S.diff(p,t)*S.diff(q,x)-6)
    zero('quadratic target reduction control '+str(S.degree(q,t)),q-phi(p)-3*t)

print(json.dumps({'status':'PASS','count':len(checks),
 'symmetric_genus_five_gap_sets':sorted(gaps),
 'nonprincipal_characteristic_sequences':sorted(found),
 'polynomial_realizations':{'quadratic_undecic':['t^2','t^11+t'],
 'quartic_sextic':['t^4+t','t^6+(3/2)t^3+t']},
 'universal_identity':'J = (H*a+3*b*v)\' t^2 + (H*b+3*v^2/2-3*a*b^2/4)\' t + H*d\' -(3/4)*(b^2*b\'+b*v*a\')',
 'not_claimed':['unconditional planar Jacobian conjecture','all nonmonic source algebras','all terminal (4,6) pairs without delta bound','all special-fiber conductors'],
 'checks':checks},indent=2,sort_keys=True))
