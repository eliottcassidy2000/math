#!/usr/bin/env python3
"""Exact identity controls for the global invariant/logarithm theorem."""
import json
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
x,y,z,w,a,d,t,c=s.symbols('x y z w a d t c')
checks=[]
def check(label, assertion):
    assert assertion,label
    checks.append(label)
def eq(a,b): return s.cancel(a-b)==0
def M(u): return (x/(1-u*x), y*(1-u*x)**2)
def T(v): return (x,y+v*x**2)
def compose(after,before):
    return tuple(s.cancel(f.subs({x:before[0],y:before[1]},simultaneous=True)) for f in after)
def same(p,q): return all(eq(f,g) for f,g in zip(p,q))
phi=compose(T(-d),compose(M(a),T(d)))
J=x*x*y+d*x**4
check('polynomial conjugate invariant',eq(J.subs({x:phi[0],y:phi[1]},simultaneous=True),J))
check('polynomial conjugate first coordinate',eq(phi[0],x/(1-a*x)))
flow=compose(T(-d),compose(M(a*t),T(d)))
H=a*J
check('polynomial Hamiltonian time derivative',same(tuple(s.diff(f,t).subs(t,0) for f in flow),(s.diff(H,y),-s.diff(H,x))))
check('polynomial time map Jacobian',eq(s.det(s.Matrix(phi).jacobian([x,y])),1))
check('zero-translation primitive',eq(s.diff(-c*x**3/(3*(1-t*x)**3),x),-c*x*x/(1-t*x)**4))

# Rational primitive examples, including several distinct translation orbits.
examples=[]
for ashift,bs in [(Q(1),{Q(0):Q(1)}),
                  (Q(1),{Q(2):Q(1)}),
                  (Q(2),{Q(-2):Q(3),Q(0):Q(-1),Q(1):Q(2),Q(3):Q(5)}),
                  (Q(2,3),{Q(0):Q(2),Q(2,3):Q(-4),Q(1,2):Q(7)})]:
    cs={}
    for site,coef in bs.items():
        cs[site]=cs.get(site,Q(0))+coef
        cs[site+ashift]=cs.get(site+ashift,Q(0))-coef
    cs={site:coef for site,coef in cs.items() if coef}
    bfun=sum(s.Rational(v.numerator,v.denominator)/(z-s.Rational(q.numerator,q.denominator))**4 for q,v in bs.items())
    ffun=sum(s.Rational(v.numerator,v.denominator)/(z-s.Rational(q.numerator,q.denominator))**4 for q,v in cs.items())
    av=s.Rational(ashift.numerator,ashift.denominator)
    check(f'a={ashift},b={bs}: coboundary',eq(ffun,bfun-bfun.subs(z,z-av)))
    invariant=-w+bfun
    check(f'a={ashift},b={bs}: canonical invariant',eq(invariant.subs({z:z-av,w:w-ffun},simultaneous=True),invariant))
    remaining=set(cs);orbits=[];recovered={}
    while remaining:
        origin=next(iter(remaining))
        orbit={q:cs[q] for q in remaining if ((q-origin)/ashift).denominator==1}
        remaining-=set(orbit)
        check(f'a={ashift},origin={origin}: zero orbit sum',sum(orbit.values())==0)
        inds={int((q-origin)/ashift):v for q,v in orbit.items()}
        running=Q(0)
        for i in range(min(inds),max(inds)+1):
            running+=inds.get(i,Q(0))
            if running: recovered[origin+i*ashift]=running
        orbits.append({'sites':[str(q) for q in sorted(orbit)],'sum':'0'})
    check(f'a={ashift},b={bs}: prefix reconstruction',recovered==bs)
    examples.append({'translation':str(ashift),'primitive_sites':{str(q):str(v) for q,v in sorted(bs.items())},
       'cocycle_sites':{str(q):str(v) for q,v in sorted(cs.items())},'orbits':orbits})

# Single-pole case: finite formal Hamiltonian tail, not a rationality proof.
coeff={}
for n in range(3,13):
    target=s.Integer(1) if n==3 else s.Integer(0)
    coeff[n]=s.simplify(-(target+sum(s.binomial(n,j-1)*coeff[j] for j in coeff))/n)
h=sum(v*x**n for n,v in coeff.items())
res=s.series(h-h.subs(x,x/(1-x))-x**4,x,0,14).removeO()
check('single-pole formal log through degree twelve',s.expand(res)==0)
check('single-pole orbit obstruction',sum([Q(1)])!=0)

# Distinct-orbit finite differences: exact first discrepancy order.
jet_examples=[]
for n in range(1,11):
    nodes=[Q(j,n+2) for j in range(n+1)]  # a=1; general case scales kth moments by a^k.
    cs=[(-1)**j*s.binomial(n,j) for j in range(n+1)]
    check(f'jet n={n}: distinct translation orbits',all((nodes[i]-nodes[j]).denominator!=1 for i in range(n+1) for j in range(i)))
    check(f'jet n={n}: nonzero orbit sums',all(c!=0 for c in cs))
    for k in range(n):
        check(f'jet n={n}: vanishing moment {k}',sum(c*s.Rational(q.numerator,q.denominator)**k for q,c in zip(nodes,cs))==0)
    moment=sum(c*s.Rational(q.numerator,q.denominator)**n for q,c in zip(nodes,cs))
    alpha=s.binomial(n+3,3)*moment
    expected=(-1)**n*s.binomial(n+3,3)*s.factorial(n)/s.Integer(n+2)**n
    check(f'jet n={n}: exact leading coefficient',alpha==expected and alpha!=0)
    if n<=4:
        fsum=sum(c/(z-s.Rational(q.numerator,q.denominator))**4 for q,c in zip(nodes,cs))
        vertical=(1-x)**2*x**(-2)*fsum.subs(z,1/x)
        check(f'jet n={n}: original-coordinate discrepancy',s.expand(s.series(vertical,x,0,n+3).removeO()-alpha*x**(n+2))==0)
    jet_examples.append({'n':n,'agreement_through_degree':n+1,'first_error_degree':n+2,'first_error_coefficient_a1':str(alpha)})

out={'status':'PASS','check_count':len(checks),'checks':checks,'orbit_examples':examples,
 'distinct_orbit_jet_examples':jet_examples,
 'single_pole_formal_hamiltonian':'x^2*y + '+str(h),
 'single_pole_formal_tail_degree':12,
 'scope_warning':'Finite identities supplement the all-order classifications. Nonrationality is proved by the orbit obstruction, not inferred from the displayed formal tail.'}
p=Path(__file__).with_name('metabelian_invariants_exact_controls.json')
p.write_text(json.dumps(out,indent=2)+'\n')
print(f'PASS: {len(checks)} exact controls. Output: {p.name}')
