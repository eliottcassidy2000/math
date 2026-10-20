#!/usr/bin/env python3
"""Reviewer-owned checks. Does not import or execute the producer checker."""
import itertools
import json
import math
import sympy as s

checks = []


def test(name, value):
    if not value:
        raise RuntimeError(name)
    checks.append(name)


def eq(name, value):
    test(name, s.expand(value) == 0)


def monoid(gens, bound):
    result = {0}
    for i in range(1, bound + 1):
        if any(i - g in result for g in gens):
            result.add(i)
    return result


# Enumerate the whole binary conductor-eight universe, not subsets of targets.
gap_sets = []
for mask in range(256):
    low = {i for i in range(8) if (mask >> i) & 1}
    if 0 not in low or 7 in low or len(low) != 4:
        continue
    if not all((i in low) != (7-i in low) for i in range(8)):
        continue
    if not all(a+b >= 8 or a+b in low for a in low for b in low):
        continue
    gap_sets.append(tuple(i for i in range(8) if i not in low))
test("all binary semigroup candidates", sorted(gap_sets) ==
     [(1, 2, 3, 7), (1, 2, 4, 7), (1, 3, 5, 7)])

# Enumerate characteristic tuples without pre-filtering to a target semigroup.
# Every strict gcd chain beginning at <=9 has <=3 nontrivial drops, so
# length <=5 is more than enough (and explicitly includes a redundant layer).
patterns = {}
for length in range(2, 6):
    for initial in itertools.product(range(1, 10), repeat=2):
        a, b = initial
        if a % b == 0 or b % a == 0:
            continue
        for tail in itertools.product(range(1, 8), repeat=length-2):
            seq = initial + tail
            ds = [seq[0]] + [math.gcd(*seq[:i]) for i in range(2, length+1)]
            ns = [ds[i] // ds[i+1] for i in range(length-1)]
            if ds[-1] != 1 or any(n <= 1 for n in ns[1:]):
                continue
            if any(ns[i-1]*seq[i] not in monoid(seq[:i], ns[i-1]*seq[i])
                   for i in range(1, length)):
                continue
            if any(seq[i] >= ns[i-2]*seq[i-1] for i in range(2, length)):
                continue
            conductor = 1-seq[0] + sum((ns[i-1]-1)*seq[i]
                                      for i in range(1, length))
            if conductor != 8:
                continue
            values = monoid(seq, 40)
            gaps = tuple(i for i in range(41) if i not in values)
            patterns.setdefault(str(gaps), []).append(seq)
patterns = {k: sorted(v) for k, v in patterns.items()}
test("all bounded characteristic sequences, independent universe",
     patterns == {
         str((1, 3, 5, 7)): [(2, 9), (6, 9, 2), (9, 2), (9, 6, 2)],
         str((1, 2, 4, 7)): [(3, 5), (5, 3)],
         str((1, 2, 3, 7)): [(4, 6, 5), (6, 4, 5)],
     })

t, x = s.symbols("t x")
a, b, c, d, e, f, A, B, L = s.symbols("a b c d e f A B L")
P = t**4+a*t**2+b*t
Q = t**6+A*t**5+c*t**4+B*t**3+e*t**2+L*t
R = s.Poly(Q**2-P**3, t)
eq("degree eleven coefficient derived", R.nth(11)-2*A)
Q = Q.subs(A, 0)
R = s.Poly(Q**2-P**3, t)
rho = R.nth(10)
eq("degree ten cancellation derived", rho-(2*c-3*a))
R = s.Poly(R.as_expr()-rho*P*Q, t)
eq("degree nine coefficient derived", R.nth(9)-(2*B-3*b))
Q = Q.subs(B, 3*b/2)
R = s.Poly(Q**2-P**3-rho*P*Q, t)
sigma = R.nth(8)
eq("degree eight cancellation derived", sigma-(a*c-c*c+2*e))
R = s.Poly(R.as_expr()-sigma*P**2, t)
eq("degree seven coefficient derived", R.nth(7)-(3*a*b/2-2*b*c+2*L))
Q = Q.subs(L, b*(c-3*a/4))
R = s.Poly(Q**2-P**3-rho*P*Q-sigma*P**2, t)
kappa = R.nth(6)
eq("degree six cancellation derived", kappa-(-a**3+a*a*c-a*e-3*b*b/4))
R = s.Poly(R.as_expr()-kappa*Q, t)
D = e-a*c+3*a*a/4
eq("degree five coefficient derived", R.nth(5)+b*D)
lam = R.nth(4)
R = s.Poly(R.as_expr()-lam*P, t)
eq("exact degree three boundary", R.as_expr().subs(e,a*c-3*a*a/4)
   -b**3*t*(3*a+2*t*t)/16)
test("b=0 is nonprimitive even locus",
     all(s.Poly(z.subs(b,0),t).nth(i)==0 for z in (P,Q) for i in (1,3,5)))

u, v = s.symbols("u v")
T = s.Matrix.hstack(s.Matrix([1,0,0,0]), s.Matrix([0,1,0,0]),
                   s.Matrix([-v,-u,1,0]), s.Matrix([0,-v,-u,1]))
test("actual power-basis change determinant", T.det()==1)
for odd in (3,5,7,9,11):
    degrees = [0,6,odd,6+odd]
    test("four residues for remainder "+str(odd),
         set(d%4 for d in degrees)=={0,1,2,3})
    gaps = [n for n in range(50)
            if not any(n>=d and (n-d)%4==0 for d in degrees)]
    test("exact residue codimension "+str(odd),
         len(gaps)==(odd+3)//2)

# Use genuine x-dependent symbolic functions, rather than a manually coded Dx.
aa, bb, cc, dd, ee, ff = [s.Function(z)(x) for z in ("a","b","c","d","e","f")]
sub = dict(zip([a,b,c,d,e,f],[aa,bb,cc,dd,ee,ff]))
p, q = (P+d).subs(sub), (Q+f).subs(sub)
J = s.Poly(s.diff(p,x)*s.diff(q,t)-s.diff(p,t)*s.diff(q,x),t)
eq("entire Jacobian constant coefficient",
   J.nth(0)-bb*((cc-3*aa/4)*s.diff(dd,x)-s.diff(ff,x)))
derivs = [s.diff(z,x) for z in (cc,dd,ee,ff)]
eqs = [J.nth(i).subs(s.diff(bb,x),0) for i in (7,5,3,2)]
M, rhs = s.linear_eq_to_matrix(eqs,derivs)
eq("forced derivative system determinant", M.det()+192*bb)
solution = list(s.linsolve((M,rhs),derivs))
expected = (3*s.diff(aa,x)/2,aa*s.diff(aa,x)/2,
            cc*s.diff(aa,x),ee*s.diff(aa,x)/2)
test("forced derivative system has unique displayed solution",
     len(solution)==1 and all(s.simplify(v-w)==0
                             for v,w in zip(solution[0],expected)))
rules = dict(zip(derivs,solution[0]))
rules[s.diff(bb,x)] = 0
whole = s.factor(J.as_expr().subs(rules))
eq("entire Jacobian after forced solution",
   whole+s.diff(aa,x)*bb*(3*aa*aa-4*aa*cc+6*bb*t+4*ee)/8)
eq("surviving linear Jacobian coefficient",
   s.Poly(whole,t).nth(1)+3*s.diff(aa,x)*bb*bb/4)

# Centering is an actual determinant-one polynomial source substitution.
y = s.symbols("y")
h = s.Function("h")(x)
eq("moving source-centering chain rule",
   (s.diff(p.subs(t,y+h),x)*s.diff(q.subs(t,y+h),y)
    -s.diff(p.subs(t,y+h),y)*s.diff(q.subs(t,y+h),x))
   -J.as_expr().subs(t,y+h))

examples = [
    (t**2,t**9+t),
    (t**3,t**5+t),
    (t**4+t,t**6+s.Rational(3,2)*t**3+t**2),
    (t**6+4*t**2,t**9+6*t**5+6*t)]
for i,(p,q) in enumerate(examples):
    test("realization immersion "+str(i),s.gcd(s.diff(p,t),s.diff(q,t))==1)
p,q = examples[2]
eq("quartic sextic remainder independently expanded",
   q*q-p**3-2*p*p+3*q/4-p-(-t**5+t**3/8-5*t*t/4-t))
p,q = examples[3]
eq("sextic nonic exact algebra recovery",q*q-p**3-8*p-4*t*t)
z = s.symbols("z")
w = z**4+6*z*z+6
test("four distinct nonzero quadratic-model node squares",
     w.subs(z,0)!=0 and s.gcd(w,s.diff(w,z))==1)

# Both orientation swaps preserve the realizations; nonprincipal status is exact.
for p,q in examples:
    m,n = s.degree(p,t),s.degree(q,t)
    test("actual nonprincipal realizing degrees "+str((m,n)),
         m%n!=0 and n%m!=0)

alpha,beta,j,U,V = s.symbols("alpha beta j U V",nonzero=True)
F = U**4-2*U**2+7
p = t*t+alpha*x+beta
q = j*t/alpha+F.subs(U,p)
eq("general quadratic Keller control",
   s.diff(p,x)*s.diff(q,t)-s.diff(p,t)*s.diff(q,x)-j)
inv_t = alpha*(V-F)/j
inv_x = (U-inv_t**2-beta)/alpha
eq("quadratic explicit inverse first coordinate",
   p.subs({x:inv_x,t:inv_t}, simultaneous=True)-U)
eq("quadratic explicit inverse second coordinate",
   q.subs({x:inv_x,t:inv_t}, simultaneous=True)-V)

print(json.dumps({
    "status":"PASS", "count":len(checks), "checks":checks,
    "gap_sets":sorted(gap_sets), "characteristic_patterns":patterns,
    "jacobian_derivative_matrix":str(M),
    "jacobian_derivative_matrix_determinant":str(s.factor(M.det())),
    "full_jacobian_after_solution":str(whole),
},indent=2,sort_keys=True))
