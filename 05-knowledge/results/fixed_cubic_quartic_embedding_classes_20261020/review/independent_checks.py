#!/usr/bin/env python3
"""Independent exact audit, with no imports of producer code.

Main independent route: the C[u]-index of C[p,q] in C[y] is the
u-degree of det(1,q,q^2,q^3) in C[u,y]/(p(y)-u).  This uses neither
producer's conductor formula nor its residue-separated basis changes.
"""
from itertools import permutations
import hashlib
import json
import math
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

checks = {}
def require(name, ok):
    if not bool(ok):
        raise RuntimeError(name)
    checks[name] = True

R,a,b,c,B,G,u = ring('a,b,c,B,G,u', QQ)
ONE = [R.one,R.zero,R.zero,R.zero]
Y = [R.zero,R.one,R.zero,R.zero]
H = [R.zero,a,R.zero,R.one]

def mul(v,w, beta):
    # Modulus y^4 + beta*y^3 + b*y^2 + (c+beta*a)*y-u.
    modulus = [-u,c+beta*a,b,beta]
    out = [R.zero]*7
    for i,vi in enumerate(v):
        for j,wj in enumerate(w):
            out[i+j] += vi*wj
    for d in range(6,3,-1):
        lead = out[d]
        for j in range(4):
            out[d-4+j] -= lead*modulus[j]
    return out[:4]

def determinant_for_generator(q, beta):
    cols=[ONE]
    for j in range(3):
        cols.append(mul(cols[-1],q,beta))
    value=R.zero
    for perm in permutations(range(4)):
        parity=sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        term=R.one
        for j in range(4):
            term *= cols[j][perm[j]]
        value += (-1 if parity % 2 else 1)*term
    return value

H2=mul(H,H,B)
Q=[H2[i]+G*H[i] for i in range(4)]
det=determinant_for_generator(Q,B)
expr=det.as_expr()
aa,bb,cc,be,ga,uu=s.symbols('a b c B G u')
require('universal_index:degree_at_most_seven', s.degree(expr,uu)==7)
require('universal_index:top_coefficient', s.expand(expr).coeff(uu,7)==9*be**2)
e0=s.expand(expr.subs(be,0))
require('beta_zero_index:degree_at_most_six', s.degree(e0,uu)==6)
require('beta_zero_index:top_coefficient', s.expand(e0.coeff(uu,6)+(2*ga-3*cc)**2)==0)
candidate=s.expand(e0.subs(ga,s.Rational(3,2)*cc))
standard=determinant_for_generator(H,R.zero).as_expr()
require('standard_index:degree_three', s.degree(standard,uu)==3)
require('standard_index:unit_leading_coefficient', s.expand(standard).coeff(uu,3)==-1)
L=8*aa**3-28*aa**2*bb+32*aa*bb**2-12*bb**3+(8*aa-12*bb)*uu-cc**2
require('candidate_index:exact_factorization', s.expand(candidate+cc**2*L**2*standard/64)==0)
require('candidate_index:degree_five_lead',s.expand(candidate.coeff(uu,5)-cc**2*(2*aa-3*bb)**2/4)==0)
require('candidate_index:c_zero_nonprimitive',s.expand(candidate.subs(cc,0))==0)
require('candidate_index:b_two_a_thirds',s.expand(candidate.subs(bb,2*aa/3)+cc**6*standard.subs(bb,2*aa/3)/64)==0)

# Test a disjoint exact rational parameter box. These specialize the universal
# module-index calculation, rather than using a derivative conductor formula.
sample_records=[]
for av in (-2,0,3):
    for bv in (-1,s.Rational(2,3)*av,4):
        for cv in (0,-2,5):
            spec={aa:av,bb:bv,cc:cv}
            value=s.expand(candidate.subs(spec))
            degree=s.degree(value,uu)
            expected = -s.oo if cv==0 else 3 if 2*av==3*bv else 5
            label=f'candidate_box:{len(sample_records)}'
            require(label, degree==expected)
            sample_records.append({'a':str(av),'b':str(bv),'c':str(cv),'index_degree':str(degree)})

# Direct polynomial identities separately reproduce the refined equality boundary.
y,x=s.symbols('y x')
h=y**3+aa*y
k=y**4+bb*y**2+cc*y
q=h*h+3*cc*h/2
rho=4*aa-3*bb
r=s.Poly(s.expand(q*q-k**3-rho*k*q+2*aa*(aa-bb)*k*k),y)
require('candidate_residue_basis:degree_at_most_seven',r.degree()==7)
require('candidate_residue_basis:degree_seven_coefficient',s.expand(r.nth(7)+cc*(2*aa-3*bb)/2)==0)
kb=k.subs(bb,2*aa/3)
identity=q*q-kb**3-2*aa*kb*q+2*aa**2*kb**2/3+(2*aa**3/27+3*cc**2/4)*q-(aa**4/9+aa*cc**2)*kb-cc**3*h/8
require('candidate_generation:universal_identity',s.expand(identity)==0)
require('candidate_generation:c_zero_even',s.expand(q.subs(cc,0).subs(y,-y)-q.subs(cc,0))==0)

# Recursive, axiom-aware sequence enumeration, rather than producer's product loop.
semigroup={3*i+4*j for i in range(30) for j in range(30)}
require('semigroup:gaps', [i for i in range(20) if i not in semigroup]==[1,2,5])
def in_generated(value,seq):
    reached={0}
    for n in range(value+1):
        if n in reached:
            reached.update(n+d for d in seq if n+d<=value)
    return value in reached
patterns=[]
def extend(seq,ds):
    d=ds[-1]
    if d==1:
        if {3,4}.issubset(seq):
            ns=[ds[i-1]//ds[i] for i in range(1,len(ds))]
            cond=1-seq[0]+sum((ns[i-1]-1)*seq[i] for i in range(1,len(seq)))
            if cond==6:patterns.append(seq)
        return
    for nxt in (3,4):
        newd=math.gcd(d,nxt)
        if newd==d:continue
        n=d//newd
        if not in_generated(n*nxt,seq):continue
        if nxt >= (ds[-2]//d)*seq[-1]:continue
        extend(seq+[nxt],ds+[newd])
for first in (3,4,6,7):
    for second in (3,4,6,7):
        if first<=second or first%second==0:continue
        d=math.gcd(first,second)
        extend([first,second],[first,d])
require('characteristic_sequences:exact_two',sorted(patterns)==[[4,3],[6,4,3]])

# Sharp witness, conormal determinant, and collision-period comparison.
h0=y**3;k0=y**4+y;z=8*y**6+12*y**3+3
require('sharp:reconstruction',s.expand(z*z-64*k0**3-9-8*h0)==0)
U,V=s.symbols('U V')
F=V**3-U*(U+1)**3
require('sharp:implicit_equation',s.expand(F.subs({U:h0,V:k0}))==0)
require('sharp:singular_gradient',F.subs({U:-1,V:0})==0 and s.diff(F,U).subs({U:-1,V:0})==0 and s.diff(F,V).subs({U:-1,V:0})==0)
D=s.Matrix([V,8*U**2+12*U+3]).jacobian([U,V]).det().subs({U:-1,V:0})
require('sharp:conormal_determinant_four',D==4)
def period_remainder(P,Q):
    return s.rem(s.integrate(s.expand(P*s.diff(Q,y)),y),y**3+1,y)
period1=period_remainder(h0,k0)
period2=period_remainder(k0,z)
require('sharp:period_one',period1==9*y/28)
require('sharp:period_two',period2==-9*y/5)
require('sharp:ratio_mismatch',s.simplify(period2/period1)==-s.Rational(28,5) and -s.Rational(28,5)!=D)

# A new entire generating cylinder, composition of two triangular maps.
P=h0+x*k0**2
Q=k0+(x*x+1)*P**2
Pr=s.expand(P);Qr=s.expand(Q)
require('cylinder:entire_inverse_k',s.expand(Qr-(x*x+1)*Pr**2-k0)==0)
require('cylinder:entire_inverse_h',s.expand(Pr-x*(Qr-(x*x+1)*Pr**2)**2-h0)==0)
remainder=period_remainder(Pr,Qr)
require('cylinder:composite_automorphism_period',s.diff(remainder,y)==s.Rational(9,28))
require('cylinder:derivative_period_zero',s.diff(remainder,x,y)==0)
# Generic generation is insufficient: exact bad family with variable period.
Pbad=x*h0;Qbad=k0
require('cylinder:generic_only_period',period_remainder(Pbad,Qbad)==9*x*y/28)
require('cylinder:generic_only_bad_fiber',Pbad.subs(x,0)==0)
# Ramification control uses arbitrary elements, not generation.
Pcus=y**3+x*y**8+x**3
Qcus=y**4+x*x*y**9+x
Jcus=s.diff(Pcus,x)*s.diff(Qcus,y)-s.diff(Pcus,y)*s.diff(Qcus,x)
require('ramification:both_y_derivatives_zero',s.diff(Pcus,y).subs(y,0)==0 and s.diff(Qcus,y).subs(y,0)==0)
require('ramification:jacobian_zero',Jcus.subs(y,0)==0)
# Cubic/quartic cusp parametrization is injective: common nonzero point's
# y value is k/h; zero in h forces y=0. Record exact reconstruction away from 0.
require('ramification:cusp_rational_inverse',s.cancel(y**4/y**3)==y)

# Exact-count refinement: eliminate the two collision endpoints using their sum.
S,T,d = s.symbols('S T d')
collision_cubic=S**3+4*aa*S/3-cc
endpoint_quadratic=y**2-S*y+S**2+aa
require('exact_count:h_common_value',s.rem(h+S**3+aa*S,endpoint_quadratic,y)==0)
# More directly, kb modulo the quadratic has y-coefficient equal to minus the cubic.
kremainder=s.Poly(s.rem(kb,endpoint_quadratic,y),y)
require('exact_count:k_collision_cubic',s.expand(kremainder.nth(1)+collision_cubic)==0)
require('exact_count:h_common_mod_cubic',s.rem(-S**3-aa*S-(-cc+aa*S/3),collision_cubic,S)==0)
require('exact_count:separation_resultant',s.factor(s.resultant(collision_cubic,-3*S**2-4*aa,S))==-27*cc**2)
# A triple root of the monic cubic must have root0 from its quadratic term.
require('exact_count:triple_root_coefficients',s.Poly(s.expand((S-T)**3-collision_cubic),S).nth(2)==-3*T)
require('exact_count:triple_root_a_zero',s.Poly(s.expand((S-T)**3-collision_cubic),S).nth(1).subs(T,0)==-4*aa/3)
require('exact_count:triple_root_c_zero',s.Poly(s.expand((S-T)**3-collision_cubic),S).nth(0).subs(T,0)==cc)
forced_det=-2*(-cc+aa*S/3)-3*cc/2
require('exact_count:forced_determinant_difference',s.expand(forced_det-forced_det.subs(S,T))==-2*aa*S/3+2*aa*T/3)
require('exact_count:standard_rescaling_h',s.expand(h.subs({aa:0,y:d*T})-d**3*T**3)==0)
require('exact_count:standard_rescaling_k',s.expand(k.subs({bb:0,cc:d**3,y:d*T})-d**4*(T**4+T))==0)
require('exact_count:second_rescaling_q',s.expand(q.subs({aa:0,cc:d**3,y:d*T})-d**6*(T**6+3*T**3/2))==0)

print(json.dumps({'status':'PASS: independent exact controls; not a substitute for the all-embedding source theorem','check_count':len(checks),'checks':checks,'characteristic_patterns':patterns,'candidate_module_index_ratio':str(s.factor(-cc**2*L**2/64)),'generic_module_determinant_terms':len(det),'generic_module_determinant_sha256':hashlib.sha256(str(expr).encode()).hexdigest(),'candidate_box':sample_records,'period_remainders':{'standard':str(period1),'second':str(period2),'composite_entire':str(remainder),'generic_only':str(period_remainder(Pbad,Qbad))}},sort_keys=True,indent=2))
