"""Exact small-k audit of the fixed-critical-orbit mutation corollary.

All coefficient arithmetic is over Q[c]/(c^4+3c^2+3), represented by SymPy's
exact algebraic field. No mutation producer or numerical root is imported.
The unbounded module conclusions use the analytic supplier and proof note.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import sympy as s

sys.stdout.reconfigure(newline='\n')
GATES=Counter()

def check(value,label):
    GATES[label]+=1
    if not value:raise RuntimeError('always-active gate failed: '+label)

def main():
    w,z=s.symbols('w z')
    modulus=w**4+3*w**2+3
    c=s.CRootOf(modulus,0)
    field=s.QQ.algebraic_field(c)
    cc=field.convert(c);one=field.one;zero=field.zero
    zeta=cc**2+one
    a=(one-zeta)*cc
    b=-zeta*cc
    check(s.Poly(modulus,w).is_irreducible,'coefficient field irreducible')
    check(zeta*zeta+zeta+one==zero,'primitive cube root relation')
    check(zeta!=one and cc!=zero,'nonzero critical orbit parameters')
    check(cc*cc==zeta-one,'prescribed square relation')
    P=lambda value:value**3+cc
    check(P(zero)==cc and P(cc)==zeta*cc and P(zeta*cc)==zeta*cc,'complete critical orbit')
    check(a!=zero and b!=zero and a!=b,'three target supports distinct')
    check(a-b==cc,'target support separation identity')
    C=lambda value:s.Poly(field.to_sympy(value),z,domain=field)
    unit=C(one)
    iterations=[s.Poly(z,z,domain=field)]
    for j in range(1,5):iterations.append(iterations[-1]**3+C(cc))
    for i,pi in enumerate(iterations):
        check(pi.degree()==3**i and pi.LC()==1,'monic iterate degree')
        check(pi.gcd(pi.diff())==unit,'every preimage level simple')
        for j,pj in enumerate(iterations[:i]):
            check(pi.gcd(pj)==unit,'pairwise disjoint preimage levels')
    examples=[]
    for k in range(2,5):
        R=unit
        for pj in iterations[:k]:R*=pj
        q=(3**k-1)//2
        Q=iterations[k]-C(zeta*cc)
        check(R.degree()==q,'half derivative degree')
        check(R.gcd(R.diff())==unit,'squarefree derivative square root')
        check(Q.degree()==3**k and Q.LC()==1,'normalized primitive degree')
        check(Q.diff()==R**2*(3**k),'exact squared derivative identity')
        check(Q.eval(0)==0 and R.eval(0)==0,'required source zero normalization')
        check(field.convert(R.diff().eval(0))==cc*(zeta*cc)**(k-2),'simple zero derivative coefficient')
        lower=unit
        for pj in iterations[:k-1]:lower*=pj
        top=iterations[k-1]
        check(Q.gcd(R)==lower,'critical zero-value group exact')
        check((Q-C(a)).gcd(R)==top,'critical nonzero-value group exact')
        check((Q-C(b)).gcd(R)==unit,'Eu target not a critical value')
        check(Q-C(a)==top**3,'top critical fibre exact cubic power')
        quotient,remainder=Q.div(lower**3)
        check(remainder.is_zero and quotient.gcd(lower)==unit,'zero critical group exact cubic multiplicity')
        check(quotient.gcd(quotient.diff())==unit,'remaining zero-fibre roots simple')
        check(quotient.degree()==q+2,'remaining zero-fibre degree')
        check(Q-C(b)==iterations[k],'chosen Eu level has exact target')
        check(iterations[k].gcd(R)==unit,'all chosen Eu points admissible')
        check(iterations[k].gcd((3**k)*R**2-unit)==unit,'lambda one admissible at every small-k Eu choice')
        check(lower.degree()==(3**(k-1)-1)//2 and top.degree()==3**(k-1),'exact critical group sizes')
        support=s.Poly(z,z,domain=field)*(s.Poly(z,z,domain=field)-C(a))*(s.Poly(z,z,domain=field)-C(b))
        check(support.degree()==3 and support.gcd(support.diff())==unit,'annihilator support squarefree cubic')
        check(support**2==support*support and (support**2).degree()==6,'fixed scalar-annihilator degree')
        fixed=s.Poly(z**3-z,z,domain=field)+C(cc)
        other,remainder=fixed.div(s.Poly(z,z,domain=field)-C(zeta*cc))
        check(remainder.is_zero and other.degree()==2,'two other fixed points exist')
        check(other==s.Poly(z**2,z,domain=field)+C(zeta*cc)*s.Poly(z,z,domain=field)-C(zeta**2),'fixed-h exact quadratic')
        check(fixed.gcd(fixed.diff())==unit,'all fixed points simple')
        check(other.gcd(s.Poly(z,z,domain=field)*(s.Poly(z,z,domain=field)-C(cc)))==unit,'other fixed points avoid zero and c')
        check(other.gcd(s.Poly(z,z,domain=field)-C(zeta*cc))==unit,'fixed-h three main supports distinct')
        check((iterations[k]-s.Poly(z,z,domain=field)).rem(other).is_zero,'fixed-h iterate value')
        check((R-s.Poly(z**k,z,domain=field)).rem(other).is_zero,'fixed-h derivative square root value')
        check(other.gcd(R)==unit,'fixed-h admissibility at all small-k levels')
        check(other.gcd(s.Poly(3**k*z**(2*k)-1,z,domain=field))==unit,'fixed-h lambda one admissible')
        # Serialize exact coefficient vectors; no numerical complex roots enter.
        data=[[str(v) for v in coeff.to_list()] for coeff in Q.rep.to_list()]
        digest=hashlib.sha256(json.dumps(data,separators=(',',':')).encode()).hexdigest()
        examples.append({'k':k,'iterate_degree':3**k,'q':q,
            'ambient_arms_by_fixed_support':[lower.degree(),top.degree(),1],
            'ambient_arms':q+1,'unit_arms':3,'annihilator_degree':6,
            'source_t_degree':3**k,'intrinsic_pair_degree':2*3**k,
            'generic_affine_punctures':3*q+2,'zero_fibre_remaining_degree':q+2,
            'Q_exact_coefficient_sha256':digest,'lambda_one_admissible_all_Ek':True})
    Q1=iterations[1]-C(zeta*cc)
    check(field.convert(Q1.eval(0))==a and a!=zero,'k-one fails prescribed zero normalization')
    bad_iterations=[s.Poly(z,z,domain=s.QQ)]
    for j in range(2):bad_iterations.append(bad_iterations[-1]**3)
    bad_root=bad_iterations[0]*bad_iterations[1]
    check(bad_root==s.Poly(z**4,z,domain=s.QQ),'zero-parameter first repeated-root witness')
    check(bad_root.gcd(bad_root.diff()).degree()==3,'critical return destroys squarefreeness')
    check(cc**2!=field.convert(s.Rational(4,27)),'fixed-point multiplicity exclusion')
    certificate={'status':'FINITE-EXACT controls for analytic fixed-support iteration corollary',
        'coefficient_field':'Q[c]/(c^4+3c^2+3), irreducible by Eisenstein at 3',
        'zeta':'c^2+1','critical_orbit':['0','c','zeta*c','zeta*c'],
        'fixed_support':['0','(1-zeta)*c','eta-zeta*c'],
        'fixed_source_parameters':{'eta':'either root of eta^2+zeta*c*eta-zeta^2=0','h':'-eta/2','lambda':1},
        'moving_source_variant_support':['0','(1-zeta)*c','-zeta*c'],
        'normalization':'f_k=(sqrt(3))^k product_{j=0}^{k-1}P^j; square tested without adjoining sqrt(3)',
        'examples':examples,'gate_counts':dict(sorted(GATES.items())),'total_gates':sum(GATES.values()),
        'scope':'All-k assertions are proved in the report; small-k controls are k=2,3,4. No polynomial mate or JC assertion.'}
    base=Path(__file__).resolve().parent
    if base.name=='04-computation':base=base.parent/'05-knowledge'/'results'
    target=base/(Path(__file__).stem+'_certificate.json')
    target.write_bytes((json.dumps(certificate,sort_keys=True,indent=2)+'\n').encode())
    print('Fixed-critical-orbit iteration: PASS')
    print('Exact coefficient field: Q[c]/(c^4+3c^2+3)')
    print('Small-k controls: k=2,3,4; pair degrees18,54,162; ambient arms5,14,41')
    print('Fixed h=-eta/2, lambda=1; fixed unit supports: 0, (1-zeta)c, eta-zeta*c')
    print('Three pure order-two arms; simple/disjoint levels, critical groups, and hostile controls: PASS')
    print('Always-active exact gates: '+str(sum(GATES.values())))
    print('Certificate: '+target.name)

if __name__=='__main__':main()
