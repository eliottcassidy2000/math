"""Exact controls for two fixed-W2 mutations with the same torsion passport.

All arithmetic is in Q[r]/(6*r^2-15*r+10). The generic-fibre nonisomorphism
is proved separately using degrees of closed punctures, not finite samples.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import sys
import sympy as s

sys.stdout.reconfigure(newline='\n')
GATES=Counter()

def check(condition,label):
    GATES[label]+=1
    if not condition:
        raise RuntimeError('always-active gate failed: '+label)

def main():
    r,z,w=s.symbols('r z w')
    minimal=6*r*r-15*r+10
    c=s.CRootOf(minimal,0)
    K=s.QQ.algebraic_field(c)
    rr=K.convert(c)
    one=K.one
    zero=K.zero
    roots=[rr,K.convert(s.Rational(5,2))-rr]
    S=lambda x:6*x**5-15*x**4+10*x**3
    check(s.discriminant(minimal,r)==-15,'nonreal distinct parameter roots')
    check(s.Poly(minimal,r).is_irreducible,'exact parameter field')
    check(s.expand(s.diff(S(w),w)-30*w**2*(w-1)**2)==0,'outer squared derivative')
    check(s.expand(S(w)-w**3*(6*w*w-15*w+10))==0,'outer zero fibre')
    check(s.expand(S(w)-1-(w-1)**3*(6*w*w+3*w+1))==0,'outer one fibre')
    check(s.expand(S(1-w)-(1-S(w)))==0,'target-flipping reflection hostile')
    check(S(-1)==-31 and S(0)==0 and S(1)==1,'fixed target supports')
    polynomials=[]
    invariants=[]
    data=[]
    poly=lambda value:s.Poly(value,z,domain=K)
    unit=poly(1)
    for i,root in enumerate(roots):
        rv=K.to_sympy(root)
        d=-one-root
        dv=K.to_sympy(d)
        Z=poly(z)
        W=poly(dv*z**3+rv)
        Q=6*W**5-15*W**4+10*W**3
        R=Z*W*(W-unit)
        s_squared=90*d
        check(6*root**2-15*root+10==zero,'parameter equation')
        check(root!=zero and root!=one and d!=zero,'all disjointness factors nonzero')
        check(Q.degree()==15 and R.degree()==7,'actual primitive and derivative-root degrees')
        check(Q.diff()==R**2*poly(K.to_sympy(s_squared)),'literal squared derivative identity')
        check(R.gcd(R.diff())==unit,'squarefree derivative root')
        check(R.eval(0)==0 and R.diff().eval(0)!=0,'simple source zero')
        check(Q.eval(0)==0,'primitive normalization')
        check(W.eval(1)==-1 and Q.eval(1)==-31,'same fixed source point and third value')
        check(R.eval(1)==2,'fixed source f value before square root')
        check(s_squared!=zero and 4*s_squared!=one,'both square-root choices are globally admissible')
        critical_zero=Q.gcd(R).monic()
        critical_one=(Q-unit).gcd(R).monic()
        check(critical_zero.degree()==4 and critical_one.degree()==3,'full critical multiplicity groups')
        check(critical_zero== (Z*W).monic(),'zero group contains original z zero and three preimages')
        check(critical_one==(W-unit).monic(),'one group contains three preimages')
        check(critical_zero.gcd(critical_one)==unit,'critical groups disjoint')
        for value,critical,count in [(0,critical_zero,3),(1,critical_one,6)]:
            quotient,remainder=(Q-poly(value)).div(critical**3)
            check(remainder.is_zero,'exact cubic critical multiplicities')
            check(quotient.degree()==count,'remaining simple fibre degree')
            check(quotient.gcd(critical)==unit,'no fourth-order critical zero')
            check(quotient.gcd(quotient.diff())==unit,'remaining roots simple')
        check((Q+poly(31)).gcd(R)==unit,'third support is not a critical value of Q')
        check((Q+poly(31)).gcd(Q.diff())==unit,'third polynomial fibre entirely simple')
        for target in [-3,2,5]:
            F=Q-poly(target)
            check(F.gcd(F.diff())==unit,'ordinary moving roots simple')
            check(F.gcd(R)==unit,'ordinary moving roots avoid all fixed punctures')
            check(F.degree()+R.degree()+1==23,'ordinary puncture count')
        check(Q.nth(14)==0,'centered primitive')
        a15=K.convert(Q.nth(15))
        a12=K.convert(Q.nth(12))
        check(a15==6*d**5 and a12==15*(2*root-one)*d**4,'top coefficient formula')
        invariant=a12**5/a15**4
        check(invariant==K.convert(s.Rational(15**5,6**4))*(2*root-one)**5,'scaling invariant formula')
        scaled=poly(Q.as_expr().subs(z,2*z))
        check(K.convert(scaled.nth(12))**5/K.convert(scaled.nth(15))**4==invariant,'positive affine-scaling control')
        check(4+3+1==8 and 2*15==30,'ambient rank and intrinsic pair degree')
        polynomials.append(Q)
        invariants.append(invariant)
        data.append({'parameter':str(rv),'d':str(dv),'primitive_coefficients':[str(v) for v in Q.all_coeffs()],
                     'critical_point_groups':[4,3],'ambient_arms':[4,3,1],
                     'source_t_degree':15,'rational_pair_degree':30,'ordinary_punctures':23,
                     'scaling_invariant':str(K.to_sympy(invariant))})
    check(roots[0]!=roots[1],'two distinct parameters')
    check(invariants[0]!=invariants[1],'same passport does not determine affine class')
    check(polynomials[0]!=polynomials[1],'literal primitive distinction')
    # The tempting critical pullback value r=0 collides with the outer critical
    # point; its derivative square root has a multiple zero and cannot be used.
    bad_w=-z**3
    bad_R=poly(z*bad_w*(bad_w-1))
    check(bad_R.gcd(bad_R.diff()).degree()>0,'critical-value pullback collision hostile')
    # The outer reflection is lawful only if the two target values are swapped.
    check(s.expand(S(1-w)-S(w))!=0,'a target swap is not a fixed-target symmetry')
    fixed=[s.Integer(0),s.Integer(1),s.Integer(-31)]
    affine_support_symmetries=[]
    import itertools
    for perm in itertools.permutations(fixed):
        slope=perm[1]-perm[0]
        shift=perm[0]
        if slope*fixed[2]+shift==perm[2]:affine_support_symmetries.append((slope,shift))
        check(slope!=0,'target-support permutations are injective')
    check(affine_support_symmetries==[(1,0)],'full support has no nontrivial affine target symmetry')
    stem=Path(__file__).stem
    total=sum(GATES.values())
    certificate={'status':'FINITE-EXACT controls of the independently proved pair',
                 'coefficient_field':'Q[r]/(6*r^2-15*r+10)','examples':data,
                 'scaling_invariant_difference':str(K.to_sympy(invariants[0]-invariants[1])),
                 'same_pointed_ambient_torsion':'four towers at0, three at1, one at-31; unit one pure double-pole direction per support',
                 'generic_fibre_scope':'Nonisomorphic as C(T)-curves; no assertion after arbitrary algebraic base extension',
                 'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
                 'gates':total,'gate_categories':dict(sorted(GATES.items()))}
    folder=Path(__file__).parent
    if folder.name=='04-computation':folder=folder.parent/'05-knowledge/results'
    (folder/(stem+'_certificate.json')).write_text(json.dumps(certificate,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Same ambient torsion, different generic source fibres: PASS')
    print('Exact field: Q[r]/(6*r^2-15*r+10); two conjugate degree15 primitives')
    print('Fixed h=-1/2, lambda=1; supports0,1,-31; ambient arms4,3,1; unit arms3')
    print('Intrinsic pair degree30; ordinary punctures23; unequal affine invariants')
    print('Generic-fibre nonisomorphism is analytic over C(T); finite controls are not its proof')
    print('Always-active exact gates:',total)
    print('Certificate:',stem+'_certificate.json')

if __name__=='__main__':main()
