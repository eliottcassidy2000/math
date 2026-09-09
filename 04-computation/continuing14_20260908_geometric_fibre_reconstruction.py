"""Exact valuation/count controls for geometric generic-fibre reconstruction."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys
import sympy as s

sys.stdout.reconfigure(newline='\n')
gates=0


def check(ok,label):
    global gates
    gates+=1
    if not bool(ok):raise RuntimeError(label)


def determinant_valuation(p,q):
    # Finite points are literal Laurent coefficient dictionaries. Infinity
    # has homogeneous representative (1,0); its determinant is a unit.
    if p is None or q is None:
        if p is None and q is None:raise RuntimeError('repeated infinity')
        return 0
    difference={k:p.get(k,0)-q.get(k,0) for k in p.keys()|q.keys()}
    nonzero=[k for k,v in difference.items() if v!=0]
    if not nonzero:raise RuntimeError('repeated point')
    return min(nonzero)


def cross_valuations(points,order):
    a,b,c,d=[points[i] for i in order]
    den=determinant_valuation(a,c)+determinant_valuation(b,d)
    lam=determinant_valuation(a,b)+determinant_valuation(c,d)-den
    # Plucker identity: 1-lambda has this numerator, up to a harmless sign.
    one_minus=determinant_valuation(a,d)+determinant_valuation(b,c)-den
    return lam,one_minus


records=[]
for q in [2,3,7]:
    n=2*q+1
    # Arbitrary distinct leading coefficients suffice for this exact local
    # cluster test. Outer subleading coefficients vary to avoid a pure-scale
    # tautology. The report separately proves the actual roots have this type.
    points=[{0:i+1} for i in range(q)]
    points += [{-1:11+2*j,0:j*j+3,1:1-j} for j in range(n)]
    points += [None]
    counts=[0]*len(points)
    degenerating=0
    quadruples=0
    for ids in combinations(range(len(points)),4):
        a,b,c,d=ids
        vals=[cross_valuations(points,order) for order in [(a,b,c,d),(a,c,b,d),(a,d,b,c)]]
        predicates=[lv!=0 or ov!=0 for lv,ov in vals]
        check(len(set(predicates))==1,'cross-ratio ordering invariance')
        degenerate=predicates[0]
        expected=sum(i<q for i in ids)==2
        check(degenerate==expected,'complete two-plus-two classification')
        if degenerate:
            degenerating+=1
            for i in ids:counts[i]+=1
        quadruples+=1
    inner=(q-1)*s.binomial(2*q+2,2)
    outer=(2*q+1)*s.binomial(q,2)
    check(inner>outer,'distinct cluster counts')
    for i,value in enumerate(counts):
        check(value==(inner if i<q else outer),'exact pointwise intrinsic count')
    check(quadruples==s.binomial(3*q+2,4),'complete declared quadruple universe')
    check(degenerating==s.binomial(q,2)*s.binomial(2*q+2,2),'complete degenerating count')
    check(s.Rational(inner,outer)==s.Rational(2*(q+1),q),'all-degree count ratio sample')
    w=s.symbols('w')
    disc=s.discriminant(w**n-1,w)
    check(disc!=0,'moving scaled-root discriminant unit model')
    records.append({'q':q,'moving_degree':n,'quadruples':quadruples,
                    'degenerating':degenerating,'inner_count':int(inner),
                    'outer_count':int(outer),'pointwise_counts':counts,
                    'scaled_reduction_discriminant':str(disc)})

# General combinatorial difference is symbolic, independently of samples.
k=s.symbols('k',integer=True,positive=True)
DI=(k-1)*(2*k+2)*(2*k+1)/2
DO=(2*k+1)*k*(k-1)/2
check(s.expand(DI-DO-(k-1)*(2*k+1)*(k+2)/2)==0,'symbolic strict count difference')

# A balanced two-scale set admits an actual nonconstant cluster swap.
t,z=s.symbols('t z',nonzero=True)
phi=1/(t*z)
for c in [s.Integer(1),s.Integer(2),s.Integer(3)]:
    check(s.cancel(phi.subs(z,c)-1/(t*c))==0,'balanced cluster swap outward')
    check(s.cancel(phi.subs(z,1/(t*c))-c)==0,'balanced cluster swap inward')
balanced=[{0:s.Integer(c)} for c in [1,2,3]]+[{-1:s.Rational(1,c)} for c in [1,2,3]]
balanced_counts=[0]*6
for ids in combinations(range(6),4):
    lv,ov=cross_valuations(balanced,ids)
    if lv!=0 or ov!=0:
        for i in ids:balanced_counts[i]+=1
check(balanced_counts==[6]*6,'balanced counts fail separation')

# Two constants alone do not force a Mobius map to be constant. This is
# only a hostile to that implication, not an example in the mutation family.
psi=z/(t*z+1-t)
check(s.cancel(psi.subs(z,0))==0,'q2 hostile fixes first constant')
check(s.cancel(psi.subs(z,1)-1)==0,'q2 hostile fixes second constant')
check(s.diff(psi,t)!=0,'q2 hostile is nonconstant over base')
check(s.cancel(s.limit(psi,z,s.oo)-1/t)==0,'q2 hostile moves infinity')

source=Path(__file__).resolve()
folder=source.parent
if folder.name=='04-computation':folder=folder.parent/'05-knowledge'/'results'
cert={
 'status':'ANALYTIC geometric reconstruction q>=3 plus independent FINITE-EXACT valuation controls',
 'gates':gates,'producer_imports':False,'records':records,
 'exact_pair_q7_counts':[720,315],
 'all_degree_scope':'Same target; algebraic closure of C(T); mutation q>=3',
 'q2_boundary':'Counts recover only two constants; geometric iff not claimed here',
 'hostiles':['equal cluster sizes can be swapped by 1/(t*z)',
             'two constant points fixed by z/(t*z+1-t) do not force constant coefficients'],
 'old_packet_status':'Strengthened, not retracted; old frozen bytes unchanged',
}
target=folder/(source.stem+'_certificate.json')
target.write_bytes((json.dumps(cert,indent=2,sort_keys=True)+'\n').encode())
print('PASS: geometric generic-fibre reconstruction valuation audit')
print('q=2,3,7: all 9255 quadruples, three cross-ratio pairings each')
print('q=7 recovered intrinsic counts: inner720, outer315')
print('The geometric iff is proved for all q>=3; q2 remains outside this proof')
print('Always-active gates: '+str(gates))
print('Certificate SHA256: '+hashlib.sha256(target.read_bytes()).hexdigest())
