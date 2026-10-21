#!/usr/bin/env python3
"""Independent audit: projective slopes and Eisenstein Euclidean ideal generators.
Does not import producer code. All checks survive python -O.
"""
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
COUNT = 0

def need(ok, message):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(message)

def norm(z):
    x,y=z
    return x*x-x*y+y*y

def multiply(z,t):
    x,y=z; u,v=t
    return x*u-y*v, x*v+y*u-y*v

def conjugate(z):
    x,y=z
    return x-y,-y

def subtract(z,t):
    return z[0]-t[0],z[1]-t[1]

def gcd_eisenstein(a,b):
    """Exact norm-Euclidean algorithm, independent of binary-form reduction."""
    while b != (0,0):
        numerator=multiply(a,conjugate(b)); den=norm(b)
        x,y=numerator[0]//den,numerator[1]//den
        options=[]
        for s in range(x-1,x+3):
            for t in range(y-1,y+3):
                r=subtract(a,multiply((s,t),b))
                options.append((norm(r),r))
        _,r=min(options)
        need(norm(r)<norm(b),'Eisenstein Euclidean norm descent')
        a,b=b,r
    return a

def slope(z,f):
    """Unique representative 1+s*omega, s even: rotate mod 2 then scale."""
    need(norm(z)%2==1,'slope input is a dyadic unit')
    for unused in range(3):
        a,b=z
        if a%2 and b%2==0:
            return (b*pow(a,-1,f))%f
        z=(-b,a-b)  # multiplication by omega
    raise RuntimeError('unit failed F4 normalization')

def operation(s,t,f):
    return ((s+t-s*t)*pow(1-s*t,-1,f))%f

def power_slope(s,n,f):
    result=0
    for unused in range(n):
        result=operation(result,s,f)
    return result

def order(s,f):
    value=0
    for n in range(1,f+1):
        value=operation(value,s,f)
        if value==0:
            return n
    raise RuntimeError('group order exceeded')

def forms(f):
    """Enumerate by positive b and explicitly add inverse when not ambiguous."""
    result=set()
    d=3*f*f
    for a in range(1,f+1):
        for b in range(0,a+1,2):
            if (b*b+d)%(4*a):
                continue
            c=(b*b+d)//(4*a)
            if c<a or gcd(a,gcd(b,c))!=1:
                continue
            result.add((a,b,c))
            if 0<b<a<c:
                result.add((a,-b,c))
    return sorted(result)

def form_generator(q,f):
    a,b,c=q
    if a%2==0:
        a,b,c=c,-b,a
    need(a%2==1,'primitive dyadic form has an odd outer coefficient')
    # I(q) = <a, (-b+f*sqrt(-3))/2> with the MIT/Cox ideal convention.
    beta=((f-b)//2,f)
    alpha=gcd_eisenstein((a,0),beta)
    need(norm(alpha)==a,'extended proper form ideal has exactly norm a')
    # Both generators divided by alpha are integral in R.
    for z in ((a,0),beta):
        numerator=multiply(z,conjugate(alpha))
        need(all(v%a==0 for v in numerator),'ideal basis lies in the gcd principal ideal')
    return alpha

def generic_quotient(f):
    units={(a,b) for a in range(f) for b in range(f) if gcd(norm((a,b)),f)==1}
    def m(z,t):
        return tuple(v%f for v in multiply(z,t))
    denominator={m((s,0),u) for s in range(f) if gcd(s,f)==1
                 for u in ((1,0),(0,1),(-1,-1))}
    canon={}; reps=[]
    for x in sorted(units):
        if x in canon:
            continue
        reps.append(x)
        for h in denominator:
            canon[m(x,h)]=x
    squares={canon[m(x,x)] for x in reps}
    norm_kernel={x for x in reps if norm(x)%8==1}
    return {'conductor':f,'units':len(units),'denominator':len(denominator),
            'classes':len(reps),'squares':len(squares),'norm8_kernel':len(norm_kernel),
            'nonsquare_norm1_witnesses':sorted(norm_kernel-squares)}

def main():
    pin_file=HERE/'FINAL_INPUT_PINS.json'
    if not pin_file.exists():
        pin_file=HERE/'INPUT_PINS.json'
    pins=json.loads(pin_file.read_text())
    snapshot=pins.get('snapshot','pinned')
    for row in pins['files']:
        actual=sha256((HERE/row['audit_copy']).read_bytes()).hexdigest()
        need(actual==row['sha256'],'frozen input hash '+row['audit_copy'])
    producer=json.loads((HERE/snapshot/'verification.json').read_text())
    claimed={level['conductor']:level for level in producer['levels']}
    data=[]; level_labels={}
    for f in (2,4,8,16,32,64,128):
        slopes=list(range(0,f,2))
        counts=Counter()
        for a,b in product(range(f),repeat=2):
            if norm((a,b))%2==0:
                continue
            s=slope((a,b),f)
            counts[s]+=1
            if f>=8:
                need(norm((a,b))%8==norm((1,s))%8,'norm descends under both unit and scalar ambiguity')
        need(set(counts)==set(slopes),'all quotient slopes present')
        need(set(counts.values())=={3*f//2},'correct scalar/unit denominator fibres')
        for s in slopes:
            for t in slopes:
                need(operation(s,t,f)==slope(multiply((1,s),(1,t)),f),'slope multiplication law')
        coords={(0,0):0} if f==2 else {
            (t,e):operation(power_slope(4%f,t,f),power_slope(2%f,e,f),f)
            for t in range(f//4) for e in range(2)}
        need(set(coords.values())==set(slopes) and len(coords)==len(slopes),'w,j coordinate bijection')
        if f>=4:
            need(order(4%f,f)==f//4 and order(2%f,f)==2,'generator orders')
        fforms=forms(f)
        q_slopes={q:slope(form_generator(q,f),f) for q in fforms}
        need(len(q_slopes)==f//2 and set(q_slopes.values())==set(slopes),'independent form/gcd bijection')
        inv_coord={s:c for c,s in coords.items()}
        labels={inv_coord[s]:q for q,s in q_slopes.items()}
        level_labels[f]=labels
        printed={tuple(r['coordinate_w_j']):tuple(r['form']) for r in claimed[f]['rows']}
        need(labels==printed,'literal producer labels via independently computed ideal gcd')
        if f>=4:
            lower_inverse={slope(form_generator(q,f//2),f//2):q for q in forms(f//2)}
            for row in claimed[f]['rows']:
                q=tuple(row['form']); alpha=form_generator(q,f)
                actual=lower_inverse[slope(alpha,f//2)]
                need(actual==tuple(row['extension_form']),'exact contraction/extension form target')
                coord=tuple(row['coordinate_w_j'])
                expected=(coord[0]%max(1,f//8),coord[1]) if f>=8 else (0,0)
                need(actual==level_labels[f//2][expected],'coordinate extension target')
        if f>=8:
            squares={operation(s,s,f) for s in slopes}
            kernel={s for s in slopes if norm((1,s))%8==1}
            reduced_kernel={s for s in slopes if s%8==0}
            need(kernel==squares==reduced_kernel,'norm kernel = square subgroup = extension-to-8 kernel')
            need({norm((1,s))%8 for s in slopes}=={1,3,5,7},'all four genus signatures')
        data.append({'conductor':f,'classes':len(slopes),'forms':len(fforms),
                     'labels':[{'coordinate':c,'form':q} for c,q in sorted(labels.items())]})
    # Exact finite valuation checks. The report separately checks the all-order induction.
    z=(1,4)
    for r in range(10):
        need(z[0]%2==1 and (z[1] & -z[1])==2**(r+2),'integer coefficient valuation')
        z=multiply(z,z)
    # Finite fibre states include genuinely distinct parent finite subsets and infinity.
    infinity=None; subsets=(0,1,2,3)
    def theta(bits):
        left,right=bits
        return (right,(left+right)%2)
    def plus(x,y):
        return (x[0]^y[0],x[1]^y[1])
    E,A,B,C=(0,0),(1,0),(0,1),(1,1)
    expected_stabilizers={('finite','finite'):{E},('infinity','finite'):{E,B},
                          ('finite','infinity'):{E,C},('infinity','infinity'):{E,A,B,C}}
    for parent in product((*subsets,infinity),repeat=2):
        choices=[(None,) if p is None else (0,1) for p in parent]
        states=list(product(*choices)); origin=states[0]
        stabilizer=set()
        for move in product((0,1),repeat=2):
            image=tuple(None if v is None else v^m for v,m in zip(origin,move))
            if image==origin:
                stabilizer.add(theta(move))
        kind=tuple('infinity' if p is None else 'finite' for p in parent)
        need(stabilizer==expected_stabilizers[kind],'actual inactive shore stabilizer')
        need(len(states)*len(stabilizer)==4,'orbit stabilizer')
        for bits in states:
            child=tuple(None if p is None else p+4*b for p,b in zip(parent,bits))
            blue=child[0]==child[1]
            if None not in parent:
                need(blue==((parent[0]==parent[1]) and bits[0]==bits[1]),'parent coincidence AND bit equality')
                t,e=theta(bits); n=pow(5,t,8)*pow(3,e,8)%8
                chi4=1 if n%4==1 else -1
                chi8=1 if n in (1,7) else -1
                need((chi8,chi4*chi8,chi4)==((-1)**bits[0],(-1)**bits[1],(-1)**(bits[0]^bits[1])),
                     'all three port/blue characters')
    need(expected_stabilizers[('infinity','finite')]!=expected_stabilizers[('finite','infinity')],
         'two size-two stabilizers differ')
    def add16(x,y): return ((x[0]+y[0])%4,(x[1]+y[1])%2)
    def project(x): return (x[0]%2,x[1])
    group=list(product(range(4),range(2)))
    liftsA=[x for x in group if project(x)==A]
    need(len(liftsA)==2 and all(add16(x,x)==(2,0) for x in liftsA),'all A lifts have nontrivial square')
    two_torsion=[x for x in group if add16(x,x)==(0,0)]
    need({project(x) for x in two_torsion}=={E,B},'characteristic two-torsion image')
    autos=[]
    for imagew,imagej in product(group,repeat=2):
        if add16(imagej,imagej)!=(0,0): continue
        image={((t*imagew[0]+e*imagej[0])%4,(t*imagew[1]+e*imagej[1])%2) for t,e in group}
        if len(image)==8: autos.append((imagew,imagej))
    need(len(autos)==8,'independent automorphism census')
    need(not any(project(u)==A and project(v)==C for u,v in autos),'no commuting shore-reflection lift')
    for t,e in product(range(2),repeat=2):
        need(((-t)%2,(-e)%2)==(t,e),'conjugation is identity at level8')
    hostile=generic_quotient(40)
    need(norm((1,8))==57 and pow(57%5,2,5)==4,'additional mod5 norm character obstructs a square at f40')
    need(hostile['classes']==24 and hostile['norm8_kernel']==6 and hostile['squares']==3,
         'non-dyadic scope hostile at conductor40')
    return {'status':('PASS for the corrected dyadic scope' if snapshot!='pinned' else 'PASS within the dyadic scope; one wording correction required'),
            'independent_checks':COUNT,'method':'even projective slopes plus Eisenstein Euclidean ideal gcd; no producer imports',
            'levels':data,'conductor40_scope_hostile':hostile,'g16_automorphisms':len(autos),
            'written_all_order_induction':'audited separately in INDEPENDENT_AUDIT.md'}

if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
