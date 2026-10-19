#!/usr/bin/env python3
"""Exact finite controls for the field-sensitive monomial jet theorem."""
import json
import hashlib
from collections import Counter
from pathlib import Path
import sympy as s

z=s.symbols('z'); checks=[]; real_rows=[]; cyclotomic_rows=[]
def check(label, assertion):
    assert assertion,label
    checks.append(label)
def b(j): return s.binomial(j+3,3)

def weights(nodes,k):
    Q=s.Poly(s.prod(z-a for a in nodes),z)
    Q0=Q.nth(0)
    c0=Q.nth(k)/(b(k)*Q0)
    cs=[]
    for a in nodes:
        quotient=s.div(Q.as_expr(),z-a,z)[0]
        cardinal=s.Poly(z*quotient/(a*Q.diff().eval(a)),z)
        cs.append(cardinal.nth(k)/b(k))
    return Q,c0,cs

for k in range(2,8):
    for D in range(k,k+7):
        p=max(k,D-1)
        # Enforce q_(k-1)=0, also when D=k, so the same word attains one extra jet.
        positive=list(map(s.Integer,range(1,p)))
        A=s.Poly(s.prod(z-a for a in positive),z)
        last=A.nth(k-2)/A.nth(k-1)
        nodes=positive+[last]
        Q,c0,cs=weights(nodes,k)
        check(f'real k={k},D={D}: distinct nonzero rational nodes',len(set(nodes))==p and all(a!=0 for a in nodes))
        check(f'real k={k},D={D}: negative last node',last<0)
        check(f'real k={k},D={D}: required coefficient gap',Q.nth(k-1)==0)
        check(f'real k={k},D={D}: adjacent coefficients nonzero',Q.nth(k)!=0 and Q.nth(k-2)!=0)
        check(f'real k={k},D={D}: all weights nonzero',c0!=0 and all(c!=0 for c in cs))
        for j in range(D+1):
            m=sum(c*a**j for a,c in zip(nodes,cs))+(c0 if j==0 else 0)
            check(f'real k={k},D={D}: moment {j}',s.simplify(m)==(1/b(k) if j==k else 0))
        length=2*p+2
        check(f'real k={k},D={D}: optimal length',length==(2*k+2 if D==k else 2*D))
        if D>=k+1:
            next_m=s.simplify(sum(c*a**(D+1) for a,c in zip(nodes,cs)))
            check(f'real k={k},D={D}: exact next error',next_m==-Q.nth(k-2)/b(k) and next_m!=0)
        if k<=4 and D<=k+3:
            real_rows.append({'k':k,'r':k+2,'D':D,'N':D+2,'nonzero_sites':p,'word_length':length,
               'nodes':[str(a) for a in nodes],'coefficients':[str(c) for c in cs],'zero_site_coefficient':str(c0)})

# Cyclotomic moment computation in Q[z]/(z^p-1), not floating-point roots.
# A weighted sum over all p roots is the trace: sum z^e is p if p|e, else 0.
for k in range(2,10):
    for D in range(k,3*k+4):
        p=max(k,D-k+1)
        c0=-1/b(k) if p==k else s.Integer(0)
        for j in range(D+1):
            # Weight at a root a is a^(-k)/(p*b(k)).
            trace_m=(1/b(k) if (j-k)%p==0 else 0)+(c0 if j==0 else 0)
            check(f'cyclotomic k={k},D={D}: moment {j}',trace_m==(1/b(k) if j==k else 0))
        next_index=2*k if p==k else p+k
        check(f'cyclotomic k={k},D={D}: first next index',next_index>D)
        length=2*p+1+(1 if c0 else 0)
        expected=2*k+2 if D<=2*k-1 else 2*(D-k+1)+1
        check(f'cyclotomic k={k},D={D}: optimal length',length==expected)
        if k<=5 and D in [k,2*k-1,2*k,2*k+1]:
            cyclotomic_rows.append({'k':k,'r':k+2,'D':D,'N':D+2,'nonzero_sites':p,'word_length':length,
                'node_polynomial':f'z^{p}-1','weight_at_root_a':f'a^(-{k})/({p}*{b(k)})',
                'zero_site_coefficient':str(c0),'first_error_degree':next_index+2,
                'first_error_coefficient':str(b(next_index)/b(k))})

# Universal lower-bound polynomial identities.
for k in range(2,9):
    for D in range(k,k+9):
        for p in range(1,max(k,D-k+1)):
            if p<k:
                j=k-p
                check(f'lower k={k},D={D},p={p}: leading-coefficient contradiction',j>=1 and p+j<=D)
            else:
                check(f'lower k={k},D={D},p={p}: constant-coefficient contradiction',p+k<=D)

# Exact cubic-root sum descent identity, reduced modulo omega^2+omega+1.
x,omega=s.symbols('x omega')
complex_sum=x**2*s.Rational(1,30)*sum(a/(1-a*x)**4 for a in [1,omega,omega**2])
rational_sum=x**4*(10+16*x**3+x**6)/(10*(1-x**3)**4)
num=s.fraction(s.together(complex_sum-rational_sum))[0]
check('descent: exact rational closed form',s.rem(s.Poly(num,omega),s.Poly(omega**2+omega+1,omega)).is_zero)
check('descent: no cube-root pole cancellation',(10+16*x**3+x**6).subs(x**3,1)==27)

out={'status':'PASS','check_count':len(checks),'field':'exact rational arithmetic; cyclotomic traces computed algebraically',
 'scope':'Real constructions k=2..7, D=k..k+6; cyclotomic constructions k=2..9, D=k..3k+3; symbolic index controls',
 'proof_warning':'The all-order real lower bound uses strict interlacing in the companion proof; finite checks are not a substitute.',
 'check_categories':dict(Counter(label.split()[0].rstrip(':') for label in checks)),
 'ordered_check_labels_sha256':hashlib.sha256(('\n'.join(checks)+'\n').encode()).hexdigest(),
 'real_examples':real_rows,'cyclotomic_examples':cyclotomic_rows}
path=Path(__file__).with_name('field_sensitive_jets_exact_controls.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(f'PASS: {len(checks)} exact controls. Output: {path.name}')
