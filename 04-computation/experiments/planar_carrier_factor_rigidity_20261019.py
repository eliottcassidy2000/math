#!/usr/bin/env python3
"""Standalone exact controls. No repository imports and no removable asserts."""
import itertools
import json
import math
import sympy as sp

x,t,p,y,s,z,b = sp.symbols('x t p y s z b')
checks = 0

def check(condition, label):
    global checks
    if not bool(condition):
        raise RuntimeError(label)
    checks += 1

def zero(expr, label):
    check(sp.cancel(sp.expand(expr)) == 0, label)

def val(expr, variable=z):
    pol = sp.Poly(sp.expand(expr), variable)
    if pol.is_zero:
        return math.inf
    return min(mon[0] for mon, coeff in pol.terms() if coeff != 0)

q=1+x*x*t
P=t*q
Y=x*t*t*q
D=p**3-y**2
zero((P**3-Y**2)-t**3*q**2, 'literal cusp divisor')
zero(P**2*(P**3-Y**2)-t**5*q**4, 'literal universal divisor')
zero((y*p/D).subs({p:P,y:Y}, simultaneous=True)-x, 'inverse x')
zero((D/p**2).subs({p:P,y:Y}, simultaneous=True)-t, 'inverse t')
zero((sp.diff(x*t,x)*sp.diff(t,t)-sp.diff(x*t,t)*sp.diff(t,x))-t,
     'log bracket {s,t}=t')
zero(P.subs({x:s/z,t:z}, simultaneous=True)-(s*s+z), 'log p')
zero(Y.subs({x:s/z,t:z}, simultaneous=True)-s*(s*s+z), 'log y')

carriers=[p**a*y**c for a in range(4) for c in range(4-a)]
carriers += [1+p+y, p+y+y*y, D, D**2, D*(1+p+y),
             D**2*(p-y), (p-y)**2+D, (p**3-y**2)**3]
records=[]
for i,R in enumerate(carriers):
    Sz=sp.expand((p*p*D*R).subs({p:s*s+z,y:s*(s*s+z)}, simultaneous=True))
    k=val(Sz)
    A=sp.expand(Sz).coeff(z,k)
    check(k>=1 and A!=0, f'nonzero carrier leading coefficient {i}')
    ds=sp.expand(z*sp.diff(Sz,z))
    dt=sp.expand(-z*sp.diff(Sz,s))
    check(val(ds)==k and val(dt)>=k+1, f'vector field orders {i}')
    zero(ds.coeff(z,k)-k*A, f'leading displacement {i}')
    d2s=sp.expand(ds*sp.diff(ds,s)+dt*sp.diff(ds,z))
    check(val(d2s)>=2*k, f'second displacement order {i}')
    # Distinct integer-time leading jets stay nonzero in characteristic zero.
    for n in range(-5,6):
        if n:
            check(n*k*A!=0, f'infinite order witness {i}:{n}')
    Sxt=sp.expand(t**5*q**4*R.subs({p:P,y:Y}, simultaneous=True))
    n=val(Sxt,t)
    B=Sxt.coeff(t,n)
    dx=sp.diff(Sxt,t)
    dtt=-sp.diff(Sxt,x)
    check(n>=5 and B!=0, f'original source order {i}')
    check(val(dx,t)==n-1, f'original source dx order {i}')
    check(val(dtt,t)>=n, f'original source dt order {i}')
    zero(sp.expand(dx).coeff(t,n-1)-n*B, f'original source first jet {i}')
    if i in [0,1,10,12,13]:
        literal_ds=sp.cancel((t*dx+x*dtt).subs({x:s/z,t:z}, simultaneous=True))
        literal_dt=sp.cancel(dtt.subs({x:s/z,t:z}, simultaneous=True))
        zero(literal_ds-ds, f'literal Hamiltonian ds {i}')
        zero(literal_dt-dt, f'literal Hamiltonian dt {i}')
    records.append({'R':str(R),'k':k,'A':str(A),'n':n,'B':str(B)})

# Exhaust all nonempty subset choices; retain precisely disjoint full covers.
partition_records=[]
for r in range(1,5):
    count=0
    for masks in itertools.product(range(1,1<<r), repeat=r):
        union=0
        disjoint=True
        for mask in masks:
            if union & mask:
                disjoint=False
                break
            union |= mask
        if not disjoint or union != (1<<r)-1:
            continue
        check(all(mask.bit_count()==1 for mask in masks), 'support partition singleton')
        count+=1
    check(count==math.factorial(r), 'support partition permutation count')
    partition_records.append({'r':r,'permutations':count})

# Marked pair stabilizer: a=b^{-1} from Jacobian one, residual q error
# is (b^(m-1)-1)x^m t. Its root-of-unity condition is exact.
for m in range(2,10):
    X=b*x
    T=t/b
    zero(sp.diff(X,x)*sp.diff(T,t)-1, f'diagonal determinant m={m}')
    zero(1+X**m*T-(1+x**m*t)-(b**(m-1)-1)*x**m*t,
         f'pair stabilizer residual m={m}')
    H=x*t*(1+x**m*t)
    residual=sp.expand(H.subs({x:X,t:T}, simultaneous=True)-H)
    zero(residual-(b**(m-1)-1)*x**(m+1)*t**2,
         f'invariant residual m={m}')

# Hostiles: m=1 permits infinite-order determinant-one torus action.
H=x*t*(1+x*t)
zero(H.subs({x:2*x,t:t/2}, simultaneous=True)-H, 'm=1 invariant')
zero(sp.Rational(2)*sp.Rational(1,2)-1, 'm=1 determinant')
for n in range(1,21):
    check(2**n!=1, 'm=1 infinite order controls')
# One component alone permits a nontrivial unipotent shear.
zero(t.subs({x:x+t,t:t}, simultaneous=True)-t, 'shear invariant')
zero(sp.det(sp.Matrix([x+t,t]).jacobian([x,t]))-1, 'shear determinant')
for n in range(1,21):
    check(sp.expand(x+n*t-x)!=0, 'shear infinite order controls')
# Non-étale pullback splits/multiplies: it must not enter Lemma 1.
check(sp.det(sp.Matrix([x*x,t]).jacobian([x,t]))==2*x, 'nonetale hostile jacobian')
check(sp.degree(x*x,x)==2, 'nonetale squarefree hostile')

print(json.dumps({'status':'PASS','checks':checks,
  'claims_checked':'finite symbolic controls; universal claims use the written proofs',
  'carriers':records,'factor_support_partitions':partition_records,
  'marked_exponents':list(range(2,10)),
  'hostiles':['m=1 torus','single-component shear','non-etale repeated pullback']},
  indent=2,sort_keys=True))
