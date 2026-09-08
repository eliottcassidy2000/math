#!/usr/bin/env python3
"""Exact 67-to-17 labelled compiler for the fixed-DG quartic proof.

This source checks finite coverage and algebraic interfaces. The individual
unrestricted-mate theorems are explicit analytic dependencies; their truth
is not inferred from filenames, a finite search, or this compiler's gates.
"""
from collections import Counter
import hashlib
import json
import sympy as S

gates = 0
def need(name, test):
    global gates
    if not bool(test):
        raise RuntimeError(name)
    gates += 1

def zero(name, expr):
    need(name, S.cancel(expr)==0)

# Each key is the FINITE partition. Infinity is retained as a separate
# marked multiplicity; no projective permutation identifies these keys.
# Tags are the complete positional conditions from boundary_exactness.
ROWS = [
    ((),8,'always',('infinity_octuple',)),
    ((1,),7,'always',('seven_one_infinity',)),
    ((3,),5,'always',('five_three_infinity',)),
    ((4,),4,'always',('four_four','four_four_transport')),
    ((3,1),4,'always',('four_three_one_infinity',)),
    ((5,),3,'always',('five_three_infinity',)),
    ((6,),2,'always',('degree_six_constant_d','const_d_translation','constant_d_shared_six')),
    ((5,1),2,'always',('five_two_one',)),
    ((4,1,1),2,'midpoint',('four_two_one_one',)),
    ((3,3),2,'always',('three_three_two',)),
    ((7,),1,'always',('seven_one_two_locations',)),
    ((8,),0,'always',('octuple_transport',)),
    ((7,1),0,'always',('seven_one_two_locations',)),
    ((6,1,1),0,'six_residue',('six_one_one',)),
    ((5,3),0,'always',('five_three',)),
    ((5,1,1,1),0,'sparse_elliptic',('five_one_one_one',)),
    ((4,3,1),0,'weighted_midpoint',('four_three_one',)),
]
routes = {part:(inf,tag,files) for part,inf,tag,files in ROWS}
need('no duplicated finite row',len(routes)==len(ROWS)==17)

# Independent enumeration by multiplicity COUNTS, not the predecessor's
# descending-first-part recursion. Every total degree from 0 to 8 occurs.
def count_vectors(m, remaining, prefix=()):
    if m==0:
        yield prefix
    else:
        for count in range(remaining//m+1):
            yield from count_vectors(m-1,remaining-count*m,prefix+(count,))

bank = []
for counts in count_vectors(8,8):
    part = tuple(m for m,cnt in zip(range(8,0,-1),counts) for _ in range(cnt))
    bank.append(part)
need('complete count-vector universe',len(bank)==67 and len(set(bank))==67)
counts_by_degree = Counter(map(sum,bank))
# A separate truncated Euler product supplies every partition count.
euler = [1]+[0]*8
for m in range(1,9):
    for n in range(m,9):
        euler[n] += euler[n-m]
need('Euler product degree counts',euler==[1,1,2,3,5,7,11,15,22])
need('independent enumerations agree',[counts_by_degree[n] for n in range(9)]==euler)

def leading_decision(part):
    n = sum(part)
    if not part: return 'always'
    if 2 in part: return 'NO: finite double residue'
    if all(m%2==0 for m in part):
        return 'always' if len(part)==1 else 'NO: rational square-field residue'
    if n==1: return 'always'
    if n==2: return 'NO: infinity simple residue'
    capacity = sum(m-2 for m in part if m>=3)
    required = n-2 if n%2 else n//2-1
    if capacity<required: return 'NO: primitive-degree capacity'
    if len(part)==1: return 'always'
    if len(part)==2 and all(m%2 for m in part): return 'always'
    positional = {(4,1,1):'midpoint',(6,1,1):'six_residue',
                  (4,3,1):'weighted_midpoint',(5,1,1,1):'sparse_elliptic'}
    need('complete residual position tag '+str(part),part in positional)
    return positional[part]

decisions = {}
for part in sorted(bank,key=lambda q:(sum(q),tuple(-x for x in q))):
    n = sum(part)
    need('degree in declared universe '+str(part),0<=n<=8)
    tag = leading_decision(part)
    decisions[part] = tag
    if tag.startswith('NO:'):
        need('rejected row has no exact route '+str(part),part not in routes)
    else:
        need('every exact row is routed '+str(part),part in routes)
        inf,expected,suppliers = routes[part]
        need('exact positional tag '+str(part),tag==expected)
        need('marked infinity degree '+str(part),n+inf==8)
        need('at least one explicit supplier '+str(part),bool(suppliers))
need('all and only seventeen exact rows',sum(not v.startswith('NO:') for v in decisions.values())==17)
need('exactly fifty uniform rejection rows',sum(v.startswith('NO:') for v in decisions.values())==50)
need('four position-sensitive rows',sum(tag!='always' for _,_,tag,_ in ROWS)==4)
need('all proposed routes accounted for',set(routes)=={p for p,v in decisions.items() if not v.startswith('NO:')})

# Homogenization proves the marked infinity rule at EVERY degree, including
# the nonzero constant and the absent-infinity degree-eight cases.
x,r,t,bD = S.symbols('x r t bD')
for n in range(9):
    aa = S.symbols('a0:'+str(n))
    NN = x**n+sum(aa[j]*x**j for j in range(n))
    transformed = S.Poly(S.expand(r**8*NN.subs(x,1/r)),r)
    need('literal infinity valuation '+str(n),min(m[0] for m in transformed.monoms())==8-n)
    zero('literal infinity unit '+str(n),transformed.coeff_monomial(r**(8-n))-1)

# Same unmarked binary partition does NOT imply the same actual entry.
need('both marked octuple placements retained',routes[()][0]==8 and routes[(8,)][0]==0)
need('all three marked 7+1 placements retained',
     [routes[p][0] for p in [(1,),(7,),(7,1)]]==[7,1,0])
need('all three marked 5+3 placements retained',
     [routes[p][0] for p in [(3,),(5,),(5,3)]]==[5,3,0])
need('4+4 finite/infinity survives',leading_decision((4,))=='always')
need('4+4 all-finite fails',leading_decision((4,4)).startswith('NO:'))
need('431 other two infinity placements fail',all(leading_decision(p).startswith('NO:') for p in [(4,1),(4,3)]))

# Literal labelled conditions, invariant under affine scaling/translation,
# never under an unweighted arbitrary projective change of the source line.
p,q,s,h = S.symbols('p q s h')
alpha = S.symbols('alpha',nonzero=True)
mid = 2*p-q-s
six = 3*(2*p-q-s)**2-4*(p-q)*(p-s)
weighted = (p-q)+3*(p-s)
for name, expr, degree in [('midpoint',mid,1),('sixfold',six,2),('weighted',weighted,1)]:
    moved = expr.subs({p:alpha*p+h,q:alpha*q+h,s:alpha*s+h},simultaneous=True)
    zero('labelled affine covariance '+name,moved-alpha**degree*expr)
zero('weighted labelled denominator',
     (3/(p-q)+1/(p-s))*(p-q)*(p-s)-weighted)

# Independent binomial jet coefficients recover the three local residue
# conditions; their nonzero unit multipliers do not change their vanishing.
d1,d2,z = S.symbols('d1 d2 z',nonzero=True)
B1 = sum(S.binomial(-S.Rational(1,2),j)*(z/d1)**j for j in range(3))
B2 = sum(S.binomial(-S.Rational(1,2),j)*(z/d2)**j for j in range(3))
zero('411 actual first residue numerator',S.expand(B1*B2).coeff(z,1)+(d1+d2)/(2*d1*d2))
zero('611 actual second residue numerator',
     S.expand(B1*B2).coeff(z,2)-(3*(d1+d2)**2-4*d1*d2)/(8*d1*d1*d2*d2))
B3 = sum(S.binomial(-S.Rational(3,2),j)*(z/d1)**j for j in range(2))
zero('431 actual first residue numerator',S.expand(B3*B2).coeff(z,1)+(3/d1+1/d2)/2)

# Positive and hostile positions, with the source root labels retained.
zero('411 positive midpoint',mid.subs({p:0,q:1,s:-1}))
need('411 hostile nonmidpoint',mid.subs({p:0,q:1,s:2})!=0)
k = S.symbols('k')
poly_k = 3*k*k+2*k+3
zero('611 positive algebraic locus',six.subs({p:0,q:1,s:k})-poly_k)
need('611 locus roots remain distinct and nonzero',S.gcd(poly_k,k*(k-1))==1)
need('611 quadratic is separable',S.discriminant(poly_k,k)!=0)
need('611 midpoint is a hostile',six.subs({p:0,q:1,s:-1})!=0)
zero('431 positive labelled positions',weighted.subs({p:0,q:3,s:-1}))
need('431 hostile wrong positions',weighted.subs({p:0,q:1,s:2})!=0)
a,d = S.symbols('a d',nonzero=True)
b,c = S.symbols('b c')
Nell = x**5*(a*x**3+b*x**2+c*x+d)
Rell = -2/(3*d*x**4)
def radical_derivative(N,R):
    return S.expand(N*S.diff(R,x)+S.diff(N,x)*R/2)
zero('5111 exact elliptic operator',radical_derivative(Nell,Rell)-1-(b*x*x+2*c*x)/(3*d))
zero('5111 sparse cubic discriminant',S.discriminant(a*x**3+d,x)+27*a*a*d*d)
zero('5111 positive primitive',radical_derivative(x**5*(x**3+1),-2/(3*x**4))-1)
need('5111 wrong-coefficient hostile remains smooth',S.discriminant(x**3+x+1,x)!=0)
need('5111 hostile primitive defect is nonzero',radical_derivative(x**5*(x**3+x+1),-2/(3*x**4))-1!=0)

# Two opposite finite residues must not be replaced by their zero sum.
zero('4+4 finite residue at zero',S.residue(1/(x*x*(x-1)**2),x,0)-2)
zero('4+4 finite residue at one',S.residue(1/(x*x*(x-1)**2),x,1)+2)
zero('4+4 finite-infinity primitive',S.diff(-1/x,x)-1/(x*x))

# Exact algebra behind the independently proved global root-gluing supplier.
NN,aa = S.symbols('NN aa',nonzero=True)
PP,QQ,MM,RR,bb,Y = S.symbols('PP QQ MM RR bb Y')
HH = NN*t*t+PP*t+QQ
FF = HH**2+MM*t+RR
def root_with_leading(f,var,V):
    coeff = S.Poly(S.expand(f),var)
    beta,gamma = coeff.coeff_monomial(var**3),coeff.coeff_monomial(var**2)
    return V*var*var+beta*var/(2*V)+(4*gamma*V*V-beta*beta)/(8*V**3)
zero('canonical approximate root',root_with_leading(FF,t,NN)-HH)
zero('affine fibre preserves canonical root',
     root_with_leading(FF.subs(t,aa*Y+bb),Y,NN*aa*aa)-HH.subs(t,aa*Y+bb))
xt = S.Matrix([[S.diff(1/r,r),S.diff(1/r,bD)],
               [S.diff(-r*r-r**4*bD,r),S.diff(-r*r-r**4*bD,bD)]]).det()
zero('actual source volume multiplier',xt-r*r)

def jac(f,g):
    return S.diff(f,x)*S.diff(g,t)-S.diff(f,t)*S.diff(g,x)
def inf(f):
    return S.expand(f.subs({x:1/r,t:-r*r-r**4*bD},simultaneous=True))

# F global alone is not a root-gluing premise. The inherited hostile has a
# polynomial canonical root on U0 and a genuine pole on the second chart.
Fbad = t**4+2*x**4*t*t
Hbad = t*t+x**4
need('root-gluing hostile F global',S.denom(S.cancel(inf(Fbad)))==1)
zero('root-gluing hostile canonical root',root_with_leading(Fbad,t,1)-Hbad)
zero('root-gluing hostile genuine boundary pole',inf(Hbad).coeff(r,-4)-1)
zero('root-gluing hostile critical first derivative',S.diff(Fbad,x).subs(t,0))
zero('root-gluing hostile critical second derivative',S.diff(Fbad,t).subs(t,0))

# The complete result is polynomial, not rational: actual global F with an
# exact rational mate and nonconstant L, from the proved shared-six table.
qq = 1+x*x*t
vv = x*qq
Hsharp = x*x*qq*qq+qq
Lsharp = vv
Gsharp = 1/(2*vv)-Hsharp
zero('global rational hostile exact Jacobian',jac(Hsharp**2+Lsharp,Gsharp)-1)
need('global rational hostile H extends',S.denom(S.cancel(inf(Hsharp)))==1)
need('global rational hostile L extends',S.denom(S.cancel(inf(Lsharp)))==1)
need('global rational hostile actual quartic',S.Poly(S.expand(Hsharp**2+Lsharp),t).degree()==4)
zero('global rational hostile affine pole',S.limit(x*Gsharp,x,0)-S.Rational(1,2))

# Constant output recombination preserves globality analytically; this is
# its exact bracket factor, with no assumption on alpha^2+beta^2.
aa,bb,cc,dd,fx,ft,gx,gt = S.symbols('aa bb cc dd fx ft gx gt')
zero('output pencil determinant',
     (aa*fx+bb*gx)*(cc*ft+dd*gt)-(aa*ft+bb*gt)*(cc*fx+dd*gx)
     -(aa*dd-bb*cc)*(fx*gt-ft*gx))

metadata = {
    'universe':[{'finite':list(p),'decision':decisions[p]} for p in sorted(decisions,key=lambda q:(sum(q),q))],
    'routes':[{'finite':list(p),'infinity':i,'condition':tag,
               'suppliers':['planar_jc48_sep08_'+s+'.md' for s in files]}
              for p,i,tag,files in ROWS],
    'conditions':[str(mid),str(six),str(weighted),'b=c=0 with a*d!=0'],
    'gates':gates,
}
digest = hashlib.sha256(json.dumps(metadata,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('Fixed-DG quartic exclusion: complete labelled compiler')
print('Exact gates:',gates)
print('Universe: 67 finite partitions; 50 uniform rejections; 17 exact rows')
print('Exact rows: 13 unrestricted-position, 4 conditional-position')
for p,i,tag,files in ROWS:
    finite = '+'.join(map(str,p)) if p else 'empty'
    print('finite='+finite+'; infinity='+str(i)+'; '+tag+' -> '+','.join(files))
print('Analytic supplier truth and audit status are separate from finite coverage.')
print('Global root-gluing interfaces and rational/pole hostiles: PASS')
print('Semantic SHA256:',digest)
