#!/usr/bin/env python3
"""Independent Session59 audit. No producer imports. All gates survive python -O."""
from itertools import product, combinations
from math import gcd, comb
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
from pathlib import Path

checks=0
def check(condition, message):
    global checks
    checks+=1
    if not condition: raise RuntimeError(message)

# Explicit fields: prime fields plus F4 and F9, independently enumerated.
class Field:
    def __init__(self,p,extension=False):
        self.p=p; self.ext=extension; self.q=p*p if extension else p
    def add(self,a,b):
        if not self.ext:return (a+b)%self.p
        return (a%self.p+b%self.p)%self.p+self.p*((a//self.p+b//self.p)%self.p)
    def neg(self,a):
        if not self.ext:return -a%self.p
        return (-a%self.p)%self.p+self.p*((-(a//self.p))%self.p)
    def mul(self,a,b):
        if not self.ext:return a*b%self.p
        a0,a1=a%self.p,a//self.p;b0,b1=b%self.p,b//self.p
        # F4: alpha²=alpha+1; F9: alpha²=-1.
        const=1 if self.p==2 else -1
        linear=1 if self.p==2 else 0
        return ((a0*b0+const*a1*b1)%self.p+self.p*((a0*b1+a1*b0+linear*a1*b1)%self.p))
    def square(self,a):return self.mul(a,a)

def projective_points(q,n):
    for first in range(n):
        for tail in product(range(q),repeat=n-first-1):
            yield (0,)*first+(1,)+tail

count_rows=[]
for field in [Field(2),Field(3),Field(5),Field(7),Field(11),Field(13),Field(2,True),Field(3,True)]:
    q=field.q; mul=field.mul; add=field.add; sq=field.square
    for r in ([1,2,3] if q<=7 else [1,2]):
        total=boundary=singular=0
        for point in projective_points(q,2*r+1):
            c=point[0]; target=add(sq(c),sq(c))
            if all(add(sq(point[2*j+1]),sq(point[2*j+2]))==target for j in range(r)):
                total+=1;boundary+=c==0
                # For r2 in odd char, Jacobian rank drops iff an entire root pair is zero.
                if r==2 and field.p!=2 and c==0 and (point[1:3]==(0,0) or point[3:5]==(0,0)):singular+=1
        if field.p==2:
            predicted=sum(q**j for j in range(r+1))
            expected_boundary=sum(q**j for j in range(r))
        else:
            epsilon=1 if field.neg(1) in {sq(a) for a in range(q)} else -1
            expected_boundary=sum(2**s*comb(r,s)*(q-1)**(s-1) for s in range(1,r+1)) if epsilon==1 else 0
            predicted=(q-epsilon)**r+expected_boundary
            if r==2:check(singular==(4 if epsilon==1 else 0),f'singularity count q={q}')
        check(total==predicted,f'B_r projective count q={q} r={r}: {total}!={predicted}')
        check(boundary==expected_boundary,f'boundary count q={q} r={r}')
        count_rows.append(dict(q=q,characteristic=field.p,r=r,total=total,boundary=boundary,singular_r2=singular if r==2 and field.p!=2 else None))

# Independently generated eight line constraints, no assumed magic parameterization.
lines=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
C=S.Matrix([[int(j in line)-int(j in lines[0]) for j in range(9)]for line in lines[1:]])
D=smith_normal_form(C,domain=ZZ)
invariants=[D[i,i] for i in range(min(D.shape)) if D[i,i]]
check(invariants==[1]*6,'magic SNF')
M=S.Matrix([[1,1,0],[1,-1,-1],[1,0,1],[1,-1,1],[1,0,0],[1,1,-1],[1,0,-1],[1,1,1],[1,-1,0]])
check(C*M==S.zeros(7,3),'magic normal form forward')
minor=C.extract([0,1,2,3,5,6],[1,3,5,6,7,8])
check(minor.det()==1,'magic integer inverse minor')
# Free coordinates are A,C,E. Unique integral solution from six independent equations.
free=[0,2,4];rest=[1,3,5,6,7,8]
A,E,Cc=S.symbols('A E C')
unknown=-minor.inv()*C.extract([0,1,2,3,5,6],free)*S.Matrix([A,Cc,E])
full=[None]*9
for j,v in zip(free,[A,Cc,E]):full[j]=v
for j,v in zip(rest,unknown):full[j]=v
check(S.Matrix(full)==M*S.Matrix([E,A-E,Cc-E]),'magic scheme inverse over Z')
check(C*S.Matrix(full)==S.zeros(7,1),'seventh relation redundant integrally')

# Universal interpolation identities and denominator-cleared cycles.
u,v,w,t=S.symbols('u v w t')
op=lambda a,b:a+b+t*a*b
check(S.expand(op(op(u,v),w)-op(u,op(v,w)))==0,'associative polynomial')
check(S.cancel(op(u,-u/(1+t*u)))==0,'inverse identity')
check(S.expand(1+t*op(u,v)-(1+t*u)*(1+t*v))==0,'multiplicative coordinate')
ring_rows=[]
for n in range(2,13):
    for parameter in range(n):
        carrier=[u for u in range(n) if gcd(1+parameter*u,n)==1]
        def law(a,b):return (a+b+parameter*a*b)%n
        def inverse(a):return (-a*pow((1+parameter*a)%n,-1,n))%n
        for a,b in product(carrier,repeat=2):
            check(law(a,b) in carrier,'group closure over finite ring')
            check(law(a,inverse(a))==0,'group inverse over finite ring')
        # For two paths with edge labels a,b and c,d, group equality is divided product relation.
        for a,b,c,d in product(carrier,repeat=4):
            divided=(a+b-c-d+parameter*(a*b-c*d))%n
            check((law(a,b)==law(c,d))==(divided==0),'divided cycle over finite ring')
        ring_rows.append(dict(modulus=n,t=parameter,carrier_size=len(carrier)))
# Hostile: t=2 mod4, all factors are units; raw multiplicative relation loses a genuine condition.
n=4;t0=2;a,b,c,d=2,0,0,0
check(((1+t0*a)*(1+t0*b)-(1+t0*c)*(1+t0*d))%n==0,'raw products hostile')
check((a+b-c-d+t0*(a*b-c*d))%n==2,'divided equation retains lost condition')

# Two connected three-edge stars intersect in three isolated leaf vertices.
# Their edge arrays are arbitrary; cross-circuit map has differences of local gaps.
graph_rows=[]
for n in (2,3,4):
    buckets={};global_edge_labels=set()
    for left,right in product(product(range(n),repeat=3),repeat=2):
        gap=[(left[j]-right[j])%n for j in range(3)]
        obstruction=((gap[1]-gap[0])%n,(gap[2]-gap[0])%n)
        buckets[obstruction]=buckets.get(obstruction,0)+1
    for potentials in product(range(n),repeat=4):
        # Gauge first hub to 0, second hub potential b; leaves p0,p1,p2.
        b,*leaves=potentials
        global_edge_labels.add(tuple(leaves)+tuple((p-b)%n for p in leaves))
    check(len(buckets)==n*n,'k3 overlap obstruction onto')
    check(set(buckets.values())=={n**4},'k3 overlap equal fibers')
    check(len(global_edge_labels)==n**4,'k3 graph rank')
    for edges in global_edge_labels:
        gap=[(edges[j]-edges[j+3])%n for j in range(3)]
        check(len(set(gap))==1,'global edges kernel obstruction')
    graph_rows.append(dict(modulus=n,middle_size=n**6,obstruction_size=n**2,kernel_size=n**4))

# Crosspolytope semigroup normal form: c^(n-|m|_1) product signed generators.
semigroup_tests=0
for r in (1,2,3):
    for n in range(7):
        for m in product(range(-n,n+1),repeat=r):
            if sum(map(abs,m))<=n:
                factors=[(1,tuple(0 for _ in range(r)))]*(n-sum(map(abs,m)))
                for j,value in enumerate(m):
                    unit=tuple((1 if value>0 else -1) if k==j else 0 for k in range(r))
                    factors.extend([(1,unit)]*abs(value))
                check(len(factors)==n,'semigroup degree')
                check(tuple(sum(g[1][j] for g in factors) for j in range(r))==m,'semigroup generation')
                semigroup_tests+=1

# Literal blueprint maps to the initial pointed-monoid blueprint retain multiplicities.
blue_counts={}
for name,size in [('B',5),('Hplus',7),('Hminus',7),('M',9)]:
    maps=[]
    for z in product((0,1),repeat=size):
        c,a,b,d,e=z[:5]
        okay=(a+b==2*c and d+e==2*c)
        if name in ('Hplus','M'):
            h,hp=z[5:7]
            okay=okay and h+c==a+d and hp+a+d==3*c
        if name=='Hminus':
            h,hp=z[5:7]
            okay=okay and h+d==c+a and hp+a==c+d
        if name=='M':
            h,hp=z[7:9]
            okay=okay and h+d==c+a and hp+a==c+d
        if okay:maps.append(z)
    check(maps==[(0,)*size,(1,)*size],f'blue maps {name}')
    blue_counts[name]=len(maps)

# Geometry identities and genuine simple branch in the target torus.
a,b,d,e,c=S.symbols('a b d e c')
iu=(a+S.I*b)/(1+S.I);iv=(a-S.I*b)/(1-S.I)
check(S.simplify(iu*iv-(a*a+b*b)/2)==0,'split coordinate identity')
point={c:1,a:0,b:S.sqrt(2),d:1,e:1}
F=a*a+d*d-c*c;J=c*c+a*a-d*d
for f in (a*a+b*b-2*c*c,d*d+e*e-2*c*c,F,J):
    check(f.subs(point)==0,'branch base equations')
check((S.diff(F,d)-S.diff(F,e)).subs(point)==2,'F branch tangent')
check((S.diff(J,d)-S.diff(J,e)).subs(point)==-2,'J branch tangent')

result=dict(status='PASS',always_active_checks=checks,finite_field_counts=count_rows,blue_assignment_counts=blue_counts,magic_snf=[int(x) for x in invariants],magic_pivot_minor=int(minor.det()),finite_ring_fibers=ring_rows,graph_gluing=graph_rows,semigroup_elements_checked=semigroup_tests,hostile={'ring':'Z/4Z','t':2,'labels':[2,0,0,0],'raw_product_difference':0,'divided_difference':2})
Path(__file__).with_name('certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
