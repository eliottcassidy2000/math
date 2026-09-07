#!/usr/bin/env python3
"""Exact next-curve geometry and an actual marked infinity A6 witness.

Universe: one displayed rational sextic; one raw labelled representation.
The representation is not an affine-complement or Keller realization.
"""
import hashlib
import json
import sympy as S

gates=0
def check(x,label):
    global gates
    gates+=1
    if not bool(x): raise RuntimeError(label)
def eq(x,y,label): check(S.cancel(x-y)==0,label)
def leading(f,v):
    n,d=map(lambda a:S.Poly(a,v),S.fraction(S.cancel(f)))
    ni=min(i[0] for i,c in n.terms());di=min(i[0] for i,c in d.terms())
    return ni-di,S.cancel(n.nth(ni)/d.nth(di))

s,t,p,q,z,a,b=S.symbols('s t p q z a b')
U=t**4+t**3+t**2
V=16*t**6+24*t**5-19*t**3-19*t**2
W=S.expand(V+19*U-19*U**2)
eq(W,-14*t**5-41*t**6-38*t**7-19*t**8,'finite cusp cancellation')
check(S.gcd(S.diff(U,t),S.diff(V,t))==t,'only critical parameter')
check(S.gcd(U,V)==t*t,'only cusp preimage')
N=S.cancel((U.subs(t,s)-U)/(s-t));M=S.cancel((V.subs(t,s)-V)/(s-t))
qpair=p*(p*p+p+1)/(2*p+1)
eq(S.rem(N.subs(t,p-s),s*s-p*s+q,s),p*(p*p+p+1)-(2*p+1)*q,'pair equation')
check((p*(p*p+p+1)).subs(p,-S.Rational(1,2))==-S.Rational(3,8),'denominator unit')
h=4*p**3+10*p*p+15*p+14
eq(S.rem(M.subs(t,p-s),s*s-p*s+q,s).subs(q,qpair),-p*p*h,'complete pair algebra')
eq(p*p-4*qpair,-p*(2*p*p+3*p+4)/(2*p+1),'pair discriminant')
check(S.discriminant(h,p)==-20972,'three distinct pair sums')
check(S.resultant(h,p*(2*p*p+3*p+4)*(2*p+1),p)==220864,'pairs distinct and off diagonal')
T=S.diff(U,t).subs(t,s)*S.diff(V,t)-S.diff(V,t).subs(t,s)*S.diff(U,t)
check(S.groebner([N,M,T],s,t,domain=S.QQ)==S.groebner([s+t**3+t*t+t,t**4],s,t,domain=S.QQ),'every off-diagonal collision transverse')
R=S.rem(V-b,U-a,t)
eq(R,-3*t**3+(16*a+5)*t*t+8*a*t-24*a-b,'triple cubic remainder')
r2=256*a*a+232*a+49
r1=128*a*a-8*a-3*b
r0=-384*a*a-16*a*b-201*a-8*b
eq(9*S.rem(U-a,R,t),r2*t*t+r1*t+r0,'triple divisibility')
check(S.groebner([r2,r1,r0],a,b,domain=S.QQ)==S.groebner([1],a,b,domain=S.QQ),'no triple fibre')

finite=[(U,W),(U,W/U),(U,W/U**2),(U**3/W,W/U**2),(U**5/W**2-S.Rational(1,196),W/U**2)]
for (x,y),order in zip(finite,[(2,5),(2,3),(2,1),(1,1),(1,1)]):
    check((leading(x,t)[0],leading(y,t)[0])==order,'finite chart order')
X=S.cancel(U.subs(t,1/z)/V.subs(t,1/z));Z=S.cancel(1/V.subs(t,1/z))
rho=S.cancel(Z/X**3-256)
tau=S.cancel(rho/X+15360)
eta=S.cancel(tau/X-753664)
inf=[(X,Z),(X,Z/X),(X,Z/X**2),(X,rho),(X,tau),(X,eta),(X/eta,eta),(X/eta**2-S.Rational(1,4848615424),eta)]
orders=[(2,6),(2,4),(2,2),(2,2),(2,2),(2,1),(1,1),(1,1)]
for (x,y),order in zip(inf,orders):
    check((leading(x,z)[0],leading(y,z)[0])==order,'infinity chart order')
check(leading(Z-256*X**3+15360*X**4-753664*X**5,z)==(11,-S.Rational(17,1024)),'first odd infinity exponent')
check(leading(eta,z)==(1,S.Integer(-17408)),'final transverse parameter')
mult=[2,2,1,1,2,2,2,2,2,1,1]
check(sum(m*(m-1)//2 for m in mult)+3==10,'rational sextic genus')
check(36-sum(m*m for m in mult)==4,'strict transform square')
check(18-2-sum(mult)==4-2*3==-2,'strict Nori criterion fails')

weights={'L':1};edges=set()
centres=[('E1',('L',)),('E2',('L','E1')),('E3',('L','E2')),('E4',('E3',)),('E5',('E4',)),('E6',('E5',)),('E7',('E5','E6'))]
for name,parents in centres:
    for v in parents:weights[v]-=1
    if len(parents)==2:edges.remove(tuple(sorted(parents)))
    weights[name]=-1
    for v in parents:edges.add(tuple(sorted((name,v))))
names=['L']+['E'+str(i) for i in range(1,8)]
check([weights[n] for n in names]==[-2,-2,-2,-2,-2,-3,-2,-1],'actual boundary weights')
matrix=S.Matrix([[-weights[x] if x==y else -int(tuple(sorted((x,y))) in edges) for y in names] for x in names])
check(matrix.det()==-1,'boundary determinant')
check(matrix*S.Matrix([-6,-4,-8,-12,-10,-8,-7,-14])==S.Matrix([0]*7+[1]),'marked meridian abelianization')

one=tuple(range(6))
def mul(p,q):return tuple(p[q[i]] for i in range(6))
def power(p,n):
    r=one
    for _ in range(n):r=mul(r,p)
    return r
def closure(gens):
    seen={one};todo=[one]
    while todo:
        p=todo.pop()
        for q in gens:
            z=mul(p,q)
            if z not in seen:seen.add(z);todo.append(z)
    return seen
l=(0,1,3,2,5,4);a=(1,2,0,4,5,3);c=one
d=(1,2,4,0,3,5);g=(2,4,3,1,0,5);e=(0,4,3,2,1,5);f=one
bperm=power(a,2);mu=(1,2,0,3,4,5)
fibres=dict(zip(names,[l,a,bperm,c,d,g,e,f]))
incident={'L':[c],'E1':[bperm],'E2':[a,c],'E3':[l,bperm,d],'E4':[c,g],'E5':[d,f],'E6':[f],'E7':[g,e,mu]}
for v in names:
    product=one
    for hperm in incident[v]:product=mul(product,hperm)
    check(power(fibres[v],-weights[v])==product,'full plumbing vertex relation')
for x,y in sorted(edges):check(mul(fibres[x],fibres[y])==mul(fibres[y],fibres[x]),'plumbing crossing commutator')
check(mul(f,mu)==mul(mu,f),'arrow commutator')
for v in list(fibres.values())+[mu]:
    check(sum(v[i]>v[j] for i in range(6) for j in range(i+1,6))%2==0,'all generators even')
check(mu==(1,2,0,3,4,5),'single three-cycle meridian')
left=closure([l,a]);right=closure([g,e]);whole=closure([l,a,e])
check(len(left)==len(right)==60 and len(whole)==360,'two order60 pieces generate A6')
check(all(any(v[0]==i for v in left) for i in range(6)),'left piece transitive on six')
check([i for i in range(6) if all(v[i]==i for v in right)]==[5],'right piece fixes exactly one label')

report={'status':'PROVED curve geometry; FINITE-EXACT marked infinity A6 witness; full affine consumer OPEN','curve':[str(U),str(V)],'finite_cusp':[2,5],'ordinary_nodes':3,'infinity_cusp':[2,11],'pair_polynomial':str(h),'Nori_margin':-2,'weights':weights,'edges':sorted(edges),'fibres':fibres,'arrow':mu,'group_orders':[len(left),len(right),len(whole)],'gates':gates}
print(json.dumps(report,sort_keys=True,indent=2))
print('semantic_sha256='+hashlib.sha256(json.dumps(report,sort_keys=True,separators=(',',':')).encode()).hexdigest())
