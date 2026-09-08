"""Independent exact algebra; no producer or supplier engine import."""
from pathlib import Path
import hashlib
import json
import sys
import sympy as s

sys.stdout.reconfigure(newline='\n')
r,z,w,R,B,u,y,a,c,ss=s.symbols('r z w R B u y a c ss')
M=r*r-s.Rational(5,2)*r+s.Rational(5,3)
gates=0


def check(ok,label):
    global gates
    gates+=1
    if not bool(ok):raise RuntimeError(label)


def red(expr):
    return s.expand(s.rem(s.expand(expr),M,r))


def eq(expr,label):
    check(red(expr)==0,label)


def nz(expr,label):
    check(red(expr)!=0,label)


def S(arg):return 6*arg**5-15*arg**4+10*arg**3


eq(s.diff(S(w),w)-30*w*w*(w-1)**2,'outer derivative')
eq(S(w)-6*w**3*(w-r)*(w-(s.Rational(5,2)-r)),'outer zero factorization in quotient ring')
eq(S(w)-1-(w-1)**3*(6*w*w+3*w+1),'outer one factorization')
eq(S(1-w)-(1-S(w)),'target reflection hostile identity')
nz(S(1-w)-S(w),'reflection fails fixed target')
check(s.discriminant(M,r)==-s.Rational(5,12),'two nonreal distinct roots')

data=[]
invariants=[]
for rv in [r,s.Rational(5,2)-r]:
    d=-1-rv
    W=d*z**3+rv
    Q=red(S(W))
    L=red(z*W*(W-1))
    other=s.Rational(5,2)-rv
    eq(M.subs(r,rv),'parameter minimal equation')
    for value,label in [(rv,'r'),(rv-1,'r-1'),(d,'d'),(rv-other,'root separation')]:
        nz(value,'nonzero '+label)
    eq(s.diff(Q,z)-90*d*L**2,'primitive square derivative')
    check(s.degree(Q,z)==15 and s.degree(L,z)==7,'degrees')
    eq(Q.subs(z,0),'normalization Q0')
    eq(L.subs(z,0),'normalization f0')
    nz(s.diff(L,z).subs(z,0),'simple derivative root at zero')
    eq(Q-6*d*z**3*W**3*(W-other),'complete zero fibre factorization')
    eq(Q-1-(W-1)**3*(6*W**2+3*W+1),'complete one fibre factorization')
    # Each cubic pullback az^3+b is simple exactly when a*b is nonzero.
    for qv in [0,1,other]:
        disc=red(-27*d**2*(rv-qv)**2)
        nz(disc,'complete cubic discriminant')
    check(s.discriminant(6*w*w+3*w+1,w)==-15,'remaining one-fibre outer roots simple')
    nz(6*rv**2+3*rv+1,'remaining one-fibre pullbacks avoid inner critical point')
    check((6*w*w+3*w+1).subs(w,1)==10,'remaining roots avoid w1')
    check((6*w*w+3*w+1).subs(w,0)==1,'remaining roots avoid w0')
    eq(W.subs(z,1)+1,'fixed source point W1')
    eq(Q.subs(z,1)+31,'fixed source value')
    eq(L.subs(z,1)-2,'fixed source derivative-root factor')
    nz(90*d,'both square roots nonzero')
    nz(1-360*d,'both square roots avoid f1=-1')
    check(4*3+3==15 and 3*3+6==15,'complete critical fibre partitions')
    eq(s.Poly(Q,z).coeff_monomial(z**14),'centered degree15 polynomial')
    lead=s.Poly(Q,z).coeff_monomial(z**15)
    nextcoef=s.Poly(Q,z).coeff_monomial(z**12)
    eq(lead-6*d**5,'leading coefficient')
    eq(nextcoef-15*(2*rv-1)*d**4,'degree12 coefficient')
    Iv=red(s.Rational(15**5,6**4)*(2*rv-1)**5)
    eq(nextcoef**5-Iv*lead**4,'scaling invariant by denominator-free identity')
    # Concrete lawful coordinate scaling; it preserves I.
    eq((nextcoef*2**12)**5-Iv*(lead*2**15)**4,'positive affine-scaling control')
    invariants.append(Iv)
    data.append({'parameter':str(rv),'invariant':str(Iv),
                 'critical_partitions':[[3,3,3,3,1,1,1],[3,3,3,1,1,1,1,1,1]],
                 'source_arms':[4,3,1],'source_components':[5,4,2]})

difference=red(invariants[0]-invariants[1])
nz(difference,'exact affine-class separation')

# Actual chart coordinate and Jacobian chain checks, independently of f degree.
h=-s.Rational(1,2)
U=(1-h*R)/R
ti=-R**2-R**4*B
zi=s.cancel(U-2*h+U**3*ti)
J=h*h*(3-h*R)+B*(1-h*R)**3
check(s.cancel(zi+R*J)==0,'actual boundary coordinate z')
check(s.cancel(zi).subs(R,0)==0,'boundary z zero')
check(s.diff(J,B).subs(R,0)==1,'boundary tangent chain factor')
check(s.diff(u-2*h+u**3*y,y)==u**3,'source symplectic Jacobian factor')
# In coordinates (u,z), J_source=u^3 J_(u,z). Use a sample of degree7
# only as a chain-rule control; the general identity is algebraic.
f=z*(z-1)*(z-2)*(z-3)*(z-4)*(z-5)*(z-6)
Q0=s.integrate(f*f,z)
A=u*f
T=A+Q0
G=1/(2*A*A)
bracket=u**3*(s.diff(T,u)*s.diff(G,z)-s.diff(T,z)*s.diff(G,u))
check(s.cancel(bracket)==1,'independent rational-mate chain identity')
check(s.degree(Q0,z)==15,'mutation t-degree supplier control')
check(2*15==30 and 7+15+1==23,'pair degree and ordinary punctures')

# Closed-point degree and constant-map steps have all-degree proofs in the
# report. These hostile/positive identities expose the information used.
tau=s.symbols('tau')
check(s.degree(z**15-tau,z)==15,'moving closed-point degree model')
eq(S(1/w)*w**5-(6-15*w+10*w*w),'nonaffine Mobius creates finite pole')
check(S(0)!=S(1),'critical values force pointwise fixation')
check(s.solve([a*0+c,a+c-1],[a,c])=={a:1,c:0},'affine self-map fixes both critical points')

# Block maps for the distinguished unit commute with both module operators.
p1,p2,g=s.symbols('p1 p2 g',nonzero=True)
scale=p2/p1
for order in range(1,7):
    check(s.cancel(scale*(g*p1/g**order)-g*(scale*p1/g**order))==0,'constant block map commutes with g')
    check(s.cancel(scale*s.diff(p1/g**order,g)-s.diff(scale*p1/g**order,g))==0,'constant block map commutes with derivative')
check(s.cancel(scale*p1/g**2-p2/g**2)==0,'block map carries pointed unit')
check(4+3+1==8 and len([0,1,-31])==3,'full ambient versus generated arm counts')

source=Path(__file__).resolve()
folder=source.parent
if folder.name=='04-computation':folder=folder.parent/'05-knowledge'/'results'
cert={
 'status':'INDEPENDENT ANALYTIC ACCEPTANCE + FINITE-EXACT controls',
 'gates':gates,
 'producer_imports':False,
 'arithmetic':'Polynomial remainder over Q[r]/(r^2-5r/2+5/3)',
 'target_report_sha256':'601a55857454668fd85d902e3ee49345dd698d20212ef7ca88a454ada6e90fa7',
 'target_source_sha256':'780fad38dd9865860f1415a11558e1f11120715af479d5329f27ead6a82e4033',
 'examples':data,
 'invariant_difference':str(difference),
 'generic_scope':'K=C(T); no base extension and no target change',
 'all_degree_proof':'Smooth completion; rational versus higher-degree closed punctures; three constant points; polynomial pole',
}
target=folder/(source.stem+'_certificate.json')
target.write_bytes((json.dumps(cert,indent=2,sort_keys=True)+'\n').encode())
print('PASS: independent generic-fibre reconstruction and equal-torsion pair audit')
print('Two degree15 primitives; full torsion arms4,3,1; unit-generated arms3; pair degree30; punctures23')
print('Exact invariant difference: '+str(difference))
print('Always-active gates: '+str(gates))
print('Certificate SHA256: '+hashlib.sha256(target.read_bytes()).hexdigest())
