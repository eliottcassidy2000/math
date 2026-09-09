"""Independent exact exhaustion of the q=2 geometric reconstruction boundary."""
from pathlib import Path
import hashlib
import json
import sys
import sympy as s

sys.stdout.reconfigure(newline='\n')
z,w,v,T,C,C1,C2,A,k,a,beta=s.symbols('z w v T C C1 C2 A k a beta')
gates=0


def check(ok,label):
    global gates
    gates+=1
    if not bool(ok):raise RuntimeError(label)


def eq(expr,label):check(s.cancel(expr)==0,label)


def S(x):return 6*x**5-15*x**4+10*x**3


f=a*z*(z-beta)
Q=s.integrate(f*f,z)
eq(Q.subs(z,beta*w)-a*a*beta**5*S(w)/30,'whole q2 normalized primitive')
eq(Q.subs(z,0),'primitive normalization')
eq(s.diff(Q,z)-f*f,'literal square derivative')
eq(S(1-w)-(1-S(w)),'target-reflection identity')
check(s.expand(S(1-w)-S(w))!=0,'reflection is not fixed-target symmetry')

def B(cc):
    return s.expand((v-1)*(cc*v**3*(v*v-5*v+10)-T*(v-1)**5))


BB=B(C)
eq(BB-s.cancel((v-1)**6*(C*S(v/(v-1))-T)),'actual puncture polynomial')
expected=[C-T,-6*(C-T),15*(C-T),20*T-10*C,-15*T,6*T,-T]
check(s.Poly(BB,v).all_coeffs()==expected,'all seven coefficients')
eq(BB.subs(v,1),'old infinity puncture')
eq(s.diff(BB,v).subs(v,1)-6*C,'old infinity simple')
disc=s.factor(s.discriminant(BB,v))
check(disc!=0,'generic six points distinct')
check(s.Poly(BB,v).degree()==6,'six outer points')
eq(BB.subs(v,0)+T,'outer points avoid first inner point')
eq(s.Poly(BB,v).LC()-(C-T),'outer points avoid second inner point')

scaling=s.Poly(s.expand(B(C2).subs(v,A*v)-k*B(C1)),v)
eq(scaling.nth(0)-T*(k-1),'scaling constant implication')
eq(scaling.nth(1)-6*T*(A-k),'scaling linear implication')
eq(scaling.nth(6).subs({k:1,A:1})-(C2-C1),'scaling leading implication')
eq(scaling.as_expr().subs({k:1,A:1,C2:C1}),'positive identical normalized curve')

inverted=s.Poly(s.expand(s.cancel(v**6*B(C2).subs(v,A/v))-k*B(C1)),v)
eq(inverted.nth(6)-(-T-k*(C1-T)),'inversion leading implication')
eq(inverted.nth(5)-(6*T*A+6*k*(C1-T)),'inversion next implication')
eq(inverted.nth(5)+6*inverted.nth(6)-6*T*(A-1),'inversion forces unit scale')
residual=s.cancel(inverted.nth(0).subs({A:1,k:-T/(C1-T)})*(C1-T))
eq(residual-(C1*C2-(C1+C2)*T),'complete inversion obstruction')
eq(s.Poly(residual,T).nth(0)-C1*C2,'nonzero constants obstruct inversion')
eq(s.Poly(residual,T).nth(1)+(C1+C2),'target coefficient of inversion obstruction')

# Declared exact nonzero parameter controls; an all-parameter identity above
# proves exhaustion, while these guard hidden sign and target inversions.
for c1 in [s.Integer(1),s.Integer(-2),s.Rational(1,3)]:
    for c2 in [s.Integer(1),s.Integer(-2),s.Rational(1,3)]:
        check(s.Poly(residual.subs({C1:c1,C2:c2}),T).as_expr()!=0,'inversion hostile parameter pair')
        eq((B(C2)-B(C1)).subs({C1:c1,C2:c2})-(c2-c1)*(v-1)*v**3*(v*v-5*v+10),'literal normalized difference')

# Independent positive affine changes in the original source coordinate.
for b1 in [-2,1,3]:
    for b2 in [-2,1,3]:
        q1=C*S(z/s.Integer(b1))
        q2=C*S(z/s.Integer(b2))
        eq(q2.subs(z,s.Rational(b2,b1)*z)-q1,'lawful beta-coordinate affine equivalence')

source=Path(__file__).resolve()
folder=source.parent
if folder.name=='04-computation':folder=folder.parent/'05-knowledge'/'results'
cert={
 'status':'Complete q2 analytic exhaustion + independent FINITE-EXACT controls',
 'gates':gates,'producer_imports':False,
 'puncture_polynomial_coefficients':[str(cc) for cc in expected],
 'puncture_discriminant':str(disc),
 'scaling_conclusion':'k=A=1, C1=C2',
 'inversion_obstruction':str(residual),
 'combined_geometric_theorem':'Entire mutation supplier q>=2, same target, algebraic closure of C(T)',
 'scope_boundary':'Not arbitrary quintic first functions or total fibrations',
 'parameter_controls':{'nonzero_C_pairs':9,'beta_coordinate_pairs':9},
}
target=folder/(source.stem+'_certificate.json')
target.write_bytes((json.dumps(cert,indent=2,sort_keys=True)+'\n').encode())
print('PASS: complete quintic geometric-reconstruction boundary')
print('Scaling forces identity; inversion contradicts C1*C2!=0')
print('Combined with q>=3: geometric reconstruction iff for entire supplier q>=2')
print('Always-active gates: '+str(gates))
print('Certificate SHA256: '+hashlib.sha256(target.read_bytes()).hexdigest())
