"""Independent exact bracket/chart audit; no primary engine import."""
from pathlib import Path
import hashlib,json,sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
x,t,R,B,a,b0,b1,b2,c=s.symbols('x t R B a b0 b1 b2 c')
gates=0
def check(ok,label):
    global gates
    gates+=1
    if not bool(ok):raise RuntimeError(label)
def eq(expr,label):check(s.cancel(expr)==0,label)
def J(F,H):return s.diff(F,x)*s.diff(H,t)-s.diff(F,t)*s.diff(H,x)
def inf(F):return s.expand(F.subs({x:1/R,t:-R**2-R**4*B},simultaneous=True))
def valuation(F):return min(int(z.as_powers_dict().get(R,0)) for z in s.expand(F).as_ordered_terms())

# All-parameter polynomial-coordinate identities, independently expanded.
BB=b2*x*x+b1*x+b0
U=2*a*t+BB+x+c;V=2*a*t+BB-x-c
H=U*V/(4*a)
eq(J(U,V)-4*a,'all-parameter coordinate Jacobian')
eq((U-V)/2-c-x,'inverse coordinate x')
eq(J(U,H)-U,'all-parameter positive eigen-coordinate')
eq(J(V,H)+V,'all-parameter negative eigen-coordinate')

records=[];critical_samples=[]
for av in [1,-2]:
    for BB in [x*x+2,-x,x,s.Integer(0),x+1,-x+1]:
        U=2*av*t+BB+x;V=2*av*t+BB-x
        H=s.expand(U*V/s.Integer(4*av))
        for m in [0,1]:
            F=s.expand(U*(H+2)**m)
            eq(J(F,H)-F,'constant-A source bracket')
            HI=inf(H);FI=inf(F)
            order=valuation(FI)
            expected_global=BB==-x or (BB==-x+1 and m==0)
            check((order>=0)==expected_global,'constant-A complete sampled globality')
            if expected_global:
                eq(s.diff(FI,R).subs(R,0),'constant-A whole-D normal critical')
                eq(s.diff(FI,B).subs(R,0),'constant-A whole-D tangent critical')
            records.append({'kind':'constant A','a':av,'B':str(BB),'P_degree':m,'order':order,'global':bool(expected_global)})
            if av==1 and BB in [-x,x*x+2] and m==1:critical_samples.append(F)

for rho in [-1,1]:
    av=-s.Rational(1,2*rho)
    for alpha in [0,1]:
        X=x-alpha
        for LL in [s.Integer(rho),s.Integer(-rho),s.Integer(0),x+1]:
            tau=t+LL
            H=s.expand(av*X*(tau*tau-rho*rho)+tau)
            hminus=-rho;hplus=rho
            for m in [0,1]:
                P=(H-hminus)*(H-hminus+1)**m
                F=s.expand((tau-rho)*(av*X*(tau-rho)+1)*(H-hminus+1)**m)
                eq(F-(tau-rho)*P/(tau+rho),'linear-A actual polynomial cancellation')
                eq(J(F,H)-F,'linear-A source bracket and orientation')
                check(hplus!=hminus and (hplus-hminus+1)!=0,'linear-A complete root conditions')
                FI=inf(F);order=valuation(FI)
                expected=2 if LL==rho else (-1 if LL==-rho else (-(m+1) if LL==0 else -3*(m+1)))
                check(order==expected,'linear-A exact boundary order')
                if LL==rho:
                    eq(s.diff(FI,R).subs(R,0),'linear-A whole-D normal critical')
                    eq(s.diff(FI,B).subs(R,0),'linear-A whole-D tangent critical')
                records.append({'kind':'linear A','rho':rho,'alpha':alpha,'L':str(LL),'P_degree':m+1,'order':order,'global':bool(LL==rho)})
                if rho==1 and alpha==0 and m==0 and LL in [s.Integer(1),x+1]:critical_samples.append(F)

check(len(records)==56,'complete declared pair count')
for F in critical_samples:
    check(list(s.groebner([s.diff(F,x),s.diff(F,t)],t,x))==[1],'independent full affine critical ideal')
check(len(critical_samples)==4,'declared critical ideal count')

# Raw eigenfunctions need not be source submersions.
H=t*t-x*x/16;U=2*t+x/2;F=U*U
eq(J(F,H)-F,'raw reciprocal-integer residue control')
eq(s.diff(F,x).subs(x,-4*t),'raw repeated component critical x partial')
eq(s.diff(F,t).subs(x,-4*t),'raw repeated component critical t partial')

source=Path(__file__).resolve();folder=source.parent
if folder.name=='04-computation':folder=folder.parent/'05-knowledge'/'results'
cert={'status':'INDEPENDENT ANALYTIC ACCEPTANCE + FINITE-EXACT controls',
      'gates':gates,'producer_imports':False,
      'target_report_sha256':'2b78aedb18e49e396c4477f71f6cd861e203f882dc58094e2fe85841962208f3',
      'target_source_sha256':'162360786a05ed8d592f635c61ee00094522055cd2e8ac2f90401ad63b02b875',
      'declared_pairs':len(records),'independent_critical_ideals':len(critical_samples),
      'global_survivors':sum(z['global'] for z in records),'all_global_survivors_boundary_critical':True,
      'records':records,'scope':'All-degree original quadratic witness with squarefree generic discriminant; repeated corridor excluded'}
target=folder/(source.stem+'_certificate.json')
target.write_bytes((json.dumps(cert,indent=2,sort_keys=True)+'\n').encode())
print('PASS: independent squarefree quadratic-witness global obstruction audit')
print('56 declared original-source pairs; four independent full critical ideals')
print('All '+str(cert['global_survivors'])+' global survivors have the whole added divisor critical')
print('Always-active gates: '+str(gates))
print('Certificate SHA256: '+hashlib.sha256(target.read_bytes()).hexdigest())
