"""Independent recurrence/Groebner/chart audit of source-linear logarithmic pairs."""
from pathlib import Path
from hashlib import sha256
from collections import Counter
import json,sys
import sympy as s

sys.stdout.reconfigure(newline='\n')
gates=Counter()
def check(ok,label):
    gates[label]+=1
    if not ok:raise RuntimeError('always-active gate failed: '+label)

def main():
    x,t,v,r,b=s.symbols('x t v r b')
    cases=[]
    for alpha in [-2,0]:
        A=x-alpha
        for B in [s.Integer(7),x,x*x+1]:
            gamma=B.subs(x,alpha)+1
            for P in [s.Integer(1),(v-gamma)*(v-gamma-1)]:
                H=A*t+B
                cases.append(('+',s.expand((x-alpha)*P.subs(v,H)),H,False))
    for beta in [0,2]:
        for C in [s.Integer(0),s.Integer(2),x+1]:
            for R in [s.Integer(1),v-7]:
                H=5+(x-beta)*(-t+C)
                F=(-t+C)*R.subs(v,H)
                global_expected=C==0 or (C==2 and R==1)
                cases.append(('-',s.expand(F),s.expand(H),global_expected))
    for alpha in [-1,0]:
        for beta in [1,2]:
            A=(x-alpha)*(x-beta)/s.Integer(alpha-beta)
            for B in [x,x+(x-alpha)*(x-beta)]:
                b0=B.subs(x,beta)
                for P in [v-b0,(v-b0)*(v-b0-2)]:
                    H=s.expand(A*t+B)
                    F=s.cancel((x-alpha)/(x-beta)*P.subs(v,H))
                    cases.append(('+-',s.expand(F),H,False))
    check(len(cases)==40,'complete declared forty-pair grid')
    rows=[]
    for tag,F,H,global_expected in cases:
        fp=s.Poly(F,t)
        A=s.diff(H,t);B=H.subs(t,0)
        residual=s.expand(s.diff(F,x)*s.diff(H,t)-s.diff(F,t)*s.diff(H,x)-F)
        check(residual==0,'literal original-source logarithmic equation')
        for j in range(fp.degree()+1):
            fj=fp.nth(j);fk=fp.nth(j+1)
            recurrence=s.expand(A*s.diff(fj,x)-j*s.diff(A,x)*fj-(j+1)*s.diff(B,x)*fk-fj)
            check(recurrence==0,'independent coefficient recurrence')
        gradient=s.groebner([s.diff(F,x),s.diff(F,t)],t,x)
        check(list(gradient)==[s.Integer(1)],'complete affine critical ideal is unit')
        chart=s.expand(F.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))
        num,den=s.fraction(s.cancel(chart))
        is_global=not den.has(r)
        check(is_global==global_expected,'complete actual-chart Laurent globality')
        if is_global:
            check(s.diff(chart,r).subs(r,0)==0 and s.diff(chart,b).subs(r,0)==0,
                  'entire added divisor is critical')
        rows.append({'family':tag,'F':str(F),'H':str(H),'global':bool(is_global)})
    # The high-degree residue hostile satisfies the equation but is critical.
    A=x*(x*x-1)/2
    Y=(x*x-1)/(x*x)
    H=A*t
    F=(x*x-1)**3*t*t/4
    check(s.cancel(s.diff(Y,x)/Y-1/A)==0,'higher-degree rational logarithmic exactness hostile')
    check(s.expand(s.diff(F,x)*s.diff(H,t)-s.diff(F,t)*s.diff(H,x)-F)==0,'higher-degree hostile bracket')
    check(s.diff(F,x).subs(t,0)==0 and s.diff(F,t).subs(t,0)==0,'higher-degree hostile fails submersion')
    # Correctly signed labelled primitive principal parts at the two fibres.
    F=t*(1+x*t);H=-x*t
    check(H.subs(t,0)==0,'regular-component primitive numerator')
    check(s.cancel(H.subs(t,-1/x))==1,'positive simple principal-part coefficient')
    check(s.expand(s.diff(F,x)*s.diff(H,t)-s.diff(F,t)*s.diff(H,x)-F)==0,'unit-one hostile identity')
    F0=-t+x*x;H0=x*F0+5
    check(s.expand(s.diff(F0,x)*s.diff(H0,t)-s.diff(F0,t)*s.diff(H0,x))==F0,'zero-unit logarithmic boundary')
    check(-s.diff(F0,t)==1,'zero-unit boundary has polynomial mate x')
    # Exact Newton recurrence verifies the logical final step independently.
    a0,a1,a2=s.symbols('a0 a1 a2')
    for degree in [2,3]:
        coeff=[a0,a1,a2][:degree]
        poly=v**degree+sum(coeff[j]*v**j for j in range(degree))
        companion=s.zeros(degree)
        for j in range(1,degree):companion[j,j-1]=1
        for j in range(degree):companion[j,degree-1]=-coeff[j]
        recovered=[s.Integer(1)]
        powers=[s.trace(companion**j) for j in range(1,degree+1)]
        for k in range(1,degree+1):
            recovered.append(s.expand(sum((-1)**(i-1)*recovered[k-i]*powers[i-1] for i in range(1,k+1))/k))
        for k in range(1,degree+1):
            check(s.expand((-1)**k*recovered[k]-s.Poly(poly,v).nth(degree-k))==0,'Newton power sums recover root polynomial')
    stem=Path(__file__).stem
    cert={'universe':'forty declared source-submersive normal-form pairs; complete critical ideals and actual-chart checks',
          'rows':rows,'global_survivors':sum(row['global'] for row in rows),
          'all_degree_claim':'independent analytic proof, not extrapolation from the finite bank',
          'gates':sum(gates.values()),'categories':dict(sorted(gates.items())),
          'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    folder=Path(__file__).parent
    if folder.name=='04-computation':folder=folder.parent/'05-knowledge/results'
    (folder/(stem+'_certificate.json')).write_text(json.dumps(cert,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Independent source-linear logarithmic witness audit: PASS')
    print('Forty original-source pairs; independent coefficient recurrences and unit critical ideals')
    print('Six globally regular survivors; every entire added divisor critical')
    print('Signed unit-one and unit-zero hostiles; degree-three residue hostile: PASS')
    print('Always-active exact gates:',sum(gates.values()))
    print('Certificate:',stem+'_certificate.json')

if __name__=='__main__':main()
