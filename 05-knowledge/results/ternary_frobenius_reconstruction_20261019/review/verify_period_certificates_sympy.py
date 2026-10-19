#!/usr/bin/env python3
"""Independent literal SymPy expansion of both displayed rational formulas."""
from sympy import symbols,Matrix,Poly,expand,cancel,diff
from collections import defaultdict
from pathlib import Path
import json
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--base',type=Path,help='directory containing the two rational coefficient certificates')
args=parser.parse_args()
HERE=Path(__file__).resolve().parent
candidates=[args.base] if args.base is not None else [HERE.parent,HERE.parent/'session6-periods']
BASE=next((q for q in candidates if (q/'rational_kernel_d12_k2_eps1.json').is_file()),None)
if BASE is None:raise FileNotFoundError('Pass --base with the certificate directory')
x,y,z,X,Y,Z,c=symbols('x y z X Y Z c');vs=(x,y,z,X,Y,Z)
A={x:-y,y:x-y,z:z,X:-Y,Y:X-Y,Z:Z}
A2={x:y-x,y:-x,z:z,X:Y-X,Y:-X,Z:Z}
cyc={x:y,y:z,z:x,X:Y,Y:Z,Z:X}
swap={x:y,y:x,z:z,X:Y,Y:X,Z:Z}
checks=[]
def check(ok,label):
    if not ok:raise RuntimeError(label)
    checks.append(label)
def alt(columns):
    return Matrix([[t**a*T**b for a,b in columns]for t,T in ((x,X),(y,Y),(z,Z))]).det()
G2=(alt([(1,0),(3,0),(-1,3)])-alt([(1,0),(1,1),(1,2)])
    -2*alt([(1,0),(2,1),(0,2)])-alt([(2,0),(0,1),(1,2)])
    -alt([(2,0),(1,1),(0,2)])+alt([(3,0),(0,1),(0,2)])
    -2*alt([(0,1),(1,1),(2,1)]))
Palt=Poly(expand(x*y*z*G2),*vs)
def H(u,s,t):
    return (175*s*u**6+(-70*s*s+315*t)*u**5+(-175*s**3+525*s*t)*u**4
            +(136*s**4-754*s*s*t+639*t*t)*u**3+(21*s**5-14*s**3*t-231*s*t*t)*u*u
            +(-34*s**6+263*s**4*t-519*s*s*t*t+108*t**3)*u
            +11*s**7-111*s**5*t+287*s**3*t*t-68*s*t**3)
Fsym=0
for u,v,w,U,V,W in ((x,y,z,X,Y,Z),(y,z,x,Y,Z,X),(z,x,y,Z,X,Y)):
    Fsym+=(v*v-w*w)**2*H(u,v+w,v*w)*(u*V-v*U)*(u*W-w*U)/u
Psym=Poly(expand(x*y*z*Fsym),*vs)
output={'status':'PASS','families':[]}
for label,P,eps,file in [('alternating',Palt,-1,'rational_kernel_d3_k3_eps-1.json'),('symmetric',Psym,1,'rational_kernel_d12_k2_eps1.json')]:
    expr=P.as_expr()
    residual=Poly(expand((x-y)*expr-x*expr.xreplace(A)+y*expr.xreplace(A2)),*vs)
    check(residual.is_zero,label+'_universal_norm')
    check(expand(expr.xreplace(cyc)-expr)==0,label+'_cyclic')
    check(expand(expr.xreplace(swap)-eps*expr)==0,label+'_symmetry')
    check(expand(x*diff(expr,X)+y*diff(expr,Y)+z*diff(expr,Z))==0,label+'_gauge')
    cert={tuple(map(int,k.split(','))):int(v)for k,v in json.loads((BASE/file).read_text())[0].items()}
    pd={tuple(k):int(v)for k,v in P.terms()}
    sign=1 if pd==cert else (-1 if pd=={k:-v for k,v in cert.items()} else 0)
    check(sign!=0,label+'_saved_certificate_up_to_unit')
    specs=[]
    for p in (3,5,7,11,13,17,19,23,29,31):
        ds=defaultdict(int)
        for ex,coeff in P.terms():
            dst=tuple(ex[i]-1+p*ex[i+3]for i in range(3))
            check(min(dst)>0,label+'_interior_p'+str(p))
            ds[dst]=(ds[dst]+int(coeff))%p
        ds={k:v for k,v in ds.items()if v}
        if label=='symmetric':
            check(bool(ds)==(p>=11),'symmetric_small_prime_boundary')
            if p==11:check(ds.get((1,6,27))==10,'symmetric_p11_coefficient')
            if p==13:check(ds.get((1,12,25))==11,'symmetric_p13_coefficient')
        elif p>=5:
            check(ds.get((1,3,3*p-1))==1,'alternating_uniform_coefficient')
        specs.append({'p':p,'terms':len(ds),'nonzero':bool(ds)})
    aa=[ex[i]-1 for ex in pd for i in range(3)]
    output['families'].append({'name':label,'P_terms':len(pd),'certificate_file_sign':sign,'rational_original_exponent_min':min(aa),'rational_original_exponent_max':max(aa),'specializations':specs})
check(Psym.coeff_monomial(x**2*y**13*Z**2)==11,'symmetric_uniform_coefficient')
# Verify the degree-eleven p=2 example agrees with the paper's xyz D^2.
D=(x+y+z)*(x+y)*(y+z)*(z+x)
Alt137=Matrix([[t,t**3,t**7]for t in (x,y,z)]).det()
check(all(int(v)%2==0 for v in Poly(expand(Alt137-x*y*z*D**2),x,y,z).coeffs()),'characteristic_two_paper_example')
output['check_count']=len(checks)
print(json.dumps(output,indent=2))
