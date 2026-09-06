"""Independent literal-fibre/Berkowitz audit of the five endpoint certificates.
Integer matrices and exact rational arithmetic; no producer imports.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb,gcd,lcm,prod
import argparse,hashlib,json,sys,time
import sympy as S
sys.stdout.reconfigure(newline='\n')
COUNT=0

def check(ok,msg):
    global COUNT
    COUNT+=1
    if not ok:raise ArithmeticError(msg)

def literal_fibre(h,x,mass):
    g=x+3*h+1;N=6*h+3;out={}
    # Eliminate n_beta from the actual charge equation, then enforce all
    # nonnegativity constraints, independent of the producer's double loop.
    for nc in range(mass+1):
        numerator=N*mass-3*g*nc
        if numerator%(2*g):continue
        nb=numerator//(2*g);na=mass-nb-nc
        if min(na,nb)<0:continue
        out[(na,nb,nc)]=F(comb(mass,nc)*comb(mass-nc,nb))
    return out

def numeric_rows(h,x):
    g=x+3*h+1;raw1=literal_fibre(h,x,g);raw2=literal_fibre(h,x,2*g)
    expect1={(x+j,3*h-3*j,1+2*j) for j in range(h+1)}
    expect2={(2*x+e,6*h-3*e,2+2*e) for e in range(-1,2*h+1)}
    check(set(raw1)==expect1 and set(raw2)==expect2,'all literal fibres, with complete carry')
    scale1=comb(g,2*h+1)
    scale2=prod(range(2*g-4*h-1,2*g+1))
    p=[raw1[(x+j,3*h-3*j,1+2*j)]/scale1 for j in range(h+1)]
    q={e:raw2[(2*x+e,6*h-3*e,2+2*e)]/scale2 for e in range(-1,2*h+1)}
    check(p[h]==1 and p[0]>0 and q[-1]>0,'monic first row and retained nonzero inverse carry')
    check(scale2==factorial(2*g)//factorial(2*g-4*h-2),'double-scale factorial length')
    carry=F(2**(h+1)*factorial(3*h),factorial(6*h+3)*factorial(2*h+1))*x*prod(2*x+2*j+1 for j in range(h))
    check(q[-1]/p[0]==carry,'independent literal inverse-carry cancellation')
    return p,q,len(raw1)+len(raw2)

def times_t(v,p):
    h=len(v);lead=v[-1]
    return [(-lead*p[0])]+[v[i-1]-lead*p[i] for i in range(1,h)]

def matrix_from_literal(p,q):
    h=len(p)-1;R=[F(0)]*h
    # Horner in the companion operator; no polynomial long-division routine.
    for exponent in range(2*h,-1,-1):
        R=times_t(R,p);R[0]+=q[exponent]
    inverse=[-p[j+1]/p[0] for j in range(h)]
    invcheck=times_t(inverse,p)
    check(invcheck==[F(1)]+[F(0)]*(h-1),'literal t inverse in the quotient')
    R=[a+q[-1]*b for a,b in zip(R,inverse)]
    cols=[R]
    for _ in range(1,h):cols.append(times_t(cols[-1],p))
    denominator=lcm(*(v.denominator for col in cols for v in col))
    A=[[int(cols[j][i]*denominator) for j in range(h)] for i in range(h)]
    # Fraction-free Berkowitz on integers; producer uses polynomial
    # Faddeev-LeVerrier, so this is a distinct characteristic algorithm.
    coefficients=[int(c) for c in S.Matrix(A).berkowitz_charpoly().all_coeffs()]
    return [F(c,denominator**k) for k,c in enumerate(coefficients)]

def horner(coefficients,y):
    value=0
    for coefficient in reversed(coefficients):value=value*y+coefficient
    return value

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--certificate',type=Path,default=Path(__file__).with_name('continuing3_20260906_laurent_endpoints_certificate.json'))
    ap.add_argument('--quick',action='store_true',help='timing diagnostic only; never frozen as full proof audit')
    args=ap.parse_args()
    raw=args.certificate.read_bytes()
    check(hashlib.sha256(raw).hexdigest()=='b55bec0ba7f1063396af1a8db9be725f279fd01a403adb94b9b026b9490d8e62','full producer certificate pin')
    data=json.loads(raw);families=data['families']
    check([f['h'] for f in families]==list(range(6,11)),'complete named family scope')
    total_coefficients=0;matrices=0;fibres=0;primitive=0
    for family in families:
        h=family['h'];packets=family['characteristic']
        check(family['endpoint']==6*h+3 and family['g_min']==3*h+2,'family endpoint typing')
        check([p['k'] for p in packets]==list(range(1,h+1)),'all characteristic coefficients present')
        for packet in packets:
            k=packet['k'];co=packet['coefficients_ascending']
            check(packet['degree']==2*h*k and len(co)==2*h*k+1,'exact degree/list length')
            check(F(packet['content'])>0 and all(isinstance(v,int) and v>0 for v in co),'every certificate coefficient strictly positive')
            check(gcd(*co)==1,'primitive integer polynomial content')
            total_coefficients+=len(co)
        samples=(1,2*h*h+1) if args.quick else range(1,2*h*h+2)
        for x in samples:
            p,q,nf=numeric_rows(h,x);cs=matrix_from_literal(p,q)
            for k,packet in enumerate(packets,1):
                expected=F(packet['content'])*horner(packet['coefficients_ascending'],x-1)
                check(cs[k]==expected,'independent numeric characteristic equals certificate evaluation')
                check(cs[k]>0,'real-parameter shifted positivity control')
            matrices+=1;fibres+=nf;primitive+=gcd(x+3*h+1,6*h+3)==1
        # Earlier support return without gcd hypothesis, and all forbidden
        # smaller masses checked by the literal single-loop charge engine.
        gx=2;gh=gx+3*h+1
        check(gh==3*h+3 and gcd(gh,6*h+3)==3,'uniform nonprimitive hostile')
        for mass in range(1,h+1):check(not literal_fibre(h,gx,mass),'gcd-hostile no earlier mass')
        earlier=literal_fibre(h,gx,h+1)
        check((1,h-1,1) in earlier,'gcd-hostile earlier return attained')
        print('AUDIT_FAMILY','h',h,'endpoint',6*h+3,'distinct_parameter_values',len(samples),
              'max_required_degree',2*h*h,'literal_Berkowitz_matches',len(samples)*h,flush=True)
    check(total_coefficients==3170,'all shifted polynomial coefficients audited')
    if not args.quick:check(matrices==665,'complete independent identity sample universe')
    print('TOTAL','matrices',matrices,'complete_fibres',fibres,'primitive_parameter_controls',primitive,'positive_coefficients',total_coefficients)
    print('METHOD: literal multinomial fibres -> exact rational quotient -> integer Berkowitz characteristic -> bounded-degree identity test.')
    print('SCOPE: five fixed endpoints, all allowed integer g; all-h positivity and other actual support shapes remain OPEN.')
    print('PASS',COUNT,'always-active independent exact gates','QUICK_DIAGNOSTIC' if args.quick else 'FULL_IDENTITY_AUDIT')

if __name__=='__main__':main()
