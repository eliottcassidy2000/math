"""Independent literal multinomial/Berkowitz audit of regular reflected rows.
No producer imports. The all-height divisor and degree proof licenses the
complete 258-point residual-polynomial identity reconstruction.
"""
from fractions import Fraction as F
from math import comb, factorial, gcd, lcm, prod
from pathlib import Path
import argparse, hashlib, json, sys
import sympy as S
sys.stdout.reconfigure(newline='\n')
N=0

def need(ok,label):
    global N
    N+=1
    if not ok:raise ArithmeticError(label)

def fibre(h,z,m):
    g=3*h+z;out={}
    # Solve the literal charge equation for n_beta, while n_gamma ranges
    # over every possible value. No factorial progression is presupposed.
    for nc in range(m+1):
        numerator=3*h*m-3*g*nc
        if numerator%g:continue
        nb=numerator//g;na=m-nb-nc
        if min(na,nb)<0:continue
        out[(na,nb,nc)]=F(comb(m,nc)*comb(m-nc,nb))
    return out

def literal_rows(h,z):
    g=3*h+z;rows=[];count=0
    for mult in (1,2):
        m=mult*g;raw=fibre(h,z,m)
        triples=tuple((mult*z+2*j,mult*3*h-3*j,j) for j in range(mult*h+1))
        need(set(raw)==set(triples),'complete literal count fibre')
        scale=comb(m,mult*h)
        row=[raw[t]/scale for t in triples]
        need(row[-1]==1 and all(v>0 for v in row),'monic positive complete regular row')
        for j,value in enumerate(row):
            direct=F(factorial(mult*h)*prod(range(mult*z+2*j+1,mult*z+2*mult*h+1)),factorial(j)*factorial(mult*3*h-3*j))
            need(value==direct,'literal multinomial equals stated rising-product coefficient')
        rows.append(row);count+=len(raw)
    return rows[0],rows[1],count

def times_t(a,p):
    n=len(a);lead=a[-1]
    return [-lead*p[0]]+[a[i-1]-lead*p[i] for i in range(1,n)]

def characteristic(p,q):
    n=len(p)-1;v=[F(0)]*n
    # Horner companion evaluation, distinct from polynomial division.
    for c in reversed(q):v=times_t(v,p);v[0]+=c
    columns=[v]
    for _ in range(1,n):columns.append(times_t(columns[-1],p))
    denominator=lcm(*(x.denominator for col in columns for x in col))
    a=S.Matrix([[int(columns[j][i]*denominator) for j in range(n)] for i in range(n)])
    cs=a.berkowitz_charpoly().all_coeffs()
    return tuple(F(int(c),denominator**k) for k,c in enumerate(cs))

def horner(co,z):
    value=F(0)
    for c in reversed(co):value=value*z+c
    return value

def divisor(h,k,z):
    return prod(((z+2*ell-1)*(z+2*ell))**max(0,k-h+ell) for ell in range(1,h+1))

def falling(a,n):return prod(a-i for i in range(n))

def boundary_controls():
    # Exact factor zero orders underlying the all-height local proof.
    for h in range(1,7):
        for m in range(1,2*h+1):
            ell=(m+1)//2
            for j in range(h+1):
                factors=tuple(-m+s for s in range(2*j+1,2*h+1))
                zeros=factors.count(0)
                need(zeros==(1 if j<ell else 0),'Phi local coefficient zero exactly simple below ell')
                if j==0:need(prod(v for v in factors if v)!=0,'Phi constant derivative nonzero')
            for j in range(2*h+1):
                factors=tuple(-2*m+s for s in range(2*j+1,4*h+1))
                need(factors.count(0)==(1 if j<m else 0),'Psi coefficient zero exactly simple below m')
            need(m>=ell,'every small regular response has at least one full delta order')
        for r in range(1,h+1):
            H=h-r;z=2*r+1
            pp,qq=([F(1)],[F(1)]) if H==0 else literal_rows(H,z)[:2]
            for j in range(h+1):
                old=F(factorial(2*h+1)*falling(h-r,h-j),factorial(3*h-3*j)*factorial(1+2*j))
                want=F(0) if j<r else pp[j-r]
                need(old==want,'complete first boundary reflection identity')
            for j in range(-1,2*h+1):
                old=F(falling(2*h-2*r,2*h-j),factorial(6*h-3*j)*factorial(2+2*j))
                want=F(0) if j<2*r else qq[j-2*r]/factorial(4*h+2)
                need(old==want,'complete double boundary reflection including vanished inverse carry')
            need((-3*H,z,6*H+3*z)==(-3*(h-r),2*r+1,6*h+3),'literal exponent reflection and reordered support')
            g=3*h-r+1
            need(g==3*H+z and comb(g,2*h+1)==comb(g,H),'reflected first normalization and mass exactly agree')
            need(F(factorial(2*g),factorial(2*H)*factorial(4*h+2))==comb(2*g,2*H),'reflected doubled normalization exactly agrees')
            need(gcd(g,6*h+3)==gcd(z,3*H) and (-r+r,3*h-3*r,1+2*r)==(0,3*H,z),'primitive gcd and source monomial reflection agree')
    # The canceled old inverse carry cannot be evaluated naively at its
    # singular zero-root block, even when the raw specialized row is regular.
    x=S.symbols('x');p0=x+1
    qm1=(2*x+2)*(2*x+1)*(2*x)/S.factorial(9)
    q0=(2*x+2)*(2*x+1)/(S.factorial(6)*S.factorial(2))
    q1=(2*x+2)/(S.factorial(3)*S.factorial(4));q2=S.Rational(1,720)
    old=S.cancel(q0-p0*q1+p0*p0*q2-qm1/p0)
    need(old.subs(x,-1)==-S.Rational(1,90720),'old generic quotient nonzero at singular zero block')
    need(qm1.subs(x,-1)==q0.subs(x,-1)==q1.subs(x,-1)==0,'naive specialized old Laurent row loses that singular response')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--certificate',type=Path,default=Path(__file__).with_name('continuing4_20260906_regular_duality_certificate.json'))
    args=ap.parse_args();raw=args.certificate.read_bytes()
    need(hashlib.sha256(raw).hexdigest()=='0d5a65f03fc4f4295f3db38bca8609375cfea8805f21499a28e9a5e0d9a1ccd4','complete producer certificate pin')
    bank=json.loads(raw);need([r['H'] for r in bank['rows']]==list(range(1,7)),'exact certified height universe')
    matrices=fibres=primitive=coefficients=matches=0
    for row in bank['rows']:
        h=row['H'];packets=row['residuals']
        need([p['k'] for p in packets]==list(range(1,h+1)),'all characteristic coefficients supplied')
        for packet in packets:
            k=packet['k'];degree=4*h*k-k*(k+1);co=tuple(map(F,packet['coefficients']))
            need(packet['degree']==degree and len(co)==degree+1,'full residual degree and coefficient inventory')
            for value in co:need(value>0,'strict sign of every rational residual coefficient')
            coefficients+=len(co)
            need(2*sum(max(0,k-h+ell) for ell in range(1,h+1))==k*(k+1),'paired divisor degree')
        number=3*h*h-h+1
        need(all(p['degree']<number for p in packets),'identity sample count exceeds every residual degree bound')
        for z in range(1,number+1):
            p,q,nf=literal_rows(h,z);cs=characteristic(p,q)
            for k,packet in enumerate(packets,1):
                co=tuple(map(F,packet['coefficients']));div=divisor(h,k,z)
                need(div>0 and cs[k]/div==horner(co,z),'independent literal Berkowitz residual equals complete polynomial certificate')
                matches+=1
            matrices+=1;fibres+=nf
            if gcd(z,3*h)==1:
                primitive+=1;g=3*h+z
                for m in range(1,g):need(not fibre(h,z,m),'literal primitive first-support exclusion')
            if z==1:
                t=S.symbols('t');phi=S.Poly(sum(S.Rational(a.numerator,a.denominator)*t**j for j,a in enumerate(p)),t)
                need(phi.count_roots(-S.oo,0)==h and S.gcd(phi,phi.diff()).degree()==0,'independent small real-root/simple control for typed supplier')
        need(not fibre(h,3*h,1) and fibre(h,3*h,2)=={(1,1,0):F(2)},'dropped gcd hypothesis gives first return2 instead of6H')
        print('HEIGHT',h,'identity_nodes',number,'residual_degrees',[p['degree'] for p in packets],flush=True)
    need(coefficients==833 and bank['positive_coefficient_count']==833,'all positive polynomial coefficients audited')
    need(matrices==258,'complete degree-bound interpolation universe')
    boundary_controls()
    # Both actual alternatives, with a rational first cancellation at H=z=1.
    raw=fibre(1,1,4)
    need(sum(raw.values())==8 and sum(c*(-1)**nc for (na,nb,nc),c in raw.items())==0,'actual massg positive and cancellation alternatives')
    for m in range(1,8):
        value=sum(c*(-1)**nc for (na,nb,nc),c in fibre(1,1,m).items())
        need(value==0,'actual canceled example no nonzero moment before2g')
    need(sum(c*(-1)**nc for (na,nb,nc),c in fibre(1,1,8).items())==-224,'actual2g detection with complete regular row')
    print('TOTAL literal_matrices',matrices,'literal_fibres',fibres,'primitive_parameters',primitive,'characteristic_matches',matches,'positive_coefficients',coefficients)
    print('ANALYTIC_AUDIT paired-divisor orders and weighted degree hold for every H; reflection is exact off the discarded zero-root block')
    print('SINGULAR_CONTROL old h1,x=-1 generic response=-1/90720 while naive raw specialized row vanishes at t0')
    print('SCOPE H1..6,z>=1,gcd(z,3H)=1: actual first momentg or2g; all-H residual positivity remains OPEN')
    print('PASS',N,'always-active independent exact gates')

if __name__=='__main__':main()
