"""Independent carried-boundary norm jets via exact rational Sylvester resultants.
No producer imports. Complete low-height norm polynomials are reconstructed
at the proved degree bound; regular norm positivity uses previously audited data.
"""
from fractions import Fraction as F
from math import factorial,prod,comb
from pathlib import Path
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
G=0

def need(ok,label):
    global G
    G+=1
    if not ok:raise ArithmeticError(label)

def falling(x,n):return prod((x-i for i in range(n)),start=F(1))
def evalpoly(p,x):
    v=F(0)
    for c in reversed(p):v=v*x+c
    return v

def add(p,q):
    r=[F(0)]*max(len(p),len(q))
    for i,c in enumerate(p):r[i]+=c
    for i,c in enumerate(q):r[i]+=c
    return r

def mul(p,q):
    r=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]+=a*b
    return r

def det(a):
    a=[list(map(F,r)) for r in a];out=F(1)
    for i in range(len(a)):
        k=next((k for k in range(i,len(a)) if a[k][i]),None)
        if k is None:return F(0)
        if k!=i:a[k],a[i]=a[i],a[k];out=-out
        v=a[i][i];out*=v
        for j in range(i+1,len(a)):
            c=a[j][i]/v
            for k in range(i+1,len(a)):a[j][k]-=c*a[i][k]
    return out

def resultant(p,q):
    m,n=len(p)-1,len(q)-1
    matrix=[]
    for j in range(n):matrix.append([F(0)]*j+list(reversed(p))+[F(0)]*(n-1-j))
    for j in range(m):matrix.append([F(0)]*j+list(reversed(q))+[F(0)]*(m-1-j))
    return det(matrix)

def rows(h,x):
    y=x+h
    p=[F(factorial(2*h+1),factorial(3*h-3*j)*factorial(1+2*j))*falling(y,h-j) for j in range(h+1)]
    q=[falling(2*y,2*h-e)/F(factorial(6*h-3*e)*factorial(2+2*e)) for e in range(-1,2*h+1)]
    need(p[-1]==1 and len(q)==2*h+2,'literal monic p and complete t*q including lower carry')
    return p,q

def oldnorm(h,x):
    p,q=rows(h,x)
    need(p[0]!=0,'positive interpolation samples avoid singular inverse')
    return resultant(p,q)/p[0]

def interpolate(values):
    # Forward differences on x=1,...,D+1; Newton falling basis.
    out=[F(0)];basis=[F(1)];layer=list(values)
    for k in range(len(values)):
        out=add(out,[layer[0]*v/F(factorial(k)) for v in basis])
        layer=[b-a for a,b in zip(layer,layer[1:])]
        basis=mul(basis,[F(-k-1),F(1)])
    return out

def regular(H,z):
    if H==0:return [F(1)],[F(1)]
    phi=[F(factorial(H),factorial(j)*factorial(3*H-3*j))*prod((F(z+2*j+k) for k in range(1,2*H-2*j+1)),start=F(1)) for j in range(H+1)]
    psi=[F(factorial(2*H),factorial(j)*factorial(6*H-3*j))*prod((F(2*z+2*j+k) for k in range(1,4*H-2*j+1)),start=F(1)) for j in range(2*H+1)]
    return phi,psi

def regnorm(H,z):
    if H==0:return F(1)
    phi,psi=regular(H,z)
    return (-1)**H*resultant(phi,psi)

def constants(h,r):
    H=h-r
    a0=F(factorial(2*h+1)*factorial(H),factorial(3*H)*factorial(2*r+1))
    A=F((-1)**(r-1)*factorial(2*h+1)*factorial(H)*factorial(r-1),factorial(3*h))
    B=F(2*factorial(2*H)*factorial(2*r),factorial(6*h+3))
    return a0,A,B

def main():
    # Independent determinant convention controls.
    need(resultant([F(-2),F(1)],[F(3),F(4)])==11,'Sylvester orientation')
    need(resultant([F(6),F(-5),F(1)],[F(1),F(1)])==12,'quadratic product convention')
    triples=[];samples=0
    for h in range(1,5):
        D=2*h*h
        values=[oldnorm(h,F(x)) for x in range(1,D+2)];samples+=len(values)
        norm=interpolate(values)
        need(len(norm)==D+1 and norm[-1]!=0,'exact complete norm degree')
        for x,v in enumerate(values,1):need(evalpoly(norm,F(x))==v,'complete interpolation identity')
        need(evalpoly(norm,F(1,2))==oldnorm(h,F(1,2)),'fresh noninteger full-resultant control')
        for r in range(1,h+1):
            H=h-r;z=2*r+1;a0,A,B=constants(h,r)
            nh=regnorm(H,z)
            expected=B**r*a0**(2*r+1)*nh/(A*factorial(4*h+2)**H)
            jet=[sum(norm[j]*comb(j,k)*F(-r)**(j-k) for j in range(k,len(norm))) for k in range(r)]
            for c in jet[:-1]:need(c==0,'all lower negative-boundary jets vanish')
            need(jet[-1]==expected and (-1)**(r-1)*expected>0,'exact coefficient and sign from fresh resultants')
            # Reflection is checked directly before any quotient specialization.
            pp,qq=rows(h,F(-r));phi,psi=regular(H,z)
            need(pp==[F(0)]*r+phi,'complete first reflection')
            need(qq==[F(0)]*(2*r+1)+[v/factorial(4*h+2) for v in psi],'complete carried reflection away from singular zero')
            need(phi[0]==a0,'regular constant normalization')
            triples.append((h,r,H))
        print('LOW_HEIGHT',h,'samples',D+1,'boundaries',h,flush=True)
    need(samples==64 and len(triples)==10,'complete bounded fresh matrix universe')
    # Direct falling-product differentiation, not factorial expression evaluation.
    for h in range(1,13):
        for r in range(1,h+1):
            H=h-r;a0,A,B=constants(h,r)
            first=F(factorial(2*h+1),factorial(3*h))*prod((F(H-j) for j in range(h) if j!=H),start=F(1))
            carry=F(2,factorial(6*h+3))*prod((F(2*H-j) for j in range(2*h+1) if j!=2*H),start=F(1))
            need(first==A and carry==B and B>0,'direct full falling derivative constants')
            phi,_=regular(H,2*r+1)
            need(phi[0]==a0 and a0>0,'general complementary constant')
    # Previously audited regular characteristic data: positivity of the norms alone.
    path=Path(__file__).with_name('continuing4_20260906_regular_duality_certificate.json')
    raw=path.read_bytes()
    need(hashlib.sha256(raw).hexdigest()=='0d5a65f03fc4f4295f3db38bca8609375cfea8805f21499a28e9a5e0d9a1ccd4','prior independently audited regular certificate pin')
    norms={}
    for row in json.loads(raw)['rows']:
        H=row['H'];entry=row['residuals'][-1]
        need(entry['k']==H,'norm, not trace or another exterior coefficient')
        coeff=list(map(F,entry['coefficients']))
        need(len(coeff)==4*H*H-H*(H+1)+1,'complete residual norm degree')
        for c in coeff:need(c>0,'all residual norm coefficients strictly positive')
        norms[H]=coeff
    need(set(norms)==set(range(1,7)),'exact seven diagonals including H0 separately')
    for h,r,H in triples:
        if not H:continue
        z=2*r+1;n=evalpoly(norms[H],F(z))
        for ell in range(1,H+1):n*=((z+2*ell-1)*(z+2*ell))**ell
        need(n==regnorm(H,z),'independent direct complementary norm versus inherited certificate')
    need(constants(1,1)[1]==1 and regnorm(0,3)==1,'H0 and r1 boundary normalization')
    print('AUDIT PASS general leading jet formula; exact order r-1 iff complementary regular norm is nonzero')
    print('COROLLARY H=h-r in0..6: seven exact boundary diagonals at every height; sign(-1)^(r-1)')
    print('INDEPENDENT 64 positive interpolation points plus4 noninteger Sylvester resultants; all10 low-height jets and reflections')
    print('SCOPE all-h coefficient positivity remains open; no global actual sign inferred from negative-parameter jets')
    print('PASS',G,'always-active exact gates; raw LF')

if __name__=='__main__':main()
