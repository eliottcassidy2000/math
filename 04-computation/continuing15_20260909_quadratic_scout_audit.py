"""Independent normalized-differential and valuation audit; no producer import."""
from pathlib import Path
from hashlib import sha256
import json,sys,math
import sympy as s
sys.stdout.reconfigure(newline='\n')
gates=0
def check(ok,label):
    global gates
    gates+=1
    if not ok:raise RuntimeError(label)
def eq(value,label):check(s.cancel(value)==0,label)
x,t,y,u,v,k,d1,d0,r,b=s.symbols('x t y u v k d1 d0 r b')
def J(F,H):return s.diff(F,x)*s.diff(H,t)-s.diff(F,t)*s.diff(H,x)
def main():
    # Independent conic vector field, no original-source coefficient expansion.
    D=k*k*x*x+d1*x+d0
    U=y+k*x+d1/(2*k);V=y-k*x-d1/(2*k)
    deltaU=y*s.diff(U,x)+s.diff(D,x)*s.diff(U,y)/2
    eq(deltaU-k*U,'symbolic conic eigenvector for arbitrary coefficients')
    eq(s.expand(U*V).subs(y*y,D)-(d0-d1*d1/(4*k*k)),'symbolic conic norm')
    eq((U-V-d1/k)/(2*k)-x,'rational conic inverse x')
    eq((U+V)/2-y,'rational conic inverse y')
    for m in range(1,16):
        if m%2:
            form=s.diff(u*u,u)/u**m
            eq(form-2*u**(1-m),'odd finite normalization')
        else:
            eq(s.diff(u,u)/u**(m//2)-u**(-m//2),'even finite normalization')
    # Generic double residues vary unless the B-square strictly wins.
    a2,b1,c0=s.symbols('a2 b1 c0',nonzero=True)
    check(s.diff(b1*b1+4*a2*(v-c0),v)!=0,'order-two A gives nonconstant residue square')
    bx,cx,cxx,bxx=s.symbols('bx cx cxx bxx')
    L=s.Matrix([[bx,0],[-cxx+bxx*cx/bx,-bx]])
    z=s.symbols('z')
    eq(L.charpoly(z).as_expr()-(z*z-bx*bx),'source-fixed-point linearization spectrum')
    spectrum=[]
    for m in range(1,11):
        H=x**(m+2)*t*t+x*t;K0=1+x**(m+1)*t
        pos=x**(m+2)*t*t if m%2 else x**((m+2)//2)*t
        neg=t**m*K0**(m+2) if m%2 else t**(m//2)*K0**((m+2)//2)
        ell=m//math.gcd(m,2);power=m+2 if m%2 else (m+2)//2
        eq(J(pos,H)-ell*pos,'positive polynomial spectrum generator')
        eq(J(neg,H)+ell*neg,'negative polynomial spectrum generator')
        eq(pos*neg-H**power,'both generators product is first integral')
        Y=1+2*x**(m+1)*t;R=s.cancel((Y-1)/(Y+1))
        eq(J(R,H)-m*R,'normalized deck ratio logarithmic rate')
        eq(x**m*H*(1-R)**2-R,'cyclic cover relation')
        spectrum.append([m,ell])
    for m in range(1,26):
        admissible=[lam for lam in range(0,3*m+1) if (2*lam)%m==0]
        ell=m//math.gcd(m,2)
        check(admissible==list(range(0,3*m+1,ell)),'valuation lattice equals exact spectrum multiples')
    for shift in [-2,1,3]:
        H=t*t-x*x/4;U=x+2*t
        for degree in [1,2,3]:
            roots=[shift+7*j for j in range(degree)]
            P=s.prod(v-z for z in roots)
            F=s.expand(U*P.subs(v,H))
            eq(J(F,H)-F,'shifted minimum quadratic witness bracket')
            eq(J(F,H/F)-1,'original rational mate')
            check(P.subs(v,0)!=0 and s.discriminant(P,v)!=0,'source-submersion exact factor conditions')
            check(s.degree(F,t)==2*degree+1,'all polynomial repair gauge degrees separated')
            pull=s.expand(F.subs({x:1/r,t:-r*r-b*r**4},simultaneous=True))
            valuation=min(term.as_powers_dict().get(r,0) for term in s.Add.make_args(pull))
            check(valuation==-(2*degree+1),'actual W2 pole survives')
            if degree==1:check(list(s.groebner([s.diff(F,x),s.diff(F,t)],t,x))==[1],'full independent affine critical ideal')
    F=x*(x*t+1);H=x*t+3*F*F
    eq(J(F,H)-F,'admissible fixed-double hostile')
    check(list(s.groebner([s.diff(F,x),s.diff(F,t)],t,x))==[1],'admissible-double source critical ideal')
    cert=dict(status='Independent analytic audit accepted; finite exact controls.',
              producer_report_sha256='e969113c9048dcc9ea30195aadaa1993e445235a4f00ed57014f988c510e04cc',
              source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),gates=gates,
              spectrum_controls=spectrum,producer_imported=False)
    dest=Path(__file__).with_name(Path(__file__).stem+'_certificate.json')
    if Path(__file__).parent.name=='04-computation':dest=Path(__file__).parents[1]/'05-knowledge/results'/dest.name
    dest.write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('INDEPENDENT_QUADRATIC_SCOUT_AUDIT PASS')
    print('TOTAL_ALWAYS_ACTIVE_GATES',gates)
    print('CERTIFICATE_SHA256',sha256(dest.read_bytes()).hexdigest())
if __name__=='__main__':main()
