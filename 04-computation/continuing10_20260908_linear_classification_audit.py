"""Independent controls for the complete global source-linear classification.

The all-degree statements use the companion analytic audit; no producer is
imported. Exact controls reconstruct the global polynomial part by division.
"""
from pathlib import Path
from hashlib import sha256
from math import comb
import json
import sys
import sympy as s

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'05-knowledge/results' if HERE.name=='04-computation' else HERE
STEM='continuing10_20260908_linear_classification_audit'
gates=0

def need(ok,label):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(label)

def equal(expr,label):
    need(s.cancel(expr)==0,label)

def main():
    x,t,r,b=s.symbols('x t r b')
    a,h=s.symbols('a h',nonzero=True)
    B0=s.symbols('B0')
    records=[]
    for m in range(1,9):
        coefficients=s.symbols('a0:'+str(2*m+1))
        A=sum(coefficients[i]*x**i for i in range(2*m+1))
        quotient,remainder=s.div(A,x**m,x)
        B=B0+quotient-quotient.subs(x,0)
        F=A*t+B
        expected=B0-sum(coefficients[i]*r**(m-i) for i in range(m+1))-b*sum(coefficients[i]*r**(2*m-i) for i in range(2*m+1))
        equal(F.subs({x:1/r,t:-r**m-r**(2*m)*b},simultaneous=True)-expected,'complete source-linear global coefficient space')
        equal(expected.subs(r,0)-(B0-coefficients[m]-coefficients[2*m]*b),'complete boundary restriction')
        equal(s.diff(expected,r).subs(r,0)+coefficients[m-1]+coefficients[2*m-1]*b,'entire boundary first normal derivative')
        equal(s.diff(expected,b).subs(r,0)+coefficients[2*m],'entire boundary tangential derivative')
        # A degree above2m leaves an uncancellable b/r polar term.
        bad=x**(2*m+1)*t+x**(m+1)
        equal(bad.subs({x:1/r,t:-r**m-r**(2*m)*b},simultaneous=True)+b/r,'degree2m+1 globality hostile')
        # Constant A is globally submersive only on W1.
        Fconst=a*t+B0
        constinf=-a*r**m-a*r**(2*m)*b+B0
        need((s.diff(constinf,r).subs(r,0)!=0)==(m==1),'constant-A boundary classification')
        equal(-s.diff(Fconst,t)*s.diff(-x/a,x)-1,'positive polynomial mate for constant A')
        need(s.denom(-1/(a*r)).has(r),'constant-A mate still has a boundary pole')
        for n in range(2,2*m+1):
            A=a*(x-h)**n
            quotient,remainder=s.div(s.expand(A),x**m,x)
            B=B0+quotient-quotient.subs(x,0)
            F=A*t+B
            derivative=s.expand(s.diff(B,x).subs(x,h))
            if n<=m:
                equal(derivative,'low powers have a critical affine line')
            else:
                target=a*(-1)**(n-m-1)*comb(n-2,m-1)*h**(n-m-1)
                equal(derivative-target,'general pure-power first-jet identity')
            Finf=s.expand(F.subs({x:1/r,t:-r**m-r**(2*m)*b},simultaneous=True))
            need(s.denom(s.cancel(Finf))==1,'pure power gives an actual global function')
            equal(s.diff(F,x).subs(x,h)-derivative,'complete affine critical-line test')
            G0=1/(a*(n-1)*(x-h)**(n-1))
            DF=lambda P:s.diff(F,x)*s.diff(P,t)-s.diff(F,t)*s.diff(P,x)
            equal(DF(G0)-1,'whole pure-power rational primitive')
            for zero_h in [False,True]:
                d=derivative.subs(h,0) if zero_h else derivative
                Fr=s.diff(Finf,r).subs(r,0)
                Fb=s.diff(Finf,b).subs(r,0)
                if zero_h:Fr=Fr.subs(h,0);Fb=Fb.subs(h,0)
                affine_ok=(d!=0)
                boundary_ok=(Fb!=0) or (s.diff(Fr,b)==0 and Fr!=0)
                actual=bool(affine_ok and boundary_ok)
                wanted=(n==2 if m==1 else (not zero_h and (m+1<=n<=2*m-2 or n==2*m)))
                need(actual==wanted,'complete pure-power global submersion classification')
                if n==2*m-1:
                    need(s.degree(Fr,b)==1 and Fb==0,'middle exponent has exactly one boundary critical point')
                if actual:
                    f0=B.subs(x,h)
                    K=s.cancel((F-f0)/(x-h))
                    need(s.denom(K)==1,'complete reduced special-factor polynomial')
                    equal(K.subs(x,h)-derivative,'special components remain distinct')
                    equal(DF(K)-a*(x-h)**(n-1)*K,'factorized exact annihilator witness')
                records.append(dict(m=m,n=n,h_zero=zero_h,globally_submersive=actual,unit_order=n-1 if actual else None))
        orders=sorted({row['unit_order'] for row in records if row['m']==m and row['globally_submersive']})
        expected_orders=[1] if m==1 else list(range(m,2*m-2))+[2*m-1]
        need(orders==expected_orders,'exact positive order spectrum for this surface')
        need(2 not in orders,'order2 absent from complete globally submersive list')
    # A globally smooth two-root source need not have a rational mate.
    hostile=(x*x-1)**2*t+x*x
    equal(s.diff(hostile,x).subs(x,1)-2,'two-root positive affine derivative at1')
    equal(s.diff(hostile,x).subs(x,-1)+2,'two-root positive affine derivative at-1')
    hi=s.expand(hostile.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))
    equal(s.diff(hi,b).subs(r,0)+1,'two-root positive boundary derivative')
    equal(s.residue(1/(x*x-1)**2,x,1)+s.Rational(1,4),'nonzero residue obstruction at1')
    equal(s.residue(1/(x*x-1)**2,x,-1)-s.Rational(1,4),'nonzero residue obstruction at-1')
    equal(s.residue(1/x,x,0)-1,'simple-root logarithmic obstruction')
    # Actual m1 ring change: the affine primitive is excluded by the added D.
    g=a*t; P=-x*t
    equal(-a*s.diff(P,x)-g,'global m1 order-one annihilator witness on original chart')
    Pinf=1+r*b; ginf=-a*r*(1+r*b); qinf=-1/(a*r)
    equal(P.subs({x:1/r,t:-r-r*r*b},simultaneous=True)-Pinf,'annihilator witness is global across both charts')
    Vr=a*r*r; Vb=-a*(1+2*r*b)
    need(s.Poly(Vr,r,b).total_degree()==2,'global Hamiltonian r coefficient is polynomial')
    need(s.Poly(Vb,r,b).total_degree()==2,'global Hamiltonian b coefficient is polynomial')
    equal(Vr*s.diff(Pinf,r)+Vb*s.diff(Pinf,b)-ginf,'global m1 annihilator identity in second chart')
    equal(Vr*s.diff(qinf,r)+Vb*s.diff(qinf,b)-1,'global m1 rational primitive in second chart')
    equal((1+r*b)-b*r-1,'two global special components are comaximal')
    need(s.gcd(r,1+r*b)==1,'primitive boundary pole is regular on the other component')
    equal(s.limit(ginf*qinf,r,0)-1,'global unit has nonzero scalar principal coefficient')
    cert=dict(status='FINITE-EXACT controls; complete analytic classification in companion audit',
        symbolic_surface_m=list(range(1,9)),pure_power_records=records,
        module='original C[x,t]/D_F(C[x,t])',positive_order_spectrum='1 and every integer>=3; order2 absent',
        m1_constant_ring_change=dict(source_unit_zero=True,global_unit_order=1),gates=gates)
    raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
    (OUT/(STEM+'_certificate.json')).write_bytes(raw)
    print('INDEPENDENT SOURCE-LINEAR AUDIT: symbolic full coefficient spaces m1..8; every allowed exponent and h0 branch')
    print('PASS rational primitive pole-degree obstruction; complete global-submersion list; exact per-surface order spectrum')
    print('HOSTILES degree2m+1 boundary pole; middle exponent boundary critical point; smooth two-root source has nonzero residues')
    print('SCOPE source-linear global functions and original polynomial response module; all-degree proof is analytic')
    print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
    print('Always-active exact gates:',gates)

if __name__=='__main__':main()
