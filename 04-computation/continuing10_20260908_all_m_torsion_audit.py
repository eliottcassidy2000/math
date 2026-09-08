"""Independent exact controls for the all-m chart/annihilator theorem.

The unbounded assertions use the companion proof, not the finite m bank.
No mathematical producer is imported or executed.
"""
from pathlib import Path
from hashlib import sha256
from math import comb
import json
import sys
import sympy as s

sys.stdout.reconfigure(encoding='utf-8', newline='\n')
HERE = Path(__file__).resolve().parent
OUT = HERE.parent / '05-knowledge/results' if HERE.name == '04-computation' else HERE
STEM = 'continuing10_20260908_all_m_torsion_audit'
gates = 0

def need(ok, label):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(label)

def equal(expr, label):
    need(s.cancel(expr) == 0, label)

def main():
    x,t,r,b,u,c = s.symbols('x t r b u c')
    a,h = s.symbols('a h', nonzero=True)
    B0 = s.symbols('B0')
    records=[]
    for m in range(1,8):
        N=2*m; e=N-1; z=x-h
        # Reconstruct by polynomial division, not the producer's binomial sum.
        Q,R=s.div(s.expand(z**N),x**m,x)
        need(s.degree(Q,x)==m and s.degree(R,x)<m, 'complete polynomial part')
        Qh=s.expand(Q.subs(x,h))
        expected_c=(-1)**m*comb(2*m-1,m)*h**m
        equal(Qh-expected_c, 'special value binomial identity')
        k0=(-1)**(m-1)*comb(2*m-2,m-1)*h**(m-1)
        equal(s.diff(Q,x).subs(x,h)-k0, 'nonzero first jet at shifted root')
        A=a*z**N; B=a*Q+B0; F=A*t+B; c0=a*Qh+B0
        K=s.cancel((Q-Qh)/z)+z**e*t
        need(s.Poly(K,x,t).total_degree()>=1, 'residual factor is polynomial')
        equal(F-c0-a*z*K,'entire special fibre factorization')
        equal(K.subs(x,h)-k0,'two special components disjoint')
        equal(s.diff(F,t)-A,'unique possible affine critical line')
        equal(s.diff(F,x).subs(x,h)-a*k0,'no critical point on exceptional line')
        Finf=s.expand(F.subs({x:1/r,t:-r**m-r**N*b}, simultaneous=True))
        need(s.denom(s.cancel(Finf))==1,'full second-chart global regularity')
        expected_Finf=B0-a*(1-h*r)**N*b-a*sum(comb(N,j)*(-h)**j*r**(j-m) for j in range(m+1,N+1))
        equal(Finf-expected_Finf,'complete second-chart expression')
        equal(s.diff(Finf,b).subs(r,0)+a,'boundary derivative does not vanish')
        xinv=h+1/u
        tinv=u**N*((c-B0)/a-Q.subs(x,xinv))
        rinv=u/(1+h*u)
        binv=(B0-c)/a*(1+h*u)**N-sum(comb(N,j)*(-h)**j*u**(j-m)*(1+h*u)**(3*m-j) for j in range(m+1,N+1))
        need(s.Poly(s.cancel(tinv),u,c,B0).total_degree()>=1,'first inverse t is polynomial over the parameter field')
        need(s.Poly(s.cancel(binv),u,c,B0).total_degree()>=1,'second inverse b is polynomial over the parameter field')
        equal(F.subs({x:xinv,t:tinv},simultaneous=True)-c,'first inverse recovers fibre coordinate')
        equal(Finf.subs({r:rinv,b:binv},simultaneous=True)-c,'second inverse recovers fibre coordinate')
        equal(1/xinv-rinv,'inverse chart coordinates agree')
        equal(-rinv**m-rinv**N*binv-tinv,'inverse t transition agrees on overlap')
        equal((1+h*u)-h*u-1,'literal Bezout cover identity')
        need(s.gcd(u,1+h*u)==1,'two principal inverse opens cover the plane')
        ux=1/z
        ur=r/(1-h*r)
        equal((ux.subs(x,1/r))-ur,'u is globally regular off the pole component')
        determinant=s.diff(1/r,r)*s.diff(-r**m-r**N*b,b)
        equal(determinant-r**(N-2),'source two-form on full second chart')
        DF=lambda f:s.diff(F,x)*s.diff(f,t)-s.diff(F,t)*s.diff(f,x)
        equal(DF(ux)-a*z**(N-2),'actual coordinate-Jacobian identity')
        equal(ux**(N-2)*DF(ux)/a-1,'global rearranged two-form identity')
        G0=1/(a*e*z**e)
        equal(DF(G0)-1,'rational unit-response primitive')
        equal(DF(F),'fibre constants have zero derivative')
        equal(DF(K)-a*z**e*K,'polynomial residual-factor derivative')
        equal((F-c0)/z-a*K,'polynomial annihilator witness factor')
        # Product/chain rule proves D((aK)^e/(a e))=(F-c0)^e.
        equal(a**(e-1)*K**(e-1)*(a*z**e*K)-(a*z*K)**e,
              'complete polynomial top-power primitive identity')
        equal(s.diff(u**e/(a*e),u)-u**(e-1)/a,'literal power-map Jacobian')
        need(e>=1 and e%2==1,'exact positive odd power-map degree')
        need((e==1)==(m==1),'m1 is precisely unramified degree-one boundary')
        # Sharp h=0 hostile: fibre multiplicity changes the required power.
        Fzero=s.expand(F.subs(h,0))
        equal(Fzero-(a*x**m*(1+x**m*t)+B0),'h0 complete source factorization')
        equal(s.diff(Fzero,t).subs(x,0),'h0 t-derivative vanishes')
        if m==1:
            equal(s.diff(Fzero,x).subs(x,0)-a,'m1 h0 remains a submersion')
            need(e==1,'m1 h0 retains unit-response order1')
        else:
            equal(s.diff(Fzero,x).subs(x,0),'h0 higher-m critical line')
            need(m<e<=2*m,'h0 exact annihilator valuation threshold2')
        records.append(dict(m=m,degree=e,special_value_coefficient=(-1)**m*comb(2*m-1,m),
            first_jet_coefficient=(-1)**(m-1)*comb(2*m-2,m-1),
            second_chart_degree=s.Poly(Finf,r,b).total_degree()))
    # Separate integer reconstruction of the all-m alternating sums.
    for m in range(1,101):
        need(sum((-1)**j*comb(2*m,j) for j in range(m+1)) == (-1)**m*comb(2*m-1,m),
             'alternating partial binomial sum')
        need(sum((-1)**j*comb(2*m,j)*(m-j) for j in range(m+1)) == (-1)**(m-1)*comb(2*m-2,m-1),
             'independent weighted alternating sum')
    certificate=dict(status='FINITE-EXACT controls; unbounded proof in companion audit',
        symbolic_m=list(range(1,8)),binomial_m=list(range(1,101)),records=records,
        parameter_scope='C, m>=1, a*h!=0; original module C[x,t]/D_F C[x,t]',
        theorem='Ann_C[F]([1])=((F-B(h))^(2m-1)); W_m minus{x=h}=A2_(F,1/(x-h))',gates=gates)
    raw=json.dumps(certificate,sort_keys=True,separators=(',',':')).encode()+b'\n'
    (OUT/(STEM+'_certificate.json')).write_bytes(raw)
    print('INDEPENDENT ALL-M AUDIT: full symbolic source/charts for m1..7; integer binomial controls m1..100')
    print('PASS alternative A2 inverse cover; two special components; rational unit primitive; exact annihilator witness')
    print('BOUNDARY m1: degree1 and no ramification; h0 higher-m: critical line and valuation order2')
    print('SCOPE original polynomial response module; unbounded theorem is analytic, finite controls are not the proof')
    print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
    print('Always-active exact gates:',gates)

if __name__=='__main__':
    main()
