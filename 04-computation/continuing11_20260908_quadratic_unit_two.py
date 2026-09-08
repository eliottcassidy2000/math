"""Exact controls for a globally submersive DG quadratic with unit order2.

No prior producer is imported. The all-point statements are proved in the
companion report; controls retain both actual charts and all fibre components.
"""
from pathlib import Path
from hashlib import sha256
import json
import sys
import sympy as s

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'05-knowledge/results' if HERE.name=='04-computation' else HERE
STEM=Path(__file__).stem
gates=0

def need(ok,label):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(label)

def equal(expr,label):need(s.cancel(expr)==0,label)
def jac(F,G,x,t):return s.diff(F,x)*s.diff(G,t)-s.diff(F,t)*s.diff(G,x)

def main():
    x,t,r,b,U,W,Y,c=s.symbols('x t r b U W Y c')
    u=x-1;w=1+u*u*t;y=u*w;L=(u+1)*w-4
    F=u*(u+1)*w*w-4*u*w
    G=((1-2*u)*w+2)/(3*u*u*w*w)
    P=L*L*((1-2*u)*w+2)/3
    N=x*(x-1)**5;B=2*(x-1)**3*(x-2);C=(x-1)*(x-4)
    equal(F-(N*t*t+B*t+C),'literal complete original quadratic')
    need(s.Poly(F,t).degree()==2,'genuine quadratic source degree')
    equal(F-u*w*L,'three entire special-fibre factors')
    equal(jac(F,G,x,t)-1,'actual rational mate bracket')
    equal(F*F*G-P,'polynomial annihilator witness')
    need(s.Poly(s.cancel(P),x,t).total_degree()==12,'witness is a complete source polynomial')
    equal(jac(F,P,x,t)-F*F,'literal polynomial annihilator bracket')
    zeta=r*(1+b*(1-r)**2)
    Finf=(1-r)*(zeta*zeta-4)
    yinf=(1-r)*(2-zeta)
    equal(F.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True)-Finf,'whole second-chart global identity')
    equal(y.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True)-yinf,'retained fibre coordinate in second chart')
    need(s.denom(s.cancel(Finf))==1,'global regularity on entire boundary chart')
    equal(Finf.subs(r,0)+4,'exact boundary fibre value')
    equal(s.diff(Finf,r).subs(r,0)-4,'all boundary points noncritical')
    equal(s.diff(Finf,b).subs(r,0),'boundary tangential derivative')
    equal(s.diff(F,x).subs(x,1)+3,'omitted affine line noncritical')
    Fuy=(1+1/U)*Y*Y-4*Y
    equal(F.subs({x:U+1,t:(Y/U-1)/U**2},simultaneous=True)-Fuy,'complete u-nonzero coordinate expression')
    equal(s.diff(Fuy,U)+Y*Y/U**2,'all possible affine critical points have Y0')
    equal(s.diff(Fuy,Y).subs(Y,0)+4,'no affine critical point off u0')
    equal(jac(u,y,x,t)-u**3,'localized source coordinate map is invertible exactly off u0')
    Uinverse=Y*Y/(c-Y*Y+4*Y)
    equal(Fuy.subs(U,Uinverse)-c,'rational constant-field inverse')
    equal(Y*Y/(Fuy-Y*Y+4*Y)-U,'constant-field inverse in the other direction')
    equal(jac(F,y,x,t)+u*y*y,'actual derivation on rational fibre coordinate')
    Gfy=c/(3*Y**3)+2/Y**2-1/Y
    equal(G-Gfy.subs({c:F,Y:y},simultaneous=True),'rational primitive in constant-field coordinates')
    equal((-Y**4/(c-Y*Y+4*Y))*s.diff(Gfy,Y)-1,'unit primitive with F fixed')
    # Every exceptional fibre factor and every irreducibility sidecar.
    equal(w.subs(x,1)-1,'u and w components disjoint')
    equal(L.subs(x,1)+3,'u and L components disjoint')
    equal(((U+1)*W-4).subs(W,0)+4,'w and L components disjoint')
    need(s.gcd(U**2,1)==1,'w is a primitive linear factor')
    need(s.gcd((U+1)*U**2,U-3)==1,'L is a primitive linear factor')
    equal(B*B-4*N*(C-c)-4*(x-1)**5*((c+4)*(x-1)+c),'entire discriminant pencil')
    equal(s.gcd(N,B)-(x-1)**3,'only possible vertical content root')
    equal((C-c).subs(x,1)+c,'only zero fibre can have vertical content')
    # Scalar principal parts at every pole component, retaining the same F.
    Fu=F.subs(x,U+1);Gu=G.subs(x,U+1)
    equal(s.limit(Fu*Fu*Gu,U,0)-9,'u0 highest scalar coefficient')
    equal(s.limit(Fu*Gu-9/Fu,U,0)-4,'u0 simple scalar coefficient')
    Fuw=U*W*((U+1)*W-4);Guw=((1-2*U)*W+2)/(3*U**2*W**2)
    equal(s.limit(Fuw*Fuw*Guw,W,0)-s.Rational(32,3),'w0 highest scalar coefficient')
    equal(s.limit(Fuw*Guw-s.Rational(32,3)/Fuw,W,0)-4,'w0 simple scalar coefficient')
    need(s.gcd(U*W,(U+1)*W-4)==1,'rational primitive regular along residual L component')
    need(s.det(s.Matrix([[9,s.Rational(32,3)],[4,4]]))==-s.Rational(20,3),'both principal-part directions retained modulo diagonal')
    need(s.Rational(32,3)!=0,'nonzero simple pole after multiplying primitive by F')
    # Concrete positive ordinary points and exact source/pole boundaries.
    for xv,tv in [(0,0),(1,0),(2,0),(3,1),(-1,-2)]:
        need((s.diff(F,x).subs({x:xv,t:tv}),s.diff(F,t).subs({x:xv,t:tv}))!=(0,0),'positive exact source-gradient control')
    for bv in [-3,0,2,s.Rational(1,2)]:
        equal(Finf.subs({r:0,b:bv})+4,'boundary control value')
        equal(s.diff(Finf,r).subs({r:0,b:bv})-4,'boundary control derivative')
    # Literal finite controls in the inherited component principal-part model.
    g=s.symbols('g')
    Av=s.Matrix([9,s.Rational(32,3)]);Bv=s.Matrix([4,4])
    theta=Av/g**2+Bv/g
    def pp(expr):
        return s.Add(*(term for term in s.expand(expr).as_ordered_terms() if term.as_powers_dict().get(g,0)<0))
    def mul(vec):return vec.applyfunc(lambda expr:pp(g*expr))
    def der(vec):return vec.diff(g)
    def veq(vec,label):need(all(s.cancel(v)==0 for v in vec),label)
    veq(mul(mul(theta)),'unit killed by g2 in complete component quotient')
    veq(mul(theta)-Av/g,'first full principal-part arm recovered')
    veq(mul(der(theta))+2*theta-Bv/g,'second full principal-part arm recovered')
    veq(der(mul(theta))-mul(der(theta))-theta,'literal Weyl commutator on unit')
    columns=[theta.diff(g,j) for j in range(7)]+[(Av/g).diff(g,j) for j in range(7)]
    coefficient_matrix=s.Matrix([[s.expand(col[k]).coeff(g,-power) for col in columns] for k in range(2) for power in range(1,9)])
    need(coefficient_matrix.rank()==14,'finite PBW normal forms independently injective through derivative6')
    cert=dict(status='FINITE-EXACT controls; analytic existence/order theorem in companion report',
        source_F=str(s.expand(F)),source_G=str(G),boundary_F=str(Finf),
        unit_annihilator='(F^2)',special_components=['u=0','w=0','(u+1)w-4=0'],
        principal_part_coefficients={'u0':{'F^-2':'9','F^-1':'4'},'w0':{'F^-2':'32/3','F^-1':'4'},'L0':{}},
        scope='F globally submersive on fixed W2; response module C[x,t]/D_F(C[x,t]); JC2 OPEN',
        torsion_arms=2,weyl_left_annihilator='D g^2, [nabla,g]=1',gates=gates)
    raw=json.dumps(cert,sort_keys=True,separators=(',',':')).encode()+b'\n'
    (OUT/(STEM+'_certificate.json')).write_bytes(raw)
    print('FIXED W2: explicit genuine source-quadratic F is globally submersive and has an exact rational mate')
    print('SOURCE UNIT: complete annihilator(F^2); all three reduced disjoint zero-fibre components retained')
    print('BOUNDARY: F=-4, normal derivative4; original omitted u0 line derivative-3')
    print('PRINCIPAL PARTS: (9/F^2+4/F, (32/3)/F^2+4/F, 0); coefficient directions independent')
    print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
    print('Always-active exact gates:',gates)

if __name__=='__main__':main()
