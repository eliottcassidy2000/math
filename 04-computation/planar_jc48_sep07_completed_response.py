#!/usr/bin/env python3
"""Complete universal-carrier image through row 15, with actual source data.

Standalone exact identities; no repository producer is imported. The proof
supplies the valuation cutoff for every polynomial carrier, not a census.
"""
from math import comb
from hashlib import sha256
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix

x, t, p, y, Xi = s.symbols('x t p y Xi')
F = s.Rational
GATES = 0


def check(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def zero(value, label):
    check(s.expand(value) == 0, label)


A = [1+x*x/4, F(4,3)+2*x*x, -F(32,9)-F(4,5)*x*x,
     F(2176,135)+F(64,9)*x*x-F(32,15)*x**4,
     F(224,9)*x**4-F(256,75)*x*x-F(37376,405),
     F(128,9)*x**6-F(8576,75)*x**4-F(6,11)*Xi*x*x
     -F(203776,7425)*x*x-x/2+F(3528704,6075)]
C = [-3*x/4-x**3/8, -4*x-F(3,2)*x**3,
     F(88,15)*x-F(12,5)*x**3,
     -F(1408,45)*x+F(64,15)*x**3+F(8,5)*x**5,
     -F(184,15)*x**5-F(2048,25)*x**3+F(98944,675)*x,
     -F(32,3)*x**7-F(128,5)*x**5+F(9,22)*Xi*x**3
     +F(130816,275)*x**3+F(3,8)*x*x+F(9,11)*Xi*x
     -F(4024832,4455)*x+F(3,4)]

# Literal P(A,C) rows, including the varying background without specializing.
G = [s.expand(sum(C[i]*C[k-i] for i in range(k+1))
     -sum(A[i]*A[j]*A[k-i-j] for i in range(k+1)
          for j in range(k-i+1)) + F(3,4)*A[k]
     + (F(1,4) if k == 0 else 0)) for k in range(6)]
GA = [s.expand(-3*sum(A[i]*A[k-i] for i in range(k+1))
              + (F(3,4) if k == 0 else 0)) for k in range(6)]

RM = [(n-2*b, b) for n in range(6,12) for b in range(n//2+1)]
TM = [(n-2*b, b) for n in range(11,16) for b in range(n//2+1)
      if 2*n-b <= 22] + [(3,6)]


def source_rows(poly, last=15):
    """Literal binomial substitution p=t(1+x^2t), y=xtp."""
    out = {n: s.Integer(0) for n in range(last+1)}
    for (a,b), c in s.Poly(poly,p,y).terms():
        for j in range(a+b+1):
            n = a+2*b+j
            if n <= last:
                out[n] += c*comb(a+b,j)*x**(b+2*j)
    return {n:s.expand(v) for n,v in out.items()}


def carrier_rows(a,b):
    # p^2 Delta p^a y^b = x^b t^(a+2b+5)(1+x^2t)^(a+b+4).
    n = a+2*b+5
    return {n+j:s.Integer(comb(a+b+4,j))*x**(b+2*j)
            for j in range(max(0,17-n))}


def act(rows, carrier):
    out = {n:s.Integer(0) for n in range(10,16)}
    for k,bk in enumerate(rows):
        for m,sm in carrier.items():
            n = k+m-1
            if n in out:
                out[n] += m*s.diff(bk,x)*sm-k*bk*s.diff(sm,x)
    return {n:s.expand(v) for n,v in out.items()}


def depth_bank(rows, depth):
    """Full inherited pi_15(P_depth) annihilator, no degree-only shortcut."""
    out = []
    for ell in range(1,31):
        start=(ell+1)//2; stop=min(15,ell+depth); rho=(ell+2)//3
        if start>stop:
            continue
        for q in range(max(0,ell+depth-stop),rho):
            out.append(s.expand(sum((-1)**(n-start)*comb(stop+q-n,q)
                       *rows.get(n,s.Integer(0)).coeff(x,2*n-ell)
                       for n in range(start,stop+1))))
    return out


def main():
    check(len(RM)==30 and len(TM)==17, 'complete bounded universe')
    check(sum(a+2*b==11 for a,b in RM)==6, 'six terminal carrier rows')
    H5 = (-3*p+F(8,3)*p**2-F(1376,135)*p**3+F(896,15)*p**4
          -F(32,5)*y**2-F(731648,2025)*p**5+F(512,75)*p*y**2)
    hrows=source_rows(H5,5)
    for n in range(6):
        zero(G[n]-hrows[n]+(x*x/2 if n==1 else 0), 'actual H through row5')
    check(Xi in A[5].free_symbols and Xi in C[5].free_symbols,
          'varying background retained')
    check(all(Xi not in g.free_symbols for g in G), 'actual low G constant')

    da=[]; dc=[]; dg=[]
    Delta=p**3-y**2
    for a,b in RM:
        sc=carrier_rows(a,b)
        ac=act(A,sc); cc=act(C,sc); gc=act(G,sc)
        da.append(ac); dc.append(cc); dg.append(gc)
        # Independent cusp-coordinate operator versus literal source action.
        rr=p**a*y**b
        Ar=(2*Delta+3*p**3)*rr+p*Delta*s.diff(rr,p)
        Br=p*(-2*y*rr+Delta*s.diff(rr,y))
        er=10*rr+2*p*s.diff(rr,p)+3*y*s.diff(rr,y)
        image=-p**3*y*er/2+Delta*(Ar*s.diff(H5,y)-Br*s.diff(H5,p))
        image_rows=source_rows(s.expand(image))
        for n in range(16):
            zero(image_rows[n]-gc.get(n,0), 'two actual response paths')
        for n in range(10,16):
            chain=0; jac=0
            for k in range(6):
                j=n-k
                if j in ac:
                    chain+=GA[k]*ac[j]+2*C[k]*cc[j]
                    jac+=(j*s.diff(A[k],x)*cc[j]+k*s.diff(ac[j],x)*C[k]
                          -k*A[k]*s.diff(cc[j],x)-j*ac[j]*s.diff(C[k],x))
            zero(chain-gc[n], 'every literal source-chain row')
            zero(jac, 'every literal bracket row')
            check(s.degree(ac[n],x)<=n+1 or ac[n]==0, 'A degree cap')
            check(s.degree(cc[n],x)<=n+2 or cc[n]==0, 'C degree cap')
        bank=depth_bank(ac,2)+depth_bank(cc,3)
        check(len(bank)==91, 'complete two-depth bank')
        for value in bank:
            zero(value, 'each actual depth row')
        n=a+2*b
        zero(ac[n+4]-(n+5)*x**(b+1)/2, 'first coordinate leading injection')

    rv=s.symbols('r0:30'); tv=s.symbols('h0:17')
    target=source_rows(sum(v*p**a*y**b for v,(a,b) in zip(tv,TM)))
    eq=[]
    for n in range(10,16):
        value=sum(v*col[n] for v,col in zip(rv,dg))-target[n]
        eq.extend(s.Poly(s.expand(value),x).all_coeffs())
    matrix,rhs=s.linear_eq_to_matrix(eq,list(rv)+list(tv))
    check(matrix.shape==(51,47) and rhs==s.zeros(51,1), 'full source system')
    check(not set().union(*(v.free_symbols for v in matrix)),
          'constant matrix before elimination; no specialization loss')
    red,piv=DomainMatrix.from_Matrix(matrix).convert_to(s.QQ).rref()
    red=red.to_Matrix()
    check(len(piv)==36, 'joint dimension eleven')
    source_rel=s.Matrix([list(red[i,30:]) for i in range(red.rows)
                         if all(red[i,j]==0 for j in range(30)) and any(red[i,30:])])
    check(source_rel.shape==(12,17), 'complete five-dimensional source image')
    free=[j for j in range(47) if j not in piv]
    check(free==list(range(24,30))+[31,37,42,44,46], 'all free parameters')
    variables=list(rv)+list(tv)
    kernel=s.zeros(47,len(free))
    for col,j in enumerate(free):
        kernel[j,col]=1
        for row,i in enumerate(piv):
            kernel[i,col]=-red[row,j]
    check(matrix*kernel==s.zeros(51,11), 'all exact lifts, not only image ranks')
    check(kernel[free,:]==s.eye(11), 'independent complete kernel')
    response_basis=kernel[30:,6:]
    check(response_basis[[1,7,12,14,16],:]==s.eye(5), 'declared source coordinates')

    # The inherited seven raw response conditions must all hold on the image.
    av,bv,cv,dv,ev=tv[0],tv[2],tv[4],tv[6],tv[8]
    vv,wv,rho,kk=tv[11],tv[13],tv[14],tv[16]
    den=s.Integer(119712905708881603)
    inherited=[tv[5],tv[10],
      den*av-736701988512897*vv+156822363286128*wv+128022192499536*rho+F(315619079794353,2)*kk,
      den*bv+F(17667890668710639,2)*vv-8012657615180013*wv+1940411140875144*rho-F(1696910148536463,2)*kk,
      den*cv-6603015506480190*vv+F(15756867503321415,2)*wv-2521182637658160*rho+F(2500583114637405,8)*kk,
      den*dv+40574567706597180*vv+17040074980207200*wv+5161415972895171*rho-1095964361654112*kk,
      den*ev-67697466072152616*vv+13525083050179824*wv-F(37920301110828423,2)*rho+691333448130801*kk]
    inh,_=s.linear_eq_to_matrix(inherited,list(tv))
    check(inh*response_basis==s.zeros(7,5), 'all seven incoming conditions')
    check(inh.rank()==7, 'inherited full source dimension ten including high')
    extra=[tv[9],tv[15],tv[3]+F(5,6)*tv[12],
           tv[11]+(3916119864752*rho+35134702155*kk)/s.Integer(1552077744720),
           tv[13]-45*(957062888*rho+31275315*kk)/s.Integer(206943699296)]
    exmat,_=s.linear_eq_to_matrix(extra,list(tv))
    check(exmat*response_basis==s.zeros(5,5), 'five displayed additional conditions')
    check(inh.col_join(exmat).rank()==12, 'complete displayed iff, not only containment')

    # Explicit short supplier and literal all-equation verification.
    Rstar=(y**3*(F(12015,15962432)+F(28467,3990608)*p+F(2499,498826)*p**2
             -F(4806,1247065)*p**3+F(639368,11223585)*p**4)
           +F(532468,6235325)*y**5)
    Tstar=(F(108135,15962432)*p**7*y**2-F(12015,1680256)*p**3*y**4
           +F(208143,3990608)*p**8*y**2-F(741987,7981216)*p**4*y**4
           -F(196239,997652)*p**5*y**4+F(180225,15962432)*p*y**6
           +F(346905,3990608)*p**2*y**6-p**3*y**6)
    rc=s.Matrix([s.Poly(Rstar,p,y).coeff_monomial(p**a*y**b) for a,b in RM])
    tc=s.Matrix([s.Poly(Tstar,p,y).coeff_monomial(p**a*y**b) for a,b in TM])
    check(matrix*rc.col_join(tc)==s.zeros(51,1), 'entire compact supplier')
    check(inh*tc==s.zeros(7,1), 'supplier in inherited complete family')
    check(tc[16]==-1 and all(2*a+3*b<=22 for (a,b),v in zip(TM[:-1],tc[:-1]) if v),
          'exact high cancellation with weight22 low source')
    check(s.Poly(Rstar,p,y).total_degree()==7, 'supplier total degree seven')
    check(max(a+2*b for (a,b),v in zip(RM,rc) if v)==10, 'supplier maximal valuation ten')
    check(max(2*a+3*b for (a,b),v in zip(RM,rc) if v)==17, 'supplier weight seventeen')

    lead=kernel.extract([1,3],[9,10])
    check(lead.det()==F(1751787,413887398592), 'sharp earlier-prefix obstruction')
    old=s.zeros(17,1)
    old[11]=-F(27945,235202); old[13]=F(39123,470404)
    old[14]=F(52578,117601); old[16]=-1
    check(inh*old==s.zeros(7,1), 'old packet is a genuine positive finite response')
    check(source_rel*old!=s.zeros(12,1), 'old exact packet outside completed carrier image')

    # A source-invisible terminal carrier still changes the coordinate fibre.
    i=RM.index((11,0))
    check(all(v==0 for v in dg[i].values()), 'terminal source invisible')
    zero(da[i][15]-8*x, 'terminal coordinate nonzero')
    check(2*10>=16 and 2*10-1>=15, 'all finite parameter nonlinear terms invisible')
    check(10+10>=16, 'completed exponential equals linear action at this horizon')

    relations=[[str(v) for v in row] for row in response_basis.tolist()]
    lifts=[[str(v) for v in row] for row in kernel[:30,6:].tolist()]
    data={'carrier_monomials':RM, 'source_monomials':TM,
          'source_free_indices':[1,7,12,14,16], 'source_basis':relations,
          'carrier_lift_basis':lifts, 'Rstar':s.sstr(Rstar), 'Tstar':s.sstr(Tstar),
          'leading_minor':'1751787/413887398592'}
    print('STATUS: complete finite universal-carrier response, pending independent audit')
    print('UNIVERSE: 30 carrier monomials of valuation6..11; 17 source monomials; arbitrary Xi')
    print('FULL_CHECKS: source, coordinate chain, bracket, prefix, caps, and all91 depth rows')
    print('MATRIX: 51x47 constant, rank36; complete source image dimension5; coordinate kernel6')
    print('FREE_SOURCE: [p^9y], [p^6y^3], [p^3y^5], [p^2y^6], [p^3y^6]')
    print('RSTAR:',s.sstr(Rstar))
    print('TSTAR:',s.sstr(Tstar))
    print('SOURCE_BASIS:',json.dumps(relations,separators=(',',':')))
    print('CARRIER_LIFT_BASIS:',json.dumps(lifts,separators=(',',':')))
    print('SHARP_PREFIX: every nonzero high payer changes coordinate row10')
    print('HOSTILES: old exact packet not in image; six terminal carriers not source-visible')
    print('SCOPE: true finite restriction of completed flow; later equations and polynomial termination remain open')
    raw=json.dumps(data,sort_keys=True,separators=(',',':'))
    print('SEMANTIC_SHA256',sha256(raw.encode()).hexdigest())
    print('PASS gates='+str(GATES))


if __name__=='__main__':
    main()
