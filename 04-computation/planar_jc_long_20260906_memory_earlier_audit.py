#!/usr/bin/env python3
"""Independent valuation11 audit: solve bracket rows before literal depth tests.

The literal depth-generator helper is inherited from the audited valuation12
independent program. No runtime repository imports or producer output input.
The primary uses180 raw coordinate coefficients; this route uses81 complete
bracket tangent parameters and a different polynomial minor certificate.
"""
from math import comb
from hashlib import sha256
from pathlib import Path
import json
import sympy as S
from sympy.polys.matrices import DomainMatrix
CHECKS=0

def check(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def literal_depth_annihilator(depth,maxrow):
    coordinates = {}
    for n in range(maxrow+1):
        for power in range(n+depth+1):
            coordinates.setdefault(2*n-power,[]).append((n,power))
    columns = {ell:set() for ell in coordinates}
    for aa in range(depth+1):
        for bb in range(depth-aa+1):
            for cc in range(maxrow+1):
                for ee in range(maxrow//2+1):
                    first = bb+cc+2*ee
                    if first>maxrow:
                        continue
                    ell = 2*cc+3*ee-aa
                    values = {coordinate:0 for coordinate in coordinates[ell]}
                    for j in range(cc+ee+1):
                        n,power = first+j,aa+2*bb+ee+2*j
                        if n<=maxrow:
                            check((n,power) in values,'literal depth generator cap')
                            values[n,power] = comb(cc+ee,j)
                    columns[ell].add(tuple(values[coordinate] for coordinate in coordinates[ell]))
    answer = []
    for ell in sorted(coordinates):
        bank = sorted(columns[ell])
        if bank:
            mat = S.Matrix(bank).T
            lefts = mat.T.nullspace()
            for left in lefts:
                check(left.T*mat==S.zeros(1,mat.cols),'literal generator annihilator')
        else:
            lefts = list(S.eye(len(coordinates[ell])).columnspace())
        for left in lefts:
            answer.append((coordinates[ell],list(left)))
    return answer


def main():
    x,Xi=S.symbols('x xi10')
    A=[1+x*x/4,S.Rational(4,3)+2*x*x,
       -S.Rational(32,9)-S.Rational(4,5)*x*x,
       S.Rational(2176,135)+S.Rational(64,9)*x*x-S.Rational(32,15)*x**4,
       S.Rational(224,9)*x**4-S.Rational(256,75)*x*x-S.Rational(37376,405),
       S.Rational(128,9)*x**6-S.Rational(8576,75)*x**4-S.Rational(6,11)*Xi*x*x
       -S.Rational(203776,7425)*x*x-x/2+S.Rational(3528704,6075)]
    C=[-3*x/4-x**3/8,-4*x-S.Rational(3,2)*x**3,
       S.Rational(88,15)*x-S.Rational(12,5)*x**3,
       -S.Rational(1408,45)*x+S.Rational(64,15)*x**3+S.Rational(8,5)*x**5,
       -S.Rational(184,15)*x**5-S.Rational(2048,25)*x**3+S.Rational(98944,675)*x,
       -S.Rational(32,3)*x**7-S.Rational(128,5)*x**5+S.Rational(9,22)*Xi*x**3
       +S.Rational(130816,275)*x**3+S.Rational(3,8)*x*x+S.Rational(9,11)*Xi*x
       -S.Rational(4024832,4455)*x+S.Rational(3,4)]
    # Independent sixth-row inheritance: THM4308 equations(6),(28),(32),(33).
    # No valuation>=7 omitted source face enters this predicted source row.
    bracket5=sum((5-i)*S.diff(A[i],x)*C[5-i]-i*A[i]*S.diff(C[5-i],x)
                 for i in range(1,5))
    check(S.expand(5*(S.diff(A[0],x)*C[5]-S.diff(C[0],x)*A[5])+bracket5)==0,
          'sixth background solves bracket5')
    bracket6=sum((6-i)*S.diff(A[i],x)*C[6-i]-i*A[i]*S.diff(C[6-i],x)
                 for i in range(1,6))
    source6=sum(C[i]*C[6-i] for i in range(1,6))-sum(
        A[i]*A[j]*A[k] for i in range(6) for j in range(6) for k in range(6) if i+j+k==6)
    U=(S.Integer(475515904)-109350*Xi)/200475
    target6=U+x+(-S.Rational(731648,405)+Xi)*x*x+(
        6*S.Rational(896,15)+3*S.Rational(512,75))*x**4+(-S.Rational(32,5)-S.Rational(1376,135))*x**6
    check(S.expand(source6+(x*x+6)*bracket6/12-target6)==0,
          'sixth background solves predicted source6 for every xi10')
    z,h,ss=S.symbols('z h s')
    actualxi=(108391820625*h+3765431711250*z+S.Integer(324974300165767168))/24542012448000
    boundary=801*h+27826*z+S.Rational(85855050266495746048,37533020625)
    slope=S.Rational(1485,269322496);intercept=S.Rational(34969998848,55604475)
    check(S.expand(actualxi-slope*boundary-intercept)==0,'recovered graph factors through actual Gm coordinate')
    check(S.diff(A[5].subs(Xi,slope*ss+intercept),ss)!=0,'sixth row retains varying background')
    da,dc,theta={},{},[]
    for n in range(10,16):
        debt=0
        for k in range(1,6):
            j=n-k
            if j in da:
                debt+=(j*S.diff(A[k],x)*dc[j]+k*S.diff(da[j],x)*C[k]
                       -k*A[k]*S.diff(dc[j],x)-j*da[j]*S.diff(C[k],x))
        f=S.expand(-debt/n)
        pa=S.Rational(4,3)*f.subs(x,0)
        pc=S.cancel(2*(f-S.Rational(3,8)*(x*x+2)*pa)/x)
        check(S.denom(pc)==1,'polynomial bracket particular')
        variables=S.symbols('theta%d_0:%d'%(n,n+1));theta.extend(variables)
        th=sum(v*x**j for j,v in enumerate(variables))
        da[n]=S.expand(pa+S.diff(A[0],x)*th)
        dc[n]=S.expand(pc+S.diff(C[0],x)*th)
        check(S.degree(da[n],x)<=n+1 and S.degree(dc[n],x)<=n+2,'coordinate caps')
        check(S.expand(n*(S.diff(A[0],x)*dc[n]-S.diff(C[0],x)*da[n])+debt)==0,'whole bracket row')
    check(len(theta)==81,'complete tangent variable count')
    # Enumerate by a rectangular lattice rather than valuation slices.
    low=sorted([(a,b) for a in range(12) for b in range(8)
                if 11<=a+2*b<=15 and 2*a+3*b<=23],key=lambda q:(q[0]+2*q[1],q[1]))
    monomials=low+[(3,6)];rs=S.symbols('r0:21')
    check(len(low)==20,'whole source universe')
    source={n:S.Integer(0) for n in range(10,16)}
    for rv,(a,b) in zip(rs,monomials):
        for j in range(a+b+1):
            n=a+2*b+j
            if n in source:source[n]+=rv*comb(a+b,j)*x**(b+2*j)
    ga=[S.expand(-3*sum(A[i]*A[k-i] for i in range(k+1))+
                 (S.Rational(3,4) if k==0 else 0)) for k in range(6)]
    eq=[]
    for n in range(10,16):
        residual=-source[n]
        for k in range(6):
            j=n-k
            if j in da:residual+=ga[k]*da[j]+2*C[k]*dc[j]
        residual=S.expand(residual)
        if n==10:check(residual==0,'first source row vanishes')
        else:
            check(S.degree(residual,x)<=n,'predicted source cap')
            eq.extend(residual.coeff(x,j) for j in range(n+1))
    check(len(eq)==70,'all source equations')
    for depth,rows in [(2,da),(3,dc)]:
        for coordinates,weights in literal_depth_annihilator(depth,15):
            eq.append(sum(w*rows.get(n,S.Integer(0)).coeff(x,power)
                          for (n,power),w in zip(coordinates,weights)))
    check(len(eq)==161,'literal depth codimension91')
    matrix,rhs=S.linear_eq_to_matrix(eq,theta+list(rs))
    check(rhs==S.zeros(161,1),'homogeneous complete system')
    M=matrix[:,:81];Z=matrix[:,81:]
    dm,piv=DomainMatrix.from_Matrix(matrix).convert_to(S.QQ.frac_field(Xi)).rref()
    red=dm.to_Matrix()
    raw=[list(red[i,81:]) for i in range(red.rows)
         if all(red[i,j]==0 for j in range(81)) and any(red[i,81:])]
    dual=S.Matrix(raw)
    check((sum(j<81 for j in piv),len(raw))==(71,7),'independent generic ranks')
    check(not set().union(*(v.free_symbols for v in dual)),'constant response relations')
    # Separate polynomial certificate for every background, using this
    # smaller tangent-coordinate matrix and no generic-denominator argument.
    tangentcols=[j for j in piv if j<81];Mc=M[:,tangentcols]
    _,rowids=DomainMatrix.from_Matrix(Mc.subs(Xi,0).T).convert_to(S.QQ).rref()
    rows=list(rowids);minor=Mc[rows,:]
    m0=minor.subs(Xi,0)
    inv0=DomainMatrix.from_Matrix(m0).convert_to(S.QQ).inv().to_Matrix()
    check(minor.diff(Xi,2)==S.zeros(71,71),'affine smaller tangent minor')
    K=inv0*minor.diff(Xi)
    check(K*K==S.zeros(71,71),'nilpotent smaller variation')
    inverse=(S.eye(71)-Xi*K)*inv0
    check((inverse*minor).applyfunc(S.expand)==S.eye(71),'polynomial inverse all backgrounds')
    check((M-Mc*inverse*M[rows,:]).applyfunc(S.expand)==S.zeros(*M.shape),'all tangent columns span uniformly')
    F=(Z-Mc*inverse*Z[rows,:]).applyfunc(S.expand)
    check(not set().union(*(v.free_symbols for v in F)),'entire response constant')
    epiv=list(dual.rref()[1])
    check((F-F[:,epiv]*dual)==S.zeros(*F.shape),'complete source cokernel')
    _,srows=DomainMatrix.from_Matrix(F[:,epiv].T).convert_to(S.QQ).rref()
    determinant=F.extract(list(srows),epiv).det(method='domain-ge')
    check(determinant!=0 and determinant.is_Rational,'nonzero constant source minor')
    for w in range(17,24):
        idx=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=w]
        ranks=(dual[:,idx].rank(),dual[:,idx].row_join(dual[:,-1]).rank())
        check(ranks=={17:(1,2),18:(3,4),19:(3,4),20:(6,7),21:(6,7),22:(7,7),23:(7,7)}[w],
              'sharp weight boundary '+str(w))
    w22=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=22]
    free=[i for i in w22 if i not in epiv]
    check(len(free)==9 and set(epiv)<=set(w22),'nine free source coordinates')
    bank=S.zeros(21,10);bank[20,0]=-1
    for col,i in enumerate(free,1):bank[i,col]=1
    for row,i in enumerate(epiv):bank[i,:]=-dual[row,:]*bank
    check(dual*bank==S.zeros(7,10),'entire affine nine source family')
    lifts=S.zeros(81,10)
    smalllifts=-(inverse*Z[rows,:]*bank).applyfunc(S.expand)
    for row,i in enumerate(tangentcols):lifts[i,:]=smalllifts[row,:]
    for v in (matrix*lifts.col_join(bank)).applyfunc(S.expand):check(v==0,'all raw equations of all ten affine columns')
    check(all(S.denom(v).is_Integer for v in lifts),'all family lifts polynomial in background')
    # Uniform ten-dimensional homogeneous tangent kernel: free columns
    # are completed by the same fixed polynomial inverse.
    freecoords=[j for j in range(81) if j not in tangentcols]
    kernel=S.zeros(81,10)
    for col,j in enumerate(freecoords):kernel[j,col]=1
    smallkernel=-(inverse*M[rows,freecoords]).applyfunc(S.expand)
    for row,i in enumerate(tangentcols):kernel[i,:]=smallkernel[row,:]
    check((M*kernel).applyfunc(S.expand)==S.zeros(161,10),'complete homogeneous fibre all backgrounds')
    for j in (5,11):
        singleton=S.zeros(1,21);singleton[0,j]=1
        check(S.Matrix.vstack(dual,singleton).rank()==7,'forced-zero hostile '+str(monomials[j]))
    old=S.zeros(21,1)
    old[13]=-S.Rational(27945,235202);old[15]=S.Rational(39123,470404)
    old[17]=S.Rational(52578,117601);old[20]=-1
    check(dual*old==S.zeros(7,1),'inherited valuation13 transport')
    print('UNIVERSE tangent81 source21 source_equations70 literal_depth91 symbolic_xi10')
    print('RANK tangent71 source7 terminal10 uniform_every_characteristic_zero_xi10')
    print('RELATIONS',dual.tolist())
    print('UNIFORM_MINOR tangent71 K_square_zero; source7 determinant',determinant)
    print('FAMILY affine9 source; ten affine columns satisfy all161 raw equations; separate affine10 coordinate fibre')
    print('HOSTILES py5,y6 forced_zero; weight21_rank6_to7; weight22_rank7_to7')
    print('RELATION_SHA256',sha256(json.dumps([[str(v) for v in row] for row in dual.tolist()],separators=(',',':')).encode()).hexdigest())
    print('CHECKS',CHECKS,'PASS')
    print('SOURCE_SHA256',sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':main()
