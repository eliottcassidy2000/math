#!/usr/bin/env python3
"""Independent tangent-row / literal-generator audit of valuation-12 memory.

No repository mathematical imports, no producer output input. The source
uses 70 complete tangent parameters after solving bracket rows, whereas
the primary computation uses all 155 raw coordinate coefficients.
"""
from math import comb
from hashlib import sha256
from pathlib import Path
import json

import sympy as S
from sympy.polys.matrices import DomainMatrix

CHECKS = 0


def check(ok,label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


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
    x = S.Symbol('x')
    A = [1+x*x/4,S.Rational(4,3)+2*x*x,
         -S.Rational(32,9)-S.Rational(4,5)*x*x,
         S.Rational(2176,135)+S.Rational(64,9)*x*x-S.Rational(32,15)*x**4,
         S.Rational(224,9)*x**4-S.Rational(256,75)*x*x-S.Rational(37376,405)]
    C = [-3*x/4-x**3/8,-4*x-S.Rational(3,2)*x**3,
         S.Rational(88,15)*x-S.Rational(12,5)*x**3,
         -S.Rational(1408,45)*x+S.Rational(64,15)*x**3+S.Rational(8,5)*x**5,
         -S.Rational(184,15)*x**5-S.Rational(2048,25)*x**3+S.Rational(98944,675)*x]
    # Independently recover the fifth background row from THM4308's row5
    # source at the boundary, before any valuation>=7 source face can enter.
    # The m5 Student operator is injective, so these two identities determine
    # A4,C4 uniquely among the degree-capped bracket solutions.
    bracket4 = sum((4-i)*S.diff(A[i],x)*C[4-i]-i*A[i]*S.diff(C[4-i],x)
                   for i in range(1,4))
    check(S.expand(4*(S.diff(A[0],x)*C[4]-S.diff(C[0],x)*A[4])+bracket4)==0,
          'inherited fifth row solves bracket4')
    bracket5 = sum((5-i)*S.diff(A[i],x)*C[5-i]-i*A[i]*S.diff(C[5-i],x)
                   for i in range(1,5))
    source5 = (sum(C[i]*C[5-i] for i in range(1,5))-
               sum(A[i]*A[j]*A[k] for i in range(5) for j in range(5) for k in range(5) if i+j+k==5))
    target5 = -S.Rational(731648,2025)+S.Rational(6144,25)*x*x-S.Rational(1952,45)*x**4
    check(S.expand(source5+(x*x+6)*bracket5/10-target5)==0,
          'inherited fifth row solves predicted source5')
    da,dc,theta = {},{},[]
    for n in range(11,16):
        # Earlier perturbation rows determine the inhomogeneous bracket debt.
        # This particular solution uses only its constant term in delta A;
        # it differs from the inherited pivot-matrix particular solution.
        debt = 0
        for k in range(1,5):
            j = n-k
            if j in da:
                debt += (j*S.diff(A[k],x)*dc[j]+k*S.diff(da[j],x)*C[k]
                         -k*A[k]*S.diff(dc[j],x)-j*da[j]*S.diff(C[k],x))
        f = S.expand(-debt/n)
        pa = S.Rational(4,3)*f.subs(x,0)
        pc = S.cancel(2*(f-S.Rational(3,8)*(x*x+2)*pa)/x)
        check(S.denom(pc)==1,'row particular has no x denominator')
        variables = S.symbols('theta%d_0:%d' % (n,n+1))
        theta.extend(variables)
        th = sum(v*x**j for j,v in enumerate(variables))
        da[n] = S.expand(pa+S.diff(A[0],x)*th)
        dc[n] = S.expand(pc+S.diff(C[0],x)*th)
        check(S.degree(da[n],x)<=n+1 and S.degree(dc[n],x)<=n+2,'complete tangent degree caps')
        check(S.expand(n*(S.diff(A[0],x)*dc[n]-S.diff(C[0],x)*da[n])+debt)==0,
              'entire bracket row solved')
    check(len(theta)==70,'complete tangent variable count')
    low = [(v-2*b,b) for v in range(12,16) for b in range(v//2+1) if 2*v-b<=23]
    monomials = low+[(3,6)]
    rs = S.symbols('r0:'+str(len(monomials)))
    check(len(low)==14 and len(monomials)==15,'complete visible source packet')
    source = {n:S.Integer(0) for n in range(11,16)}
    for rv,(a,b) in zip(rs,monomials):
        for j in range(a+b+1):
            n=a+2*b+j
            if n in source:
                source[n] += rv*comb(a+b,j)*x**(b+2*j)
    ga = [S.expand(-3*sum(A[i]*A[k-i] for i in range(k+1))+
                   (S.Rational(3,4) if k==0 else 0)) for k in range(5)]
    eq=[]
    for n in range(11,16):
        residual = -source[n]
        for k in range(5):
            j=n-k
            if j in da:
                residual += ga[k]*da[j]+2*C[k]*dc[j]
        residual=S.expand(residual)
        if n==11:
            check(residual==0,'first source row identically unchanged')
        else:
            check(S.degree(residual,x)<=n,'actual predicted-source row cap')
            eq.extend(residual.coeff(x,j) for j in range(n+1))
    check(len(eq)==58,'full source compatibilities after bracket solution')
    depth_count=0
    for depth,rows in [(2,da),(3,dc)]:
        for coordinates,weights in literal_depth_annihilator(depth,15):
            eq.append(sum(weight*rows.get(n,S.Integer(0)).coeff(x,power)
                          for (n,power),weight in zip(coordinates,weights)))
            depth_count+=1
    check(depth_count==91,'full literal-generator depth codimension')
    matrix,rhs = S.linear_eq_to_matrix(eq,theta+list(rs))
    check(rhs==S.zeros(len(eq),1),'homogeneous entire joint system')
    dm,piv = DomainMatrix.from_Matrix(matrix).convert_to(S.QQ).rref()
    reduced = dm.to_Matrix()
    tangent_rank = sum(p<len(theta) for p in piv)
    raw = [list(reduced[i,len(theta):]) for i in range(reduced.rows)
           if all(reduced[i,j]==0 for j in range(len(theta))) and any(reduced[i,len(theta):])]
    dual = S.Matrix(raw)
    expected = S.zeros(5,15)
    expected[0,1]=1;expected[0,11]=S.Rational(85731129,383445497);expected[0,14]=S.Rational(38329362,383445497)
    expected[1,3]=1;expected[1,11]=-S.Rational(200039301,766890994);expected[1,14]=-S.Rational(44717589,383445497)
    expected[2,5]=1
    expected[3,7]=1;expected[3,11]=-S.Rational(112896720,383445497);expected[3,14]=-S.Rational(192065985,766890994)
    expected[4,9]=1;expected[4,11]=-S.Rational(217326816,383445497);expected[4,14]=-S.Rational(261093861,1533781988)
    check(dual==expected,'all five source relation rows independently recovered')
    check((tangent_rank,len(piv),len(theta)-tangent_rank)==(60,65,10),'independent complete ranks')
    print('UNIVERSE tangent70 source15 source_rows58 literal_depth91')
    print('RANK tangent60 source5 joint65 terminal10')
    print('MONOMIALS',monomials)
    print('RELATIONS',dual.tolist())
    for w in range(18,24):
        idx=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=w]
        before=dual[:,idx].rank();after=dual[:,idx].row_join(dual[:,-1]).rank()
        check((before,after)=={18:(1,2),19:(1,2),20:(3,4),21:(3,4),22:(5,5),23:(5,5)}[w],
              'sharp weight boundary '+str(w))
        print('WEIGHT',w,before,after)
    idx22=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=22]
    check(len(idx22)-dual[:,idx22].rank()==5,'complete replacement family is affine five')
    odd22=[i for i in idx22 if monomials[i][1]%2]
    check(len(odd22)==4 and dual[:,odd22]==S.zeros(5,4),'four genuinely neutral odd parameters')
    # Six source columns: one affine base replacement and all five independent
    # directions. A single full raw substitution verifies the entire family.
    bank=S.zeros(15,6)
    bank[14,0]=-1
    bank[1,0]=S.Rational(38329362,383445497)
    bank[3,0]=-S.Rational(44717589,383445497)
    bank[7,0]=-S.Rational(192065985,766890994)
    bank[9,0]=-S.Rational(261093861,1533781988)
    for column,index in enumerate(odd22,1):bank[index,column]=1
    bank[11,5]=1
    bank[1,5]=-S.Rational(85731129,383445497)
    bank[3,5]=S.Rational(200039301,766890994)
    bank[7,5]=S.Rational(112896720,383445497)
    bank[9,5]=S.Rational(217326816,383445497)
    check(dual*bank==S.zeros(5,6),'every affine family column retains all source relations')
    inherited=S.zeros(15,1)
    inherited[7]=-S.Rational(27945,235202)
    inherited[9]=S.Rational(39123,470404)
    inherited[11]=S.Rational(52578,117601)
    inherited[14]=-1
    check(bank[:,0]+S.Rational(52578,117601)*bank[:,5]==inherited,
          'valuation13 replacement recovered on the exact even slice')
    lifts=S.zeros(70,6)
    for row,pivot in enumerate(piv):
        if pivot<70:
            lifts[pivot,:]=-reduced[row,70:]*bank
    for value in matrix*lifts.col_join(bank):
        check(value==0,'every raw compatibility for entire affine family')
    # The terminal homogeneous fibre consists of ten independent complete
    # tangent responses, not a count inferred from the source quotient alone.
    tangent_matrix=matrix[:,:70]
    kernel=tangent_matrix.nullspace()
    check(len(kernel)==10,'complete terminal tangent kernel')
    for vector in kernel:
        check(tangent_matrix*vector==S.zeros(matrix.rows,1),'all terminal kernel equations')
    print('FAMILY affine5=(four_odd,one_even); every fibre affine10')
    print('PASS checks='+str(CHECKS))
    semantic={'monomials':monomials,'dual':[[str(v) for v in row] for row in dual.tolist()],
              'ranks':[60,65,10],'equations':[58,91]}
    print('SEMANTIC_SHA256',sha256(json.dumps(semantic,sort_keys=True,separators=(',',':')).encode()).hexdigest())
    print('SOURCE_SHA256',sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':main()
