#!/usr/bin/env python3
"""Complete earlier-memory probe, valuation>=11 through row15.

FINITE-EXACT exploration until the result and its global boundary consumer
are independently audited. No repository mathematical imports.
"""
import sympy as S
from sympy.polys.matrices import DomainMatrix
from math import comb
from hashlib import sha256
import json

CHECKS=0
def check(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def build_system(last=15):
    x=S.Symbol('x'); Xi=S.Symbol('xi10')
    A=[1+x*x/4,S.Rational(4,3)+2*x*x,
       -S.Rational(32,9)-S.Rational(4,5)*x*x,
       S.Rational(2176,135)+S.Rational(64,9)*x*x-S.Rational(32,15)*x**4,
       S.Rational(224,9)*x**4-S.Rational(256,75)*x*x-S.Rational(37376,405)]
    C=[-3*x/4-x**3/8,-4*x-S.Rational(3,2)*x**3,
       S.Rational(88,15)*x-S.Rational(12,5)*x**3,
       -S.Rational(1408,45)*x+S.Rational(64,15)*x**3+S.Rational(8,5)*x**5,
       -S.Rational(184,15)*x**5-S.Rational(2048,25)*x**3+S.Rational(98944,675)*x]
    A.append(S.Rational(128,9)*x**6-S.Rational(8576,75)*x**4
             -S.Rational(6,11)*Xi*x*x-S.Rational(203776,7425)*x*x
             -x/2+S.Rational(3528704,6075))
    C.append(-S.Rational(32,3)*x**7-S.Rational(128,5)*x**5
             +S.Rational(9,22)*Xi*x**3+S.Rational(130816,275)*x**3
             +S.Rational(3,8)*x*x+S.Rational(9,11)*Xi*x
             -S.Rational(4024832,4455)*x+S.Rational(3,4))
    da,dc,unknown={},{},[]
    for n in range(10,last+1):
        va=S.symbols(f'a{n}_0:{n+2}');vc=S.symbols(f'c{n}_0:{n+3}')
        da[n]=sum(v*x**j for j,v in enumerate(va))
        dc[n]=sum(v*x**j for j,v in enumerate(vc))
        unknown.extend(va+vc)
    low=[(v-2*b,b) for v in range(11,16) for b in range(v//2+1)
         if 2*v-b<=23]
    monomials=low+[(3,6)]
    rs=S.symbols(f'r0:{len(monomials)}')
    source={n:S.Integer(0) for n in range(10,last+1)}
    for r,(a,b) in zip(rs,monomials):
        for j in range(a+b+1):
            n=a+2*b+j
            if n in source: source[n]+=r*comb(a+b,j)*x**(b+2*j)
    ga=[S.expand(-3*sum(A[i]*A[k-i] for i in range(k+1))
                +(S.Rational(3,4) if k==0 else 0)) for k in range(6)]
    eq=[]
    for n in range(10,last+1):
        pe=-source[n];je=0
        for k in range(6):
            j=n-k
            if j not in da:continue
            pe+=ga[k]*da[j]+2*C[k]*dc[j]
            je+=(j*S.diff(A[k],x)*dc[j]+k*S.diff(da[j],x)*C[k]
                 -k*A[k]*S.diff(dc[j],x)-j*da[j]*S.diff(C[k],x))
        eq.extend(S.Poly(S.expand(pe),x).all_coeffs())
        eq.extend(S.Poly(S.expand(je),x).all_coeffs())
    for d,v in [(2,da),(3,dc)]:
        for ell in range(1,31):
            start=(ell+1)//2;stop=min(last,ell+d);rho=(ell+2)//3
            if start>stop:continue
            for q in range(max(0,ell+d-stop),rho):
                eq.append(sum((-1)**(n-start)*comb(stop+q-n,q)
                              *S.expand(v[n]).coeff(x,2*n-ell)
                              for n in range(start,stop+1) if n in v))
    matrix,rhs=S.linear_eq_to_matrix(eq,unknown+list(rs))
    return x,Xi,A,C,da,dc,unknown,low,monomials,rs,eq,matrix,rhs

def main():
    x,Xi,A,C,da,dc,unknown,low,monomials,rs,eq,matrix,rhs=build_system()
    print('BACKGROUND A5',A[5],flush=True)
    print('BACKGROUND C5',C[5],flush=True)
    check(rhs==S.zeros(len(eq),1),'homogeneous literal system')
    print("MATRIX",matrix.shape,flush=True)
    dm,piv=DomainMatrix.from_Matrix(matrix).convert_to(S.QQ.frac_field(Xi)).rref()
    red=dm.to_Matrix();nu=len(unknown)
    raw=[list(red[i,nu:]) for i in range(red.rows)
         if all(red[i,j]==0 for j in range(nu)) and any(red[i,nu:])]
    dual=S.Matrix(raw)
    print('VALUATION11 RAW',len(eq),'UNKNOWN',nu,'SOURCE',len(rs))
    print('MONOMIALS',[(a,b,a+2*b,2*a+3*b) for a,b in monomials])
    print('TANGENT_RANK',sum(p<nu for p in piv),'SOURCE_RANK',len(raw))
    print('RELATIONS',dual.tolist(),flush=True)
    print('RELATION_SYMBOLS',set().union(*(v.free_symbols for v in dual)),flush=True)
    check((len(eq),nu,len(rs))==(301,180,21),'complete literal universe')
    check(sum(p<nu for p in piv)==170 and dual.rows==7,'full ranks')
    check(not set().union(*(v.free_symbols for v in dual)),'background-independent source relations')
    check(low==[(a,b) for v in range(11,16) for b in range(v//2+1)
                for a in [v-2*b] if 2*a+3*b<=23],'complete low source universe')
    # All-specialization certificate: invert a fixed tangent minor over
    # Q[xi10], without dividing by a polynomial in the background.
    M=matrix[:,:nu]; Z=matrix[:,nu:]
    tangent_cols=[p for p in piv if p<nu]
    Mc=M[:,tangent_cols]
    Mc0=Mc.subs(Xi,0)
    _,row_ids=DomainMatrix.from_Matrix(Mc0.T).convert_to(S.QQ).rref()
    minor0=Mc0[list(row_ids),:]
    inverse0=DomainMatrix.from_Matrix(minor0).convert_to(S.QQ).inv().to_Matrix()
    check(inverse0*minor0==S.eye(len(tangent_cols)),'rational tangent inverse')
    variation=Mc.diff(Xi)[list(row_ids),:]
    check(Mc.diff(Xi,2)==S.zeros(*Mc.shape),'affine background matrix')
    K=inverse0*variation
    check(K*K==S.zeros(*K.shape),'nilpotent background correction')
    inverse=(S.eye(K.rows)-Xi*K)*inverse0
    check((inverse*Mc[list(row_ids),:]).applyfunc(S.expand)==S.eye(K.rows),'polynomial tangent inverse')
    check((M-Mc*inverse*M[list(row_ids),:]).applyfunc(S.expand)==S.zeros(*M.shape),
          'full tangent column span over polynomial background ring')
    F=(Z-Mc*inverse*Z[list(row_ids),:]).applyfunc(S.expand)
    check(all(Xi not in v.free_symbols for v in F),'whole residual response constant')
    epiv=list(dual.rref()[1])
    check((F-F[:,epiv]*dual).applyfunc(S.expand)==S.zeros(*F.shape),'constant full source row-space containment')
    _,source_rows=DomainMatrix.from_Matrix(F.subs(Xi,0)[:,epiv].T).convert_to(S.QQ).rref()
    small=F.extract(list(source_rows),epiv)
    small_det=S.factor(small.det(method='domain-ge'))
    check(Xi not in small_det.free_symbols and small_det!=0,'constant nonzero source minor')
    print('ALL_SPECIALIZATIONS tangent_minor',minor0.shape,'K_square_zero',True,
          'source_minor',small.shape,'source_det',small_det,flush=True)
    print('BACKGROUND_DEPENDENCE raw_response_degree',max(S.degree(v,Xi) if v else 0 for v in F),flush=True)

    for w in range(17,24):
        idx=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=w]
        lower=dual[:,idx]
        before=lower.rank();after=lower.row_join(dual[:,-1]).rank()
        print('WEIGHT',w,'rank',before,'with_high',after,flush=True)
        check((before,after)=={17:(1,2),18:(3,4),19:(3,4),20:(6,7),21:(6,7),22:(7,7),23:(7,7)}[w],
              'sharp filtered response weight '+str(w))
        if before==after:
            ans=S.linsolve((lower,dual[:,-1]))
            print('REPLACEMENTS r_high=-1',[(monomials[i]) for i in idx],ans)
    # Positive actual normalized transport inherited from the later filter.
    rvec=S.zeros(len(rs),1)
    rvec[13]=-S.Rational(27945,235202);rvec[15]=S.Rational(39123,470404)
    rvec[17]=S.Rational(52578,117601);rvec[20]=-1
    check(dual*rvec==S.zeros(7,1),'inherited actual transport retained')
    solved=-(inverse*Z[list(row_ids),:]*rvec).applyfunc(S.expand)
    allcoords=S.zeros(nu,1)
    for i,column in enumerate(tangent_cols):allcoords[column]=solved[i]
    residual=(matrix*S.Matrix(list(allcoords)+list(rvec))).applyfunc(S.expand)
    for row,value in enumerate(residual):check(value==0,'literal raw coefficient '+str(row))
    check(all(S.denom(v).is_Integer for v in allcoords),'polynomial lift no background denominator')
    odd=S.zeros(1,len(rs));odd[0,5]=1
    check(S.Matrix.vstack(dual,odd).rank()==7,'p*y5 forced zero')
    check(dual[:,5]!=S.zeros(7,1),'odd parity alone does not imply neutrality')
    for last in (13,14):
        *_,rr,ee,mm,bb=build_system(last)
        n_unknown=mm.cols-len(rs)
        reduced,_=DomainMatrix.from_Matrix(mm).convert_to(S.QQ).rref()
        reduced=reduced.to_Matrix()
        rows=[list(reduced[i,n_unknown:]) for i in range(reduced.rows)
              if all(reduced[i,j]==0 for j in range(n_unknown)) and any(reduced[i,n_unknown:])]
        response=S.Matrix(rows)
        if last==13:
            check(response[:,5]==S.zeros(response.rows,1),'odd py5 is fully neutral through row13')
        else:
            check(S.Matrix.vstack(response,odd).rank()==response.rank(),'odd py5 forced zero by row14')
        print('ODD_ONSET row',last,'response_rank',response.rank(),'py5_column_zero',response[:,5]==S.zeros(response.rows,1))
    even=S.zeros(1,len(rs));even[0,11]=1
    check(S.Matrix.vstack(dual,even).rank()==7,'old y6 depth condition retained')
    # Five even solved coordinates and two zero coordinates. The remaining
    # nine weight<=22 coordinates are free, for each fixed high coefficient.
    w22=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=22]
    check(len(w22)==16 and len(w22)-dual[:,w22].rank()==9,'complete weight22 affine nine source fibre')
    semantic={'relations':[[str(v) for v in row] for row in dual.tolist()],
              'monomials':monomials,'minor_rows':list(row_ids),
              'source_rows':list(source_rows),'source_det':str(small_det)}
    print('SEMANTIC_SHA256',sha256(json.dumps(semantic,sort_keys=True,separators=(',',':')).encode()).hexdigest())
    print('CHECKS',CHECKS,'PASS')
    print('SCOPE all_characteristic_zero_background_xi10;valuation>=11;rows<=15;minweight22;not_lower_valuations_or_termination_or_JC2')

if __name__=='__main__':main()
