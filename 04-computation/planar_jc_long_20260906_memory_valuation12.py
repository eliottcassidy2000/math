#!/usr/bin/env python3
"""Complete earlier-memory certificate, valuation>=12 through row15.

Exact full-coefficient path, with live controls and an independent tangent
reconstruction. Universal scope is fixed in the companion proof report.
No repository mathematical imports.
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

def main():
    x=S.Symbol('x')
    A=[1+x*x/4,S.Rational(4,3)+2*x*x,
       -S.Rational(32,9)-S.Rational(4,5)*x*x,
       S.Rational(2176,135)+S.Rational(64,9)*x*x-S.Rational(32,15)*x**4,
       S.Rational(224,9)*x**4-S.Rational(256,75)*x*x-S.Rational(37376,405)]
    C=[-3*x/4-x**3/8,-4*x-S.Rational(3,2)*x**3,
       S.Rational(88,15)*x-S.Rational(12,5)*x**3,
       -S.Rational(1408,45)*x+S.Rational(64,15)*x**3+S.Rational(8,5)*x**5,
       -S.Rational(184,15)*x**5-S.Rational(2048,25)*x**3+S.Rational(98944,675)*x]
    da,dc,unknown={},{},[]
    for n in range(11,16):
        va=S.symbols(f'a{n}_0:{n+2}');vc=S.symbols(f'c{n}_0:{n+3}')
        da[n]=sum(v*x**j for j,v in enumerate(va))
        dc[n]=sum(v*x**j for j,v in enumerate(vc))
        unknown.extend(va+vc)
    low=[(v-2*b,b) for v in range(12,16) for b in range(v//2+1)
         if 2*v-b<=23]
    monomials=low+[(3,6)]
    rs=S.symbols(f'r0:{len(monomials)}')
    source={n:S.Integer(0) for n in range(11,16)}
    for r,(a,b) in zip(rs,monomials):
        for j in range(a+b+1):
            n=a+2*b+j
            if n in source: source[n]+=r*comb(a+b,j)*x**(b+2*j)
    ga=[S.expand(-3*sum(A[i]*A[k-i] for i in range(k+1))
                +(S.Rational(3,4) if k==0 else 0)) for k in range(5)]
    eq=[]
    for n in range(11,16):
        pe=-source[n];je=0
        for k in range(5):
            j=n-k
            if j not in da:continue
            pe+=ga[k]*da[j]+2*C[k]*dc[j]
            je+=(j*S.diff(A[k],x)*dc[j]+k*S.diff(da[j],x)*C[k]
                 -k*A[k]*S.diff(dc[j],x)-j*da[j]*S.diff(C[k],x))
        eq.extend(S.Poly(S.expand(pe),x).all_coeffs())
        eq.extend(S.Poly(S.expand(je),x).all_coeffs())
    for d,v in [(2,da),(3,dc)]:
        for ell in range(1,31):
            start=(ell+1)//2;stop=min(15,ell+d);rho=(ell+2)//3
            if start>stop:continue
            for q in range(max(0,ell+d-stop),rho):
                eq.append(sum((-1)**(n-start)*comb(stop+q-n,q)
                              *S.expand(v[n]).coeff(x,2*n-ell)
                              for n in range(start,stop+1) if n in v))
    matrix,rhs=S.linear_eq_to_matrix(eq,unknown+list(rs))
    check(rhs==S.zeros(len(eq),1),'homogeneous equations')
    dm,piv=DomainMatrix.from_Matrix(matrix).convert_to(S.QQ).rref()
    red=dm.to_Matrix();nu=len(unknown)
    raw=[list(red[i,nu:]) for i in range(red.rows)
         if all(red[i,j]==0 for j in range(nu)) and any(red[i,nu:])]
    dual=S.Matrix(raw)
    check((len(eq),nu,len(rs))==(271,155,15),'complete raw universe')
    check(sum(p<nu for p in piv)==145 and len(raw)==5,'exact complete ranks')
    expected=S.zeros(5,15)
    expected[0,1]=1;expected[0,11]=S.Rational(85731129,383445497);expected[0,14]=S.Rational(38329362,383445497)
    expected[1,3]=1;expected[1,11]=-S.Rational(200039301,766890994);expected[1,14]=-S.Rational(44717589,383445497)
    expected[2,5]=1
    expected[3,7]=1;expected[3,11]=-S.Rational(112896720,383445497);expected[3,14]=-S.Rational(192065985,766890994)
    expected[4,9]=1;expected[4,11]=-S.Rational(217326816,383445497);expected[4,14]=-S.Rational(261093861,1533781988)
    check(dual==expected,'five full response relations')
    print('VALUATION12 RAW',len(eq),'UNKNOWN',nu,'SOURCE',len(rs))
    print('MONOMIALS',[(a,b,a+2*b,2*a+3*b) for a,b in monomials])
    print('TANGENT_RANK',sum(p<nu for p in piv),'SOURCE_RANK',len(raw))
    print('RELATIONS',dual.tolist())
    for w in range(18,24):
        idx=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=w]
        lower=dual[:,idx]
        before=lower.rank();after=lower.row_join(dual[:,-1]).rank()
        check((before,after)=={18:(1,2),19:(1,2),20:(3,4),21:(3,4),22:(5,5),23:(5,5)}[w],
              'exact sharp weight threshold')
        print('WEIGHT',w,'rank',before,'with_high',after)
        if before==after:
            ans=S.linsolve((lower,dual[:,-1]))
            print('REPLACEMENTS r_high=-1',[(monomials[i]) for i in idx],ans)
    # Complete affine basis at weight22, then substitute original equations.
    allowed=[i for i,(a,b) in enumerate(low) if 2*a+3*b<=22]
    coeff=dual[:,allowed]
    space=next(iter(S.linsolve((coeff,dual[:,-1]))))
    free=sorted(set().union(*(v.free_symbols for v in space)),key=str)
    check(len(free)==5,'five-dimensional complete source fibre')
    origin=[v.subs({p:0 for p in free}) for v in space]
    columns=[origin]+[[S.diff(v,p) for v in space] for p in free]
    for column_number,column in enumerate(columns):
        values=[S.Integer(0)]*len(rs)
        for i,value in zip(allowed,column):values[i]=value
        values[-1]=-1 if column_number==0 else 0
        check(dual*S.Matrix(values)==S.zeros(5,1),'affine source basis in full kernel')
        solution=[S.Integer(0)]*nu
        for row,pivot in enumerate(piv):
            if pivot<nu:
                solution[pivot]=-sum(red[row,nu+j]*values[j] for j in range(len(rs)))
        residual=matrix*S.Matrix(solution+values)
        for value in residual:check(value==0,'literal raw residual on affine basis')
    # The recovered old transport remains a section of the larger fibre.
    old=S.zeros(15,1)
    old[7]=-S.Rational(27945,235202)
    old[9]=S.Rational(39123,470404)
    old[11]=S.Rational(52578,117601)
    old[14]=-1
    check(dual*old==S.zeros(5,1),'earlier audited transport embeds')
    bad=old.copy();bad[5]=1
    check(dual*bad!=S.zeros(5,1),'valuation12 y6 coefficient is forced zero')
    bad=old.copy();bad[1]+=1
    check(dual*bad!=S.zeros(5,1),'p8y2 off-graph hostile')
    print('CONTROLS six affine columns vanish in all271 original equations')
    print('HOSTILES weight21_fails y6_forced_zero p8y2_off_graph_fails')
    payload={'monomials':monomials,'relations':[[str(v) for v in row] for row in dual.tolist()]}
    print('SEMANTIC_SHA256',sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest())
    print('CHECKS',CHECKS,'PASS')
    print('SCOPE exact_row15_valuation_at_least12_firstfive_backgroundrows_fixed_JC2_OPEN')

if __name__=='__main__':main()
