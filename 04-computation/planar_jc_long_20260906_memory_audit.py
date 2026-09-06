#!/usr/bin/env python3
"""Independent literal P/J and full diagonal-depth audit of the response.

No repository mathematical imports. Starts with all 128 row coefficients,
not the producer's 58 tangent variables or its row-solver routines.
"""
import sympy as S
from sympy.polys.matrices import DomainMatrix
from math import comb
from hashlib import sha256
import json

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def main():
    x = S.Symbol('x')
    A = [1+x*x/4, S.Rational(4,3)+2*x*x,
         -S.Rational(32,9)-S.Rational(4,5)*x*x,
         S.Rational(2176,135)+S.Rational(64,9)*x*x-S.Rational(32,15)*x**4]
    C = [-3*x/4-x**3/8, -4*x-S.Rational(3,2)*x**3,
         S.Rational(88,15)*x-S.Rational(12,5)*x**3,
         -S.Rational(1408,45)*x+S.Rational(64,15)*x**3+S.Rational(8,5)*x**5]
    da, dc, unknown = {}, {}, []
    for n in range(12,16):
        va = S.symbols('a%d_0:%d' % (n,n+2))
        vc = S.symbols('c%d_0:%d' % (n,n+3))
        da[n] = sum(v*x**j for j,v in enumerate(va))
        dc[n] = sum(v*x**j for j,v in enumerate(vc))
        unknown.extend(va+vc)
    check(len(unknown)==128, 'full coefficient universe')
    monomials = [(7,3),(5,4),(3,5),(1,6),(4,5),(2,6),(0,7),(1,7),(3,6)]
    rs = S.symbols('r0:9')
    source = {n:S.Integer(0) for n in range(12,16)}
    for r,(a,b) in zip(rs,monomials):
        # p^a y^b = x^b t^(a+2b) (1+x^2t)^(a+b)
        for j in range(a+b+1):
            n = a+2*b+j
            if n in source:
                source[n] += r*comb(a+b,j)*x**(b+2*j)
    ga = [S.expand(-3*sum(A[i]*A[k-i] for i in range(k+1))
                   +(S.Rational(3,4) if k==0 else 0)) for k in range(4)]
    eq = []
    p_count=j_count=0
    for n in range(12,16):
        pe = -source[n]
        je = 0
        for k in range(4):
            j = n-k
            if j not in da:
                continue
            pe += ga[k]*da[j]+2*C[k]*dc[j]
            je += (j*S.diff(A[k],x)*dc[j] + k*S.diff(da[j],x)*C[k]
                   -k*A[k]*S.diff(dc[j],x)-j*da[j]*S.diff(C[k],x))
        pv=S.Poly(S.expand(pe),x).all_coeffs()
        jv=S.Poly(S.expand(je),x).all_coeffs()
        eq.extend(pv+jv)
        p_count+=len(pv); j_count+=len(jv)
    depth_count=0
    for d,v in [(2,da),(3,dc)]:
        # Complete source-depth annihilator, derived by polynomial division.
        for ell in range(1,31):
            start=(ell+1)//2
            stop=min(15,ell+d)
            rho=(ell+2)//3
            if start>stop:
                continue
            for q in range(max(0,ell+d-stop),rho):
                e=S.Integer(0)
                for n in range(start,stop+1):
                    if n in v:
                        coeff=S.expand(v[n]).coeff(x,2*n-ell)
                        e+=(-1)**(n-start)*comb(stop+q-n,q)*coeff
                eq.append(e)
                depth_count+=1
    check(depth_count==91, 'complete depth codimension')
    matrix, rhs=S.linear_eq_to_matrix(eq,unknown+list(rs))
    check(rhs==S.zeros(len(eq),1),'homogeneous equations')
    dm,pivots=DomainMatrix.from_Matrix(matrix).convert_to(S.QQ).rref()
    reduced=dm.to_Matrix()
    raw=[]
    for i in range(reduced.rows):
        if all(reduced[i,j]==0 for j in range(128)):
            row=list(reduced[i,128:])
            if any(row): raw.append(row)
    dual=S.Matrix(raw).rref()[0]
    expected=S.zeros(3,9)
    expected[0,1]=1;expected[0,8]=-S.Rational(27945,235202)
    expected[1,3]=1;expected[1,8]=S.Rational(39123,470404)
    expected[2,5]=1;expected[2,8]=S.Rational(52578,117601)
    check(dual==expected,'independently reconstructed full response')
    tangent_rank=sum(p<128 for p in pivots)
    check(tangent_rank==118 and len(pivots)==121,'full linear ranks')
    print('INDEPENDENT MEMORY AUDIT; no repository mathematical imports')
    print('RAW rows P',p_count,'J',j_count,'depth',depth_count,'unknown128 source9')
    print('RANK tangent118 source3 joint121 terminal10')
    print('FULL_SOURCE_RELATIONS',dual.tolist())
    for w in (20,21,22,23):
        idx=[i for i,(a,b) in enumerate(monomials[:-1]) if 2*a+3*b<=w]
        before=dual[:,idx].rank()
        after=dual[:,idx].row_join(dual[:,-1]).rank()
        check((before,after)=={20:(1,2),21:(1,2),22:(3,3),23:(3,3)}[w],
              'independent minimality test')
        print('WEIGHT',w,before,after)
    # Concrete coefficient lift obtained by a fresh RREF of the full system.
    params={r:S.Integer(0) for r in rs}
    params[rs[1]]=-S.Rational(27945,235202)
    params[rs[3]]=S.Rational(39123,470404)
    params[rs[5]]=S.Rational(52578,117601)
    params[rs[8]]=-1
    solution={v:S.Integer(0) for v in unknown}
    for row,pivot in enumerate(pivots):
        if pivot<128:
            solution[unknown[pivot]]=-sum(reduced[row,128+j]*params[rs[j]] for j in range(9))
    original_residual=matrix*S.Matrix([solution[v] for v in unknown]+[params[r] for r in rs])
    for value in original_residual:
        check(value==0,'every original literal equation after independent lift')
    record={'dual':[[str(v) for v in row] for row in dual.tolist()],
            'rows':len(eq),'ranks':[tangent_rank,len(pivots)]}
    print('SEMANTIC_SHA256',sha256(json.dumps(record,sort_keys=True).encode()).hexdigest())
    print('CHECKS',CHECKS,'PASS')
    print('SCOPE valuation>=13 source weights<=23 plus designated high term; rows<=15 only')


if __name__=='__main__':
    main()
