#!/usr/bin/env python3
"""Exact rational certificates plus explicitly non-rigorous high-precision diagnostics.
No network, repository access, Git mutation, or external data is needed.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
import mpmath as mp

HERE = Path(__file__).resolve().parent
mp.mp.dps = 100

def a(m):
    return F((m+1)**m, factorial(m+1))

def dec(q, digits=80, upper=False):
    scale = 10**digits
    n = q.numerator*scale // q.denominator
    if upper and F(n,scale) != q:
        n += 1
    sign = '-' if n < 0 else ''
    n = abs(n)
    return f'{sign}{n//scale}.{n%scale:0{digits}d}'

def mpq(q):
    return mp.mpf(q.numerator)/q.denominator

def determinant(mat):
    mat = [list(row) for row in mat]
    d = F(1)
    n = len(mat)
    for k in range(n):
        row = next((i for i in range(k,n) if mat[i][k]),None)
        if row is None:
            return F(0)
        if row != k:
            mat[k],mat[row] = mat[row],mat[k]
            d = -d
        pivot = mat[k][k]
        d *= pivot
        for i in range(k+1,n):
            c = mat[i][k]/pivot
            for j in range(k+1,n):
                mat[i][j] -= c*mat[k][j]
    return d

def solve(mat,rhs):
    aug = [list(row)+[rhs[i]] for i,row in enumerate(mat)]
    n = len(aug)
    for k in range(n):
        row = next(i for i in range(k,n) if aug[i][k])
        aug[k],aug[row] = aug[row],aug[k]
        p = aug[k][k]
        aug[k] = [v/p for v in aug[k]]
        for i in range(n):
            if i != k:
                c = aug[i][k]
                aug[i] = [u-c*v for u,v in zip(aug[i],aug[k])]
    return [row[-1] for row in aug]

def pade_gauss(mom,n):
    """[n-1/n] approximant with exact match through degree 2n-1."""
    q = [F(1)] + solve([[mom[k-j] for j in range(1,n+1)]
                       for k in range(n,2*n)],[-mom[k] for k in range(n,2*n)])
    p = [sum((q[j]*mom[k-j] for j in range(k+1)),F(0)) for k in range(n)]
    for k in range(2*n):
        lhs = sum((q[j]*mom[k-j] for j in range(min(k,n)+1)),F(0))
        assert lhs == (p[k] if k < n else 0)
    return p,q

def poly(co,z):
    out = F(0)
    for c in reversed(co):
        out = out*z+c
    return out

def val(pq,z):
    p,q = pq
    assert poly(q,z)>0
    return poly(p,z)/poly(q,z)

def main():
    cert = {'status':'PASS', 'proof_dependency':'TREE_MOMENT_BRIDGE.md',
            'certification_scope':'Rational calculations are exact. Analytic enclosure directions use the proved moment and Gauss theorems. mpmath tests are diagnostics, not interval certificates.'}
    cert['initial_coefficients'] = [str(a(m)) for m in range(12)]
    hankel = []
    for shift in range(5):
        for n in range(1,9):
            d = determinant([[a(shift+i+j) for j in range(n)] for i in range(n)])
            assert d>0
            hankel.append({'shift':shift,'size':n,'determinant':str(d)})
    cert['exact_hankel_tests'] = hankel
    for m in range(101):
        assert a(m+1)/a(m) == F(m+2,m+1)**m
        if m:
            assert a(m-1)*a(m+1)>a(m)**2
    cert['exact_ratio_and_log_convexity_range'] = 'm=0..100 ratios; m=1..100 strict log-convexity'

    # e < sum_{k=0}^{K} 1/k! + 1/(K K!), by a geometric tail bound.
    K=10
    e_upper = sum((F(1,factorial(k)) for k in range(K+1)),F(0))+F(1,K*factorial(K))
    b=F(2719,1000)
    assert e_upper<b<3
    cert['rational_support_enclosure']={'e_upper':str(e_upper),'b':str(b),'proof':'Geometric bound on exp(1) Taylor tail at K=10; e<e_upper<b.'}
    rows=[]
    z=F(1,3)
    actual=-mp.lambertw(-mpq(z))/mpq(z)
    for n in [1,2,4,8,12,16,24,32]:
        moments=[a(j) for j in range(2*n+1)]
        transformed=[b*moments[j]-moments[j+1] for j in range(2*n)]
        gauss=pade_gauss(moments,n)
        gauss_b=pade_gauss(transformed,n)
        lower=val(gauss,z)
        upper=(1-z*val(gauss_b,z))/(1-b*z)
        assert lower<upper
        # Diagnostic only: enclosure's proof does not depend on Lambert W numerics.
        assert mpq(lower)<actual<mpq(upper)
        raw=sum((a(j)*z**j for j in range(2*n)),F(0))
        raw_width=a(2*n)*z**(2*n)/(1-e_upper*z)
        assert raw<=lower and upper<raw+raw_width
        rows.append({'N':n,'z':str(z),'lower_exact':str(lower),'upper_exact':str(upper),
                     'lower_decimal_outward':dec(lower),'upper_decimal_outward':dec(upper,upper=True),
                     'width_exact':str(upper-lower),'width_approx':mp.nstr(mpq(upper-lower),15),
                     'raw_series_degree':2*n-1,'raw_series_certified_width_approx':mp.nstr(mpq(raw_width),15),
                     'width_improvement_factor_approx':mp.nstr(mpq(raw_width/(upper-lower)),15),
                     'P':[str(v) for v in gauss[0]],'Q':[str(v) for v in gauss[1]],
                     'P_b':[str(v) for v in gauss_b[0]],'Q_b':[str(v) for v in gauss_b[1]]})
    cert['gauss_rational_enclosures_R_one_third']=rows

    # Positive-density integrals; the quadrature is not claimed to be rigorous.
    def x(v):
        if v == 0:
            return mp.e
        if v == mp.pi:
            return mp.mpf(0)
        return mp.exp(v*mp.cot(v))*mp.sin(v)/v

    def D(v):
        if v == 0:
            return mp.mpf(0)
        return v*v+(1-v*mp.cot(v))**2

    mesh=[0,mp.pi/4,mp.pi/2,3*mp.pi/4,mp.pi]

    def moment_formula(s):
        p=s+1
        return mp.exp((p-2)*mp.log(p))/mp.gamma(p)

    def moment_integral(s):
        p=s+1
        return mp.quad(lambda v:x(v)**p,mesh)/(mp.pi*p)

    checks=[]
    for s in [mp.mpf('-0.875'),mp.mpf('-.5'),0,1,2,7,mp.mpf('.2')+mp.mpf('.3')*1j]:
        q=moment_integral(s)
        target=moment_formula(s)
        err=abs(q-target)
        assert err<mp.mpf('1e-85')
        checks.append({'s':str(s),'value':mp.nstr(target,65),'abs_error':mp.nstr(err,6)})
    cert['fractional_moment_diagnostics']=checks
    checks=[]
    for zc in [mp.mpf(0),mp.mpf('-2'),mp.mpf('0.1'),mp.mpf(1)/3,
               mp.mpf('.2')+mp.mpf('.1')*1j,mp.mpf('1')+mp.mpf('.25')*1j]:
        q=mp.quad(lambda v:x(v)*D(v)/(1-zc*x(v)),mesh)/mp.pi
        target=-mp.lambertw(-zc)/zc if zc else mp.mpf(1)
        err=abs(q-target)
        assert err<mp.mpf('1e-80')
        checks.append({'z':str(zc),'value':mp.nstr(target,65),'abs_error':mp.nstr(err,6)})
    cert['stieltjes_integral_diagnostics']=checks
    checks=[]
    for v in [mp.mpf('.3'),mp.mpf('1'),mp.mpf('2.5')]:
        xx=x(v)
        eps=xx*mp.mpf('1e-40')
        target=v*mp.cot(v)-1j*v
        f=-mp.lambertw(-1/(xx+1j*eps))
        err=abs(f-target)
        assert f.imag<0 and err<mp.mpf('1e-37')
        checks.append({'v':str(v),'x':mp.nstr(xx,40),'upper_boundary_abs_error':mp.nstr(err,6)})
    cert['boundary_sign_diagnostics']=checks
    checks=[]
    for w in [mp.mpf(1),mp.mpf(2),mp.mpf(3),mp.mpf(8),mp.mpf(20),1+1j,2+3j]:
        p=mp.mpf(1)/w
        bridge=w*w/(mp.gamma(p)*moment_integral(p-1))
        target=mp.exp(mp.log(w)/w)
        err=abs(bridge-target)
        assert err<mp.mpf('1e-75')
        checks.append({'w':str(w),'target':mp.nstr(target,55),'integral_bridge_abs_error':mp.nstr(err,6)})
    cert['full_MRB_kernel_moment_bridge_diagnostics']=checks
    checks=[]
    for p in [mp.mpf(1)/8,mp.mpf(1)/4,mp.mpf('.7'),mp.mpf('.2')+mp.mpf('.1')*1j]:
        A=2/mp.pi*mp.exp(p*mp.log(p))*mp.sin(mp.pi*p)
        via_gamma=2*mp.exp(p*mp.log(p))/(mp.gamma(p)*mp.gamma(1-p))
        via_moment=2*p*p*moment_integral(p-1)/mp.gamma(1-p)
        err=max(abs(A-via_gamma),abs(A-via_moment))
        assert err<mp.mpf('1e-85')
        checks.append({'p':str(p),'A_one_over_p':mp.nstr(A,55),'abs_error':mp.nstr(err,6)})
    cert['branch_kernel_gamma_moment_diagnostics']=checks

    (HERE/'verification.json').write_text(json.dumps(cert,indent=2)+'\n')
    print('PASS: 40 exact Hankel determinants; 101 exact ratio identities; 100 strict log-convexity tests.')
    print('PASS: exact rational Gauss and modified-Gauss bounds at z=1/3.')
    for row in rows:
        print(f"  N={row['N']:2d}: width={row['width_approx']}; raw degree {row['raw_series_degree']:2d} width={row['raw_series_certified_width_approx']}; factor={row['width_improvement_factor_approx']}")
    print('PASS: all high-precision diagnostics (separate from exact certification).')
    print('R(1/3), diagnostic:',mp.nstr(actual,82))

if __name__ == '__main__':
    main()
