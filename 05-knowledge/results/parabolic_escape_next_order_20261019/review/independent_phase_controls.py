#!/usr/bin/env python3
"""Independent audit controls. Symbolic identities and one interval certificate
are distinguished from finite high-precision numerical checks.
Dependencies: Python 3, sympy, mpmath. No network or source-packet imports.
"""
from fractions import Fraction
import sympy as sp
import mpmath as mp

mp.mp.dps = 80
mp.iv.dps = 80


def require(ok, why):
    if not ok:
        raise AssertionError(why)


def symbolic_checks():
    y, e, M, C, z = sp.symbols('y e M C z', real=True)
    w = 1 + y
    # Derivative of the improved exact clock defect, with A_y=-1/(w^2+e).
    E_y = -(1+e/3)/(w*w+e)+w/(w*w+e)-y
    stated = (-2*y*y-y**3-e*y-e/3)/(w*w+e)
    require(sp.cancel(E_y-stated) == 0, 'Improved derivative mismatch')
    E0 = 1/w + sp.log(w) - 1 - y*y/2
    require(sp.series(E0,y,0,3).removeO() == 0, 'Cubic vanishing failed')
    E_e_at_zero = 1/(3*w)-1/(3*w**3)+1/(2*w*w)-sp.Rational(1,2)
    require(E_e_at_zero.subs(y,0) == 0, 'Parameter cancellation failed')
    t = M-sp.pi**2/(3*M)
    require(sp.series((sp.pi**2/t**2).subs(M,1/z),z,0,5).removeO()
            == sp.pi**2*z*z+2*sp.pi**4*z**4/3, 'Inversion coefficient failed')
    width = sp.pi**2/(1/z-C)**2-sp.pi**2/(1/z+1-C)**2
    require(sp.expand(sp.series(width,z,0,5).removeO()
            - (2*sp.pi**2*z**3+(6*C-3)*sp.pi**2*z**4)) == 0,
            'Adjacent coefficient failed')
    # Independently construct a high-order unperturbed normalization.
    D0 = 1/w+sp.log(w)-1
    P = sp.S.Zero
    for k in range(1,9):
        c = sp.symbols('c')
        residual = D0+P.subs(y,y+y*y)-P+c*((y+y*y)**k-y**k)
        equation = sp.series(residual,y,0,k+2).removeO().expand().coeff(y,k+1)
        P += sp.solve(equation,c)[0]*y**k
    Rprime = sp.cancel(y/w**2+sp.diff(P,y).subs(y,y+y*y)*(1+2*y)-sp.diff(P,y))
    Q = sp.cancel(Rprime*w**2/y**9)
    require(sp.denom(Q) == 1 or sp.denom(Q).is_Integer, 'Nonpolynomial Q')
    require(sp.Poly(Q,y).degree() == 8, 'Residual quotient degree mismatch')
    require(sp.cancel(Rprime-y**9*Q/w**2) == 0, 'Residual derivative failed')
    print('EXACT SYMBOLIC: defect derivative, cubic/parameter cancellations, inversion, adjacent width passed')
    print('EXACT SYMBOLIC: P8 =',P)
    print('EXACT SYMBOLIC: Rprime =',sp.factor(Rprime))
    return [P.coeff(y,k) for k in range(1,9)], [sp.expand(Q).coeff(y,k) for k in range(9)]


def ivrat(r):
    num, den = r.as_numer_denom()
    return mp.iv.mpf(int(num))/int(den)


def phase_enclosure(b_string, pcoeff, qcoeff, m=512):
    iv=mp.iv
    u=iv.mpf('-0.5');v=iv.mpf(b_string)
    for _ in range(m):
        u=u+u*u
        v=2*v/(1+iv.sqrt(1+4*v))
    eta=iv.mpf('0.01')
    require(u.a > -eta.a and u.b < 0 and v.a > 0 and v.b < eta.a,
            'Endpoint cutoff too large for certified residual bound')
    # R'(y)=y^9 Q(y)/(1+y)^2, hence |R(y)|<=K |y|^10/10.
    # Sum |y|^10 <= endpoint^8 sum y^2 = endpoint^9 on either tail.
    K=sum(abs(ivrat(c))*eta**k for k,c in enumerate(qcoeff))/(1-eta)**2
    tail=K*((-u)**9+v**9)/10
    def H(x,positive):
        return -1/x+iv.log(x if positive else -x)+sum(ivrat(c)*x**k for k,c in enumerate(pcoeff,1))
    estimate=H(v,True)-H(u,False)+2*m
    enclosure=estimate+iv.mpf([-1,1])*tail
    print('DIRECTED INTERVAL: b =',b_string,'m =',m,'C_R in',enclosure)
    print('  Absolute analytic tail bound <=',tail.b)
    lo=mp.mpf(enclosure.a._mpi_[0]);hi=mp.mpf(enclosure.b._mpi_[1])
    require(lo < hi, 'Empty phase interval')
    if b_string=='1.5':
        require(enclosure.a > iv.mpf('-2.463836664').b
                and enclosure.b < iv.mpf('-2.463836643').a,
                'Disagrees with initial phase certificate')
        if m == 512:
            require(enclosure.a > iv.mpf('-2.463836653392077114885602').b
                    and enclosure.b < iv.mpf('-2.463836653392077114885581').a,
                    'Rounded published C2 interval failed')
    return (lo+hi)/2,(hi-lo)/2,enclosure


def shoot(e,n,b):
    y=mp.mpf('-0.5');d=mp.mpf('0')
    for _ in range(n):
        d=(1+2*y)*d+1
        y=y+y*y+e
        if y>b+1:
            return None,None
    return y-b,d


def threshold_newton(n,b,c):
    lo=mp.mpf(0);hi=mp.mpf(1)
    # hi=1 escapes all selected b and n; certify this rather than infer it.
    f,_=shoot(hi,n,b)
    require(f is None or f>0, 'Upper root bracket failed')
    e=mp.pi**2/(n-c)**2
    for it in range(200):
        f,d=shoot(e,n,b)
        if f is None:
            hi=e;candidate=(lo+hi)/2
        else:
            if abs(f)<mp.mpf('1e-65'):
                return e,it+1
            if f>0:hi=e
            else:lo=e
            candidate=e-f/d
            if not lo<candidate<hi:
                candidate=(lo+hi)/2
        e=candidate
    raise AssertionError('Safeguarded Newton did not converge')


def independent_shooting_controls(phases):
    print('NUMERICAL: independent safeguarded Newton; residual is n-pi/delta-C_R-pi*delta/3')
    print('b | n | delta | residual/delta^2 | residual/(delta^2*log(1/delta)) | iterations')
    for bstr,(c,cerr,_) in phases.items():
        b=mp.mpf(bstr)
        roots={}
        for n in [30,100,300,301,1000]:
            e,it=threshold_newton(n,b,c);roots[n]=e;delta=mp.sqrt(e)
            residual=n-mp.pi/delta-c-mp.pi*delta/3
            require(cerr/e < mp.mpf('1e-12'), 'Phase certificate too coarse for this sample')
            print(bstr,n,mp.nstr(delta,13),mp.nstr(residual/e,13),
                  mp.nstr(residual/(e*mp.log(1/delta)),13),it,sep=' | ')
            # Direct check of the improved identity at this computed orbit.
            def E(y):
                w=1+y
                return (1+e/3)*mp.atan(delta/w)/delta+mp.log(w*w+e)/2-1-(y*y+e)/2
            def clock(y):
                return (1+e/3)*mp.atan(y/delta)/delta+mp.log(y*y+e)/2-y/2
            y=mp.mpf('-.5');total=mp.mpf(0)
            for k in range(n):
                total+=E(y)
                y+=y*y+e
            require(abs(clock(y)-clock(mp.mpf('-.5'))-n-total)<mp.mpf('1e-65'),
                    'Summed improved clock identity failed')
        n=300
        actual=roots[n]-roots[n+1]
        prediction=2*mp.pi**2/n**3+(6*c-3)*mp.pi**2/n**4
        print('NUMERICAL: width b =',bstr,'n = 300; remainder*n^5/log(n) =',
              mp.nstr((actual-prediction)*n**5/mp.log(n),14))


def strict_boundary_control():
    # Rational threshold for R=3/4: e_2=1/4, and e_1=1/2.
    R=Fraction(3,4)
    def escape(e):
        x=Fraction(0)
        for n in range(10):
            if x>R:return n
            x=x*x+Fraction(1,4)+e
        raise AssertionError('No rational escape')
    require(escape(Fraction(1,4))==3, 'Equality incorrectly counted as escape')
    require(escape(Fraction(1,2))==2, 'Upper band equality convention failed')
    require(escape(Fraction(3,10))==2, 'Interior band convention failed')
    print('EXACT RATIONAL: strict escape and both threshold-band endpoint conventions passed')


if __name__=='__main__':
    p,q=symbolic_checks()
    phases={b:phase_enclosure(b,p,q) for b in ['0.05','0.25','1.5','10']}
    c2, r2, iv2 = phases['1.5']
    print('CERTIFICATE: -2.463836653392077114885602 < C_2 < -2.463836653392077114885581')
    c2_short, r2_short, iv2_short = phase_enclosure('1.5',p,q,m=256)
    require(iv2.a < iv2_short.b and iv2_short.a < iv2.b, 'Independent cutoffs disagree')
    c_shift, r_shift, iv_shift = phase_enclosure('3.75',p,q,m=512)
    shift_interval = iv_shift-iv2
    require(shift_interval.a <= 1 <= shift_interval.b, 'Exact Abel shift disagrees with phase certificates')
    print('DIRECTED INTERVAL: cutoff overlap and C_(b+b^2)-C_b = 1 consistency passed')
    strict_boundary_control()
    independent_shooting_controls(phases)
    print('PASS: exact symbolic checks, directed phase certificates, rational conventions, and numerical controls')
