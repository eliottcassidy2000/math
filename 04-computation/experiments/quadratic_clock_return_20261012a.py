"""Exact quadratic cotangent-return certificates; no planar JC conclusion.

The target coefficient ring is fixed. Source field embeddings vary with depth.
Public compilers take exact rational inputs and retain parameter, depth and target.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as Q
import sympy as s

x, p, X, Y = s.symbols('x p X Y')
CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def rational(value):
    need(type(value) in (int, Q), 'exact rational input')
    return Q(value)


def positive_depth(n):
    need(type(n) is int and n >= 1, 'positive exact depth')
    return n


def symq(value):
    value = rational(value)
    return s.Rational(value.numerator, value.denominator)


def twoadic(value):
    value = rational(value)
    need(value != 0, 'nonzero valuation input')
    a, b = abs(value.numerator), value.denominator
    return (a & -a).bit_length() - (b & -b).bit_length()


def orbit(c, n):
    c, n = rational(c), positive_depth(n)
    out, z = [], Q(0)
    for _ in range(n):
        z = z*z + c
        out.append(z)
    return tuple(out)


def iterate(c, n):
    c, n = symq(c), positive_depth(n)
    f = s.Poly(x, x, domain=s.QQ)
    for _ in range(n):
        f = f*f + s.Poly(c, x, domain=s.QQ)
    return f


def postcritical(c, n):
    """The distinct-value theorem is guarded by v_2(c)<0."""
    c, n = rational(c), positive_depth(n)
    need(c != 0 and twoadic(c) < 0, 'negative two-adic parameter valuation')
    values = orbit(c, n)
    delta = s.Poly(1, X, domain=s.QQ)
    for v in values:
        delta *= s.Poly(X-symq(v), X, domain=s.QQ)
    return values, delta


def exact_poly(poly, variables):
    need(type(poly) is s.Poly and poly.gens == variables,
         'polynomial in the declared ordered variables')
    need(poly.domain in (s.QQ, s.ZZ), 'exact rational polynomial coefficients')
    return poly.set_domain(s.QQ)


def coefficients(target):
    target = exact_poly(target, (X, Y))
    out = [s.Poly(0, X, domain=s.QQ) for _ in range(max(0, target.degree(Y))+1)]
    for (i, j), c in target.terms():
        out[j] += s.Poly(c*X**i, X, domain=s.QQ)
    return tuple(out)


def admissible(c, n, target):
    _, delta = postcritical(c, n)
    cs = coefficients(target)
    return all(a.rem(delta**((j+1)//2)).is_zero for j, a in enumerate(cs))


@dataclass(frozen=True)
class ReturnReceipt:
    parameter: Q
    depth: int
    target: s.Poly
    extension: s.Poly


def compile_return(c, n, target):
    """Return the unique polynomial extension of the guarded rational pullback."""
    c, n = rational(c), positive_depth(n)
    target = exact_poly(target, (X, Y))
    need(admissible(c, n, target), 'unpaid critical-value divisor')
    f = iterate(c, n)
    derivative = f.diff()
    result = s.Poly(0, x, p, domain=s.QQ)
    for j, a in enumerate(coefficients(target)):
        numerator = s.Poly.from_dict(a.as_dict(), x, domain=s.QQ).compose(f)
        quotient, remainder = numerator.div(derivative**j)
        need(remainder.is_zero, 'independent pullback division failed')
        result += s.Poly(quotient.as_expr()*p**j, x, p, domain=s.QQ)
    return ReturnReceipt(c, n, target, result)


def audit_receipt(receipt):
    need(type(receipt) is ReturnReceipt, 'typed return receipt')
    rational(receipt.parameter)
    positive_depth(receipt.depth)
    exact_poly(receipt.target, (X, Y))
    exact_poly(receipt.extension, (x, p))
    want = compile_return(receipt.parameter, receipt.depth, receipt.target)
    need(receipt == want, 'source-bound extension mismatch')
    return receipt


def evaluate_lift(c, n, source_x, source_p):
    source_x, source_p = symq(source_x), symq(source_p)
    f = iterate(c, n)
    d = f.diff().eval(source_x)
    need(d != 0, 'removed ramification divisor')
    return Q(f.eval(source_x)), Q(source_p/d)


def direct_pullback(c, n, target):
    """Independent rational substitution; no use of divisibility criterion."""
    target = exact_poly(target, (X, Y))
    f = iterate(c, n)
    df = f.diff()
    result = 0
    for j,a in enumerate(coefficients(target)):
        if a.is_zero:
            continue
        numerator = s.Poly.from_dict(a.as_dict(),x,domain=s.QQ).compose(f)
        denominator = df**j
        common = s.gcd(numerator,denominator)
        result += numerator.exquo(common).as_expr()*p**j/denominator.exquo(common).as_expr()
    return s.cancel(result)


def rejects(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def main():
    params = (Q(-7,4), Q(-29,16), Q(3,2))
    compiled = 0
    for c in params:
        for n in range(1, 5):
            values, delta = postcritical(c, n)
            f = iterate(c, n)
            df = f.diff()
            check(len(set(values)) == n, 'distinct critical values')
            check(all(twoadic(v) == 2**(i-1)*twoadic(c)
                      for i,v in enumerate(values,1)), 'exact critical valuation growth')
            check(s.gcd(df, df.diff()).degree() == 0, 'simple critical points')
            check(df.degree() == 2**n-1, 'exponential critical-point count')
            product = s.Poly(2**n, x, domain=s.QQ)
            for j in range(n):
                product *= s.Poly(x,x,domain=s.QQ) if j == 0 else iterate(c,j)
            check(product == df, 'chain-rule factorization')
            for i,v in enumerate(values,1):
                fibre = f - s.Poly(symq(v),x,domain=s.QQ)
                ram = s.gcd(fibre,df)
                check(ram.degree() == 2**(n-i), 'critical-value multiplicity')
                check(fibre.rem(ram**2).is_zero, 'double ramification')
                companion = fibre.exquo(ram**2)
                check(companion.degree() == 2**n-2**(n-i+1), 'regular companion count')
                check(s.gcd(companion,df).degree() == 0, 'companions unramified')
            d_source = s.Poly.from_dict(delta.as_dict(),x,domain=s.QQ).compose(f)
            h, rem = d_source.div(df**2)
            check(rem.is_zero and s.gcd(h,df).degree() == 0, 'repeated-value square divisor')
            for j in range(5 if n <= 2 else 3):
                target = s.Poly(delta.as_expr()**((j+1)//2)*(X+1)*Y**j, X,Y,domain=s.QQ)
                packet = compile_return(c,n,target)
                check(audit_receipt(packet) == packet, 'authenticated return')
                check(s.cancel(packet.extension.as_expr()-direct_pullback(c,n,target)) == 0,
                      'independent pullback')
                compiled += 1
                if j:
                    unpaid = s.Poly(delta.as_expr()**(((j+1)//2)-1)*Y**j,X,Y,domain=s.QQ)
                    check(not admissible(c,n,unpaid), 'one missing full divisor rejected')
                    check(s.denom(direct_pullback(c,n,unpaid)).has(x), 'actual pole witness')
            # The two generators have the claimed exact nodal relation.
            aa = compile_return(c,n,s.Poly(delta.as_expr()*Y,X,Y,domain=s.QQ)).extension.as_expr()
            bb = compile_return(c,n,s.Poly(delta.as_expr()*Y**2,X,Y,domain=s.QQ)).extension.as_expr()
            check(s.expand(aa*aa-d_source.as_expr()*bb) == 0, 'nodal relation')
            check(s.gcd(delta,delta.diff()).degree() == 0, 'exactly n distinct nodes')
            for gen in (f.as_expr(),aa,bb):
                for variable in (x,p):
                    restricted = s.Poly(s.diff(gen,variable).subs(p,0),x,domain=s.QQ)
                    check(restricted.rem(df).is_zero, 'all observable generators have zero differential at critical zero section')
            h1 = f.as_expr()+aa**2+bb
            h2 = aa+f.as_expr()*bb
            jac0 = s.diff(h1,x).subs(p,0)*s.diff(h2,p).subs(p,0)-s.diff(h1,p).subs(p,0)*s.diff(h2,x).subs(p,0)
            check(s.Poly(jac0,x,domain=s.QQ).rem(df).is_zero, 'nonlinear reencoding retains zero Jacobian')
            if n <= 3:
                target = s.Poly(delta.as_expr()*Y,X,Y,domain=s.QQ)
                for t in range(1,n+1):
                    check(admissible(c,t,target), 'all preceding depths survive')
                check(not admissible(c,n+1,target), 'sharp next-depth failure')
            # Positive evaluation controls retain the actual open-source guard.
            for z in (-3,1,3):
                endpoint = evaluate_lift(c,n,z,2)
                for j in (1,2):
                    target = s.Poly(delta.as_expr()*Y**j,X,Y,domain=s.QQ)
                    packet = compile_return(c,n,target)
                    lhs = packet.extension.eval({x:z,p:2})
                    rhs = target.eval({X:symq(endpoint[0]),Y:symq(endpoint[1])})
                    check(lhs == rhs, 'source-endpoint evaluation')
            check(rejects(lambda:evaluate_lift(c,n,0,1)), 'critical source excluded')

    # Exhaustive monomial coefficient tests, independent rational denominators.
    monomial_controls = 0
    for c in params[:2]:
        for n in range(1,4):
            values, delta = postcritical(c,n)
            for j in range(1,5):
                for power in range(4):
                    for extra in (s.Integer(1),X,X**2-1):
                        target = s.Poly(delta.as_expr()**power*extra*Y**j,X,Y,domain=s.QQ)
                        denom = s.denom(direct_pullback(c,n,target))
                        check(admissible(c,n,target) == (not denom.has(x)), 'complete bounded coefficient test')
                        monomial_controls += 1

    # Higher ramification and finite postcritical orbits are distinct hostiles.
    z0 = direct_pullback(0,2,s.Poly(X**2*Y**3,X,Y,domain=s.QQ))
    check(z0 == p**3/(64*x), 'c0 ramification-index-four hostile')
    for e in range(2,7):
        for j in range(1,7):
            power = (j*(e-1)+e-1)//e
            remaining = e*power-j*(e-1)
            check(remaining >= 0, 'general local index pays polynomiality')
            check(remaining-e < 0, 'general local index sharpness')
            observable = x**remaining*p**j/s.Integer(e)**j
            check(all(s.diff(observable,z).subs({x:0,p:0}) == 0 for z in (x,p)),
                  'higher ramification all-observable zero differential')
    check(s.diff(2*x,x)*s.diff(p/2,p) == 1, 'linear boundary permits polynomial symplectic map')
    for n in range(1,6):
        tame = direct_pullback(-2,n,s.Poly((X**2-4)*Y**2,X,Y,domain=s.QQ))
        check(not s.denom(tame).has(x), 'finite critical orbit survives all tested depths')
    check(orbit(Q(-7,4),3) == (Q(-7,4),Q(21,16),Q(-7,256)), '63/64 critical-value identity')
    check((-7)*((-7)+4)**2+4**3 == 1, '64-minus-63 norm bracket')

    # Parabolic algebraic cycle and nonparabolic rational cycle.
    c = Q(-7,4)
    f = iterate(c,3)
    P = s.Poly(x**3+x**2/2-s.Rational(9,4)*x-s.Rational(1,8),x,domain=s.QQ)
    fixed = iterate(c,1)-s.Poly(x,x,domain=s.QQ)
    check(f-s.Poly(x,x,domain=s.QQ) == fixed*P**2, 'dynatomic double factor')
    check(s.discriminant(P.as_expr(),x) == 49, 'three distinct real roots')
    check(s.gcd(P,fixed).degree() == 0, 'exact period three')
    check((f.diff()-s.Poly(1,x,domain=s.QQ)).rem(P).is_zero, 'parabolic multiplier one')
    check(s.gcd(f.diff(),P).degree() == 0, 'cycle avoids ramification')
    vals = (Q(-7,4), Q(5,4), Q(-1,4))
    multiplier = Q(1)
    for i,z in enumerate(vals):
        check(z*z-Q(29,16) == vals[(i+1)%3], 'rational cycle')
        multiplier *= 2*z
    check(multiplier == Q(35,8), 'rational multiplier distinct from secant multiplier')
    for z in vals:
        for momentum in (0,1,-3):
            end = evaluate_lift(Q(-29,16),3,z,momentum)
            check(end == (z,Q(8,35)*momentum), 'reciprocal momentum multiplier')
            check((end == (z,Q(momentum))) == (momentum == 0), 'rational cycle fibre return iff')
    # The degree-two lift is symplectic but noninjective on its actual domain.
    q0 = x*x-s.Rational(7,4)
    q1 = p/(2*x)
    check(s.cancel(s.diff(q0,x)*s.diff(q1,p)-s.diff(q0,p)*s.diff(q1,x)) == 1,
          'rational Jacobian one')
    check(evaluate_lift(Q(-7,4),1,1,2) == evaluate_lift(Q(-7,4),1,-1,-2),
          'retained sheet collision')
    check(evaluate_lift(Q(-7,4),1,1,2) == (Q(-3,4),Q(1)), 'explicit collision endpoint')

    # Independent expression substitution, without the coefficient reader.
    for c in params[:2]:
        for n in (1,2):
            _, delta = postcritical(c,n)
            f = iterate(c,n).as_expr()
            for j in (1,2):
                target = s.Poly(delta.as_expr()*Y**j,X,Y,domain=s.QQ)
                literal = s.cancel(target.as_expr().subs({X:f,Y:p/s.diff(f,x)},simultaneous=True))
                check(literal == compile_return(c,n,target).extension.as_expr(),
                      'independent simultaneous expression substitution')

    good = compile_return(Q(-7,4),1,s.Poly((X+s.Rational(7,4))*Y,X,Y,domain=s.QQ))
    hostiles = (
        lambda:orbit(True,1), lambda:orbit(Q(-7,4),True), lambda:iterate(0.5,1),
        lambda:postcritical(0,1), lambda:postcritical(-2,1), lambda:postcritical(Q(-7,4),0),
        lambda:compile_return(Q(-7,4),1,s.Poly(Y,X,Y,domain=s.QQ)),
        lambda:compile_return(Q(-7,4),1,s.Poly(x,x,domain=s.QQ)),
        lambda:compile_return(Q(-7,4),1,s.Poly(0.1*X,X,Y)),
        lambda:audit_receipt(replace(good,depth=True)),
        lambda:audit_receipt(replace(good,extension=good.extension+1)),
        lambda:audit_receipt(replace(good,parameter=Q(-29,16))),
        lambda:evaluate_lift(Q(-7,4),1,True,1),
    )
    for fn in hostiles:
        check(rejects(fn), 'malformed/forged/unguarded input rejected')
    print('quadratic_clock_return_20261012a')
    print('Parameters: -7/4, -29/16, 3/2; full structural controls at depths 1..4.')
    print('Critical points at depth n: 2^n-1; distinct critical values and nodes: n.')
    print('Return criterion: coefficient_j divisible by Delta_n^ceil(j/2).')
    print('Every pair of full-field polynomial observables has zero Jacobian at each critical zero-section point.')
    print(f'Compiled positive receipts: {compiled}; independent coefficient controls: {monomial_controls}.')
    print('Cycle multipliers: -7/4 parabolic cubic = 1; -29/16 rational triple = 35/8.')
    print('Hostiles: repeated critical values retained; c=0 higher ramification; c=-2 finite postcritical orbit.')
    print(f'Exact malformed/forged API controls: {len(hostiles)}.')
    print('No polynomial Keller map, birational inverse, source grounding, or JC(2) conclusion.')
    print(f'Always-active exact checks: {CHECKS}')


if __name__ == '__main__':
    main()
