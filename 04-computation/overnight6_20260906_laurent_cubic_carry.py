#!/usr/bin/env python3
"""Exact four-channel family certificate, no repository imports.

Support (-21,2g-21,3g-21), integer g>=11, gcd(g,21)=1.
Symbolic characteristic positivity proves the unbounded statement.
Named literal Laurent controls verify the source/consumer map separately.
Requires SymPy; no assertions are used for verification gates.
"""
import hashlib
import json
from math import factorial, gcd
import sympy as S

G, T, Z = S.symbols('g tau s')
CHECKS = 0


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def falling(x, k):
    return S.prod(x-j for j in range(k))


def raw_channels(a, b, c, mass):
    """Original charge equation, solving for the positive c count."""
    answer = []
    for x in range(mass+1):
        numerator = a*x-b*(mass-x)
        if numerator % (c-b):
            continue
        z = numerator//(c-b)
        y = mass-x-z
        if min(y,z) >= 0:
            answer.append(((x,y,z), factorial(mass)//
                          (factorial(x)*factorial(y)*factorial(z))))
    return sorted(answer, key=lambda item: item[0][2])


def literal_rows(a, b, c, maximum):
    """Repeated Laurent multiplication retaining the gamma exponent.

    At fixed mass, charge and gamma exponent determine the other counts.
    Thus this recovers the complete coefficient polynomial, not just CT at
    one coefficient specialization.
    """
    state = {(0,0): 1}
    answer = []
    for _ in range(maximum):
        following = {}
        for (charge,z), multiplicity in state.items():
            for dc,dz in ((-a,0),(b,0),(c,1)):
                key = charge+dc,z+dz
                following[key] = following.get(key,0)+multiplicity
        state = following
        answer.append({z:value for (charge,z),value in state.items() if charge == 0})
    return answer


def coefficient_record(poly):
    """Canonical positive-denominator polynomial record in descending order."""
    numerator, denominator = S.fraction(S.cancel(poly))
    require(not denominator.has(G), 'constant coefficient denominator')
    if denominator < 0:
        numerator,denominator = -numerator,-denominator
    coefficients = S.Poly(S.expand(numerator.subs(G,Z+11)),Z).all_coeffs()
    require(all(c > 0 for c in coefficients), 'every shifted coefficient positive')
    return {'denominator':str(denominator),
            'shift_g':11, 'coefficients_descending':[str(c) for c in coefficients]}


def main():
    p_integer = (72*T**3+504*(G-7)*T**2+84*(G-7)*(G-8)*T+
                 (G-7)*(G-8)*(G-9))
    domain = S.QQ.frac_field(G)
    p = S.Poly(p_integer/72,T,domain=domain)
    first_content = falling(G,7)/factorial(9)
    first = S.Poly(sum(falling(G,10-j)*T**j/
                         (factorial(9-3*j)*factorial(1+2*j)) for j in range(4)),T)
    require(S.cancel(first.as_expr()-first_content*p_integer) == 0,
            'complete first row and its positive content')
    kfactor = falling(2*G,14)
    qbar = S.Poly(sum(falling(2*G-14,7-j)*T**j/
                       (factorial(21-3*j)*factorial(2*j)) for j in range(8)),T,domain=domain)
    qraw = S.Poly(sum(falling(2*G,21-j)*T**j/
                       (factorial(21-3*j)*factorial(2*j)) for j in range(8)),T)
    require(S.cancel(qraw.as_expr()-kfactor*qbar.as_expr()) == 0,
            'complete eight-channel second row')
    inverse = S.invert(S.Poly(T,T,domain=domain),p)
    require((S.Poly(T,T,domain=domain)*inverse).rem(p).as_expr() == 1,
            'inverse lower carry')
    response = (qbar*inverse).rem(p)
    response = S.Poly(sum(S.cancel(response.nth(j))*T**j for j in range(3)),T)
    for j in range(3):
        require(S.Poly(response.nth(j),G).degree() <= 6-j,
                'response coefficient weighted degree')
    # The apparent inverse denominator cancels in the sole j=0 term.
    require(S.cancel(qbar.nth(0)/p.nth(0)-
            1152*(G-10)*(2*G-15)*(2*G-17)*(2*G-19)/factorial(21)) == 0,
            'polynomial cancellation of inverse denominator')
    companion = S.Matrix([[0,0,-p.nth(0)],[1,0,-p.nth(1)],[0,1,-p.nth(2)]])
    multiplication = (response.nth(0)*S.eye(3)+response.nth(1)*companion+
                      response.nth(2)*companion**2).applyfunc(S.cancel)
    trace = S.factor(S.trace(multiplication))
    e2 = S.factor((trace**2-S.trace(multiplication**2))/2)
    norm = S.factor(multiplication.det(method='domain-ge'))
    characteristic = [S.factor(-trace),e2,S.factor(-norm)]
    # Independently expanding the determinant certifies all three signs belong
    # to the characteristic polynomial, rather than to unrelated statistics.
    x = S.symbols('x')
    detpoly = S.Poly((x*S.eye(3)-multiplication).det(method='domain-ge'),x)
    require(detpoly.nth(3) == 1, 'monic quotient characteristic polynomial')
    records = {}
    for i,value in enumerate(characteristic,1):
        require(S.cancel(value-detpoly.nth(3-i)) == 0,
                'quotient characteristic coefficient')
        require(S.Poly(value,G).degree() <= 6*i, 'characteristic coefficient degree')
        records[str(i)] = coefficient_record(value)
    discriminant = S.factor(S.discriminant(p.as_expr(),T))
    disc_core = 74863*G**3-1610102*G**2+11532575*G-27511000
    require(S.cancel(discriminant-(G-8)*(G-7)**2*disc_core/1728) == 0,
            'first cubic discriminant factorization')
    require(S.Poly(disc_core.subs(G,Z+11),Z).all_coeffs() ==
            [74863,860377,3285600,4167636], 'positive shifted discriminant core')
    for j in range(4):
        require(all(c > 0 for c in S.Poly(p_integer.coeff(T,j).subs(G,Z+11),Z).all_coeffs()),
                'positive first cubic coefficients')
    print('support (-21,2g-21,3g-21), g>=11 integral, gcd(g,21)=1')
    print('first_polynomial',p_integer)
    print('first_content',S.factor(first_content))
    print('second_positive_factor',S.factor(kfactor))
    print('first_discriminant',discriminant)
    print('scaled_anchored_response',S.factor(response.as_expr()))
    print('characteristic_positive_shift_certificates',json.dumps(records,sort_keys=True))

    manifest = []
    for g in (11,13,16,17,19,23,25,29):
        a,b,c = 21,2*g-21,3*g-21
        require(gcd(g,21) == 1 and gcd(a,gcd(b,c)) == 1 and 0<b<c,
                'primitive one-negative source domain')
        native = literal_rows(a,b,c,2*g)
        raw = [raw_channels(a,b,c,mass) for mass in (g,2*g)]
        require([v for v,_ in raw[0]] == [(g-10+j,9-3*j,1+2*j) for j in range(4)],
                'complete four-channel first fiber')
        require([v for v,_ in raw[1]] == [(2*g-21+j,21-3*j,2*j) for j in range(8)],
                'complete eight-channel second fiber')
        require(all(not native[m-1] for m in range(1,2*g) if m != g),
                'no omitted first or intermediate support return')
        for mass,row,formula in zip((g,2*g),raw,(first,qraw)):
            require(native[mass-1] == {v[2]:weight for v,weight in row},
                    'literal Laurent multiplication versus raw source equation')
            require(S.Poly(sum(weight*T**j for j,(_,weight) in enumerate(row)),T) ==
                    S.Poly(formula.as_expr().subs(G,g),T), 'full symbolic row specializes')
        pf,qf = [S.Poly(sum(weight*T**j for j,(_,weight) in enumerate(row)),T)
                 for row in raw]
        rf = S.Poly(S.rem(qf.as_expr()*S.invert(T,pf.as_expr(),T),pf.as_expr(),T),T)
        require(S.cancel(rf.as_expr()-kfactor.subs(G,g)*response.as_expr().subs(G,g)) == 0,
                'independent specialized anchored quotient')
        require(pf.count_roots(-S.oo,0) == 3 and S.gcd(pf,pf.diff()).degree() == 0,
                'three simple attainable negative roots')
        require(S.gcd(pf,qf).degree() == 0, 'complete two rows coprime')
        # Complete weighted Hermite trace form is an independent sign control.
        tf = S.zeros(3)
        for column in range(3):
            remainder = S.Poly(S.rem(T**(column+1),pf.as_expr(),T),T)
            for row in range(3):
                tf[row,column] = remainder.nth(row)
        vf = sum((rf.nth(j)*tf**j for j in range(3)),S.zeros(3))
        hermite = S.Matrix(3,3,lambda i,j:S.trace(vf*tf**(i+j)))
        minors = [(-1)**k*hermite[:k,:k].det() for k in range(1,4)]
        require(all(value > 0 for value in minors), 'full corrected Hermite form negative definite')
        canonical = sum((qf.nth(j)*tf**j for j in range(8)),S.zeros(3))
        wrong_form = S.Matrix(3,3,lambda i,j:S.trace(canonical*tf**(i+j)))
        require(all(wrong_form[:k,:k].det() > 0 for k in range(1,4)),
                'omitting lower carry reverses the sign at every first root')
        manifest.append({'g':g,'support':[-a,b,c],
                         'first_coefficients':[int(pf.nth(j)) for j in range(4)],
                         'signed_Hermite_minors':[str(v) for v in minors]})
    # Scope/adversarial controls, separate from admissible family parameters.
    bad_parameter = 12
    bad_rows = literal_rows(21,2*bad_parameter-21,3*bad_parameter-21,bad_parameter)
    require(next(j+1 for j,row in enumerate(bad_rows) if row) == 4,
            'dropping gcd gives earlier mass four at g=12')
    p11 = 3*T**3+84*T**2+42*T+1
    require(p11.subs(T,-28) == -1175 and p11.subs(T,-27) == 1054,
            'attainable cancellation outside real-rooted ordinary cubic-core sector')
    require(p11.count_ops() > 0 and S.Rational(-27) < -S.Rational(4,27),
            'far negative phase lies outside cubic SOS sector')
    require(S.expand((1+T**2)**2).coeff(T,2) == 2,
            'ordinary duplication requires real-rooted input')
    require(S.expand((T**2)**2).coeff(T,2) == 0,
            'ordinary duplication requires interior coefficient')
    print('named_source_controls',json.dumps(manifest,sort_keys=True))
    print('hostiles gcd_dropped_g12_first_mass=4; canonical_Q_has_positive_root_values; '
          'g11_root_in(-28,-27)_outside_real_core_SOS')
    print('PASS',CHECKS,'explicit gates; no assertions; symbolic all-height certificate')
    digest = hashlib.sha256(json.dumps({'certificates':records,'controls':manifest},
                           sort_keys=True,separators=(',',':')).encode()).hexdigest()
    print('semantic_sha256',digest)


if __name__ == '__main__':
    main()
