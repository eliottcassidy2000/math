"""Exact controls for the four-target refinement packet. No removable asserts.

Run from anywhere with Python 3 and SymPy; stdout is deterministic JSON.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd, comb, prod
import json
import sympy as sp


def check(ok, message):
    if not ok:
        raise RuntimeError(message)


def sparse_mul(a, b, n):
    out = Counter()
    for i, x in a.items():
        for j, y in b.items():
            out[(i + j) % n] += x * y
    return {i: x for i, x in out.items() if x}


def shift_sum(a, b):
    n = len(a)
    out = [0] * n
    for j, y in b.items():
        for i, x in enumerate(a):
            out[(i + j) % n] += x * y
    return out


def coeff0(a, b):
    return sum(x * b[-i % len(a)] for i, x in enumerate(a))


def packets(n, g):
    h, v = [], 1
    while v not in h:
        h.append(v)
        v = v * g % n
    check(v == 1 and len(h) % 2 == 1, "odd unit subgroup")
    check(all(gcd(v - 1, n) == 1 for v in h if v != 1), "fixed-point-free guard")
    left, out = set(range(1, n)), []
    while left:
        orbit = {min(left) * v % n for v in h}
        negative = {-v % n for v in orbit}
        check(not orbit & negative, "converse guard")
        out.append({**dict.fromkeys(orbit, 1), **dict.fromkeys(negative, -1)})
        left -= orbit | negative
    return out


def fourth_of_symbol(x):
    # Independent autocorrelation route, not the incremental convolution state.
    n = len(x)
    corr = [sum(x[j] * x[(j + d) % n] for j in range(n)) for d in range(n)]
    return n * sum(v * v for v in corr)


def enumerate_potential(n, ps, prefix):
    values = []
    for tail in product((-1, 1), repeat=len(ps) - len(prefix)):
        x = [0] * n
        for eps, b in zip(tuple(prefix) + tail, ps):
            for j, v in b.items():
                x[j] = eps * v
        values.append(fourth_of_symbol(x))
    return Fraction(sum(values), len(values))


def greedy(n, g, pins=()):
    ps = packets(n, g)
    bs2 = [sparse_mul(b, b, n) for b in ps]
    bs3 = [sparse_mul(b2, b, n) for b2, b in zip(bs2, ps)]
    r = [0] * n
    for b2 in bs2:
        for j, v in b2.items():
            r[j] += v
    check(r == [-(n - 1)] + [1] * (n - 1), "exact covariance")
    sum_b4 = n * sum(sum(v * b2.get(-j % n, 0) for j, v in b2.items()) for b2 in bs2)
    potential = 3 * n * n * (n - 1) - 2 * sum_b4
    initial = potential
    x, x2, signs, audit = [0] * n, [0] * n, [], []
    pin_potential = initial
    for k, (b, b2, b3) in enumerate(zip(ps, bs2, bs3)):
        xb = shift_sum(x, b)
        score = coeff0(x2, xb) + 3 * coeff0(xb, r)
        score -= 2 * sum(x[-j % n] * v for j, v in b3.items())
        eps = pins[k] if k < len(pins) else (-1 if score > 0 else 1)
        check(eps in (-1, 1), "pin type")
        old = potential
        potential += 4 * n * eps * score
        x2 = [a + 2 * eps * v for a, v in zip(x2, xb)]
        for j, v in b2.items():
            x2[j] += v
            r[j] -= v
        for j, v in b.items():
            x[j] += eps * v
        signs.append(eps)
        if k < len(pins):
            pin_potential = potential
        else:
            check(potential <= old, "conditional potential must not increase")
        if n <= 19:
            check(enumerate_potential(n, ps, signs) == potential, "exhaustive conditional expectation")
            audit.append(potential)
    check(not any(r), "no unassigned packet")
    check(all(x[-i % n] == -v for i, v in enumerate(x)), "final skewness")
    check(all(abs(v) == 1 for v in x[1:]) and x[0] == 0, "tournament symbol")
    independent = fourth_of_symbol(x)
    check(independent == potential == n * coeff0(x2, x2), "final fourth moment independent replay")
    check(potential <= pin_potential, "pinned certificate")
    defect = potential - n*n*(n-1)
    parity_floor = 4*n*(n-1) if n % 4 == 1 else 0
    check(defect >= parity_floor, "parity floor for every regular tournament")
    return {"n": n, "generator": g, "packets": len(ps), "pins": list(pins),
            "unconditional_fourth_budget": initial, "pinned_fourth_budget": pin_potential,
            "final_trace_fourth": potential, "conditional_checks": len(audit),
            "fourth_defect": defect, "parity_defect_floor": parity_floor,
            "attains_regular_tournament_floor": defect == parity_floor,
            "sign_word": ''.join('+' if v == 1 else '-' for v in signs)}


def moment_tests():
    cases = [(7, 2), (13, 3), (19, 7), (43, 36), (133, 11), (1807, 1764)]
    out = [greedy(n, g) for n, g in cases]
    out += [greedy(9,1), greedy(11,3), greedy(121,3)]
    for n, g in cases[:3]:
        m = len(packets(n, g))
        for pins in product((-1, 1), repeat=min(2, m)):
            out.append(greedy(n, g, pins))
    try:
        packets(21, 16)
    except RuntimeError:
        pass
    else:
        raise RuntimeError("non-free hostile accepted")
    return out


def add_scaled(a, b, scale=1):
    out = Counter(a)
    for i, v in b.items():
        out[i] += scale*v
    return {i: v for i, v in out.items() if v}


def powers(b, degree, n):
    out = [{0: 1}]
    for _ in range(degree):
        out.append(sparse_mul(out[-1], b, n))
    return out


def high_moment_greedy(n, g, m):
    ps = packets(n, g)
    degree = 2*m
    bp = [powers(b, degree, n) for b in ps]
    moments = [{0: 1}] + [{} for _ in range(degree)]
    for bs in bp:
        new = []
        for j in range(degree+1):
            v = {}
            for ell in range(0, j+1, 2):
                v = add_scaled(v, sparse_mul(bs[ell], moments[j-ell], n), comb(j, ell))
            new.append(v)
        moments = new
    initial = n * (-1)**m * moments[degree].get(0, 0)
    bound = (n-1) * prod(range(1, degree, 2)) * n**m
    check(initial <= bound, "Gaussian even-moment bound")
    old, x, signs = initial, {}, []
    for b, bs in zip(ps, bp):
        remaining = []
        for j in range(degree+1):
            v = dict(moments[j])
            for ell in range(2, j+1, 2):
                v = add_scaled(v, sparse_mul(bs[ell], remaining[j-ell], n), -comb(j, ell))
            remaining.append(v)
        candidates = []
        for eps in (-1, 1):
            y = add_scaled(x, b, eps)
            yp = powers(y, degree, n)
            value = 0
            for j in range(degree+1):
                value += comb(degree, j) * sparse_mul(yp[degree-j], remaining[j], n).get(0, 0)
            candidates.append((n * (-1)**m * value, eps, y))
        check(candidates[0][0] + candidates[1][0] == 2*old, "exact conditional average")
        value, eps, x = min(candidates, key=lambda row: row[0])
        check(0 <= value <= old, "positive decreasing moment")
        old, moments = value, remaining
        signs.append(eps)
    final = n * (-1)**m * powers(x, degree, n)[-1].get(0, 0)
    check(final == old, "direct final high moment")
    if n <= 13:
        matrix = sp.Matrix(n, n, lambda i,j: x.get((j-i) % n, 0))
        check(int(sp.trace((-matrix*matrix)**m)) == final, "independent dense matrix moment")
    return {"n": n, "m": m, "initial": initial, "gaussian_bound": bound,
            "final": final, "sign_word": ''.join('+' if v == 1 else '-' for v in signs)}


def jacobian_tests():
    x, p = sp.symbols('x p')
    def bracket(a, b):
        return sp.expand(sp.diff(a, x) * sp.diff(b, p) - sp.diff(a, p) * sp.diff(b, x))
    def in_ideal(a, e):
        return all(i >= e or j >= e or (i >= e - 1 and j >= e - 1)
                   for (i, j), c in sp.Poly(a, x, p).terms() if c)
    pairs = 0
    for e in range(2, 9):
        monomials = [(i, j) for i in range(13) for j in range(13) if (i - j) % e == 0]
        for i, j in monomials:
            for a, b in monomials:
                coeff = i * b - j * a
                if coeff:
                    u, v = i + a - 1, j + b - 1
                    check(u >= e or v >= e or (u >= e-1 and v >= e-1), "invariant Poisson ideal")
                pairs += 1
        U, V, W = x**e, x*p, p**e
        check(bracket(U, V) == e*U, "UV bracket")
        check(bracket(V, W) == e*W, "VW bracket")
        check(bracket(U, W) == e*e*V**(e-1), "UW bracket")
        check(in_ideal(bracket(U, W), e), "generator ideal")
    examples = [(x**2, 0), (x**4, 0), (x**3+x**4, 0), ((x*x-1)**2, 1),
                (x**6*(x-1)**2, 0), ((x*x-1)**3, 1)]
    rows = []
    u = sp.symbols('u')
    for F, c in examples:
        fp = sp.diff(F, x)
        resultant = sp.Poly(sp.resultant(F-u, fp, x), u).sqf_part().monic().as_expr()
        completion = sp.cancel(resultant.subs(u, F) * p / fp)
        check(sp.denom(completion) == 1 or sp.denom(completion).is_number, "global polynomial completion")
        jac = bracket(F, completion)
        check(sp.expand(jac - resultant.subs(u, F)) == 0, "sharp bracket")
        e = min(i for (i,), v in sp.Poly(sp.expand(F.subs(x,x+c)-F.subs(x,c)), x).terms() if v)
        order = min(i for (i,), v in sp.Poly(sp.expand(jac.subs(x,x+c)), x).terms() if v)
        check(e == order, "sharp local ramification order")
        rows.append({"F": str(F), "critical_point": c, "ramification_index": e,
                     "critical_value_polynomial": str(resultant), "sharp_order": order})
    # Degree one is the genuine boundary; a non-invariant source observable also escapes.
    check(bracket(2*x, p/2) == 1 and bracket(x, p) == 1, "scope hostiles")
    return {"invariant_monomial_pairs": pairs, "sharp_examples": rows}


def shared_sheet_tests():
    x, p, u = sp.symbols('x p u')
    rows = []
    for F, c in [(x**6*(x-1)**2, 0), (x**6*(x-1)**2, 1),
                 (x**3*(x-1)**2, 1), (x*x*(x-1), 1)]:
        fp = sp.diff(F,x)
        b = F.subs(x,c)
        E = max(mult for _,mult in sp.factor_list(F-b,x)[1])
        def order_at_c(f):
            terms = sp.Poly(sp.expand(f.subs(x,x+c)),x).terms()
            return min(i for (i,),a in terms if a)
        e = order_at_c(F-b)
        # Only the three local generators are needed. Global H_E is compiled
        # directly from all (rational in this test universe) critical values.
        critical = sp.roots(fp,x)
        check(sum(critical.values()) == sp.degree(fp,x), "complete critical-point list")
        indices = {}
        for root, mult in critical.items():
            value = sp.simplify(F.subs(x,root))
            indices[value] = max(indices.get(value,1), mult+1)
        def coefficient(j):
            return sp.prod((u-value)**((j*(index-1)+index-1)//index)
                           for value,index in indices.items())
        H = []
        for j in (1,E):
            val = sp.cancel(coefficient(j).subs(u,F)/fp**j)
            check(sp.denom(val).is_number, "global coefficient polynomial")
            H.append(val)
        # {F,H1*p}, {H1*p,HE*p^E}, {F,HE*p^E}: homogeneous momentum degrees.
        coefficients = [sp.expand(fp*H[0]),
                        sp.expand(E*sp.diff(H[0],x)*H[1]-H[0]*sp.diff(H[1],x)),
                        sp.expand(E*fp*H[1])]
        measured = [order_at_c(v) for v in coefficients]
        check(measured == [e,E-e,E-1], "three sharp local bracket generators")
        # Finite generator compiler, all bracket pairs through the justified L.
        L = int(sp.ilcm(*indices.values())) if len(indices)>1 else next(iter(indices.values()))
        coeffs = [F] + [sp.cancel(coefficient(j).subs(u,F)/fp**j) for j in range(1,L+1)]
        checks = 0
        for j in range(L+1):
            for k in range(j+1,L+1):
                val = sp.cancel(k*sp.diff(coeffs[j],x)*coeffs[k]-j*coeffs[j]*sp.diff(coeffs[k],x))
                if val == 0:
                    continue
                degree = j+k-1
                order = order_at_c(val)
                check(order>=e or (order>=E-1 and degree>=E-1)
                      or (order>=E-e and degree>=E), "all global generator brackets in local ideal")
                checks += 1
        rows.append({"F": str(F), "sheet": int(c), "local_index": int(e),
                     "maximum_same_value_index": int(E), "period_L": int(L),
                     "sharp_x_orders": [int(v) for v in measured], "generator_pairs": checks})
    return rows


def global_ideal_test():
    # Independent ideal arithmetic: build source brackets, then compare their
    # reduced Groebner basis with the product of the predicted comaximal ideals.
    x,p,u = sp.symbols('x p u')
    F = x*x*(x-1)
    fp = sp.diff(F,x)
    D = u*(u+sp.Rational(4,27))
    H = [F, sp.cancel(D.subs(u,F)*p/fp), sp.cancel(D.subs(u,F)*p*p/fp**2)]
    brackets = [sp.diff(a,x)*sp.diff(b,p)-sp.diff(a,p)*sp.diff(b,x)
                for i,a in enumerate(H) for b in H[i+1:]]
    expected = [(x-1)*(x+sp.Rational(1,3))*a*b for a,b in product(
                [x*x,x*p,p*p], [(x-sp.Rational(2,3))**2,(x-sp.Rational(2,3))*p,p*p])]
    left, right = sp.groebner(brackets,p,x), sp.groebner(expected,p,x)
    check(left == right, "global source bracket ideal via independent Groebner basis")
    return {"F": str(F), "source_generators": 3, "predicted_product_generators": 9,
            "reduced_basis": [str(poly) for poly in left]}


def crt_tests():
    banks = [[2], [2, 3], [2, 3, 5], [4, 2, 5, 3, 7, 13], [7, 3, 6, 4, 8, 14],
             [-1, 2, 8, -2], [2, 3, 5, 7, 11]]
    rows = []
    for bank in banks:
        primes = set()
        for a in bank:
            primes.update(p for p in sp.factorint(abs(a)) if p != 2)
        q = 8
        for v in primes:
            q *= int(v)
        counts = Counter()
        for a in range(3, q, 4):
            if gcd(a, q) == 1:
                counts[tuple(int(sp.kronecker_symbol(b, a)) for b in bank)] += 1
        check(len(set(counts.values())) == 1, "uniform feasible residue fibres")
        check(len(counts) & (len(counts)-1) == 0, "affine binary rank")
        check(sum(counts.values()) == int(sp.totient(q))//2, "source residue total")
        rows.append({"bank": bank, "Q": q, "rank": len(counts).bit_length()-1,
                     "feasible_colours": len(counts), "classes_per_colour": next(iter(counts.values()))})
    return rows


def residual_tests():
    z = Fraction(2, 3)
    weights = [Fraction(1,6), Fraction(1,3), Fraction(1,2)]
    count, moves, improved = 0, 0, 0
    best_gain = Fraction(0)
    for targets in product(range(4), repeat=3):
        # 0 is terminal; 1,2,3 are the original active sources.
        life = []
        for source in range(1,4):
            seen, y, steps = set(), source, 0
            while y and y not in seen:
                seen.add(y)
                y = targets[y-1]
                steps += 1
            life.append(steps if y == 0 else None)
        az = [1-z**t if t is not None else Fraction(1) for t in life]
        low = [1-z**min(t,2) if t is not None else 1-z*z for t in life]
        actual = sum(w*a for w,a in zip(weights, az))
        def residual(v):
            out = [w-b for w,b in zip(weights,v)]
            for b, target in zip(v, targets):
                if target:
                    out[target-1] += z*b
            return out
        def upper(v, e):
            return (1-z)*sum(v) + sum(max(a,0) + min(a,0)*l for a,l in zip(e,low))
        for v in product((Fraction(-1,4), Fraction(0), Fraction(1,2)), repeat=3):
            e = residual(v)
            exact = (1-z)*sum(v)+sum(a*b for a,b in zip(e,az))
            check(exact == actual, "signed residual identity")
            bound = upper(v, e)
            lower = (1-z)*sum(v)+sum(max(a,0)*l+min(a,0) for a,l in zip(e,low))
            check(lower <= actual <= bound, "two-sided residual certificate")
            coarse = (1-z)*sum(v)+sum(max(a,0) for a in e)
            check(bound <= coarse, "retain negative defect credit")
            improved += bound < coarse
            best_gain = max(best_gain, coarse-bound)
            for i, a in enumerate(e):
                if a > 0:
                    new = list(v)
                    new[i] += a
                    check(upper(new,residual(new)) <= bound, "monotone native refinement")
                    moves += 1
            count += 1
    # Real accelerated Collatz clocks: prefix occupation keeps exact original sources.
    rows = []
    for source in (3, 27, 75, 127):
        path, n, valuations = [], source, []
        while n != 1:
            check(n not in path and len(path) < 1000, "finite control only")
            path.append(n)
            v = 3*n+1
            a = (v & -v).bit_length()-1
            valuations.append(a)
            n = v >> a
        for k in range(len(path)+1):
            occupation = sum(z**j for j in range(k))
            residual_value = z**k*(1-z**(len(path)-k))
            check((1-z)*occupation+residual_value == 1-z**len(path), "integer-clock refinement")
        rows.append({"source": source, "odd_steps": len(path), "valuation_sum": sum(valuations)})
    return {"finite_maps": 64, "signed_trials": count, "paid_moves": moves,
            "strict_negative_credit_trials": improved, "largest_credit": str(best_gain),
            "native_collatz_controls": rows}


if __name__ == '__main__':
    print(json.dumps({"moments": moment_tests(), "jacobian": jacobian_tests(),
                      "high_moments": [high_moment_greedy(7,2,3), high_moment_greedy(13,3,4),
                                       high_moment_greedy(19,7,5), high_moment_greedy(43,36,6)],
                      "crt": crt_tests(), "residual": residual_tests(),
                      "shared_sheets": shared_sheet_tests(),
                      "global_ideal": global_ideal_test()}, indent=2, sort_keys=True))
