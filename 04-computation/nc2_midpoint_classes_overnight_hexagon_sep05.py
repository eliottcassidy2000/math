#!/usr/bin/env python3
"""Exact full cut classes. Class root SIGNS are FINITE-EXACT, not universal.

No repository producer is imported. Endpoint/path coefficient sums are replayed
against complete auxiliary-row convolution, reversal, and actual multinomials.
Root intervals are cached per first row, never replaced by floating estimates.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial, gcd
from random import Random

import sympy as sp


t = sp.symbols("t")
GATES = 0
TRACE = sha256()


def require(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def record(*row):
    TRACE.update((repr(row)+"\n").encode())


def bc(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def path_row(p, q, U, V):
    lower, upper = -(V//p), U//q
    def coefficient(j):
        X, Y = U-q*j, V+p*j
        return bc(X+Y, Y) if min(X, Y) >= 0 else 0
    return lower, upper, coefficient


def residue_row(mass, A, residue):
    return 0, (mass-residue)//A, lambda j: bc(mass, residue+A*j) if j >= 0 else 0


def convolution(left, right, j):
    lo = max(left[0], j-right[1])
    hi = min(left[1], j-right[0])
    return sum(left[2](i)*right[2](j-i) for i in range(lo, hi+1))


def poly(values, epsilon):
    return sp.Poly(sum(c*t**(j+epsilon) for j, c in values.items() if c), t, domain=sp.QQ)


def class_rows(data, check=True):
    A, B, h, r, z, x = data
    p, y = B-A, B*h+r
    mass, C, epsilon = x+y+z, p*y+B*x, 2*z//A
    alpha = residue_row(mass, A, z)
    beta = path_row(p, B, y, x)
    beta2 = path_row(p, B, 2*y, 2*x)
    P = poly({j: alpha[2](j)*beta[2](j) for j in range(h+1)}, 0)
    hit = {j: convolution(alpha, alpha, j) for j in range(-1, 2*h+2)}
    virtual, actual = {}, {}
    for j in range(-1, 2*h+2):
        virtual[j] = hit[j]*convolution(beta, beta, j)
        u, v, w = 2*x+p*j, 2*y-B*j, 2*z+A*j
        actual[j] = factorial(2*mass)//(factorial(u)*factorial(v)*factorial(w)) if min(u, v, w) >= 0 else 0
        if check:
            require(actual[j] == bc(2*mass, 2*z+A*j)*beta2[2](j), "actual multinomial equals two unrestricted path counts")
    classes = {}
    for residue in range(A):
        if residue == z:
            continue
        partner = (2*z-residue) % A
        ell = (residue+partner-2*z)//A
        left, right = residue_row(mass, A, residue), residue_row(mass, A, partner)
        values = {}
        for j in actual:
            lower = max(residue, 2*z+A*j-mass)
            lower += (residue-lower) % A
            direct = sum(bc(mass, k)*bc(mass, 2*z+A*j-k)
                         for k in range(lower, min(mass, 2*z+A*j)+1, A))
            factored = convolution(left, right, j-ell)
            if check:
                require(direct == factored, "residue midpoint sum equals full residue-factor convolution")
            values[j] = direct*beta2[2](j)
        classes[("R", residue)] = poly(values, epsilon)
    for name, dx, dy, step in (("E", 1, 0, p), ("N", 0, 1, B)):
        for delay in range(1, step):
            values = {}
            if gcd(p, B) == 1:
                u = (delay*pow(p, -1, B)) % B
                v = (p*u-delay)//B
                prefix = path_row(p, B, y-u, x+v)
                suffix = path_row(p, B, y+u-dx, x-v-dy)
                normalized_suffix = path_row(p, B, y-B+u-dx, x+p-v-dy)
                if check:
                    require(1 <= u < B and 0 <= v < p, "canonical crossing witness residue ranges")
                    require(min(y-u, x+v, y-B+u-dx, x+p-v-dy) >= 0, "both paths of the j=1 crossing witness exist")
                    require(suffix[0] == normalized_suffix[0]+1 and suffix[1] == normalized_suffix[1]+1,
                            "full suffix support translates by exactly one index")
                    require(all(suffix[2](k) == normalized_suffix[2](k-1)
                                for k in range(suffix[0], suffix[1]+1)),
                            "unshifted suffix is exactly t times the complete nonnegative-endpoint suffix")
            for j in actual:
                X, Y = 2*y-B*j, 2*x+p*j
                direct = 0
                if min(X, Y) >= 0:
                    for i in range(X-dx+1):
                        remaining = C-delay-p*i
                        if remaining % B:
                            continue
                        k = remaining//B
                        if 0 <= k <= Y-dy:
                            direct += bc(i+k, i)*bc(X-dx-i+Y-dy-k, X-dx-i)
                if check and gcd(p, B) == 1:
                    require(direct == convolution(prefix, suffix, j), "crossing-edge path sum equals complete prefix/suffix product")
                    require(direct == convolution(prefix, normalized_suffix, j-1),
                            "displayed factor t and j=l+k+1 agree with literal crossing-edge paths")
                values[j] = direct*hit[j]
            classes[(name, delay)] = poly(values, epsilon)
    Q, V = poly(actual, epsilon), poly(virtual, epsilon)
    if check:
        require(sum(classes.values(), sp.Poly(0, t)) == Q-V, "all residue/edge classes partition the actual joint defect")
        require(len(classes) == 2*B-3, "complete labelled class count")
        for label, response in classes.items():
            require(response.nth(1+epsilon) > 0, "every class has a strictly positive coefficient at raw j=1")
            name, value = label
            partner = ("R", (2*z-value) % A) if name == "R" else (name, (p if name == "E" else B)-value)
            require(response == classes[partner], "whole polynomial reversal identity")
    return P, Q, V, epsilon, classes


def interval_value(poly, a, b):
    lo = hi = F(0)
    for c in poly.all_coeffs():
        values = (lo*a, lo*b, hi*a, hi*b)
        lo, hi = min(values)+F(c), max(values)+F(c)
    return lo, hi


def sign_classes(data, P, epsilon, classes):
    intervals = []
    for (a, b), mult in P.intervals():
        require(mult == 1 and b <= 0 and P.eval(0) != 0, "first row has all simple negative roots")
        intervals.append((F(a), F(b)))
    require(len(intervals) == P.degree(), "all first roots are isolated")
    A, B = data[:2]
    z, p = data[4], B-A
    representatives, fixed = [], 0
    for label in classes:
        name, value = label
        partner = ("R", (2*z-value) % A) if name == "R" else (name, (p if name == "E" else B)-value)
        if label == partner:
            fixed += 1
        if label <= partner:
            representatives.append(label)
    require(fixed == 1 and len(representatives) == B-1, "one fixed class and exactly B-1 reversal orbits")
    for label in representatives:
        remainder = classes[label].rem(P)
        if epsilon % 2:
            remainder = -remainder
        certificates = []
        for index, (a, b) in enumerate(intervals):
            for _ in range(240):
                lo, hi = interval_value(remainder, a, b)
                if hi < 0:
                    break
                if lo > 0 or a == b:
                    raise RuntimeError(("class root sign hostile", data, label, a, b, lo, hi))
                a, b = map(F, P.refine_root(a, b, eps=sp.Rational(b-a)/4))
            else:
                raise RuntimeError(("class sign refinement exhausted", data, label))
            intervals[index] = a, b
            certificates.append(tuple(map(str, (a, b, lo, hi))))
        require(len(certificates) == P.degree(), "every reversal-orbit class has strict raw negative values at every first root")
        record("FINITE-EXACT class sign", data, label, certificates)
    return len(representatives)


def eligible(row):
    A, B, h, r, z, x = row
    return gcd(A, B) == 1 and A*x-(B-A)*z > 0 and gcd(A*(B*h+r)+B*z, x+B*h+r+z) == 1


def banks():
    head = [(A, B, h, r, z, x) for B in range(2, 7) for A in range(1, B)
            for h in (1, 2, 3) for r in range(B) for z in range(A) for x in (1, 2, 3)]
    rng, wide = Random(4440), []
    for _ in range(300):
        B = rng.randrange(3, 13)
        A = rng.randrange(1, B)
        h, r, z = rng.randrange(1, 9), rng.randrange(B), rng.randrange(A)
        x = rng.choice((1, 2, 3, 5, 7, 11, 23, 101, 997))
        wide.append((A, B, h, r, z, x))
    return [(name, [row for row in bank if eligible(row)]) for name, bank in (("head", head), ("wide", wide))]


def finite_controls():
    distinct = set()
    for name, bank in banks():
        labels = orbits = root_values = 0
        for row in bank:
            P, Q, V, epsilon, classes = class_rows(row)
            orbit_count = sign_classes(row, P, epsilon, classes)
            labels += len(classes)
            orbits += orbit_count
            root_values += len(classes)*P.degree()
            distinct.add(row)
        print(f"{name}: {len(bank)} indexed rows; {labels} nonzero labelled classes; {orbits} reversal orbits; {root_values} labelled first-root values negative")
    print(f"distinct tuples across banks: {len(distinct)}; no universal class-sign inference")


def structural_boundary_controls():
    cases = 0
    for B in range(2, 6):
        for A in range(1, B):
            if gcd(A, B) != 1:
                continue
            for h in (1, 2):
                for r in (0, B-1):
                    for z in (0, A-1):
                        for x in (0, 1):
                            P, Q, V, epsilon, classes = class_rows((A, B, h, r, z, x))
                            for response in classes.values():
                                stripped = response.exquo(sp.Poly(t**min(mon[0] for mon, coefficient in response.terms() if coefficient), t))
                                intervals = stripped.intervals()
                                require(sum(mult for _, mult in intervals) == stripped.degree(), "bounded independent real-root count of a full class response")
                                require(all(mult == 1 and b <= 0 for (_, b), mult in intervals), "nonzero class roots are simple negative, positive constants allowed")
                            cases += 1
    print(f"structural boundary bank: {cases} indexed coprime cases, x=0 and corner residues retained")


def hostiles():
    P, Q, V, epsilon, classes = class_rows((2, 4, 1, 0, 0, 0), check=False)
    missing = Q-V-sum(classes.values(), sp.Poly(0, t))
    require(Q-V == sp.Poly(396*t+32*t*t, t, domain=sp.QQ), "noncoprime exact actual defect")
    require(missing == sp.Poly(108*t, t, domain=sp.QQ), "omitted unselected beta midpoint coset is not an edge crossing")
    s = sp.symbols("s")
    H = sp.Poly((s+10)**4*(s-1)*(s-8), s)
    K = sp.Poly((s+10)**4*(s-7)*(s-100), s)
    require(H.nth(2) == -21200 and K.nth(2) == 2000 and (H*K).nth(4) == 1454400000, "common negative factor plus interlacing positive factors does not preserve mixed signed duplication")
    require(1 < 7 < 8 < 100, "strict positive-root interlacing in the mixed hostile")
    print("hostiles: noncoprime cut omits108t; shared negative binomial factor/interlacing/middle signs still allow mixed coefficient1454400000>0")
    record("hostiles", missing.all_coeffs(), H.all_coeffs(), K.all_coeffs())


def main():
    print("WHOLE MIDPOINT RESIDUE AND CROSSING CLASSES")
    print("scope=PROVED exact decomposition/reversal/nonvacuity and cited class root geometry; rootwise signs FINITE-EXACT only")
    structural_boundary_controls()
    finite_controls()
    hostiles()
    print("OPEN all-height class signs and actual joint weighted pencil; no trinomial closure claimed")
    print("trace_sha256="+TRACE.hexdigest())
    print(f"PASS explicit_gates={GATES}")


if __name__ == "__main__":
    main()
