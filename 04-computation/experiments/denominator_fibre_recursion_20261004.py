"""Exact controls for golden periodic lifts, phase towers, and signed splices.

All mathematical quantifiers are proved in the companion note. These bounded
controls are independent lattice/phase/series/replay checks, never coverage
assumptions for the signed Collatz problem. Checks survive python -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import gcd, lcm
from pathlib import Path
import json
import runpy

ROOT = Path(__file__).resolve().parents[2]
PHI, BETA = (F(0), F(1)), (F(-1), F(1))
ZERO, ONE = (F(0), F(0)), (F(1), F(0))
ROOTS = {1: ((0, 1, 2), (1, 0, 0)),
         -1: ((1, 0, 1), (1, 0)),
         -5: ((-1, 7, 11), (1, 0, 1, 0, 0)),
         -17: ((9, 41, 76), tuple(map(int, "101010100101010000")))}


def need(ok, why):
    if not ok:
        raise ValueError(why)


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a*c+b*d, a*d+b*c+b*d


def inv(x):
    a, b = x
    norm = a*a+a*b-b*b
    need(norm != 0, "zero field inverse")
    return (a+b)/norm, -b/norm


def power(x, n):
    if n < 0:
        return power(inv(x), -n)
    out = ONE
    while n:
        if n & 1:
            out = mul(out, x)
        x = mul(x, x)
        n //= 2
    return out


def normalized(x):
    q = lcm(x[0].denominator, x[1].denominator)
    return int(x[0]*q), int(x[1]*q), q


def sign(x):
    a, b, _ = normalized(x)
    return sign_integer_pair(a, b)


def sign_integer_pair(a, b):
    u, v = 2*a+b, b
    if not v:
        return (u > 0)-(u < 0)
    if not u or u*v >= 0:
        return 1 if v > 0 else -1
    return (1 if u > 0 else -1)*(1 if u*u > 5*v*v else -1)


def conjugate(x):
    return x[0]+x[1], -x[1]


def periodic_value(bits):
    acc = ZERO
    for bit in reversed(bits):
        acc = mul(BETA, add((F(bit), F(0)), acc))
    return mul(acc, inv(add(ONE, neg(power(BETA, len(bits))))))


def prefixed_value(bits, endpoint):
    acc = endpoint
    for bit in reversed(bits):
        acc = mul(BETA, add((F(bit), F(0)), acc))
    return acc


def phase_step(v, q):
    return v[1] % q, (v[0]+v[1]) % q


def jordan2(q):
    value, left, prime = q*q, q, 2
    while prime*prime <= left:
        if left % prime == 0:
            value = value//(prime*prime)*(prime*prime-1)
            while left % prime == 0:
                left //= prime
        prime += 1
    if left > 1:
        value = value//(left*left)*(left*left-1)
    return value


@lru_cache(None)
def periodic_lift(q, phase):
    """Canonical exact periodic representative of a rational torus phase.

    The zero phase is assigned 0; retain an extra boundary tag if 1 or beta
    was intended. Termination follows from the proved conjugate contraction.
    """
    need(q >= 1, "positive denominator")
    target = tuple(v % q for v in phase)
    if target == (0, 0):
        return ZERO
    a, b = target
    while sign_integer_pair(a-q, b) >= 0:
        a -= q
    state, path, positions = (a, b), [], {}
    while state not in positions:
        positions[state] = len(path)
        path.append(state)
        a, b = state
        d = int(sign_integer_pair(b-q, a+b) > 0)
        state = b-d*q, a+b
    cycle = path[positions[state]:]
    matching = [(a, b) for a, b in cycle if (a % q, b % q) == target]
    need(len(matching) == 1, "unique nonzero periodic representative")
    a, b = matching[0]
    return F(a, q), F(b, q)


@lru_cache(None)
def phase_cycles(q):
    points = {(a, b) for a in range(q) for b in range(q) if gcd(a, b, q) == 1}
    cycles = []
    while points:
        first = min(points)
        cycle, v = [], first
        while v not in cycle:
            cycle.append(v)
            points.remove(v)
            v = phase_step(v, q)
        need(v == first, "invertible phase permutation")
        cycles.append(tuple(cycle))
    return tuple(cycles)


def phase_orbit(x):
    a, b, q = normalized(x)
    start, cycle = (a % q, b % q), []
    v = start
    while v not in cycle:
        cycle.append(v)
        v = phase_step(v, q)
    need(v == start, "phase orbit closes")
    return frozenset(cycle)


@lru_cache(None)
def annihilator(q, a, b):
    # The denominator ideal is the inverse image of this kernel modulo q.
    return tuple((u, v) for u in range(q) for v in range(q)
                 if (u*a+v*b) % q == 0 and (u*b+v*(a+b)) % q == 0)


def ideal(x):
    a, b, q = normalized(x)
    return q, annihilator(q, a % q, b % q)


def lattice_census(q, inherited):
    states = inherited["trap"](q)
    cycles = inherited["cycles_in"](states, q)
    projected, seen, integer_cycles = [], set(), []
    for cycle in cycles:
        bits = [inherited["beta"](v, q)[1] for v in cycle]
        x = periodic_value(bits)
        need(x == (F(cycle[0][0], q), F(cycle[0][1], q)), "independent periodic-series lift")
        phases = tuple((a % q, b % q) for a, b in cycle)
        if q > 1:
            need(len(set(phases)) == len(phases), "golden period equals phase period")
            need(not seen.intersection(phases), "distinct periodic lifts")
            projected.append(frozenset(phases))
            seen.update(phases)
            need(periodic_lift(q, phases[0]) == x, "constructive lift from independent initial representative")
        for a, b in cycle:
            point = F(a, q), F(b, q)
            other = conjugate(point)
            need(sign(add(other, PHI)) >= 0 and sign(add(ONE, neg(other))) >= 0, "periodic conjugate window")
            if sign(add(point, neg(BETA))) > 0:
                need(sign(add(other, BETA)) >= 0, "second rectangle branch")
        decoded = inherited["decode_cycle"](cycle, q)
        if all(n.denominator == 1 for n in decoded["orbit"]):
            integer_cycles.append(sorted(int(n) for n in decoded["orbit"]))
    if q == 1:
        need({v for c in cycles for v in c} == {(0, 0), (1, 0), (-1, 1)}, "zero-phase exception")
    else:
        direct = phase_cycles(q)
        need(set(projected) == {frozenset(c) for c in direct}, "independent entire phase universe")
        need(len(seen) == sum(gcd(a, b, q) == 1 for a in range(q) for b in range(q)), "Jordan count")
        need(len(seen) == jordan2(q), "independent prime-product count")
    return dict(q=q, trapped_states=len(states), cycles=len(cycles),
                periodic_points=sum(map(len, cycles)),
                periods=dict(sorted(Counter(map(len, cycles)).items())),
                integer_cycles=integer_cycles)


def lift_phase_children(prime, q, parent):
    """Explicit child labels over q=prime^a, prime in {2,3}, a>=1."""
    need(prime in (2, 3), "only audited 2/3 towers")
    need(isinstance(q, int) and q >= prime, "positive prime-power parent candidate")
    reduced = q
    while reduced % prime == 0:
        reduced //= prime
    need(reduced == 1, "prime-power parent")
    L, v = len(parent), min(parent)
    a, b = power(PHI, L)
    mv = a*v[0]+b*v[1], b*v[0]+(a+b)*v[1]
    need(all((mv[i]-v[i]) % q == 0 for i in (0, 1)), "parent monodromy")
    c = tuple(int((mv[i]-v[i])//q) % prime for i in (0, 1))
    need(c != (0, 0), "nonzero return translation")
    children = phase_cycles(prime*q)
    lookup = {point: i for i, child in enumerate(children) for point in child}
    labels = defaultdict(set)
    for w in product(range(prime), repeat=2):
        z = v[0]+q*w[0], v[1]+q*w[1]
        label = (c[0]*w[1]-c[1]*w[0]) % prime
        labels[label].add(lookup[z])
        returned = z
        for _ in range(L):
            returned = phase_step(returned, prime*q)
        need(returned == tuple(v[i]+q*((w[i]+c[i]) % prime) for i in (0, 1)), "affine return-carry law")
    need(set(labels) == set(range(prime)) and all(len(ids) == 1 for ids in labels.values()), "determinant child selector")
    selected = {next(iter(ids)) for ids in labels.values()}
    need(len(selected) == prime, "exact child number")
    for i in selected:
        child = children[i]
        need(len(child) == prime*L, "child covers parent prime times")
        need({(x % q, y % q) for x, y in child} == set(parent), "child reduction")
    return dict(parent=list(v), translation=list(c),
                children=[dict(label=j, anchor=list(min(children[next(iter(labels[j]))])))
                          for j in range(prime)])


def word_summary(word):
    A = S = 0
    for a in word:
        S = 3*S+2**A
        A += a
    return len(word), A, S


def odd_step(n):
    z, a = 3*n+1, 0
    need(z != 0, "integer odd source cannot hit zero")
    while z % 2 == 0:
        z //= 2
        a += 1
    return z, a


def splice_family(word, root):
    p, A, S = word_summary(word)
    modulus, T, C = 3**(p+1), 2*3**p, 2**A+3*S
    choices = [a for a in range(1, T+1) if (pow(2, A+a, modulus)*root-C) % modulus == 0]
    need(len(choices) == 1, "unique full unit-clock phase")
    a = choices[0]
    while True:
        n = (2**(A+a)*root-C)//modulus
        node, okay = n, True
        for val in tuple(word)+(a,):
            okay &= abs(node) > 91 and (node > 0) == (root > 0)
            node, actual = odd_step(node)
            need(actual == val, "reverse-integrality source word")
        need(node == root, "compiled root")
        if okay:
            break
        a += T
    return dict(word=tuple(word), root=root, p=p, A=A, S=S, C=C,
                modulus=modulus, T=T, a=a, source=n)


def splice_controls():
    heads = [word for p in range(4) for word in product(range(1, 4), repeat=p)]
    roots = {r: (F(a, q), F(b, q)) for r, ((a, b, q), _) in ROOTS.items()}
    for r, (_, bits) in ROOTS.items():
        need(periodic_value(bits) == roots[r], "independent root periodic sum")
    checks, examples, phase_periods = 0, [], []
    for word in heads:
        for root in ROOTS:
            row = splice_family(word, root)
            p, A, S, C, T = (row[key] for key in ("p", "A", "S", "C", "T"))
            R = 2**T
            shift = (R-1)*C//row["modulus"]
            need((R-1)*C % row["modulus"] == 0, "integer recurrence shift")
            prefix_bits = tuple(bit for val in word for bit in (1,)+(0,)*val)+(1,)
            eta = prefixed_value(prefix_bits, ZERO)
            predicted_phase_period = len(phase_orbit(roots[root]))//gcd(len(phase_orbit(roots[root])), T)
            ra, rb, rq = ROOTS[root][0]
            phase, first, phases = (ra % rq, rb % rq), (ra % rq, rb % rq), []
            clock = power(PHI, -T)
            while phase not in phases:
                phases.append(phase)
                shifted = mul(clock, tuple(map(F, phase)))
                phase = tuple(int(v) % rq for v in shifted)
            need(phase == first and len(phases) == predicted_phase_period, "exact sibling phase period")
            if word == (1,)*p:
                phase_periods.append(dict(p=p, root=root, sibling_phase_period=predicted_phase_period))
            old_n = old_x = None
            for k in range(5):
                a = row["a"]+T*k
                n = (2**(A+a)*root-C)//row["modulus"]
                valuations = tuple(word)+(a,)
                bits = tuple(bit for val in valuations for bit in (1,)+(0,)*val)
                x = prefixed_value(bits, roots[root])
                L = A+p+a+1
                need(x == add(eta, mul(power(BETA, L), roots[root])), "golden shared shadow")
                need(ideal(x) == ideal(roots[root]), "full denominator ideal preserved")
                need(phase_orbit(x) == phase_orbit(roots[root]), "phase-orbit invariant")
                if old_n is not None:
                    need(n == R*old_n+shift, "arithmetic common recursion")
                    need(x == add(eta, mul(power(BETA, T), add(old_x, neg(eta)))), "golden common recursion")
                raw = n
                for digit in bits:
                    need(raw % 2 == digit, "independent signed ordinary parity")
                    raw = 3*raw+1 if digit else raw//2
                need(raw == root, "literal certified-home replay")
                a0, b0, q = normalized(x)
                phi_phase = power(PHI, -L)
                tail_numerator = (F(ROOTS[root][0][0]), F(ROOTS[root][0][1]))
                expected = mul(phi_phase, tail_numerator)
                need((a0 % q, b0 % q) == tuple(int(v) % q for v in expected), "clock predicts numerator phase")
                old_n, old_x = n, x
                checks += 1
            if word in ((), (1,), (1, 2)):
                examples.append({k: (list(v) if isinstance(v, tuple) else v) for k, v in row.items()})
    fake = F(1, 11), F(4, 11)
    need(ideal(fake) == ideal(roots[-5]), "same-ideal rational hostile")
    need(phase_orbit(fake).isdisjoint(phase_orbit(roots[-5])), "phase separates rational hostile")
    return dict(head_lengths="0..3", letters="1..3", heads=len(heads), roots=list(ROOTS),
                lifts="0..4", instances=checks, examples=examples,
                sibling_phase_periods=phase_periods,
                hostile="1/13 and -5 have the same denominator ideal, distinct phase orbits")


def main():
    inherited = runpy.run_path(str(ROOT / "04-computation/experiments/collatz_golden_carriers_20261004.py"))
    moduli = list(range(1, 41))+[64, 76, 81, 105]
    censuses = [lattice_census(q, inherited) for q in moduli]
    rejected_api_inputs = []
    for prime, q in ((2,0),(2,-1),(2,1),(2,6),(5,5)):
        try:
            lift_phase_children(prime, q, ((0,1),))
        except ValueError:
            rejected_api_inputs.append((prime,q))
        else:
            raise ValueError("unsupported child-selector input accepted")
    half, phi_half = (F(1,2),F(0)), (F(0),F(1,2))
    need(phase_orbit(half) == phase_orbit(phi_half), "same nonzero cycle for addition hostile")
    need(normalized(add(half,half))[2] == 1 and normalized(add(half,phi_half))[2] == 2,
         "cycle quotient loses well-defined addition")
    q105_classes = Counter()
    for cycle in phase_cycles(105):
        kind = "eigenline_at5" if (cycle[0][1]-3*cycle[0][0]) % 5 == 0 else "unit_at5"
        expected = 16 if kind == "eigenline_at5" else 80
        need(len(cycle) == expected, "105 CRT period")
        for a, b in cycle:
            need(((b-3*a) % 5 == 0) == (kind == "eigenline_at5"), "105 local support invariant")
            need((gcd(a*a+a*b-b*b,105) == 1) == (kind == "unit_at5"), "105 unit distinction")
        q105_classes[(kind, len(cycle))] += 1
    need(q105_classes == Counter({("eigenline_at5",16):96,("unit_at5",80):96}), "105 cycle split")
    primorials, product_q = [], 1
    for prime in (2,3,5,7,11,13):
        product_q *= prime
        primorials.append(dict(last_prime=prime, q=product_q,
                               primitive_phase_density=str(F(jordan2(product_q),product_q**2))))
    q4_odd_roots = []
    for cycle in inherited["cycles_in"](inherited["trap"](4), 4):
        orbit = inherited["decode_cycle"](cycle, 4)["orbit"]
        q4_odd_roots.append(sorted(str(n) for n in orbit if n.numerator % 2))
    need(sorted(q4_odd_roots) == [["1/29"], ["11/7", "5/7"]], "integer-realization loss on first binary lift")
    squares11 = {x*x % 11 for x in range(1, 11)}
    root11 = (F(-1, 11), F(7, 11))
    fake11 = (F(1, 11), F(4, 11))
    need({(a+4*b) % 11 for a, b in phase_orbit(root11)} == squares11, "minus5 phase squares")
    need({(a+4*b) % 11 for a, b in phase_orbit(fake11)} == set(range(1,11))-squares11, "rational hostile nonsquares")
    towers = []
    for prime, last in ((2, 6), (3, 4)):
        for exponent in range(1, last+1):
            q = prime**exponent
            cs = phase_cycles(q)
            expected_period = (3*2**(exponent-1) if prime == 2 else 8*3**(exponent-1))
            need(len(cs) == prime**(exponent-1), "tower cycle count")
            need(all(len(c) == expected_period for c in cs), "tower exact periods")
            if exponent < last:
                children = [lift_phase_children(prime, q, c) for c in cs]
                parent_for_child = [frozenset((a % q, b % q) for a, b in c) for c in phase_cycles(prime*q)]
                need(Counter(parent_for_child) == Counter({frozenset(c):prime for c in cs}), "each parent has prime children")
            else:
                children = []
            towers.append(dict(prime=prime, exponent=exponent, modulus=q,
                               cycles=len(cs), period=expected_period,
                               parent_child_checks=len(children), first_child_example=children[:1]))
    # The overlap list in the uniqueness proof is checked independently.
    overlaps = []
    for a, b in product(range(-4, 5), repeat=2):
        x, xp = (F(a), F(b)), (F(a+b), F(-b))
        if sign(add(ONE, x)) >= 0 and sign(add(ONE, neg(x))) >= 0 and \
           sign(add(power(PHI, 2), xp)) >= 0 and sign(add(power(PHI, 2), neg(xp))) >= 0:
            overlaps.append((a,b))
    need(set(overlaps) == {(0,0),(1,0),(-1,0),(-1,1),(1,-1),(2,-1),(-2,1)}, "complete overlap window")
    report = dict(status="PROVED all-q periodic lift and 2/3 towers; FINITE-EXACT controls; signed Collatz coverage OPEN",
                  lattice_phase_universe=moduli, exact_censuses=censuses,
                  overlap_vectors=overlaps, towers=towers,
                  rejected_child_selector_inputs=rejected_api_inputs,
                  cycle_addition_hostile="mod2: 1 and phi share an orbit; 1+1=0 but 1+phi=phi^2 is nonzero",
                  modulus105=dict(primitive_phases=9216, period16_phases=1536, period80_phases=7680,
                                  period16_cycles=96, period80_cycles=96,
                                  primitive_density=str(F(9216,105**2)), unit_phases=7680),
                  primorial_density_controls=primorials,
                  integer_lift_hostile=dict(parent_q=2, parent_odd_root=1, child_q=4, child_odd_roots=q4_odd_roots),
                  ideal_phase_hostile=dict(q=11, minus5_phase_values=sorted(squares11),
                                           one_thirteenth_phase_values=sorted(set(range(1,11))-squares11)),
                  signed_splices=splice_controls())
    print(json.dumps(report, indent=2))
    print("PASS: all checks active under -O")


if __name__ == "__main__":
    main()
