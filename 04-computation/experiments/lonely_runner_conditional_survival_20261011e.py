"""Exact common-time survival, conditioned residue profiles, and closed walls.

Finite source-bound certificates only. No product estimate proves LRC here.
Uses Fraction arithmetic, and all checks remain active under python -O.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd, lcm

CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def integer(n, lower=0):
    need(type(n) is int and n >= lower, 'exact integer in range')


def domain(speeds, delta):
    need(type(speeds) is tuple and all(type(v) is int and v > 0 for v in speeds)
         and len(set(speeds)) == len(speeds), 'distinct positive integer speeds')
    need(type(delta) is F and 0 < delta < F(1, 2), 'exact radius in (0,1/2)')


def distance(x):
    need(type(x) is F, 'exact rational phase')
    r = x % 1
    return min(r, 1-r)


def safe(speeds, delta, t):
    domain(speeds, delta)
    need(type(t) is F and 0 <= t <= 1, 'rational time in the closed period')
    return all(distance(v*t) >= delta for v in speeds)


def merged(intervals):
    answer = []
    for a, b in sorted(intervals):
        if answer and a <= answer[-1][1]:
            answer[-1] = (answer[-1][0], max(b, answer[-1][1]))
        else:
            answer.append((a, b))
    return tuple(answer)


def intersection(left, right):
    answer = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
        if a <= b:
            answer.append((a, b))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return merged(answer)


@dataclass(frozen=True)
class Survival:
    speeds: tuple
    delta: F
    intervals: tuple


def survivor(speeds, delta):
    domain(speeds, delta)
    intervals = ((F(0), F(1)),)
    for v in speeds:
        allowed = tuple(((j+delta)/v, (j+1-delta)/v) for j in range(v))
        intervals = intersection(intervals, allowed)
    return Survival(speeds, delta, intervals)


def audit_survivor(packet):
    need(type(packet) is Survival, 'typed survivor certificate')
    domain(packet.speeds, packet.delta)
    need(type(packet.intervals) is tuple and all(
        type(ab) is tuple and len(ab) == 2 and all(type(x) is F for x in ab)
        and 0 <= ab[0] <= ab[1] <= 1 for ab in packet.intervals), 'typed closed components')
    need(packet == survivor(packet.speeds, packet.delta), 'complete source-bound survivor set')
    return packet


def mass(packet):
    audit_survivor(packet)
    return sum((b-a for a, b in packet.intervals), F(0))


def witness(packet):
    audit_survivor(packet)
    if not packet.intervals:
        return None
    a, b = packet.intervals[0]
    result = (a+b)/2
    need(safe(packet.speeds, packet.delta, result), 'authenticated source witness')
    return result


def atom_engine(speeds, delta):
    """Independent complete endpoint/midpoint engine, retaining isolated walls."""
    domain(speeds, delta)
    endpoints = {F(0), F(1)}
    for v in speeds:
        endpoints.update((j+delta)/v for j in range(v))
        endpoints.update((j+1-delta)/v for j in range(v))
    points = tuple(sorted(endpoints))
    pieces = [(x, x) for x in points if safe(speeds, delta, x)]
    pieces += [(a, b) for a, b in zip(points, points[1:])
               if safe(speeds, delta, (a+b)/2)]
    return merged(pieces)


def survival_chain(speeds, delta):
    domain(speeds, delta)
    prior = F(1)
    rows = []
    for i in range(1, len(speeds)+1):
        packet = survivor(speeds[:i], delta)
        current = mass(packet)
        rows.append((current, current/prior if prior else None, witness(packet)))
        prior = current
    return tuple(rows)


def primitive(x, delta):
    """Periodic primitive of 1_[delta,1-delta] - (1-2delta)."""
    domain((), delta)
    need(type(x) is F, 'exact primitive argument')
    z, s = x % 1, 1-2*delta
    if z <= delta:
        return -s*z
    if z <= 1-delta:
        return -s*delta + 2*delta*(z-delta)
    return s*(1-z)


def append_mass(packet, speed):
    """Exact positioned boundary bill and conservative component lower bound."""
    audit_survivor(packet)
    integer(speed, 1)
    need(speed not in packet.speeds, 'new source speed')
    mu, count = mass(packet), sum(a < b for a, b in packet.intervals)
    s = 1-2*packet.delta
    bill = sum((primitive(speed*b, packet.delta)-primitive(speed*a, packet.delta)
                for a, b in packet.intervals), F(0))/speed
    lower = max(F(0), s*mu - count*2*packet.delta*s/speed)
    return s*mu+bill, bill, lower


def pair_danger(a, b, delta):
    integer(a, 1)
    integer(b, 1)
    domain((), delta)
    if a > b:
        a, b = b, a
    g = gcd(a, b)
    width, step, cap = delta*(a+b)/(a*b), F(g, a*b), 2*delta/b
    stop = int(width/step)
    return g*sum((min(max(F(0), width-abs(j)*step), cap)
                  for j in range(-stop, stop+1)), F(0))


def residue_set(modulus, values):
    integer(modulus, 1)
    need(type(values) is frozenset and all(type(x) is int and 0 <= x < modulus
         for x in values), 'exact canonical residue set')


@dataclass(frozen=True)
class Lift:
    modulus: int
    residues: frozenset
    shared_modulus: int
    left_profile: tuple
    right_profile: tuple
    conditional: object


def conditional_lift(Q, A, q, B):
    residue_set(Q, A)
    residue_set(q, B)
    d, L = gcd(Q, q), lcm(Q, q)
    left = tuple(sum(a % d == r for a in A) for r in range(d))
    right = tuple(sum(b % d == r for b in B) for r in range(d))
    results = set()
    for a in A:
        for b in B:
            if (b-a) % d:
                continue
            quotient = q//d
            shift = 0 if quotient == 1 else ((b-a)//d)*pow(Q//d, -1, quotient) % quotient
            results.add(a+Q*shift)
    size = sum(a*b for a, b in zip(left, right))
    need(len(results) == size, 'CRT fibre count')
    factor = F(size, len(A)*(q//d)) if A else None
    return Lift(L, frozenset(results), d, left, right, factor)


def clock_survival(speeds, clock, delta):
    """Residue k means actual source time k/clock; it is never k/period."""
    domain(speeds, delta)
    integer(clock, 1)
    period, alive = 1, frozenset({0})
    rows = []
    for v in speeds:
        g = gcd(clock, v)
        q, u = clock//g, v//g
        allowed = frozenset(k for k in range(q) if distance(F(u*k, q)) >= delta)
        step = conditional_lift(period, alive, q, allowed)
        rows.append(step)
        period, alive = step.modulus, step.residues
    return period, alive, tuple(rows)


def complete_clock(speeds, delta):
    """Contains every source boundary; finite but not a uniform small clock."""
    domain(speeds, delta)
    return delta.denominator*lcm(*speeds) if speeds else 1


def rejects(fn):
    try:
        fn()
    except (ValueError, TypeError):
        check(True, 'malformed or forged input rejected')
    else:
        check(False, 'malformed or forged input accepted')


def main():
    cases = 0
    radii = (F(1, 6), F(1, 7), F(1, 12), F(1, 14))
    for r in range(4):
        for speeds in combinations(range(1, 11), r):
            for delta in radii:
                packet = survivor(speeds, delta)
                check(packet.intervals == atom_engine(speeds, delta), 'two independent continuous engines')
                t = witness(packet)
                check(t is None or safe(speeds, delta, t), 'actual witness')
                product = F(1)
                for mu, rho, _ in survival_chain(speeds, delta):
                    if rho is not None:
                        product *= rho
                    check(product == mu, 'telescoping actual conditional survival')
                for v in (11, 17):
                    exact, bill, lower = append_mass(packet, v)
                    check(exact == mass(survivor(speeds+(v,), delta)), 'positioned primitive identity')
                    check(exact >= lower, 'component discrepancy sufficient bound')
                cases += 1
    for a in range(1, 31):
        for b in range(a+1, 32):
            for delta in (F(1, 12), F(1, 14)):
                check(mass(survivor((a, b), delta)) == 1-4*delta+pair_danger(a, b, delta),
                      'all-wrap pair formula against direct intervals')
    check(pair_danger(50, 51, F(1, 14)) == F(73, 3570), 'corrected inherited pair hostile')
    check(pair_danger(13, 15, F(1, 14)) == F(1, 49)-F(6, 49*13*15), 'negative pair covariance')
    check(mass(survivor((1, 2), F(1, 14))) > F(6, 7)**2, 'positive pair covariance')
    ap = []
    for k in range(3, 21):
        speeds, delta = tuple(range(1, k)), F(1, k)
        packet = survivor(speeds, delta)
        expected = tuple((F(a, k), F(a, k)) for a in range(1, k) if gcd(a, k) == 1)
        check(packet.intervals == expected and mass(packet) == 0, 'tight AP has only unit endpoint witnesses')
        ap.append((k, len(expected)))
        for eps in (delta/10, delta/100):
            check(mass(survivor(speeds, delta-eps)) >= 2*eps/max(speeds), 'Lipschitz closed-wall recovery')
    for speeds in permutations((1, 2, 5)):
        chain = survival_chain(speeds, F(1, 7))
        check(chain[-1][0] == mass(survivor((1, 2, 5), F(1, 7))), 'order changes factors not final event')
    complete_cases = 0
    for r in range(1, 4):
        for speeds in combinations(range(1, 7), r):
            for delta in (F(1, 7), F(1, 3), F(2, 5)):
                clock = complete_clock(speeds, delta)
                _, alive, _ = clock_survival(speeds, clock, delta)
                check(bool(alive) == bool(survivor(speeds, delta).intervals), 'complete source-clock weak equivalence')
                if alive:
                    check(safe(speeds, delta, F(min(alive), clock)), 'complete clock actual witness')
                complete_cases += 1
    lifts = 0
    for Q in range(1, 6):
        for q in range(1, 6):
            for left in range(1 << Q):
                A = frozenset(i for i in range(Q) if left >> i & 1)
                for right in range(1 << q):
                    B = frozenset(i for i in range(q) if right >> i & 1)
                    result = conditional_lift(Q, A, q, B)
                    brute = frozenset(i for i in range(lcm(Q, q)) if i % Q in A and i % q in B)
                    check(result.residues == brute, 'complete conditional CRT versus literal residues')
                    if gcd(Q, q) == 1 and A:
                        check(result.conditional == F(len(B), q), 'coprime native independence')
                    lifts += 1
    clock_cases = 0
    for speeds in ((), (1,), (1, 2), (2, 3, 5), (1, 3, 4), tuple(range(1, 14))):
        for clock in range(1, 49):
            period, alive, rows = clock_survival(speeds, clock, F(1, 14))
            brute = frozenset(k for k in range(clock) if safe(speeds, F(1, 14), F(k, clock)))
            lifted = frozenset(k for k in range(clock) if k % period in alive)
            check(brute == lifted, 'native common-time modular source retained')
            prob = F(1)
            for row in rows:
                if row.conditional is not None:
                    prob *= row.conditional
            check(prob == F(len(brute), clock), 'native finite survival product')
            clock_cases += 1
    coprime = []
    for moduli in ((2, 3, 5), (5, 7, 11), (2, 17, 19)):
        N = lcm(*moduli)
        speeds = tuple(N//q for q in moduli)
        period, alive, _ = clock_survival(speeds, N, F(1, 14))
        expected = F(1)
        for q in moduli:
            count = sum(distance(F(a, q)) >= F(1, 14) for a in range(q))
            expected *= F(count, q)
        check(F(len(alive), period) == expected and alive, 'cofactor family actual product and witness')
        coprime.append((moduli, expected, F(min(alive), N)))
    check(F(len(clock_survival((5, 3), 30, F(1, 14))[1]), 30) == F(23, 30)
          and F(23, 30) != F(5, 6)*F(9, 10), 'native shared-factor product hostile')
    for k, expected in ((12, F(1, 3)), (14, F(3, 7))):
        check(F(len(clock_survival(tuple(range(1, k)), k, F(1, k))[1]), k) == expected,
              'boundary counting measure remains positive')
    # Equal marginal counts can give survival zero or one after conditioning.
    yes = conditional_lift(6, frozenset({0, 2, 4}), 2, frozenset({0}))
    no = conditional_lift(6, frozenset({1, 3, 5}), 2, frozenset({0}))
    check(yes.conditional == 1 and no.conditional == 0, 'same counts, opposite conditional profile')
    # A valid continuous witness can disappear on a chosen clock grid.
    packet = survivor((1, 2), F(1, 3))
    check(witness(packet) == F(1, 3), 'closed witness exists')
    check(not clock_survival((1, 2), 2, F(1, 3))[1], 'failed grid is not failed LRC')
    period, alive, _ = clock_survival((2,), 4, F(1, 14))
    check(period == 2 and alive == frozenset({1}) and safe((2,), F(1, 14), F(1, 4))
          and not safe((2,), F(1, 14), F(1, period)), 'residue period cannot replace source clock')
    rejects(lambda: survivor((True,), F(1, 14)))
    rejects(lambda: survivor((1, 1), F(1, 14)))
    rejects(lambda: survivor((1,), 1/14))
    rejects(lambda: audit_survivor(replace(packet, intervals=())))
    rejects(lambda: audit_survivor(replace(packet, intervals=((F(1, 3), 0.333333),))))
    rejects(lambda: conditional_lift(True, frozenset({0}), 2, frozenset({1})))
    rejects(lambda: conditional_lift(2, frozenset({False}), 2, frozenset({1})))
    rejects(lambda: conditional_lift(2, frozenset({2}), 2, frozenset({1})))
    rejects(lambda: clock_survival((1,), 0, F(1, 14)))
    print('Lonely runner: conditional survival, residue profiles, and closed walls')
    print('Independent continuous cases:', cases, '(all subsets of1..10 of size0..3; four radii)')
    print('Two-speed all-wrap controls:930 plus corrected50/51 and signed-covariance hostiles')
    print('Tight AP (total runners, isolated safe points):', ap)
    print('AP12 actual conditional mass factors:', tuple(row[1] for row in survival_chain(tuple(range(1, 12)), F(1, 12))))
    print('AP14 actual conditional mass factors:', tuple(row[1] for row in survival_chain(tuple(range(1, 14)), F(1, 14))))
    print('Complete source-clock equivalence cases:', complete_cases, '(subsets1..6,size1..3; radii1/7,1/3,2/5)')
    print('Complete conditional CRT set pairs:', lifts, '(both moduli1..5, all subsets)')
    print('Native common-clock cases:', clock_cases, '(six speed tuples, clocks1..48)')
    print('Coprime cofactor families (moduli,exact survival,witness):', coprime)
    print('Same marginal profile hostile: conditional survival1 versus0; bad grid retains no global verdict')
    print('Always-active exact checks:', CHECKS)
    print('Scope: finite source certificates and conditional theorems; no product proof of LRC14')


if __name__ == '__main__':
    main()
