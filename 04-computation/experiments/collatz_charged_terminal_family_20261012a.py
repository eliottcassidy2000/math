"""Typed fixed-head ROOT completion and a charged exit for the pending source75.

Production functions use algebra, modular powers and supplied words. The finite
control-only orbit discovery is explicit, capped, and never called by them.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json


def need(ok, why):
    if not ok:
        raise ValueError(why)


def natural(n):
    need(type(n) is int and n >= 0, "exact nonnegative integer")
    return n


def odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, "exact positive odd source")
    return n


def word(w):
    need(type(w) is tuple and all(type(a) is int and a > 0 for a in w),
         "tuple of exact positive valuation letters")
    return w


def carrier(w):
    word(w)
    p, q, b = 1, 1, 0
    for a in w:
        p, q, b = 3*p, (1 << a)*q, 3*b+q
    return p, q, b


def valuation(n):
    need(type(n) is int and n > 0, "positive valuation argument")
    return (n & -n).bit_length()-1


def replay(n, w):
    odd(n), word(w)
    for a in w:
        need(n != 1, "no ROOT padding")
        z = 3*n+1
        need(valuation(z) == a, "actual valuation")
        n = z >> a
    return n


def terminal_phase(w):
    """The unique k modulo3^len(w) allowing a final valuation2k."""
    p, q, b = carrier(w)
    phase, period, modulus = 0, 1, 3
    for _ in w:
        modulus *= 3
        target = (1+3*b*pow(q, -1, modulus)) % modulus
        choices = tuple(phase+j*period for j in range(3)
                        if pow(4, phase+j*period, modulus) == target)
        need(len(choices) == 1, "unique ternary phase lift")
        phase = choices[0]
        period *= 3
    need(period == p, "phase period")
    return phase, period


def first_exponent(w):
    """Conservative cut: k>=2 and C/(Q4^k)<=1/2, without large expansion."""
    p, q, b = carrier(w)
    c = q+3*b-p
    minimum, power = 2, 16*q
    while power < 2*c:
        minimum += 1
        power *= 4
    phase, period = terminal_phase(w)
    return phase + max(0, (minimum-phase+period-1)//period)*period


@dataclass(frozen=True)
class Plan:
    head: tuple
    index: int


def audit_plan(plan):
    need(type(plan) is Plan, "exact Plan")
    word(plan.head), natural(plan.index)


def exponent(plan):
    audit_plan(plan)
    return first_exponent(plan.head)+3**len(plan.head)*plan.index


def source_mod(plan, modulus):
    audit_plan(plan)
    need(type(modulus) is int and modulus > 0, "exact positive modulus")
    p, q, b = carrier(plan.head)
    denominator = 3*p
    numerator = (q*pow(4, exponent(plan), denominator*modulus)-q-3*b) % (denominator*modulus)
    need(numerator % denominator == 0, "native division precision")
    return numerator//denominator


def source(plan, bit_cap=16384):
    audit_plan(plan), natural(bit_cap)
    p, q, b = carrier(plan.head)
    k = exponent(plan)
    need(2*k+q.bit_length() <= bit_cap, "declared materialization cap")
    numerator = (q << (2*k))-q-3*b
    need(numerator > 0 and numerator % (3*p) == 0, "positive native inverse")
    return odd(numerator//(3*p))


def compile_word(plan):
    audit_plan(plan)
    return plan.head+(2*exponent(plan),)


def decode_source(w, n):
    word(w), odd(n)
    p, q, b = carrier(w)
    numerator = 3*p*n+q+3*b
    if numerator % q:
        return None
    z = numerator//q
    if z <= 0 or z & (z-1):
        return None
    a = z.bit_length()-1
    if a % 2:
        return None
    k, first = a//2, first_exponent(w)
    if k < first or (k-first) % p:
        return None
    return Plan(w, (k-first)//p)


def normalized_score(plan):
    # This exact rational reader intentionally respects the same expansion cap.
    source(plan)
    p, q, b = carrier(plan.head)
    k = exponent(plan)
    o = sum(a % 2 for a in plan.head)
    return Q(1, 1 << o)*(1-Q(q+3*b-p, q << (2*k)))


def fixed_score(plan):
    source(plan)
    p, q, _ = carrier(plan.head)
    o = sum(a % 2 for a in plan.head)
    return Q(4*p, (q << o) << (2*exponent(plan)))


def uniform_floor(w):
    word(w)
    return Q(1, 1 << (1+sum(a % 2 for a in w)))


def localization_depth(eta):
    need(type(eta) is Q and 0 < eta <= 1, "positive exact normalized-score floor")
    d, ratio = 1, 1
    while 2688*eta*ratio**8 < 9:
        d += 1
        ratio = (4**d-1)//3
    return d


def coarse75(n):
    """Return a paid actual first descent on75mod128, not an assumed ROOT word."""
    odd(n)
    if n % 128 != 75:
        return None
    a = valuation((27*n+23)//8)
    w = (1, 2, a)
    h = (27*n+23) >> (a+3)
    need(a >= 4 and 4*(h+1) <= n+1, "two-bit source payment")
    need(replay(n, w) == h and h < n, "authenticated first descent")
    return h, w, 1+a % 2


def coarse123(n):
    """The next frozen guard: a rate-three first descent on123mod512."""
    odd(n)
    if n % 512 != 123:
        return None
    a = valuation((243*n+319)//64)
    w = (1, 2, 1, 2, a)
    h = (243*n+319) >> (a+6)
    need(a >= 3 and 2*(h+1) <= n+1, "one-bit source payment")
    need(replay(n, w) == h and h < n, "authenticated next first descent")
    return h, w, 2+a % 2


def charge_deadline(n, budget):
    odd(n), natural(budget)
    if n == 1:
        return 0
    left, right, t = (n-1) << budget, 2, 0
    while left >= right:
        left *= 3
        right *= 4
        t += 1
    return budget+t


def band_invoice(n, rate=3):
    """Total finite test of one source-owned quadratic-log band invoice.

    A false result is a failed invoice, not a nonconvergence claim.
    A true result supplies a strict smaller child, not its ROOT proof.
    """
    odd(n)
    need(n > 1 and type(rate) is int and rate > 0, "nonroot source and positive integer rate")
    e = (n+1).bit_length()-1
    budget = rate*(e*e-1)
    deadline = charge_deadline(n, budget)
    h, charge, letters = n, 0, []
    while 2*(h+1) > n+1 and charge <= budget and len(letters) < deadline:
        a = valuation(3*h+1)
        letters.append(a)
        h = (3*h+1) >> a
        charge += a % 2
    if 2*(h+1) <= n+1:
        f = (h+1).bit_length()-1
        return charge <= rate*(e*e-f*f), h, tuple(letters), deadline
    need(charge > budget, "proved ROOT-or-budget deadline")
    return False, h, tuple(letters), deadline


def compositions(length, budget, prefix=()):
    if length == 0:
        yield prefix
    else:
        for a in range(1, budget-length+2):
            yield from compositions(length-1, budget-a, prefix+(a,))


def finite_grammar(limit, include75, include123=False):
    """Declared depth8 version of inherited grammar, same15 supplied seeds.

    This returns charge proof-ledger metadata, not new orbit discoveries.
    Recognition never calls the control-only orbit routine in main().
    """
    natural(limit)
    need(type(include75) is bool and type(include123) is bool, "exact policy flags")
    data_path = Path(__file__).resolve().parents[2]/'05-knowledge/results/charge_ports_eleven_20261011.json'
    seed_data = json.loads(data_path.read_text(encoding='utf-8'))['controllers']['bank_words']
    bank = {}
    for key, letters in seed_data.items():
        n, w = int(key), tuple(letters)
        need(replay(n, w) == 1, "supplied strict seed word")
        bank[n] = sum(a % 2 for a in w)
    need(set(bank) == set(range(1, 30, 2)), "frozen seed universe")
    tables = {}
    templates = 0
    for length in range(1, 9):
        p = 3**length
        for w in compositions(length, (p-1).bit_length()-1):
            _, q, b = carrier(w)
            k = 0
            while q*4**k < p:
                c, d = q*4**k, q*(4**k-1)//3-b
                residue = (-d*pow(c, -1, p)) % p
                tables.setdefault((length, residue), []).append((c, d, -sum(a % 2 for a in w)))
                templates += 1
                k += 1
    need(templates == 1045, "complete depth8 template count")
    ledger = {n: (None, None, 'seed', charge) for n, charge in bank.items()}
    for n in range(31, limit+1, 2):
        choices = []
        if n % 4 == 1:
            a = valuation(3*n+1)
            choices.append(((3*n+1) >> a, a % 2, 'forward'))
        if n % 8 == 5:
            choices.append(((n-1)//4, 0, 'sibling'))
        k = valuation(n+1)
        if k >= 2:
            t = (n+1) >> k
            r = valuation(3**k*t-1)
            if r >= 2:
                choices.append(((n-1)//2, 1, 'reset'))
            else:
                y = (3**k*t-1)//2
                a = valuation(3*y+1)
                y = (3*y+1) >> a
                if a == 6:
                    choices.append(((n-1)//2, -1, 'sporadic26'))
                if a == 2 and valuation(3*y+1) == 10:
                    choices.append(((n-1)//2, -1, 'sporadic2210'))
        for length in range(1, 9):
            p = 3**length
            for c, d, charge in tables.get((length, n % p), ()):
                h = (c*n+d)//p
                if 0 < h < n:
                    choices.append((h, charge, 'inverse'))
        if n % 1024 == 43:
            choices.append(((81*n+101)//512, 1, '43-cylinder'))
        if n % 1024 == 39:
            choices.append(((243*n+251)//512, 3, '39-cylinder'))
        if include75 and n % 128 == 75:
            h, _, charge = coarse75(n)
            choices.append((h, charge, '75-coarse'))
        if include123 and n % 512 == 123:
            h, _, charge = coarse123(n)
            choices.append((h, charge, '123-coarse'))
        for h, delta, _ in choices:
            need(type(h) is int and 0 < h < n and h % 2,
                 "strictly smaller actual dependency")
            need(delta <= 3*((n+1).bit_length()-(h+1).bit_length()), "inherited rank invoice")
        grounded = [r for r in choices if r[0] in ledger]
        if grounded:
            h, delta, name = min(grounded, key=lambda r:(ledger[r[0]][3]+r[1], r))
            ledger[n] = (h, delta, name, ledger[h][3]+delta)
    return ledger


def main():
    checks = 0

    def check(ok, why):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(why)

    def reject(fn):
        try:
            fn()
        except ValueError:
            check(True, 'rejected')
        else:
            check(False, 'invalid typed input accepted')

    heads = tuple(w for length in range(5) for w in product(range(1, 5), repeat=length))
    for w in heads:
        p, q, b = carrier(w)
        c = q+3*b-p
        check(c >= 0 and (c == 0) == (w == ()), 'positive carry recurrence')
        for a in range(1, 5):
            p2, q2, b2 = carrier(w+(a,))
            check(q2+3*b2-p2 == 3*c+(1 << a)*q, 'carry recursion')
        phase, period = terminal_phase(w)
        check(period == p and pow(4, phase, 3*p) == (1+3*b*pow(q, -1, 3*p)) % (3*p),
              'exact terminal phase')
        for index in range(3):
            plan = Plan(w, index)
            n, root_word = source(plan), compile_word(plan)
            check(replay(n, root_word) == 1, 'strict ROOT completion')
            check(decode_source(w, n) == plan, 'actual-source membership round trip')
            charge = sum(a % 2 for a in root_word)
            K = sum((a-1)//2 for a in root_word)
            f = Q(3, 4)**len(w)/4**K
            check(charge == sum(a % 2 for a in w), 'terminal has zero odd charge')
            check(f == fixed_score(plan), 'independent exact score')
            check(Q(3*n+1, 4)*f == normalized_score(plan) >= uniform_floor(w), 'normalized formula/floor')
            for modulus in (1, 2, 3, 9, 19, 128, 729):
                check(source_mod(plan, modulus) == n % modulus, 'full-denominator modular reader')
        check(decode_source(w, source(Plan(w, 0))+2) is None, 'adjacent-source guard')
        d = localization_depth(uniform_floor(w))
        check(2688*uniform_floor(w)*((4**d-1)//3)**8 >= 9, 'uniform localization depth')
        if d > 1:
            check(2688*uniform_floor(w)*((4**(d-1)-1)//3)**8 < 9, 'least retained threshold depth')
    eta = Q(1017, 2048)
    check(terminal_phase((1, 2)) == (4, 9), 'pending75 phase')
    check(normalized_score(Plan((1, 2), 0)) == eta, 'sharp family minimum')
    for index in range(17):
        plan = Plan((1, 2), index)
        n, k = source(plan), exponent(plan)
        check(k == 4+9*index and n == ((1 << (18*index+11))-23)//27, 'one-charge formula')
        check(coarse75(n)[0] == 1, 'family grounds the coarse75 exit')
        f = fixed_score(plan)
        m, q0 = (n-1)//2, (3*n+1)//2
        nu = Q(1, (m+1)**6)-Q(1, (m+2)**6)
        for d in (1, 2, 4):
            qd = q0*(4**d-1)//3
            alias = nu+Q(1, 4**d)*(Q(1, (m+qd+1)**6)-Q(1, (m+qd+2)**6))
            error = Q(9, 448*qd**8)
            check((alias*f-error)/alias > Q(2365, 2373)*f, 'all-depth retained floor control')
    for t in range(257):
        n = 75+128*t
        h, w, delta = coarse75(n)
        check(replay(n, w) == h and 4*(h+1) <= n+1, 'all-height paid guard samples')
        check(delta <= ((n+1).bit_length()-(h+1).bit_length()), 'rate-one source invoice')
        n = 123+512*t
        h, w, delta = coarse123(n)
        check(replay(n, w) == h and 2*(h+1) <= n+1, 'next all-height paid guard samples')
        check(delta <= 3*((n+1).bit_length()-(h+1).bit_length()), 'rate-three next invoice')
    huge = Plan((1, 2), 10**60)
    check(source_mod(huge, 128) == 75, 'symbolic huge index without expansion')
    reject(lambda: source(huge))

    def control_orbit(n, cap=1000):
        # Explicit validation/discovery only, outside every production API.
        w = []
        while n != 1 and len(w) < cap:
            a = valuation(3*n+1)
            w.append(a)
            n = (3*n+1) >> a
        check(n == 1, 'declared finite replay cap completed')
        return tuple(w)

    old, new = finite_grammar(8191, False), finite_grammar(8191, True)
    third = finite_grammar(8191, True, True)
    check(set(old) <= set(new) and 75 not in old and 75 in new, 'same-seed proper finite gain')
    check(set(new) <= set(third) and 123 not in new and 123 in third, 'separate next-rule policy gain')
    for ledger in (old, new, third):
        for n, row in ledger.items():
            o = sum(a % 2 for a in control_orbit(n))
            check(o == row[3] <= 3*((n+1).bit_length()-1)+16, 'independent compiled charge replay')
    failures = [n for n in range(1, 328, 2)
                if sum(a % 2 for a in control_orbit(n)) > ((n+1).bit_length()-1)+24]
    check(failures == [327], 'frozen universal charge-potential hostile')
    root63 = control_orbit(63)
    check(replay(63, root63) == 1, 'literal63 separate from quadratic and order claims')

    # Cheap typed bridges: inherited dynamics identities, not a conjugacy to U.
    cycle = (Q(-7, 4), Q(5, 4), Q(-1, 4))
    multiplier = Q(1)
    for i, x in enumerate(cycle):
        check(x*x-Q(29, 16) == cycle[(i+1) % 3], 'rational quadratic3 cycle')
        multiplier *= 2*x
    check(multiplier == Q(35, 8), 'quadratic cycle derivative product')
    x, denominators = Q(0), []
    for _ in range(5):
        x = x*x-Q(7, 4)
        denominators.append(x.denominator)
    check(denominators == [4, 16, 256, 65536, 4294967296], 'nonintegral critical denominator growth')
    for modulus in (9, 21, 63):
        check(pow(2, 6, modulus) == 1 and all(pow(2, d, modulus) != 1 for d in (1, 2, 3)),
              'finite multiplicative order, not Collatz clock')
    c1, c2 = carrier((1, 1, 1, 3)), carrier((1, 1, 2, 2))
    check(c1 == (81, 64, 65) and c2 == (81, 64, 73), 'same slopes lose ordered carry')
    check(Q(c1[2], c1[1]-c1[0]) == Q(-65, 17)
          and Q(c2[2], c2[1]-c2[0]) == Q(-73, 17), 'different affine fixed points')

    # Control-only discovery for a falsifiable band-exit invoice, no production use.
    max_ratio, band_witness, failures1, max_steps, early_max = Q(0), None, [], 0, 0
    for n in range(3, 1000002, 2):
        h, charge, steps = n, 0, 0
        while 2*(h+1) > n+1 and steps < 10000:
            a = valuation(3*h+1)
            h = (3*h+1) >> a
            charge += a % 2
            steps += 1
        check(steps < 10000 and 2*(h+1) <= n+1, 'declared finite band-exit cap')
        e, f = (n+1).bit_length()-1, (h+1).bit_length()-1
        ratio = Q(charge, e*e-f*f)
        if ratio > max_ratio:
            max_ratio, band_witness = ratio, (n, h, charge, e, f, steps)
        if ratio > 1 and n <= 100001:
            failures1.append(n)
        max_steps = max(max_steps, steps)
        if n <= 100001:
            early_max = max(early_max, steps)
    check(max_ratio == Q(26, 9) and band_witness == (55, 23, 26, 5, 4, 37), 'sharp finite band invoice')
    check(len(failures1) == 136 and failures1[0] == 27 and early_max == 96 and max_steps == 135,
          'finite C1 hostile and extended C3 scope')
    for n in range(3, 256, 2):
        ok, h, w, deadline = band_invoice(n)
        check(ok and replay(n, w) == h and len(w) <= deadline, 'total source-owned band receipt')
    check(not band_invoice(55, 2)[0] and band_invoice(55, 3)[0], 'band coefficient boundary')
    check(not band_invoice(27, 1)[0], 'source-owned rejection without convergence assumption')

    from dataclasses import replace
    good = Plan((1, 2), 0)
    exponent(good)  # warm the legitimate path before typed aliases.
    for fn in (lambda: exponent(replace(good, index=False)), lambda: exponent(Plan((1, 2.0), 0)),
               lambda: carrier([1, 2]), lambda: source_mod(good, True),
               lambda: decode_source((1, 2), 75.0), lambda: coarse75(True),
               lambda: replay(1, (2,)), lambda: localization_depth(Q(0)),
               lambda: localization_depth(0.5), lambda: finite_grammar(31, 1),
               lambda: band_invoice(1), lambda: band_invoice(55, True),
               lambda: charge_deadline(3, 1.0), lambda: finite_grammar(31, True, 1),
               lambda: coarse123(123.0)):
        reject(fn)
    gains = sorted(set(new)-set(old))
    pending = next(n for n in range(1, 8192, 2) if n not in new)
    print('Charged fixed-head ROOT completion; no all-source coverage assertion')
    print(f'Head universe: {len(heads)} words, alphabet1..4, length0..4; three completions each')
    print('Generic terminal phase, exact source decoder, symbolic modular reader: PASS')
    print('75 family: k=4+9j; n=(2^(18j+11)-23)/27; ROOT word(1,2,18j+8)')
    print('Odd charge1; minimum Z=1017/2048; every native refinement retains >2365/2373 of f')
    print('New guarded dependency: n75mod128; exact(1,2,a),a>=4;4(h+1)<=n+1')
    print(f'Same15seeds, depth8+39/43 baseline, odd1..8191: {len(old)} -> {len(new)}, gain{len(gains)}')
    print('First newly grounded requests:', gains[:16])
    print('First still pending:', pending)
    print(f'Additional123mod512 policy: {len(new)} -> {len(third)}, gain{len(third)-len(new)}; first pending{next(n for n in range(1,8192,2) if n not in third)}')
    print('Least failure of proposed universal o<=floor(log2(n+1))+24:327, charge35')
    print(f'Dynamics roles: quadratic cycle multiplier35/8; critical denominator growth; Collatz63 ROOT rank{len(root63)}')
    print('Band-exit control odd3..1000001, cap10000: max ratio26/9 at55->23; C3 passes; max steps135')
    print('Initial band universe odd3..100001: C1 failures136, least27; max steps96')
    print('Source-owned band invoice: total ROOT-or-charge deadline; no assumed child ROOT')
    print(f'checks={checks}')


if __name__ == '__main__':
    main()
