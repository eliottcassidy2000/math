"""Exact numeric, prime-predecessor tube, and source-preserving join controls.

The accepted Poisson--Dirichlet and Partition Principle papers are premises,
not tested here. All executable claims are elementary finite identities.
"""
from fractions import Fraction as F
from itertools import product


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(n, minimum=0):
    need(type(n) is int and n >= minimum, "exact integer in domain required")
    return n


def valuation(n):
    integer(n, 1)
    return (n & -n).bit_length() - 1


def checked_word(word):
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word),
         "tuple of exact positive valuations required")
    return word


def step(n):
    integer(n, 1)
    need(n % 2 == 1, "odd source required")
    a = valuation(3 * n + 1)
    return (3 * n + 1) >> a, a


def replay(n, word, root=False):
    integer(n, 1)
    need(n % 2 == 1, "odd source required")
    checked_word(word)
    for a in word:
        if root:
            need(n != 1, "earlier ROOT forbidden")
        n, actual = step(n)
        need(actual == a, "wrong actual valuation")
    if root:
        need(n == 1, "ROOT endpoint required")
    return n


def carrier(word):
    checked_word(word)
    p, q, b = 1, 1, 0
    for a in word:
        p, q, b = 3 * p, q * (1 << a), 3 * b + q
    return p, q, b


def common_fiber(entries):
    """Return (H, endpoint period, (source, word, source period) rows)."""
    need(type(entries) is tuple and len(entries) >= 1, "nonempty tuple of entries")
    rows = []
    endpoints = []
    for entry in entries:
        need(type(entry) is tuple and len(entry) == 2, "source/word pair")
        source, word = entry
        endpoints.append(replay(source, word))
        p, q, b = carrier(word)
        rows.append((source, word, q, len(word)))
    need(len(set(endpoints)) == 1, "common endpoint required")
    r = max(row[3] for row in rows)
    return (endpoints[0], 2 * 3 ** r,
            tuple((n, w, 2 * q * 3 ** (r - length)) for n, w, q, length in rows))


def select_source(entries, index, source):
    """Recognize a supplied source; never replace it by another residue lift."""
    h, period, rows = common_fiber(entries)
    integer(index)
    need(index < len(rows), "row index")
    integer(source, 1)
    base, word, modulus = rows[index]
    difference = source - base
    need(difference >= 0 and difference % modulus == 0, "source outside forward fiber")
    t = difference // modulus
    selected = tuple(n + modulus * t for n, _, modulus in rows)
    endpoint = h + period * t
    for n, (_, word, _) in zip(selected, rows):
        need(replay(n, word) == endpoint, "fiber receipt mismatch")
    return t, selected, endpoint


def transport_root(source, source_word, child, child_word, child_root_word):
    """Conditional ROOT splice with an authenticated supplied child proof."""
    replay(child, child_root_word, root=True)
    checked_word(child_word)
    need(child_root_word[:len(child_word)] == child_word, "child proof lacks join prefix")
    need(replay(source, source_word) == replay(child, child_word), "wrong join endpoint")
    result = source_word + child_root_word[len(child_word):]
    replay(source, result, root=True)
    return result


def root_word_for_control(n, cap=10000):
    """Bounded literal search used only to supply independent finite controls."""
    integer(n, 1)
    integer(cap)
    need(n % 2 == 1, "odd source")
    result = []
    while n != 1 and len(result) < cap:
        n, a = step(n)
        result.append(a)
    need(n == 1, "control search exhausted")
    return tuple(result)


def prime(n):
    integer(n)
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def factors(n):
    integer(n, 1)
    result, d = [], 2
    while d * d <= n:
        while n % d == 0:
            result.append(d)
            n //= d
        d += 1
    if n > 1:
        result.append(n)
    return tuple(sorted(result, reverse=True))


def kernel6(n):
    integer(n, 1)
    for p in (2, 3):
        while n % p == 0:
            n //= p
    return n


def predecessor_tube(p):
    integer(p, 3)
    need(prime(p), "odd prime sampling label required")
    a = valuation(p - 1)
    t = (p - 1) >> a
    sources = tuple(3 ** j * (1 << (a - j)) * t - 1 for j in range(a))
    for x, y in zip(sources, sources[1:]):
        need(step(x) == (y, 1), "one-run identity")
    return a, sources


W223 = (1, 1, 1, 1, 3)
W233 = (2, 1, 1, 1, 2, 3, 1, 1, 2, 1)
W161_TO_233 = (2, 2, 1, 2, 1, 1)
PAIR = ((223, W223), (233, W233))
TRIPLE = ((161, W161_TO_233 + W233),) + PAIR


def main():
    checks = 0

    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)

    check(carrier((2, 3, 1, 1)) == (81, 128, 223), "exact223 carry")
    check(carrier((5, 2, 1)) == (27, 256, 233), "exact233 carry")
    check(carrier((2, 5, 1)) == (27, 256, 149), "ordered-carry hostile")
    carry_words = 0
    for length in range(1, 6):
        for word in product(range(1, 5), repeat=length):
            check(carrier(word)[2] % 2 == 1, "every nonempty carry odd")
            carry_words += 1
    fib = [0, 1]
    for _ in range(12):
        fib.append(fib[-1] + fib[-2])
    check(fib[13] == 233 and fib[11] + fib[13] == 322, "F13/L12 recovery")
    check(factors(222) == (37, 3, 2) and factors(232) == (29, 2, 2, 2), "prime predecessors")
    check(factors(322) == (23, 7, 2) and not prime(323), "322 is not a prime predecessor")
    check(replay(161, W161_TO_233) == 233, "322 even edge to161 then233")
    check(replay(223, W223) == replay(233, W233) == 425, "literal common future")
    check((len(W223), sum(W223), len(W233), sum(W233),
           len(W161_TO_233 + W233), 1 + sum(W161_TO_233 + W233))
          == (5, 7, 10, 15, 16, 25), "distinct odd/T clocks")

    h, endpoint_period, rows = common_fiber(PAIR)
    check((h, endpoint_period, rows[0][2], rows[1][2]) == (425, 118098, 62208, 65536),
          "primitive two-way fiber")
    h3, p3, rows3 = common_fiber(TRIPLE)
    check((h3, p3, tuple(r[2] for r in rows3))
          == (425, 86093442, (33554432, 45349632, 47775744)), "primitive three-way fiber")
    splices = 0
    for t in range(64):
        n, child = 233 + 65536 * t, 223 + 62208 * t
        selected_t, selected, endpoint = select_source(PAIR, 1, n)
        check(selected_t == t and selected == (child, n), "supplied-source selection")
        check(child == F(243 * n + 469, 256) and child < n, "paid dependency")
        child_proof = root_word_for_control(child)
        received = transport_root(n, W233, child, W223, child_proof)
        check(received == root_word_for_control(n), "independent literal ROOT word")
        splices += 1
        t3, selected3, endpoint3 = select_source(TRIPLE, 0, 161 + 33554432 * t)
        check(t3 == t and endpoint3 == 425 + 86093442 * t, "three-way endpoint")
        check(2 * selected3[0] == 322 + 67108864 * t, "even source sidecar")
    # Necessity of the endpoint residue, independently using all carrier equations.
    fiber_tests = 0
    for endpoint in range(1, 120000, 2):
        compatible = True
        for source, word in PAIR:
            p, q, b = carrier(word)
            numerator = q * endpoint - b
            compatible &= numerator > 0 and numerator % p == 0
        check(compatible == (endpoint >= 425 and (endpoint - 425) % 118098 == 0),
              "all positive endpoints in declared interval")
        fiber_tests += 1

    prime_count = tube_points = 0
    for p in range(3, 5001, 2):
        if not prime(p):
            continue
        prime_count += 1
        a, tube = predecessor_tube(p)
        for x in tube:
            check(kernel6(x + 1) == kernel6(p - 1), "adaptive tube kernel invariant")
            tube_points += 1
        full, rough = factors(p - 1), factors(kernel6(p - 1))
        for j in range(8):
            x = full[j] if j < len(full) else 1
            y = rough[j] if j < len(rough) else 1
            check(x == y or (y == 1 and x <= 3), "uniform normalized-log error premise")
    check(kernel6(58) == kernel6(232) == 29, "same rough kernel")
    check(step(57)[0] < 57 and step(231)[0] > 231, "opposite supplied-source drift")
    a, tube = predecessor_tube(233)
    check(tube == (231, 347, 521), "prime233 predecessor tube")
    check(kernel6(step(tube[-1])[0] + 1) == 49 != kernel6(232), "reset destroys invariant")

    proofs = {n: root_word_for_control(n) for n in (223, 233)}
    for n, proof in proofs.items():
        check(replay(n, proof, root=True) == 1, "valid bank record")
    reversed_injection = {223: (233, proofs[233]), 233: (223, proofs[223])}
    check(len({record[0] for record in reversed_injection.values()}) == 2,
          "cardinal reverse map is injective")
    check(all(reversed_injection[n][0] != n for n in reversed_injection),
          "reverse injection is not a section")
    hostiles = [lambda: predecessor_tube(323), lambda: predecessor_tube(True),
                lambda: predecessor_tube(233.0), lambda: carrier((1, True)),
                lambda: replay(True, (), root=True), lambda: replay(1.0, (), root=True),
                lambda: select_source(PAIR, 1, 235), lambda: select_source(PAIR, True, 233),
                lambda: transport_root(233, W233, 223, W223, proofs[233]),
                lambda: common_fiber(((223, W223), (233, (2,))))]
    for action in hostiles:
        try:
            action()
        except ValueError:
            check(True, "typed/source hostile")
        else:
            check(False, "hostile accepted")

    print("PROVED elementary transfers; accepted paper theorems are premises, not audited")
    print("numeric_roles: carries223/233; F13=233; L12=322; even322->161->233")
    print("carry_parity_universe: lengths1..5, exponents1..4; words=" + str(carry_words))
    print("pair_fiber: sources223+62208t,233+65536t; endpoint425+118098t")
    print("triple_fiber: odd161+33554432t,223+45349632t,233+47775744t")
    print("triple_endpoint=425+86093442t; even_source=322+67108864t")
    print("source_preserving_ROOT_splices_t0..63=" + str(splices))
    print("independent_endpoint_iff_tests=" + str(fiber_tests))
    print("odd_primes3..5000=" + str(prime_count) + "; initial_tube_points=" + str(tube_points))
    print("hostile_prime_labels59/233: same_kernel29, opposite_first_drift")
    print("hostile_reverse_injection: two_valid_ROOT_records, both_wrong_source_fibers")
    print("typed_and_source_hostiles=" + str(len(hostiles)))
    print("checks=" + str(checks))


if __name__ == "__main__":
    main()
