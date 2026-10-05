"""Finite exact W-superlevel ROOT receipts; no forward orbit discovery.

Public interface: threshold_receipts(Fraction) and certificate_weight(n, word).
An absent source has W(n)<epsilon; this is not a nonconvergence verdict.
"""
from collections import deque
from fractions import Fraction
from itertools import product
from math import comb


def require(ok, message):
    if not ok:
        raise ValueError(message)


def threshold(epsilon):
    require(type(epsilon) is Fraction and 0 < epsilon <= 1,
            "epsilon must be an exact Fraction in (0,1]")
    return epsilon


def odd_source(n):
    require(type(n) is int and n > 0 and n % 2 == 1,
            "source must be a positive odd exact integer")
    return n


def counter_weight(l, k):
    require(type(l) is int and l >= 0 and type(k) is int and k >= 0,
            "nonnegative exact counter integers")
    return Fraction(2, (l+k+2)*comb(l+k+1, k))


def certificate_weight(source, word):
    """Authenticate the source and strict first-hit word, then evaluate W."""
    odd_source(source)
    require(type(word) is tuple and all(type(a) is int and a > 0 for a in word),
            "word must be a tuple of positive exact integers")
    current = source
    for a in word:
        require(current != 1, "root padding is not a first-hit receipt")
        raw = 3*current+1
        actual = (raw & -raw).bit_length()-1
        require(actual == a, "false valuation")
        current = raw >> a
    require(current == 1, "word does not reach ROOT")
    if source == 1:
        require(not word, "ROOT has the empty word")
        return Fraction(1)
    require(len(word) >= 1 and word[-1] >= 4 and word[-1] % 2 == 0,
            "nonroot last edge is an even valuation at least four")
    l = len(word)-1
    k = sum((a-1)//2 for a in word)
    require(k >= 1, "nonroot ROOT receipt has positive K")
    return counter_weight(l, k)


def counter_bounds(epsilon):
    threshold(epsilon)
    kmax = 0
    while counter_weight(0, kmax+1) >= epsilon:
        kmax += 1
    lmax = -1
    while counter_weight(lmax+1, 1) >= epsilon:
        lmax += 1
    return lmax, kmax


def compile_threshold(epsilon, rectangular=False):
    """Return the exact threshold bank, or the comparison counter rectangle."""
    threshold(epsilon)
    require(type(rectangular) is bool, "rectangular flag must be bool")
    lmax, kmax = counter_bounds(epsilon)
    records = {1: ((), 0, 0)}
    queue = deque([1])
    attempted = 0
    while queue:
        parent = queue.popleft()
        suffix, l, k = records[parent]
        allowance = kmax-k if rectangular or parent == 1 else kmax-l-k-1
        for a in range(1, 2*allowance+3):
            attempted += 1
            numerator = (1 << a)*parent-1
            if numerator % 3:
                continue
            child = numerator//3
            if child <= 1:
                continue
            child_l = 0 if parent == 1 else l+1
            child_k = k+(a-1)//2
            if child_l > lmax or child_k > kmax:
                continue
            if not rectangular and child_l+child_k > kmax:
                continue
            if not rectangular and counter_weight(child_l, child_k) < epsilon:
                continue
            require(child not in records, "duplicate strict ROOT inverse address")
            records[child] = ((a,)+suffix, child_l, child_k)
            queue.append(child)
    words = {n: item[0] for n, item in records.items()}
    return words, {"Lmax": lmax, "Kmax": kmax, "inverse_candidates": attempted}


def threshold_receipts(epsilon):
    """Exact {source: first-hit word} for all W(source)>=epsilon."""
    return compile_threshold(epsilon)[0]


def source_receipt(source, epsilon):
    """None means W(source)<epsilon, not that source never reaches ROOT."""
    odd_source(source)
    return threshold_receipts(epsilon).get(source)


def weak_compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for tail in weak_compositions(total-first, length-1):
                yield (first,)+tail


def independent_word_bank(epsilon):
    """Small audit only: enumerate whole words, then solve their affine source."""
    lmax, kmax = counter_bounds(epsilon)
    result = {1: ()}
    words_tested = 0
    for l in range(lmax+1):
        length = l+1
        for k in range(1, kmax+1):
            if counter_weight(l, k) < epsilon:
                continue
            for depths in weak_compositions(k, length):
                for bits in product((1, 2), repeat=length):
                    word = tuple(2*d+b for d, b in zip(depths, bits))
                    words_tested += 1
                    p, q, carry = 1, 1, 0
                    for a in word:
                        p, q, carry = 3*p, q*(1 << a), 3*carry+q
                    if (q-carry) <= 0 or (q-carry) % p:
                        continue
                    source = (q-carry)//p
                    try:
                        value = certificate_weight(source, word)
                    except ValueError:
                        continue
                    require(value >= epsilon, "independent source weight")
                    require(source not in result, "independent duplicate")
                    result[source] = word
    return result, words_tested


def refinement_row(bank, source, modulus):
    odd_source(source)
    require(type(modulus) is int and modulus > 0, "positive exact modulus")
    candidates = [(certificate_weight(n, word), n)
                  for n, word in bank.items() if n % modulus == source % modulus]
    if not candidates:
        return None
    return max(candidates)


def main():
    from collatz_floor_transport_deadlines_20261005 import threshold_receipt
    checks = 0
    def check(ok, label):
        nonlocal checks
        require(ok, label)
        checks += 1

    levels = (Fraction(1, 3), Fraction(1, 6), Fraction(1, 10),
              Fraction(1, 20), Fraction(1, 100), Fraction(1, 1000))
    previous = {}
    banks = {}
    print("W-superlevel receipts: PROVED finite compiler; FINITE-EXACT controls")
    for epsilon in levels:
        bank, stats = compile_threshold(epsilon)
        banks[epsilon] = bank
        check(previous.keys() <= bank.keys(), "persistent receipts under refinement")
        for n, word in bank.items():
            check(certificate_weight(n, word) >= epsilon, "authenticated threshold member")
            if n != 1:
                l, k = len(word)-1, sum((a-1)//2 for a in word)
                check(l <= stats["Lmax"] and k <= stats["Kmax"], "counter bounds")
                check(l+k <= stats["Kmax"], "sharper shared counter bound")
                check(sum(word) <= 2*(k+l+1), "halving-cost bound")
                check(n <= ((1 << (2*(stats["Kmax"]+1)))-1)//3,
                      "sharp source-height bound")
        check(max(bank) == ((1 << (2*(stats["Kmax"]+1)))-1)//3,
              "root sibling attains the maximum source")
        check(len(bank)*epsilon <= Fraction(16, 3), "inherited mass cardinality bound")
        for source in range(1, 256, 2):
            decision = threshold_receipt(source, epsilon)
            check((decision["status"] == "met") == (source in bank),
                  "independent bounded forward decision versus inverse whole bank")
        print("epsilon", epsilon, "sources", len(bank), "candidates", stats["inverse_candidates"],
              "Lmax", stats["Lmax"], "Kmax", stats["Kmax"],
              "max_source", max(bank),
              "mass", sum((certificate_weight(n, word) for n, word in bank.items()), Fraction(0)))
        if epsilon >= Fraction(1, 20):
            rectangle, rstats = compile_threshold(epsilon, rectangular=True)
            filtered = {n: word for n, word in rectangle.items()
                        if certificate_weight(n, word) >= epsilon}
            independent, candidates = independent_word_bank(epsilon)
            check(filtered == bank == independent, "three independent finite bank paths")
            print("  rectangle", len(rectangle), "candidates", rstats["inverse_candidates"],
                  "; independent whole words", candidates)
        previous = bank
    check(banks[Fraction(1, 3)] == {1: (), 5: (4,)}, "equality at threshold1/3")
    check(set(banks[Fraction(1, 6)]) == {1, 3, 5, 21}, "equality at threshold1/6")
    bank = banks[Fraction(1, 1000)]
    for source in (7, 27, 53, 113, 155, 703):
        rows = []
        for j in range(1, 7):
            result = refinement_row(bank, source, 5**j)
            rows.append((5**j, None if result is None else (str(result[0]), result[1])))
        print("refinement source", source, rows)
    check([refinement_row(bank, 27, 5**j) for j in range(1, 5)]
          == [(Fraction(1, 30), 17), (Fraction(2, 105), 227),
              (Fraction(1, 140), 277), None], "changing source witnesses")
    check(refinement_row(bank, 7, 25) == (Fraction(1, 84), 7), "grounded persistent floor")
    check(27 not in bank, "absence is only a weight upper bound")
    # A supplied receipt, not a search input, independently verifies a smaller floor.
    word27 = (1, 2, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 3, 1, 1,
              2, 1, 2, 1, 1, 1, 1, 1, 3, 1, 1, 1, 4, 2, 2, 4, 3, 1, 1, 5, 4)
    value27 = certificate_weight(27, word27)
    check(0 < value27 < Fraction(1, 1000), "27 is rooted below current threshold")
    print("supplied independent27 receipt weight", value27, "; not an enumeration seed")
    # Exact outputs of the independently grounded polynomial-dual package.
    # These are readout controls, never used as inverse-tree seeds.
    floors = ((3, Fraction(80890082949936335617, 537701516570945126400)),
              (9, Fraction(14270710919159706722440741432397,
                           2522492030901519845183309102448640)))
    for source, epsilon in floors:
        receipt = source_receipt(source, epsilon)
        check(receipt is not None, "external true atom floor discharges to a receipt")
        value = certificate_weight(source, receipt)
        check(value >= epsilon, "independent authentication of external floor")
        print("polynomial-dual floor source", source, "word", receipt, "actual_weight", value)
    bad = (lambda: threshold_receipts(0), lambda: threshold_receipts(0.1),
           lambda: threshold_receipts(True), lambda: threshold_receipts(Fraction(0)),
           lambda: certificate_weight(True, ()), lambda: certificate_weight(1.0, ()),
           lambda: certificate_weight(2, ()), lambda: certificate_weight(1, (2,)),
           lambda: certificate_weight(5, (4, 2)), lambda: certificate_weight(5, (True,)),
           lambda: certificate_weight(3, (2,)), lambda: source_receipt(False, Fraction(1, 3)))
    for job in bad:
        try:
            job()
        except ValueError:
            check(True, "typed/first-hit hostile rejected")
        else:
            raise ValueError("invalid receipt or threshold accepted")
    print("checks", checks)
    print("bounded-forward cross-component controls: 128 sources at each of six thresholds")
    print("No source coverage theorem: a floor must be supplied or generated by a verified receipt.")


if __name__ == "__main__":
    main()
