"""Exact whole-block Collatz routing; finite controls do not imply coverage."""
from dataclasses import dataclass
from fractions import Fraction
from itertools import product


CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


def word_guard(word):
    if not isinstance(word, tuple) or any(type(a) is not int or a < 1 for a in word):
        raise ValueError("word must be a tuple of positive exact integers")


def carrier(word):
    word_guard(word)
    b = a = 0
    for e in word:
        b = 3 * b + (1 << a)
        a += e
    return len(word), a, b


def independent_carry(word):
    return sum(3 ** (len(word) - i - 1) * (1 << sum(word[:i]))
               for i in range(len(word)))


def source_class(word):
    length, cost, b = carrier(word)
    modulus = 1 << (cost + 1)
    return ((1 << cost) - b) * pow(3 ** length, -1, modulus) % modulus, modulus


def replay(n, word):
    word_guard(word)
    if type(n) is not int or n < 1 or n % 2 == 0:
        raise ValueError("positive odd source required")
    states = [n]
    for exponent in word:
        q = 3 * n + 1
        actual = (q & -q).bit_length() - 1
        if actual != exponent:
            raise ValueError("source does not obey the word guard")
        n = q >> actual
        states.append(n)
    return tuple(states)


def swap_shift(word, i, j):
    word_guard(word)
    if type(i) is not int or type(j) is not int or not 0 <= i < j < len(word):
        raise ValueError("two distinct ordered exact indices required")
    prefix, middle = word[:i], word[i + 1:j]
    _, p, _ = carrier(prefix)
    _, c, b = carrier(middle)
    return Fraction((1 << p) * ((1 << word[i]) - (1 << word[j]))
                    * (3 * b + (1 << c)), 3 ** (j + 1))


@dataclass(frozen=True)
class BlockRoute:
    source_word: tuple
    child_word: tuple
    decrement: int
    residue: int
    modulus: int


def compile_route(source_word, child_word):
    left = carrier(source_word)
    right = carrier(child_word)
    if sorted(source_word) != sorted(child_word):
        raise ValueError("this compiler requires the same valuation multiset")
    length, cost, b = left
    _, _, B = right
    if B <= b or (B - b) % (3 ** length):
        raise ValueError("no integral strictly smaller translated child")
    decrement = (B - b) // (3 ** length)
    residue, modulus = source_class(source_word)
    return BlockRoute(source_word, child_word, decrement, residue, modulus)


def route(n, receipt):
    if not isinstance(receipt, BlockRoute):
        raise ValueError("typed compiled receipt required")
    if any(type(v) is not int for v in (receipt.decrement, receipt.residue, receipt.modulus)):
        raise ValueError("receipt metadata must be exact integers")
    canonical = compile_route(receipt.source_word, receipt.child_word)
    if receipt != canonical:
        raise ValueError("receipt metadata was altered")
    if type(n) is not int or n <= receipt.decrement or n % receipt.modulus != receipt.residue:
        raise ValueError("positive source outside this exact family")
    m = n - receipt.decrement
    left, right = replay(n, receipt.source_word), replay(m, receipt.child_word)
    if left[-1] != right[-1]:
        raise ValueError("common endpoint mismatch")
    return m, left[-1]


def root_word(n, limit=10000):
    if type(n) is not int or n < 1 or n % 2 == 0:
        raise ValueError("positive odd source required")
    out = []
    while n != 1 and len(out) < limit:
        q = 3 * n + 1
        e = (q & -q).bit_length() - 1
        out.append(e)
        n = q >> e
    if n != 1:
        raise ValueError("finite control exceeded its search limit")
    return tuple(out)


def splice_root(n, receipt, supplied_child_word):
    m, _ = route(n, receipt)
    child_states = replay(m, supplied_child_word)
    if child_states[-1] != 1 or 1 in child_states[:-1]:
        raise ValueError("supplied child certificate is not first-hit ROOT")
    k = len(receipt.child_word)
    if supplied_child_word[:k] != receipt.child_word:
        raise ValueError("supplied child certificate does not contain the join")
    out = receipt.source_word + supplied_child_word[k:]
    states = replay(n, out)
    if states[-1] != 1 or 1 in states[:-1]:
        raise ValueError("spliced certificate is not first-hit ROOT")
    return out


def expect_reject(fn):
    try:
        fn()
    except ValueError:
        check(True, "hostile rejected")
    else:
        check(False, "hostile accepted")


def census():
    """Complete declared universe, not a global minimality claim."""
    rows = []
    for length in range(1, 13):
        modulus = 3 ** length
        buckets = {}
        growing = []
        for word in product((1, 2, 3), repeat=length):
            b = cost = c2 = c3 = 0
            pow3 = 1
            prefix_growth = True
            for e in word:
                b = 3 * b + (1 << cost)
                cost += e
                c2 += e == 2
                c3 += e == 3
                pow3 *= 3
                if (1 << cost) > pow3:
                    prefix_growth = False
            key = (c2, c3, b % modulus)
            prior = buckets.get(key)
            if prior is None or b > prior[0]:
                buckets[key] = (b, word)
            if prefix_growth:
                growing.append((key, b, word))
        hits = []
        for key, b, word in growing:
            B, partner = buckets[key]
            if B > b:
                hits.append((word, partner, (B - b) // modulus))
        rows.append((length, 3 ** length, len(growing), len(hits)))
        if length < 12:
            check(not hits, "no earlier growing source in this finite universe")
        else:
            check(hits == [(W, V, 4)], "unique first growing source and maximal-carry partner")
    return rows


W = (1, 1, 1, 1, 1, 2, 2, 2, 3, 1, 1, 3)
V = (1, 2, 1, 3, 1, 3, 1, 2, 1, 1, 2, 1)


def main():
    transpositions = 0
    for length in range(1, 6):
        for word in product(range(1, 6), repeat=length):
            _, _, b = carrier(word)
            check(b == independent_carry(word), "independent carry formula")
            for i in range(length):
                for j in range(i + 1, length):
                    other = list(word)
                    other[i], other[j] = other[j], other[i]
                    B = carrier(tuple(other))[2]
                    shift = swap_shift(word, i, j)
                    check(shift == Fraction(b - B, 3 ** length), "nonadjacent defect formula")
                    integer = (word[i] - word[j]) % (2 * 3 ** j) == 0
                    check((shift.denominator == 1) == integer, "exact swap divisibility")
                    transpositions += 1
    # Include large differences that really do pass the divisibility guard.
    for j in range(1, 6):
        word = (1,) * j + (1 + 2 * 3 ** j,)
        check(swap_shift(word, 0, j).denominator == 1, "positive legal individual swap")

    receipt = compile_route(W, V)
    check((receipt.residue, receipt.modulus, receipt.decrement)
          == (257727, 1048576, 4), "canonical family")
    check(carrier(W) == (12, 19, 923953), "source carrier")
    check(carrier(V) == (12, 19, 3049717), "partner carrier")
    for j in range(1, len(W) + 1):
        check((1 << sum(W[:j])) <= 3 ** j, "every prefix has noncontracting slope")
    for word in (W, V):
        for i in range(len(word)):
            for j in range(i + 1, len(word)):
                if word[i] != word[j]:
                    check(swap_shift(word, i, j).denominator != 1, "individual move is illegal")
    for t in range(128):
        n = receipt.residue + receipt.modulus * t
        m, endpoint = route(n, receipt)
        check((m, endpoint) == (n - 4, 261245 + 1062882 * t), "exact lifted endpoint")
        check(min(replay(n, W)[1:]) > n, "literal prefixes grow")
        supplied = root_word(m)
        obtained = splice_root(n, receipt, supplied)
        check(obtained == root_word(n), "independent literal ROOT splice")
    for n in range(1, 8192, 2):
        # A short sibling pair is a separate direct residue/trajectory control.
        word = (2, 3, 1, 3)
        residue, modulus = source_class(word)
        try:
            replay(n, word)
            legal = True
        except ValueError:
            legal = False
        check(legal == (n % modulus == residue), "source cylinder iff direct replay")

    rows = census()
    malformed = [(), (0,), (1, False), (1, 2.0)]
    for w in malformed:
        expect_reject(lambda w=w: compile_route(w, V))
    for n in (True, 257727.0, 257729, -1, 1):
        expect_reject(lambda n=n: route(n, receipt))
    expect_reject(lambda: swap_shift(W, True, 2))
    expect_reject(lambda: route(257727, BlockRoute(W, V, 2, 257727, 1048576)))
    expect_reject(lambda: splice_root(257727, receipt, root_word(233)))

    print("SWAPS", transpositions, "; independent carry, exact defect and divisibility checked")
    print("BLOCK", W, "=>", V)
    print("FAMILY n=257727+1048576t; child=n-4; endpoint=261245+1062882t; t>=0")
    print("ALL12 source prefixes grow; every nontrivial single transposition is nonintegral")
    print("FINITE ROOT controls:128 parameters; authentic child words reproduce literal source words")
    print("CENSUS length,total_words,prefix_slope_growing_sources,with_smaller_permutation")
    for row in rows:
        print(*row)
    print("FINITE universe total words", sum(row[1] for row in rows))
    print("SCOPE: direct12 coverage is extended; inherited adaptive selector already covers t=0")
    print("CHECKS", CHECKS)


if __name__ == "__main__":
    main()
