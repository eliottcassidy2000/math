"""Exact marked-cycle capacity and ordered-carry gates; no orbit-convergence oracle.

Run with python -B [-O] this_file.py.  Output is deterministic JSON.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import combinations, product
from math import prod
import json


def integer(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError("exact integer outside domain")
    return x


def odd(x):
    integer(x, 1)
    if x % 2 != 1:
        raise ValueError("positive odd source required")
    return x


def word_type(w):
    if type(w) is not tuple:
        raise ValueError("word must be an exact tuple")
    for a in w:
        integer(a, 1)
    return w


def step(n):
    odd(n)
    z = 3 * n + 1
    a = (z & -z).bit_length() - 1
    return z >> a, a


def carrier(w):
    word_type(w)
    p, q, b = 1, 1, 0
    for a in w:
        p, q, b = 3 * p, (1 << a) * q, 3 * b + q
    return p, q, b


def decode_carry(length, cost, carry):
    """Complete membership test for the specified positive valuation-word carrier."""
    integer(length, 1)
    integer(cost, 1)
    integer(carry, 1)
    r, a, b = length, cost, carry
    if a < r:
        return None
    w = []
    while r > 1:
        z = b - 3 ** (r - 1)
        if z <= 0:
            return None
        first = (z & -z).bit_length() - 1
        if first < 1 or a - first < r - 1:
            return None
        w.append(first)
        a, b, r = a - first, z >> first, r - 1
    if b != 1 or a < 1:
        return None
    return tuple(w + [a])


def primitive_cycle_word(source, length, cost):
    """Exact first-return cycle receipt, or None; bounded by the requested length."""
    odd(source)
    integer(length, 1)
    integer(cost, 1)
    delta = (1 << cost) - 3 ** length
    if delta <= 0:
        return None
    w = decode_carry(length, cost, source * delta)
    if w is None:
        return None
    seen, n = set(), source
    for a in w:
        if n in seen:
            return None
        seen.add(n)
        n, actual = step(n)
        if actual != a:
            return None
    return w if n == source else None


@dataclass(frozen=True)
class Prefix:
    source: int
    word: tuple
    states: tuple


def mark_prefix(source, word):
    """Authenticate a possible cycle-minimum prefix; do not infer a future cycle."""
    odd(source)
    word_type(word)
    if source == 1 or source % 3 == 0:
        raise ValueError("nontrivial cycle minimum must be an odd 3-unit above ROOT")
    states = [source]
    for a in word:
        n, actual = step(states[-1])
        if a != actual:
            raise ValueError("non-native valuation")
        if n < source or n in states or n % 3 == 0:
            raise ValueError("prefix already violates the primitive minimum hypothesis")
        states.append(n)
    return Prefix(source, word, tuple(states))


def audit_prefix(p):
    if type(p) is not Prefix or type(p.states) is not tuple:
        raise ValueError("typed prefix required")
    for n in p.states:
        odd(n)
    expected = mark_prefix(p.source, p.word)
    if p != expected:
        raise ValueError("source prefix does not replay")
    return expected


def observe(source, edges):
    integer(edges)
    w, n = [], source
    odd(n)
    for _ in range(edges):
        n, a = step(n)
        w.append(a)
    return mark_prefix(source, tuple(w))


def occupancy_capacity(minimum, marked, length):
    """Maximum over distinct odd 3-units >= minimum containing every mark."""
    odd(minimum)
    if minimum % 3 == 0 or type(marked) is not tuple:
        raise ValueError("odd 3-unit minimum and tuple marks required")
    integer(length, 1)
    for x in marked:
        odd(x)
        if x < minimum or x % 3 == 0:
            raise ValueError("mark violates the minimum/3-unit condition")
    if len(set(marked)) != len(marked) or len(marked) > length:
        raise ValueError("distinct marked states must fit the cycle")
    chosen = set(marked)
    x = minimum
    while len(chosen) < length:
        if x % 3 and x not in chosen:
            chosen.add(x)
        x += 2
    ordered = tuple(sorted(chosen))
    value = prod((Fraction(3*x + 1, 3*x) for x in ordered), start=Fraction(1))
    return value, ordered


def floor_log2(value):
    if type(value) is not Fraction or value <= 0:
        raise ValueError("positive exact Fraction required")
    n, d = value.numerator, value.denominator
    k = n.bit_length() - d.bit_length()
    if (n < (d << k)) if k >= 0 else ((n << -k) < d):
        k -= 1
    return k


@dataclass(frozen=True)
class ClockGate:
    prefix: Prefix
    length: int
    capacity: Fraction
    extremizer: tuple
    candidates: tuple


def clock_gate(prefix, length):
    p = audit_prefix(prefix)
    integer(length, len(p.states))
    cap, extremizer = occupancy_capacity(p.source, p.states, length)
    lower = max((3**length).bit_length(), sum(p.word) + length - len(p.word))
    upper = floor_log2(cap * 3**length)
    return ClockGate(p, length, cap, extremizer, tuple(range(lower, upper + 1)))


def audit_gate(g):
    if type(g) is not ClockGate or type(g.capacity) is not Fraction:
        raise ValueError("typed gate required")
    if type(g.extremizer) is not tuple or type(g.candidates) is not tuple:
        raise ValueError("exact tuple rows required")
    for x in g.extremizer:
        odd(x)
    for a in g.candidates:
        integer(a, 1)
    expected = clock_gate(g.prefix, g.length)
    if g != expected:
        raise ValueError("noncanonical gate")
    return expected


def log_interval(x, terms):
    """Strict rational bounds on ln(x), x>1, via the atanh series."""
    if type(x) is not Fraction or x <= 1:
        raise ValueError("exact Fraction above one required")
    integer(terms, 1)
    z = (x - 1) / (x + 1)
    lower = 2 * sum((z ** (2*j+1) / (2*j+1) for j in range(terms)), Fraction(0))
    tail = 2 * z**(2*terms+1) / ((2*terms+1) * (1-z*z))
    return lower, lower + tail


def alpha_interval(terms=64):
    lo2, hi2 = log_interval(Fraction(2), terms)
    lo3, hi3 = log_interval(Fraction(3), terms)
    return lo3 / hi2, hi3 / lo2


def certified_cf(terms=64, count=15):
    integer(count, 1)
    lo, hi = alpha_interval(terms)
    digits = []
    for _ in range(count):
        a = lo.numerator // lo.denominator
        if a != hi.numerator // hi.denominator:
            raise ValueError("insufficient rational precision")
        digits.append(a)
        lo, hi = lo - a, hi - a
        if lo <= 0:
            raise ValueError("interval touched a continued-fraction boundary")
        lo, hi = 1 / hi, 1 / lo
    return tuple(digits)


def rhin_minimum_possible(minimum, length):
    """Necessary only, conditional on the cited Simons--de Weger Lemma12."""
    odd(minimum)
    integer(length, 1)
    return (3 * minimum)**10 < 512**10 * length**143


def compositions(total, length):
    for cuts in combinations(range(1, total), length - 1):
        marks = (0,) + cuts + (total,)
        yield tuple(marks[i+1] - marks[i] for i in range(length))


CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def rejects(call):
    try:
        call()
    except (ValueError, TypeError):
        check(True, "guard rejected")
        return
    check(False, "guard accepted malformed input")


def main():
    # Independent composition and exact decoder controls, all 4-ary words <=5.
    words = 0
    for length in range(1, 6):
        for w in product(range(1, 5), repeat=length):
            p, q, b = carrier(w)
            prefix = 0
            independent = 0
            for i, a in enumerate(w):
                independent += 3**(length-1-i) * 2**prefix
                prefix += a
            check((p, q, b) == (3**length, 2**sum(w), independent), "carry sum")
            check(decode_carry(length, sum(w), b) == w, "exact word round trip")
            words += 1
    # Complete carry-membership comparison on every clock length<=5, cost<=10.
    for length in range(1, 6):
        for cost in range(length, 11):
            image = {carrier(w)[2]: w for w in compositions(cost, length)}
            check(len(image) == len(tuple(compositions(cost, length))), "carry injectivity")
            candidates = set(image)
            candidates.update(range(1, 301, 2))
            for b in sorted(candidates):
                check(decode_carry(length, cost, b) == image.get(b), "carry membership")
    # Brute subset maximization, independent from greedy capacity filling.
    subsets = 0
    for minimum in (1, 5, 7, 11):
        pool = tuple(x for x in range(minimum, minimum+40, 2) if x % 3)[:10]
        for length in range(1, 6):
            marked = (pool[-1],) if length > 1 else ()
            cap, _ = occupancy_capacity(minimum, marked, length)
            values = []
            for xs in combinations(pool, length):
                if not set(marked).issubset(xs):
                    continue
                value = prod((Fraction(3*x+1, 3*x) for x in xs), start=Fraction(1))
                check(value <= cap, "capacity upper bound")
                values.append(value)
                subsets += 1
            check(max(values) == cap, "capacity attained")
    cf = certified_cf()
    check(cf == (1,1,1,2,2,3,1,5,2,23,2,2,1,1,55), "certified CF prefix")
    a0, a1 = alpha_interval(32)
    b0, b1 = alpha_interval(64)
    check(a0 < b0 < b1 < a1, "strict nested log bounds")
    # Exact resonance is not a cycle receipt.
    carries = [carrier(w)[2] for w in compositions(8, 5)]
    check(len(carries) == 35 and all(b % 13 for b in carries), "8/5 carry hostile")
    check(primitive_cycle_word(1, 1, 2) == (2,), "ROOT cycle")
    check(primitive_cycle_word(1, 2, 4) is None, "primitive/ROOT padding")
    check(Fraction(carrier((1,3))[2], 2**4-3**2) == Fraction(5,7), "rational cycle")
    check((3*Fraction(5,7)+1)/2 == Fraction(11,7), "rational first edge")
    check((3*Fraction(11,7)+1)/8 == Fraction(5,7), "rational return")
    # Arbitrarily long real prefixes whose clock ratio is 1, not close to alpha.
    for k in range(2, 65):
        n = 2**k - 1
        ratio = Fraction(1)
        for _ in range(k-1):
            ratio *= Fraction(3*n+1, 3*n)
            n, a = step(n)
            check(a == 1, "Mersenne native prefix")
        check(n == 2*3**(k-1)-1, "Mersenne endpoint")
        check(Fraction(2**(k-1), 3**(k-1)) == Fraction(2**k-1, n)*ratio,
              "endpoint-corrected product identity")
    # The safe constant 512 follows exactly from two terms of ln2's series.
    lo2, _ = log_interval(Fraction(2), 2)
    check(lo2 == Fraction(56,81) and Fraction(133,10)*Fraction(46057,100000) < 9*lo2,
          "exact weakening of cited Rhin constant")
    check(rhin_minimum_possible(91, 29), "Rhin does not exclude the marked survivor")
    check(not rhin_minimum_possible(10**40+1, 2), "Rhin minimum-height gate")
    table = []
    for source in (31, 47, 71, 91):
        p = observe(source, 11)
        counts = [0, 0, 0]
        first = [None, None, None]
        for length in range(12, 257):
            gate = clock_gate(p, length)
            check(audit_gate(gate) == gate, "canonical gate")
            c0 = Fraction(3*source+1, 3*source)**length
            c1, _ = occupancy_capacity(source, (), length)
            check(gate.capacity <= c1 <= c0, "monotone refinement")
            lower = max((3**length).bit_length(), sum(p.word)+length-len(p.word))
            for j, cap in enumerate((c0, c1, gate.capacity)):
                viable = (2**lower) * cap.denominator <= 3**length * cap.numerator
                counts[j] += viable
                if viable and first[j] is None:
                    first[j] = length
                if j == 2:
                    check(viable == bool(gate.candidates), "independent cross multiplication")
            # Refinement must retain the original source and only remove clocks.
            previous = clock_gate(mark_prefix(source, ()), length)
            for depth in range(1, len(p.word)+1):
                refined = clock_gate(mark_prefix(source, p.word[:depth]), length)
                check(refined.capacity <= previous.capacity, "marked capacity refinement")
                check(set(refined.candidates).issubset(previous.candidates), "clock refinement")
                previous = refined
        table.append({"source": source, "word": p.word, "states": p.states,
                      "candidate_length_counts": counts, "least_lengths": first})
    p = observe(91, 11)
    g17, g29 = clock_gate(p, 17), clock_gate(p, 29)
    check(g17.candidates == () and g29.candidates == (46,), "marked improvement and survivor")
    check(primitive_cycle_word(91, 29, 46) is None, "capacity survivor still needs carry")
    for bad in (True, 91.0, 0, 2, Fraction(91)):
        rejects(lambda bad=bad: mark_prefix(bad, ()))
    rejects(lambda: mark_prefix(27, ()))
    rejects(lambda: mark_prefix(1, ()))
    rejects(lambda: mark_prefix(91, (2,)))
    rejects(lambda: mark_prefix(5, (4,)))
    rejects(lambda: audit_prefix(replace(p, states=(91.0,) + p.states[1:])))
    rejects(lambda: audit_prefix(replace(p, states=p.states[:-1] + (253,))))
    rejects(lambda: audit_gate(replace(g29, candidates=(46.0,))))
    rejects(lambda: audit_gate(replace(g29, capacity=1.0)))
    rejects(lambda: audit_gate(replace(g29, candidates=(45,46))))
    rejects(lambda: decode_carry(True, 2, 1))
    rejects(lambda: decode_carry(1, 2.0, 1))
    rejects(lambda: primitive_cycle_word(Fraction(5,7), 2, 4))
    report = {
        "status": "FINITE-EXACT; necessary cycle gates, not convergence",
        "checks": CHECKS,
        "universes": {"word_alphabet": [1,2,3,4], "max_word_length": 5,
                      "word_count": words, "subset_controls": subsets,
                      "sources": [31,47,71,91], "observed_edges": 11,
                      "candidate_lengths_inclusive": [12,256]},
        "certified_cf_prefix": cf,
        "capacity_modes": ["repeated minimum", "distinct odd 3-units", "marked actual prefix"],
        "capacity_table": table,
        "source91_clock27_over17_capacity": [str(g17.capacity.numerator), str(g17.capacity.denominator)],
        "source91_clock46_over29": "capacity passes; ordered carry decoder rejects",
        "clock8_over5": {"words": 35, "gap": 13, "integer_cycles": 0},
        "rational_cycle_hostile": {"word": [1,3], "states": ["5/7", "11/7"]},
        "mersenne_prefix_hostile": {"exponents_checked": [2,64], "word": "1^(K-1)",
                                    "endpoint": "2*3^(K-1)-1", "clock_ratio": 1},
        "scope": "No assertion mu(log_2(3))=2 is assumed; no unbounded trajectory is excluded."
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
