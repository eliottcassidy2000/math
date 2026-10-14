"""Exact power fibres, signed Collatz observers and false squared-anchor joins.

No orbit discovery or universal ROOT conclusion. All guards keep their source.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt
import json

CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise ValueError(message)


def integer(value, lower=None):
    need(type(value) is int and (lower is None or value >= lower), "exact integer")
    return value


def odd(value):
    integer(value)
    need(value % 2 == 1, "nonzero signed odd source")
    return value


def valuations(word, reduced=False):
    need(type(word) is tuple and bool(word), "nonempty exact valuation tuple")
    need(all(type(a) is int and a >= 1 for a in word), "positive exact valuations")
    need(not reduced or word[0] >= 2, "reduced head")
    return word


def prime(p):
    integer(p, 3)
    need(p % 2 == 1 and all(p % d for d in range(3, isqrt(p)+1, 2)), "odd prime")
    return p


def pushforward(modulus, exponent, weights):
    integer(modulus, 2)
    integer(exponent, 1)
    need(type(weights) is tuple and len(weights) == modulus, "one weight per source residue")
    need(all(type(w) is int and w >= 0 for w in weights), "nonnegative integer weights")
    out = [0]*modulus
    for x, weight in enumerate(weights):
        out[pow(x, exponent, modulus)] += weight
    return tuple(out)


def character(p, x):
    prime(p)
    integer(x)
    y = pow(x % p, (p-1)//2, p)
    return -1 if y == p-1 else y


def correct_nonresidue(p, x, multiplier):
    """Choose a quadratic-residue representative; this changes the source."""
    prime(p)
    integer(x)
    integer(multiplier)
    need(character(p, multiplier) == -1, "marked nonresidue multiplier")
    return (multiplier*x if character(p, x) == -1 else x) % p


@dataclass(frozen=True)
class MarkedSquare:
    square: int
    sign: int


def encode(source):
    odd(source)
    return MarkedSquare(source*source, 1 if source > 0 else -1)


def decode(record):
    need(type(record) is MarkedSquare, "exact marked-square record")
    integer(record.square, 1)
    need(type(record.sign) is int and record.sign in (-1, 1), "exact sign marker")
    magnitude = isqrt(record.square)
    need(magnitude*magnitude == record.square and magnitude % 2 == 1,
         "positive odd perfect square")
    return record.sign*magnitude


def step(source):
    odd(source)
    z = 3*source+1
    a = (abs(z) & -abs(z)).bit_length()-1
    return z//(1 << a), a


def marked_step(record, expected_source):
    odd(expected_source)
    source = decode(record)
    need(source == expected_source, "unchanged supplied source")
    target, a = step(source)
    return encode(target), a


def replay(source, word, terminal=None):
    odd(source)
    valuations(word)
    x = source
    for a in word:
        need(terminal is None or x != terminal, "no supplied-terminal padding")
        x, actual = step(x)
        need(a == actual, "actual source valuation")
    return x


def root_receipt(record, source, word):
    integer(source, 1)
    need(source % 2 == 1 and decode(record) == source, "positive ROOT source binding")
    if source == 1:
        need(type(word) is tuple and word == (), "ROOT has the empty word")
    else:
        need(replay(source, word, 1) == 1, "supplied strict ROOT word")
    return {"source": source, "word": word, "terminal": 1}


def carrier(word):
    valuations(word)
    p = q = 1
    b = 0
    for a in word:
        p, q, b = 3*p, q*(1 << a), 3*b+q
    return p, q, b


def anchor(word):
    valuations(word, True)
    p, q, b = carrier(word)
    return F(b-p, q)


def anchor_mark(word):
    """Total, absolute odd numerator, sign; not a decoder of the word."""
    value = anchor(word)
    total = sum(word)
    k = value*(1 << (total-1))
    need(k.denominator == 1 and k.numerator % 2 == 1, "odd anchor numerator")
    return total, abs(k.numerator), 1 if k > 0 else -1


def false_square_families(total):
    integer(total, 4)
    return ((total,), (2, total-2)), ((3, total-3), (2, 1, total-3))


def native_residue(word):
    p, q, b = carrier(word)
    modulus = 2*q
    return ((q-b)*pow(p, -1, modulus)) % modulus, modulus


def false_switch(total, run, parameter):
    """An authenticated source and a deliberately invalid squared-anchor partner."""
    integer(total, 4)
    integer(run, 1)
    integer(parameter, 1)
    left = (1,)*run+(total,)
    right = (1,)*(run-1)+(3, total-3)
    r, modulus = native_residue(left)
    source = r+modulus*parameter
    child = (source+1)//2-1
    endpoint = replay(source, left, 1)
    p, q, b = carrier(right)
    fake_endpoint = F(p*child+b, q)
    need(fake_endpoint-endpoint == F(4, 1 << total), "parallel affine gap")
    need(fake_endpoint.denominator > 1, "partner is not a native integer endpoint")
    return {"source": source, "child": child, "left": left, "right": right,
            "endpoint": endpoint, "fake_endpoint": str(fake_endpoint)}


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for a in range(1, total-length+2):
            for tail in compositions(total-a, length-1):
                yield (a,)+tail


def reject(call):
    try:
        call()
    except ValueError:
        return
    raise ValueError("malformed or hostile input accepted")


def main():
    primes = [p for p in range(3, 102, 2) if all(p % d for d in range(3, isqrt(p)+1, 2))]
    for p in primes:
        squares = pushforward(p, 2, (1,)*p)
        need(squares[0] == 1, "zero is fixed, not sent to one")
        for y in range(1, p):
            need(squares[y] == 1+character(p, y), "exact square fibre count")
        nonresidue = next(x for x in range(1, p) if character(p, x) == -1)
        for x in range(p):
            y = correct_nonresidue(p, x, nonresidue)
            need(character(p, y) >= 0 and correct_nonresidue(p, y, nonresidue) == y,
                 "nonresidue correction is idempotent")
            if p % 4 == 3:
                need(correct_nonresidue(p, x, -1) == pow(x*x % p, (p+1)//4, p),
                     "canonical sign representative from square")
            else:
                need(character(p, -x) == character(p, x), "negation cannot repair this prime")
        for j in range(1, 7):
            d = gcd(1 << j, p-1)
            counts = pushforward(p, 1 << j, (0,)+(1,)*(p-1))
            need(set(counts)-{0} == {d} and sum(c != 0 for c in counts) == (p-1)//d,
                 "iterated-square pushforward fibres")
    for modulus in (3, 4, 8, 12, 24):
        weights = tuple((x+1)*(x+2) for x in range(modulus))
        need(sum(pushforward(modulus, 2, weights)) == sum(weights), "source mass retained")
    for n in range(-1001, 1002):
        need((n*n % 24 == 1) == (gcd(n, 6) == 1), "square mod24 exact source predicate")
    need(3*3 % 3 == 0 and 3*3 % 24 == 9 and 35*35 % 24 == 1,
         "prime3 exception and composite survivor")
    for n in range(1, 2000, 2):
        plus, ap = step(n)
        minus, am = step(-n)
        need((plus*plus == minus*minus) == (n == 1), "no deterministic square dynamics")
        need(plus % 3 == pow(-1, ap, 3) and minus % 3 == pow(-1, am, 3),
             "signed endpoint retains valuation parity")
        for source in (n, -n):
            target, a = marked_step(encode(source), source)
            need(decode(encode(source)) == source and decode(target) == step(source)[0],
                 "marked source step exactly decoded")
    for E in range(3, 65):
        low, high = (1 << E)-1, (1 << E)+1
        need((low*low-1) % (1 << (E+1)) == 0 and
             (high*high-1) % (1 << (E+1)) == 0, "square forgets anchor branch")
        need(step(low)[1] == 1 and step(high)[1] == 2, "lost anchor has different reset")
    exact_false = 0
    reduced_count = 0
    for total in range(2, 17):
        fibres = {}
        for length in range(1, total):
            for w in compositions(total, length):
                if w[0] < 2:
                    continue
                reduced_count += 1
                a = anchor(w)
                mark = anchor_mark(w)
                fibres.setdefault(a*a, []).append((w, a, mark))
        for fibre in fibres.values():
            for w, a, ma in fibre:
                for v, b, mb in fibre:
                    need((a == b) == (ma == mb), "one sign marker repairs equality")
                    if a != b:
                        neg, pos = false_square_families(total)
                        need((w in neg and v in pos) or (v in neg and w in pos),
                             "complete false squared-anchor classification")
                        exact_false += 1
        if total >= 4:
            neg, pos = false_square_families(total)
            need(all(anchor(w) == -F(2, 1 << total) for w in neg) and
                 all(anchor(w) == F(2, 1 << total) for w in pos), "all four listed words exist")
    for A, r, t in product(range(4, 13), range(1, 6), (1, 2, 7)):
        false_switch(A, r, t)
    witness = false_switch(4, 1, 1)
    need(witness["source"] == 67 and witness["child"] == 33 and witness["endpoint"] == 19,
         "small unchanged-source hostile")
    reject(lambda: replay(33, (3, 1)))
    need(replay(63, (1,)*5+(4,)) == 91, "actual Mersenne squared-switch hostile")
    p, q, b = carrier((1,)*4+(3, 1))
    need(F(p*31+b, q) == F(365, 4), "Mersenne partner retains the nonzero gap")
    reject(lambda: replay(31, (1,)*4+(3, 1)))
    need(root_receipt(encode(3), 3, (1, 4))["terminal"] == 1, "positive supplied ROOT control")
    need(root_receipt(encode(1), 1, ())["word"] == (), "ROOT empty-word boundary")
    reject(lambda: root_receipt(encode(1), 1, (2,)))
    reject(lambda: root_receipt(encode(-3), -3, (3,)))
    reject(lambda: root_receipt(encode(3), 5, (4,)))
    reject(lambda: root_receipt(encode(3), 3, (1, 4, 2)))
    for bad in (MarkedSquare(True, 1), MarkedSquare(9, True), MarkedSquare(8, 1),
                MarkedSquare(4, 1), MarkedSquare(9, 0)):
        reject(lambda bad=bad: decode(bad))
    reject(lambda: marked_step(encode(3), -3))
    reject(lambda: anchor((1, 4)))
    reject(lambda: anchor((True,)))
    reject(lambda: correct_nonresidue(5, 2, -1))
    reject(lambda: pushforward(3, 2, (1, True, 1)))
    tables = {m: pushforward(m, 2, (1,)*m) for m in (3, 4, 24)}
    signed = {"source": 3, "U_plus": step(3)[0],
              "negative_source": -3, "U_negative": step(-3)[0]}
    print(json.dumps({"status": "PROVED scoped power/guard interfaces; universal ROOT OPEN",
                      "checks": CHECKS, "odd_primes_through": 101,
                      "square_pushforward_mod3": tables[3],
                      "square_pushforward_mod4": tables[4],
                      "square_pushforward_mod24": tables[24],
                      "reduced_words_cost_le16": reduced_count,
                      "ordered_false_squared_anchor_pairs": exact_false,
                      "actual_guard_hostile": witness,
                      "signed_hostile": signed,
                      "scope": "No new source coverage or ROOT discovery; prime-race counts in separate companion"},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
