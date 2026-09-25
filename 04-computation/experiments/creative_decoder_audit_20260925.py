"""Independent finite audit of the guarded common-future proof compiler.

Source note: 05-knowledge/results/creative_decoder_20260925.md.
Subject: creative_decoder_20260925.py in this directory.
This audits finite universes and verifier guards, not universal termination.
The native implementation uses repeated integer division, not affine carries.
"""

from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def native_step(n, sign):
    check(type(n) is int and n > 0 and n % 2 == 1, "native source domain")
    n = 3 * n + sign
    exponent = 0
    while n % 2 == 0:
        n //= 2
        exponent += 1
    return n, exponent


def native_word(n, length, sign):
    word = []
    for _ in range(length):
        n, exponent = native_step(n, sign)
        word.append(exponent)
    return n, tuple(word)


def native_inverse(n, word, sign):
    for exponent in reversed(word):
        numerator = (2 ** exponent) * n - sign
        if numerator % 3:
            return None
        n = numerator // 3
        if n <= 0 or n % 2 == 0:
            return None
    return n


def main():
    path = Path(__file__).with_name("creative_decoder_20260925.py")
    spec = spec_from_file_location("audited_creative_decoder", path)
    decoder = module_from_spec(spec)
    spec.loader.exec_module(decoder)

    cases = forward_accepts = inverse_accepts = 0
    for sign in (-1, 1):
        for length in range(1, 5):
            for exponents in product(range(1, 5), repeat=length):
                word = decoder.Word.make(exponents)
                for n in range(1, 256, 2):
                    endpoint, actual = native_word(n, length, sign)
                    forward = word.forward(n, sign)
                    check(forward == (endpoint if actual == exponents else None),
                          "forward guard iff exact native exponent word")
                    inverse = word.inverse(n, sign)
                    check(inverse == native_inverse(n, exponents, sign),
                          "inverse guard iff independently reversed edges")
                    if forward is not None:
                        check(word.inverse(forward, sign) == n, "forward/inverse round trip")
                        forward_accepts += 1
                    if inverse is not None:
                        check(native_word(inverse, length, sign) == (n, exponents),
                              "inverse source has exact native exponent word")
                        inverse_accepts += 1
                    cases += 1
    check(cases == 87040, "declared finite universe")
    print("WORD GUARDS: 87040 cases; signs=-1,+1; lengths=1..4; exponents=1..4; odd inputs=1..255")
    print("Accepted forward / inverse:", forward_accepts, inverse_accepts)

    for sign in (-1, 1):
        joins = 0
        for n in range(1, 10000, 2):
            for join in decoder.base_joins(n, sign):
                check(0 < join.y < join.x == n and join.y % 2 == 1,
                      "base join integer decrease")
                left = native_word(join.x, join.a, sign)[0]
                right = native_word(join.y, join.b, sign)[0]
                check(left == right, "base join independent common future")
                join.verify(sign)
                joins += 1
        print("BASE JOINS: sign", sign, "; odd inputs<10000; verified joins", joins)

    hostiles = [
        ("wrong sign", lambda: decoder.Join(5, 1, 0, 1, "hostile").verify()),
        ("circular rank", lambda: decoder.Join(7, 7, 0, 0, "hostile").verify()),
        ("lost endpoint parity", lambda: decoder.Join(31, 26, 0, 2, "hostile").verify()),
        ("negative clock", lambda: decoder.Join(5, 3, -1, 1, "hostile").verify()),
        ("even inverse target", lambda: decoder.Word.make((1,)).inverse(2)),
        ("forged carry metadata", lambda: decoder.Word((1,), 2, 1)),
        ("invalid forward sign", lambda: decoder.Word.make((1,)).forward(1, 0)),
        ("invalid inverse sign", lambda: decoder.Word.make((1,)).inverse(1, 0)),
    ]
    for label, operation in hostiles:
        try:
            operation()
        except RuntimeError:
            pass
        else:
            raise RuntimeError("hostile verifier accepted: " + label)
    print("VERIFIER HOSTILES: rejected", len(hostiles), "of", len(hostiles))

    for seed, period in ((1, 1), (5, 2), (17, 7)):
        check(native_word(seed, period, -1)[0] == seed, "minus cycle closure")
        check(not list(decoder.base_joins(seed, -1)), "minus cycle minimum base obstruction")
    print("MINUS CONTROLS: distinct cycle minima 1,5,17 retain no strict base joins")
    print("ALL FINITE AUDIT CHECKS PASS; universal root coverage remains OPEN")


if __name__ == "__main__":
    main()
