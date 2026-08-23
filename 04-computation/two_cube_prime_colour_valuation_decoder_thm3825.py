#!/usr/bin/env python3
"""Exact assertion-free probe for THM-3825's prime-colour decoder."""

from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path

from sympy import factorint, primerange


CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def factor(n: int) -> dict[int, int]:
    return {int(p): int(e) for p, e in factorint(n).items()}


def inert_supported(n: int) -> bool:
    return n >= 1 and all(p % 3 == 2 for p in factor(n))


def admissible_shell(s: int) -> bool:
    fs = factor(s)
    return s >= 2 and all(p % 3 == 2 and e <= 2 for p, e in fs.items())


def decode(n: int) -> tuple[int, int, int, int, int] | None:
    """Return (g,p,q,s,Q) on the restricted image, otherwise None."""
    if n <= 0:
        return None
    fn = factor(n)
    g = 1
    s = 1
    for ell, exponent in fn.items():
        if ell % 3 == 2:
            quotient, remainder = divmod(exponent, 3)
            g *= ell**quotient
            s *= ell**remainder
    denominator = g**3 * s
    if denominator == 0 or n % denominator:
        return None
    quadratic = n // denominator

    # In a valid primitive shell, the quadratic factor has only 1 mod 3
    # prime colours.  This redundant filter supplies a cheap hostile gate.
    if any(ell % 3 != 1 for ell in factor(quadratic)):
        return None

    numerator = 4 * quadratic - s * s
    if numerator <= 0 or numerator % 3:
        return None
    discriminant = numerator // 3
    root = math.isqrt(discriminant)
    if root * root != discriminant or (s - root) % 2:
        return None
    p = (s - root) // 2
    q = (s + root) // 2
    if not (1 <= p < q and math.gcd(p, q) == 1):
        return None
    if not admissible_shell(s) or not inert_supported(g):
        return None
    if quadratic != p * p - p * q + q * q:
        return None
    if n != (g * p) ** 3 + (g * q) ** 3:
        return None
    return g, p, q, s, quadratic


def primes_through(limit: int) -> list[int]:
    return [int(p) for p in primerange(2, limit + 1)]


ORDINARY_PRIMES = primes_through(2000)
INERT_PRIMES = [p for p in primes_through(10000) if p % 3 == 2]
ORDINARY_INDEX = {p: i for i, p in enumerate(ORDINARY_PRIMES)}
INERT_INDEX = {p: i for i, p in enumerate(INERT_PRIMES)}


def encode_tag(tag: int) -> int:
    if tag < 1:
        raise ValueError("tags are positive integers")
    encoded = 1
    for p, exponent in factor(tag).items():
        encoded *= INERT_PRIMES[ORDINARY_INDEX[p]] ** exponent
    return encoded


def decode_tag(encoded: int) -> int | None:
    tag = 1
    for p, exponent in factor(encoded).items():
        index = INERT_INDEX.get(p)
        if index is None:
            return None
        tag *= ORDINARY_PRIMES[index] ** exponent
    return tag


def main() -> None:
    ratios: list[tuple[int, int, int]] = []
    for s in range(3, 357):
        if not admissible_shell(s):
            continue
        for p in range(1, (s + 1) // 2):
            q = s - p
            if p < q and math.gcd(p, q) == 1:
                ratios.append((p, q, s))

    require(len(ratios) == 5855, "admissible ratio census changed")
    require(len({s for _, _, s in ratios}) == 94,
            "admissible shell census changed")

    decoder_rows: list[dict[str, int]] = []
    for g in (1, 2, 5, 25):
        require(inert_supported(g), "scale control left inert support")
        for p, q, s in ratios:
            n = (g * p) ** 3 + (g * q) ** 3
            quadratic = p * p - p * q + q * q
            require(decode(n) == (g, p, q, s, quadratic),
                    "prime-colour decoder failed on admissible state")
        decoder_rows.append({"g": g, "states": len(ratios)})

    # The same inert prime may occur in both scale and shell: 5^7=5^(3*2+1).
    shared = 65 * 25**3
    require(factor(shared)[5] == 7, "shared-prime valuation control changed")
    require(decode(shared) == (25, 1, 4, 5, 13),
            "shared scale/shell quotient-remainder failed")

    # Sharp scope hostiles.
    require(decode(1729) is None, "split-shell collision was accepted")
    require(decode(515375) is None, "exponent-three shell collision was accepted")
    split_scale = 7**3 + 28**3
    require(split_scale == 65 * 7**3, "split-scale hostile changed")
    require(decode(split_scale) is None, "split scale was silently decoded")

    # Direct restricted-image test through a hostile finite universe.
    accepted = 0
    for n in range(2, 20001):
        packet = decode(n)
        if packet is None:
            continue
        g, p, q, s, quadratic = packet
        require(n == (g * p) ** 3 + (g * q) ** 3,
                "accepted value failed coordinate reconstruction")
        require(s == p + q and quadratic == p * p - p * q + q * q,
                "accepted value failed shell reconstruction")
        accepted += 1

    # Every positive integer tag is transported into inert scale.  The fixed
    # primitive shell (1,4) then gives N(tag)=65*Enc(tag)^3.
    tagged = 0
    for tag in range(1, 501):
        encoded = encode_tag(tag)
        require(inert_supported(encoded), "tag encoding lost prime colour")
        require(decode_tag(encoded) == tag, "tag inverse failed")
        address = 65 * encoded**3
        require(decode(address) == (encoded, 1, 4, 5, 13),
                "tagged two-cube address failed")
        require(address == encoded**3 + (4 * encoded) ** 3,
                "tagged address lost two-cube form")
        tagged += 1

    source = Path(__file__).read_text(encoding="utf-8")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source))),
            "assert statement found")

    semantic = {
        "universe": "p<q,gcd=1,p+q<=356,inert shell exponents<=2",
        "ratios": len(ratios),
        "shells": len({s for _, _, s in ratios}),
        "scales": [row["g"] for row in decoder_rows],
        "decoder": "v_inert=3*v_scale+v_shell;Delta=(4Q-s^2)/3",
        "tag": "N(m)=65*Enc(m)^3, ordinary-prime index -> inert-prime index",
        "hostiles": [1729, 515375, split_scale],
    }
    semantic_hash = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("object=THM-3825-prime-colour-valuation-two-cube-decoder")
    print("universe=p<q;gcd=1;s=p+q<=356;inert-shell-exponents<=2;inert-scale")
    print(f"admissible_shells={len({s for _, _, s in ratios})}")
    print(f"admissible_ratios={len(ratios)}")
    print(f"all_scale_replays={len(ratios) * len(decoder_rows)}")
    print("decoder=v_ell(M)=3*v_ell(g)+v_ell(s);quotient=g;remainder=s")
    print("square_test=Delta=(4Q-s^2)/3;pair=(s-sqrtDelta,s+sqrtDelta)/2")
    print(f"accepted_values_through_20000={accepted}")
    print(f"tagged_addresses={tagged}")
    print("shared_prime_control=5^7=>scale_5^2,shell_5")
    print("hostiles=1729:split-shell;515375:exponent-three;22295:split-scale")
    print(f"CHECKS={CHECKS}")
    print(f"semantic_sha256={semantic_hash}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
