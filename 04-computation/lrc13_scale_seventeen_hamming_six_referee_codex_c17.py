#!/usr/bin/env python3
"""Independent algebraic-CRT referee for the scale-seventeen H6 face.

The proof-bearing carrier is the labelled owner-local reachable-mask family.
Scalar capacity first leaves the two quadratic color classes in F_13^*, then
every one of their owner projections misses a sheet.  Tournament rankings are
reported only as deliberately lossy telemetry.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm, prod


P = 13
C = 17
FULL = (1 << C) - 1
LABELS = tuple(range(1, P))
DIVISORS = (1, 17)
STATES = ((1, 0),) + tuple((17, unit) for unit in range(1, 17))
STATE_INDICES = {
    divisor: tuple(
        state for state, (order, _unit) in enumerate(STATES)
        if order == divisor
    )
    for divisor in DIVISORS
}
QR = (1, 3, 4, 9, 10, 12)
NQR = (2, 5, 6, 7, 8, 11)
PROJECTIVE_CYCLE = (
    frozenset((1, 12)),
    frozenset((2, 11)),
    frozenset((4, 9)),
    frozenset((5, 8)),
    frozenset((3, 10)),
    frozenset((6, 7)),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def algebraic_crt_base(label: int, state: int) -> int:
    """Solve x=D*label (mod 13), x=e (mod D) without CRT search."""
    order, unit = STATES[state]
    coefficient = (label - unit * pow(order, -1, P)) % P
    answer = unit + order * coefficient
    require(answer % P == order * label % P, "mod-thirteen CRT mismatch")
    require(answer % order == unit % order, "effective-order CRT mismatch")
    require(0 <= answer < P * order, "CRT representative outside range")
    return answer


def sheet_mask(label: int, state: int, owner: int) -> int:
    order, _unit = STATES[state]
    base = algebraic_crt_base(label, state)
    owner_inverse = pow(owner, -1, P)
    return sum(
        1 << sheet
        for sheet in range(C)
        if -order
        < centered(base * (owner_inverse + P * sheet), P * order)
        <= order
    )


MASK = {
    (provider, state, owner): sheet_mask(provider, state, owner)
    for provider in LABELS
    for state in range(len(STATES))
    for owner in LABELS
}
MASK_DIGEST = sha256(
    b"".join(
        MASK[provider, state, owner].to_bytes(3, "little")
        for provider in LABELS
        for state in range(len(STATES))
        for owner in LABELS
    )
).hexdigest()


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, P) % P


CARDINALITY = {}
for divisor in DIVISORS:
    for provider_owner_ratio in LABELS:
        sizes = {
            MASK[provider_owner_ratio, state, 1].bit_count()
            for state in STATE_INDICES[divisor]
        }
        require(len(sizes) == 1, "mask cardinality depends on unit")
        CARDINALITY[divisor, provider_owner_ratio] = sizes.pop()

for provider in LABELS:
    for owner in LABELS:
        for divisor in DIVISORS:
            sizes = {
                MASK[provider, state, owner].bit_count()
                for state in STATE_INDICES[divisor]
            }
            require(
                sizes == {CARDINALITY[divisor, ratio(provider, owner)]},
                "provider-owner ratio reduction failed",
            )


ORDER_WORDS = tuple(
    word
    for word in product(DIVISORS, repeat=6)
    if all(
        lcm(*(word[:omitted] + word[omitted + 1 :])) == C
        for omitted in range(6)
    )
)
STATE_WORDS_PER_SUPPORT = sum(
    prod(len(STATE_INDICES[divisor]) for divisor in word)
    for word in ORDER_WORDS
)
ORDER_DIGEST = sha256(
    b"".join(bytes(0 if divisor == 1 else 1 for divisor in word)
             for word in ORDER_WORDS)
).hexdigest()


def literal_state_digest() -> tuple[int, str]:
    count = 0
    digest = sha256()
    buffer = bytearray()
    for word in ORDER_WORDS:
        fibres = tuple(STATE_INDICES[divisor] for divisor in word)
        for states in product(*fibres):
            count += 1
            buffer.extend(states)
            if len(buffer) >= 1 << 20:
                digest.update(buffer)
                buffer.clear()
    digest.update(buffer)
    return count, digest.hexdigest()


def scalar_capacity(labels: tuple[int, ...], orders: tuple[int, ...]) -> bool:
    return all(
        sum(
            CARDINALITY[divisor, ratio(provider, owner)]
            for provider, divisor in zip(labels, orders)
        )
        >= C
        for owner in labels
    )


def owner_reachable(
    labels: tuple[int, ...], orders: tuple[int, ...], owner: int
) -> set[int]:
    reachable = {0}
    for provider, divisor in zip(labels, orders):
        options = {
            MASK[provider, state, owner]
            for state in STATE_INDICES[divisor]
        }
        reachable = {
            partial | option
            for partial in reachable
            for option in options
        }
    return reachable


def multiply_support(
    labels: tuple[int, ...], multiplier: int
) -> tuple[int, ...]:
    return tuple(sorted(multiplier * label % P for label in labels))


def projective_multiplicities(labels: tuple[int, ...]) -> tuple[int, ...]:
    support = set(labels)
    return tuple(len(pair & support) for pair in PROJECTIVE_CYCLE)


def histogram_text(counter: Counter) -> str:
    return " ".join(f"{key}:{counter[key]}" for key in sorted(counter))


def main() -> None:
    require(len(STATES) == 17, "literal state alphabet mismatch")
    require(STATE_INDICES[1] == (0,), "order-one state mismatch")
    require(
        {STATES[state][1] for state in STATE_INDICES[17]}
        == set(range(1, 17)),
        "order-seventeen unit grammar mismatch",
    )
    require(len(ORDER_WORDS) == 57, "hereditary order-word census mismatch")
    require(
        STATE_WORDS_PER_SUPPORT == 24_137_472,
        "weighted state-word census mismatch",
    )
    literal_count, state_digest = literal_state_digest()
    require(
        literal_count == STATE_WORDS_PER_SUPPORT,
        "literal and weighted state censuses disagree",
    )
    require(
        MASK_DIGEST
        == "469e954c1cffe3827414701864f3fb0c944883c085883c27fa0488737ee809dd",
        "algebraic CRT mask-table digest mismatch",
    )
    require(
        ORDER_DIGEST
        == "c78fafe41717a9194f5f329036eee7689a3913d98d78b010df57cfc3af3e3144",
        "order-grammar digest mismatch",
    )
    require(
        state_digest
        == "9c0a8805ccee6f09b21e4979d176f89ee9100c337053b95b4522b52e70057dd4",
        "literal-state-grammar digest mismatch",
    )

    supports = tuple(combinations(LABELS, 6))
    require(len(supports) == 924, "support census mismatch")
    scalar_bank = tuple(
        (labels, orders)
        for labels in supports
        for orders in ORDER_WORDS
        if scalar_capacity(labels, orders)
    )
    expected_bank = ((QR, (17,) * 6), (NQR, (17,) * 6))
    require(scalar_bank == expected_bank, "scalar QR/NQR classification mismatch")
    require(
        projective_multiplicities(QR) == (2, 0, 2, 0, 2, 0)
        and projective_multiplicities(NQR) == (0, 2, 0, 2, 0, 2),
        "projective alternating multiplicity classification mismatch",
    )
    require(
        {multiply_support(QR, multiplier) for multiplier in LABELS}
        == {QR, NQR},
        "multiplication orbit mismatch",
    )

    scalar_digest = sha256(
        "\n".join(
            f"{','.join(map(str, labels))}|{','.join(map(str, orders))}"
            for labels, orders in scalar_bank
        ).encode()
    ).hexdigest()
    require(
        scalar_digest
        == "74410e072b3341df3545ec0f023d3098fd822b40b52eb7fb65bb20c4bf352826",
        "scalar-bank digest mismatch",
    )

    feasible_contexts = Counter()
    maximum_union = Counter()
    reachable_count = Counter()
    maximum_mask_count = Counter()
    missing_sheet = Counter()
    owner_digest = sha256()
    for labels, orders in scalar_bank:
        feasible_owners = 0
        summaries = []
        for owner in labels:
            reachable = owner_reachable(labels, orders, owner)
            maximum = max(mask.bit_count() for mask in reachable)
            maxima = sorted(
                mask for mask in reachable if mask.bit_count() == maximum
            )
            feasible = FULL in reachable
            require(feasible == (maximum == C), "owner threshold mismatch")
            require(maximum == 16, "owner maximum is not exactly sixteen")
            require(len(reachable) == 38_540, "reachable-mask census mismatch")
            require(len(maxima) == 16, "maximum-mask census mismatch")
            for mask in maxima:
                missing = FULL ^ mask
                require(missing.bit_count() == 1, "maximum mask deficit mismatch")
                missing_sheet[missing.bit_length() - 1] += 1
            feasible_owners += feasible
            summaries.append((feasible, maximum))
            maximum_union[maximum] += 1
            reachable_count[len(reachable)] += 1
            maximum_mask_count[len(maxima)] += 1
            owner_digest.update(bytes(labels + (owner, maximum)))
            for mask in sorted(reachable):
                owner_digest.update(mask.to_bytes(3, "little"))
        feasible_contexts[feasible_owners] += 1
        require(len(set(summaries)) == 1, "owner summaries unexpectedly differ")

    require(feasible_contexts == Counter({0: 2}), "feasibility census mismatch")
    require(maximum_union == Counter({16: 12}), "maximum-union census mismatch")
    require(reachable_count == Counter({38_540: 12}), "reachability census mismatch")
    require(maximum_mask_count == Counter({16: 12}), "maximum-mask census mismatch")
    require(
        tuple(missing_sheet[sheet] for sheet in range(C))
        == (12, 11, 11, 11, 12, 11, 11, 11, 12, 11, 11, 11, 12, 11, 11, 11, 12),
        "missing-sheet profile mismatch",
    )
    owner_digest_hex = owner_digest.hexdigest()
    require(
        owner_digest_hex
        == "7dc7010e6e75741a54b3e3722901d60fa6a3bc526bc4c81c55fe5d48c53f6040",
        "owner-reachable-bank digest mismatch",
    )

    print("scale-seventeen independent algebraic-CRT Python referee")
    print("divisor grammar 1,17; literal states 17")
    print(f"supports {len(supports)}; hereditary order words {len(ORDER_WORDS)}; "
          f"labelled order contexts {len(supports) * len(ORDER_WORDS)}")
    print(f"state words/support {STATE_WORDS_PER_SUPPORT}; "
          f"raw labelled states {len(supports) * STATE_WORDS_PER_SUPPORT}")
    print(f"mask-table SHA256 {MASK_DIGEST}")
    print(f"order-grammar SHA256 {ORDER_DIGEST}")
    print(f"literal-state-grammar SHA256 {state_digest}")
    print("scalar contexts 2; all-D17; supports "
          f"QR {QR} and NQR {NQR}")
    print("projective multiplicities 2,0,2,0,2,0 and 0,2,0,2,0,2; "
          "one multiplication orbit of size 2")
    print(f"scalar-bank SHA256 {scalar_digest}")
    print("owner-local rows 12; feasible 0; maximum-union histogram 16:12")
    print("reachable-mask count histogram " + histogram_text(reachable_count))
    print("maximum-mask count histogram " + histogram_text(maximum_mask_count))
    print("missing-sheet counts among maximum masks " + histogram_text(missing_sheet))
    print(f"owner-reachable-bank SHA256 {owner_digest_hex}")
    print("global literal unit fibres 0; common-sheet survivors 0")
    print("tournament pair observable exact owner summaries; all fifteen pairs tie; "
          "coordinate Hamiltonian path gives two transitive completions with scores "
          "0,1,2,3,4,5, cycles 0, SCCs 6, paths 1")
    print("challenged vertices: signed projective classes classify scalar support; "
          "labelled owner maxima prove the terminal deficit; the oriented tournament "
          "and provider/residue/sheet quotients lose the absolute threshold")


if __name__ == "__main__":
    main()
