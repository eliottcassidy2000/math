"""Exact checks for THM-2041 and the HYP-8800 LRC transfer audit.

The mathematical theorem is proved in THM-2041.  This script independently
checks its finite cyclic group-algebra identities, records the bad-prime
failure, and compares the proof carriers used by the LRC transfer.

Tournament Analysis
-------------------
Vertices are proof carriers, not runners.  The pairwise observable is the
lexicographic tuple

    (certifies_strict_LRC_predicate, retains_endpoint_owner,
     preserves_nonzero_layer, preserves_whole_tied_layer,
     Frobenius_stable, linear_projector, compression).

The higher tuple wins; declaration order is the tie Hamiltonian path.  This
quotient preserves the proof obligations relevant to seed-and-amplify and
destroys the actual runner locations, arc endpoints, and safe-set geometry.
The challenged assumption is that a Frobenius-stable arithmetic scalar is
automatically a witness: it is not, unless the strict LRC predicate or its
endpoint-labelled dual certificate is retained.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import gcd


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def mobius(n: int) -> int:
    if n == 1:
        return 1
    parity = 0
    q = n
    prime = 2
    while prime * prime <= q:
        if q % prime == 0:
            q //= prime
            parity += 1
            if q % prime == 0:
                return 0
            while q % prime == 0:
                q //= prime
        prime += 1
    if q > 1:
        parity += 1
    return -1 if parity % 2 else 1


def ramanujan(d: int, n: int) -> int:
    return sum(e * mobius(d // e) for e in divisors(gcd(d, n)))


def primes_through(limit: int) -> list[int]:
    ans = []
    for n in range(2, limit + 1):
        if all(n % q for q in range(2, int(n**0.5) + 1)):
            ans.append(n)
    return ans


def cyclic_convolution(a: list[int], b: list[int], prime: int) -> list[int]:
    d = len(a)
    out = [0] * d
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[(i + j) % d] = (out[(i + j) % d] + ai * bj) % prime
    return out


def exact_period_idempotent(period: int, order: int, prime: int) -> list[int]:
    """Project C_period onto the characters of exact order ``order``."""
    assert period % order == 0 and gcd(period, prime) == 1
    inv_period = pow(period, -1, prime)
    return [
        (inv_period * ramanujan(order, n)) % prime for n in range(period)
    ]


def frobenius_group_algebra(a: list[int], prime: int) -> list[int]:
    """p-th power in F_p[C_d], written in the standard group basis."""
    d = len(a)
    out = [0] * d
    for exponent, coefficient in enumerate(a):
        out[(prime * exponent) % d] = (
            out[(prime * exponent) % d] + pow(coefficient, prime, prime)
        ) % prime
    return out


def integer_energy(f: list[int], order: int) -> int:
    period = len(f)
    return sum(
        f[a] * f[b] * ramanujan(order, a - b)
        for a in range(period)
        for b in range(period)
    )


def unit_dilate(f: list[int], unit: int) -> list[int]:
    d = len(f)
    assert gcd(unit, d) == 1
    out = [0] * d
    for index, value in enumerate(f):
        out[(unit * index) % d] = value
    return out


def exact_checks() -> dict[str, int]:
    counts: Counter[str] = Counter()
    for d in range(1, 31):
        for order in divisors(d):
            counts["period_order_pairs"] += 1
            for prime in primes_through(47):
                if gcd(d, prime) != 1:
                    continue
                e = exact_period_idempotent(d, order, prime)
                assert cyclic_convolution(e, e, prime) == e
                assert frobenius_group_algebra(e, prime) == e
                counts["good_period_order_prime_triples"] += 1

                tests = [
                    [(i * i + 3 * i + 1) % prime for i in range(d)],
                    [((i + 2) * (d - i + 1)) % prime for i in range(d)],
                    [1 if gcd(i, d) == 1 else 0 for i in range(d)],
                ]
                for f in tests:
                    projected = cyclic_convolution(e, f, prime)
                    frobenius_projected = frobenius_group_algebra(projected, prime)
                    projected_frobenius = cyclic_convolution(
                        e, frobenius_group_algebra(f, prime), prime
                    )
                    assert frobenius_projected == projected_frobenius
                    assert (any(projected) == any(frobenius_projected))
                    assert integer_energy(
                        unit_dilate(f, prime), order
                    ) == integer_energy(f, order)
                    counts["projector_commutation_tests"] += 1
                    counts["nonzero_preservation_tests"] += 1
                    counts["ramanujan_energy_dilation_tests"] += 1

    for prime in (2, 3, 5, 7, 11):
        # In F_p[C_p], h=u-1 is nonzero but h^p=u^p-1=0.
        h = [(-1) % prime, 1] + [0] * (prime - 2)
        assert any(h)
        assert not any(frobenius_group_algebra(h, prime))
        counts["bad_prime_noninjective_witnesses"] += 1
    return dict(counts)


CARRIERS = {
    # strict, endpoint, nonzero, whole layer, Frobenius, linear, compression
    "endpoint_labelled_primitive_packet": (1, 1, 1, 1, 1, 1, 0),
    "primitive_safe_residue_count": (1, 0, 1, 0, 1, 0, 1),
    "primitive_period_idempotent": (0, 0, 1, 1, 1, 1, 1),
    "ramanujan_shell_energy": (0, 0, 0, 1, 1, 0, 1),
    "raw_ramanujan_trace": (0, 0, 0, 0, 1, 0, 1),
    "raw_denominator_blockedness": (0, 0, 0, 0, 0, 0, 1),
}


def tournament_fingerprint() -> dict[str, object]:
    names = list(CARRIERS)
    wins = {name: set() for name in names}
    for left, right in combinations(names, 2):
        if CARRIERS[left] >= CARRIERS[right]:
            wins[left].add(right)
        else:
            wins[right].add(left)
    score_hist = dict(sorted(Counter(map(len, wins.values())).items()))
    cycles = 0
    for a, b, c in combinations(names, 3):
        if (
            (b in wins[a] and c in wins[b] and a in wins[c])
            or (c in wins[a] and b in wins[c] and a in wins[b])
        ):
            cycles += 1
    order = sorted(names, key=lambda name: (CARRIERS[name], -names.index(name)), reverse=True)
    assert all(order[j + 1] in wins[order[j]] for j in range(len(order) - 1))
    return {
        "score_hist": score_hist,
        "directed_3cycles": cycles,
        "scc_sizes": [1] * len(names),
        "hamiltonian_path_count": 1,
        "tie_hamiltonian_path": " > ".join(order),
    }


def main() -> None:
    print("THM-2041 EXACT-PERIOD FROBENIUS AUDIT")
    for key, value in exact_checks().items():
        print(f"{key}={value}")
    print("good_prime_verdict=PASS")
    print("bad_prime_boundary=Frobenius is noninjective when p divides the period")
    print("scope=LRC seed preservation only; no seed-existence or lonely-time theorem")
    print("TOURNAMENT ANALYSIS")
    for key, value in tournament_fingerprint().items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
