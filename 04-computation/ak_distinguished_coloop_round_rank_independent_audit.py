#!/usr/bin/env python3
"""Independent SymPy audit of the four THM-TBD arithmetic-Kakeya profiles.

The frozen certificate constructors are imported from the canonical audit,
but every rank, firing/coloop test, and topology defect below is recomputed
with SymPy.  In particular, delta is maximized over every row subset rather
than inferred from a witness of the requested size.
"""

from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sys

from sympy import Matrix


ROOT = Path(__file__).resolve().parents[1]
COMPUTATION = ROOT / "04-computation"
CANONICAL = COMPUTATION / "ak_distinguished_coloop_round_rank_audit.py"
EXPECTED_CANONICAL_SHA256 = (
    "b9212c6f27d839b591ab675149f3c2168ceb51d10cd9d828a670eeb26ad3b90c"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


payload = CANONICAL.read_bytes().replace(b"\r\n", b"\n")
require(
    sha256(payload).hexdigest() == EXPECTED_CANONICAL_SHA256,
    "canonical round-rank audit changed",
)
sys.path.insert(0, str(COMPUTATION))

import ak_distinguished_coloop_round_rank_audit as canonical  # noqa: E402


def rank(rows, ncols):
    return Matrix(rows).rank() if rows and ncols else 0


def validate_paid_graph(cert):
    require(
        len(cert.generators) <= cert.cost,
        f"{cert.name}: more generator occurrences than paid slots",
    )
    for row_index, sparse in enumerate(cert.generators):
        labels = tuple(sparse.values())
        require(
            len(labels) in (1, 2),
            f"{cert.name}: row {row_index} is neither seed nor edge",
        )
        require(
            all(label != (0, 0) and sum(label) != 0 for label in labels),
            f"{cert.name}: row {row_index} has an illegal label",
        )
        if len(labels) == 2:
            require(
                tuple(labels[0][i] + labels[1][i] for i in range(2))
                == (0, 0),
                f"{cert.name}: row {row_index} is not a signed edge",
            )


def coloop_firing_set(order, arows):
    ncols = 2 * len(order)
    full_rank = rank(arows, ncols)
    fired = []
    for vertex_index, vertex in enumerate(order):
        deleted = 2 * vertex_index + 1
        columns = tuple(j for j in range(ncols) if j != deleted)
        reduced_rows = [[row[j] for j in columns] for row in arows]
        if full_rank - rank(reduced_rows, ncols - 1) == 1:
            fired.append(vertex)
    return tuple(fired)


def exact_delta(brows, u):
    best = None
    row_count = len(brows)
    for size in range(row_count + 1):
        for subset in combinations(range(row_count), size):
            value = size - 2 * rank([brows[i] for i in subset], u)
            candidate = (value, -size, tuple(-i for i in subset))
            if best is None or candidate > best[0]:
                best = (candidate, subset)
    value = best[0][0]
    subset = best[1]
    return value, subset, rank([brows[i] for i in subset], u)


def bridge_count(brows, u):
    full_rank = rank(brows, u)
    return sum(
        rank(brows[:i] + brows[i + 1 :], u) < full_rank
        for i in range(len(brows))
    )


def audit(cert):
    validate_paid_graph(cert)
    live = set(range(cert.n)) - set(cert.t0)
    rounds = []
    while live:
        order, brows, arows = cert.matrices(live)
        u = len(order)
        e = len(arows)
        rb = rank(brows, u)
        r = rank(arows, 2 * u)
        sigma = e - r
        c = e - u
        z = bridge_count(brows, u)
        delta, witness, witness_rank = exact_delta(brows, u)
        fired = coloop_firing_set(order, arows)

        require(rb == u, f"{cert.name}: incidence block not pinned")
        require(delta == sigma, f"{cert.name}: delta != sigma")
        require(
            sigma >= max(0, 2 * c - (e - z)),
            f"{cert.name}: bicirculation/bridge bound failed",
        )
        require(fired, f"{cert.name}: forcing stalled")
        require(
            fired == canonical.fired_vertices(order, arows),
            f"{cert.name}: SymPy and custom RREF firing sets differ",
        )

        next_live = live - set(fired)
        _, _, next_arows = cert.matrices(next_live)
        r_next = rank(next_arows, 2 * len(next_live))
        p = r - len(fired) - r_next
        require(0 <= p <= len(fired), f"{cert.name}: invalid p")
        rounds.append(
            (
                u,
                r,
                fired,
                p,
                sigma,
                len(witness),
                witness_rank,
                z,
            )
        )
        live = next_live
    return tuple(rounds)


def main():
    strict = canonical.strict_13_over_7()
    records = (
        canonical.QuotientCertificate.from_strict("strict-13/7", strict),
        canonical.QuotientCertificate.from_mode3(
            "per-suffix-7/4", canonical.per_suffix_7_over_4()
        ),
        canonical.QuotientCertificate.from_mode3(
            "quotient-12/7", canonical.quotient_12_over_7()
        ),
        canonical.QuotientCertificate.from_mode3(
            "quotient-9/5", canonical.quotient_9_over_5()
        ),
    )
    expected = {
        "strict-13/7": ((7, 13, (3, 4), 1, 0), (5, 10, (0, 2, 5, 6, 7), 5, 1)),
        "per-suffix-7/4": (
            (8, 14, (5,), 0, 0),
            (7, 13, (0, 1, 4, 6), 3, 1),
            (3, 6, (2, 3, 7), 3, 1),
        ),
        "quotient-12/7": (
            (7, 12, (4,), 0, 0),
            (6, 11, (0, 3, 5), 2, 1),
            (3, 6, (1, 2, 6), 3, 1),
        ),
        "quotient-9/5": (
            (5, 9, (2, 3), 1, 0),
            (3, 6, (0, 1, 4), 3, 2),
        ),
    }
    for cert in records:
        rounds = audit(cert)
        projection = tuple(round_data[:5] for round_data in rounds)
        require(projection == expected[cert.name], f"{cert.name}: profile changed")
        print(f"{cert.name}: {rounds}")
    print("INDEPENDENT_SYMPY_RANK_DEFECT_AUDIT_OK")


if __name__ == "__main__":
    main()
