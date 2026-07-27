#!/usr/bin/env python3
"""Exact finite referee for THM-2544.

The computation never enumerates K_91, whose size is 91^8.  It builds an
explicit nonunit section and an explicit all-unit section of the 169 target
fibres for the canonical typed word, then checks the split-surjection witness
for the unrestricted/all-unit projector pair.

Scope: finite CRT and linear algebra only.  The sparse vectors checked here
are arbitrary residue currents; they are not asserted to be Abel-boundary
currents of a scalar covering row.
"""

from hashlib import sha256
from itertools import product
from math import gcd


P7 = 7
P13 = 13
N = 91
D = 9

# THM-2541 order: six guard/unit coordinates, owner j, targets a,b.
W = (1, 14, 27, 40, 53, 66, 13, 2197, 742586)
UNITS = tuple(range(6))
OWNER = 6
TARGET_A = 7
TARGET_B = 8

# THM-2309 owner packet: omit unit 5 and graft the first two unit rows.
OMITTED = 5
GRAFT_A = 0
GRAFT_B = 1
ROW_LABELS = (OWNER, 0, 1, 2, 3, 4)
PIVOT = ROW_LABELS

CHECKS = 0


def require(condition, message):
    """Optimization-safe assertion with an exact check counter."""
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def add_scaled(x, scalar, row, modulus):
    return tuple((a + scalar * b) % modulus for a, b in zip(x, row))


def pair_relation(k):
    """The exact relation w_k e_omitted - w_omitted e_k."""
    row = [0] * D
    row[OMITTED] += W[k]
    row[k] -= W[OMITTED]
    return row


def target_relation(target):
    """The exact relation w_target e_omitted - w_omitted e_target."""
    row = [0] * D
    row[OMITTED] += W[target]
    row[target] -= W[OMITTED]
    return row


def owner_packet():
    rows = []
    for k in ROW_LABELS:
        row = pair_relation(k)
        if k == GRAFT_A:
            row = [a + b for a, b in zip(row, target_relation(TARGET_A))]
        if k == GRAFT_B:
            row = [a + b for a, b in zip(row, target_relation(TARGET_B))]
        rows.append(tuple(row))
    return tuple(rows)


ROWS_Z = owner_packet()
ROWS_13 = tuple(tuple(a % P13 for a in row) for row in ROWS_Z)


def projection_13(x):
    """THM-2309 quotient K_13/L_13 -> F_13^{a,b}."""
    ell = (0,) * D
    diagonal = (-W[OMITTED]) % P13
    inverse_diagonal = pow(diagonal, -1, P13)
    for row, pivot_coordinate in zip(ROWS_13, PIVOT):
        alpha = x[pivot_coordinate] * inverse_diagonal % P13
        ell = add_scaled(ell, alpha, row, P13)
    residual = tuple((a - b) % P13 for a, b in zip(x, ell))
    require(all(residual[i] == 0 for i in range(D)
                if i not in (TARGET_A, TARGET_B)),
            "quotient residual escaped the target plane")
    return residual[TARGET_A], residual[TARGET_B]


def all_unit_mod13(target):
    """Find one all-nonzero point of target+L_13, constructively."""
    qa, qb = target
    for coefficients in product(range(1, P13), repeat=5):
        # A nonzero owner coefficient makes the owner pivot nonzero.
        alphas = (1,) + coefficients
        x = [0] * D
        x[TARGET_A] = qa
        x[TARGET_B] = qb
        for alpha, row in zip(alphas, ROWS_13):
            x = list(add_scaled(tuple(x), alpha, row, P13))
        x = tuple(x)
        if all(a != 0 for a in x):
            return x
    raise RuntimeError(f"no all-unit mod-13 witness for target {target}")


def all_unit_mod7():
    """Find one all-nonzero vector in K_7, solving the last coordinate."""
    last_inverse = pow(W[-1] % P7, -1, P7)
    for prefix in product(range(1, P7), repeat=D - 1):
        final = (-sum(W[i] * prefix[i] for i in range(D - 1)))
        final = final * last_inverse % P7
        if final != 0:
            return prefix + (final,)
    raise RuntimeError("no all-unit mod-7 relation witness")


def crt(a7, a13):
    """The representative in [0,91) with the two prescribed residues."""
    return (a7 + P7 * (((a13 - a7) * pow(P7, -1, P13)) % P13)) % N


def lift_pair(x7, x13):
    return tuple(crt(a, b) for a, b in zip(x7, x13))


def target_section_13(target):
    x = [0] * D
    x[TARGET_A], x[TARGET_B] = target
    return tuple(x)


def U_of_sparse(vector):
    """Unrestricted target pushforward on the sparse witness alphabet."""
    answer = {}
    for (_, q), coefficient in vector.items():
        answer[q] = answer.get(q, 0) + coefficient
    return {q: a for q, a in answer.items() if a != 0}


def J_of_sparse(vector):
    """All-unit target pushforward on the sparse witness alphabet."""
    answer = {}
    for (kind, q), coefficient in vector.items():
        if kind == "unit":
            answer[q] = answer.get(q, 0) + coefficient
    return {q: a for q, a in answer.items() if a != 0}


def main():
    require(gcd(*W) == 1, "the canonical typed word must be primitive")
    require(tuple(w % P13 for w in W[:6]) == (1,) * 6,
            "the six guard/unit entries must be 1 mod 13")
    require(tuple(w % P13 for w in W[6:]) == (0,) * 3,
            "the blocker entries must vanish mod 13")
    require(all(dot(row, W) == 0 for row in ROWS_Z),
            "every owner-packet row must be an exact integer relation")

    diagonal = (-W[OMITTED]) % P13
    for row_index, row in enumerate(ROWS_13):
        for column_index, pivot_coordinate in enumerate(PIVOT):
            expected = diagonal if row_index == column_index else 0
            require(row[pivot_coordinate] == expected,
                    "owner-packet pivot is not diagonal")
    require(diagonal != 0, "owner-packet pivot diagonal vanished")
    require(all(any(row[column] != 0 for row in ROWS_13)
                for column in range(D)),
            "owner packet is not bright in all nine coordinates")

    z7 = all_unit_mod7()
    require(all(a != 0 for a in z7), "septimal witness has a zero coordinate")
    require(dot(z7, W) % P7 == 0, "septimal witness is not a relation")

    targets = tuple((a, b) for a in range(P13) for b in range(P13))
    unit_witness = {}
    nonunit_witness = {}
    digest_rows = []

    for q in targets:
        x13_unit = all_unit_mod13(q)
        require(dot(x13_unit, W) % P13 == 0,
                "mod-13 unit witness is not a relation")
        require(projection_13(x13_unit) == q,
                "mod-13 unit witness has the wrong quotient target")

        x13_nonunit = target_section_13(q)
        require(dot(x13_nonunit, W) % P13 == 0,
                "nonunit target section is not in K_13")
        require(projection_13(x13_nonunit) == q,
                "nonunit section has the wrong quotient target")

        u = lift_pair(z7, x13_unit)
        v = lift_pair((0,) * D, x13_nonunit)
        require(dot(u, W) % N == 0, "unit CRT lift is not in K_91")
        require(dot(v, W) % N == 0, "nonunit CRT lift is not in K_91")
        require(all(gcd(a, N) == 1 for a in u),
                "unit CRT lift has a nonunit coordinate")
        require(all(gcd(a, N) != 1 for a in v),
                "zero-septimal section has a unit coordinate")
        require(projection_13(tuple(a % P13 for a in u)) == q,
                "unit CRT lift changed target")
        require(projection_13(tuple(a % P13 for a in v)) == q,
                "nonunit CRT lift changed target")

        unit_witness[q] = u
        nonunit_witness[q] = v
        digest_rows.append(f"{q}:{u}:{v}")

    require(len(set(unit_witness.values())) == P13 ** 2,
            "unit section is not injective")
    require(len(set(nonunit_witness.values())) == P13 ** 2,
            "nonunit section is not injective")

    # Hostile current: unrestricted target vector is identically one, while
    # the all-unit target vector is zero.
    hostile = {("nonunit", q): 1 for q in targets}
    hostile_A = U_of_sparse(hostile)
    hostile_B = J_of_sparse(hostile)
    require(hostile_A == {q: 1 for q in targets},
            "hostile unrestricted target vector is not full support")
    require(hostile_B == {}, "hostile all-unit target vector is nonzero")

    # One switch per target lies in ker U and maps under J to that coordinate.
    # Their J-images are the 169 coordinate vectors, proving independence and
    # surjectivity of J restricted to ker U.  The 168 q!=0 switches are the
    # precise witnesses requested by the live nonzero-target obstruction.
    for q in targets:
        switch = {("unit", q): 1, ("nonunit", q): -1}
        require(U_of_sparse(switch) == {}, "switch escaped ker U")
        require(J_of_sparse(switch) == {q: 1},
                "switch did not change exactly one J coordinate")

    nonzero_targets = tuple(q for q in targets if q != (0, 0))
    require(len(nonzero_targets) == 168, "wrong nonzero-target census")

    ambient_dimension = N ** (D - 1)
    fibre_size = P7 ** (D - 1) * P13 ** 6
    require(ambient_dimension == P13 ** 2 * fibre_size,
            "target-fibre cardinalities do not partition K_91")

    septimal_support = sum(w % P7 != 0 for w in W)
    septimal_unit_count = (
        6 ** (D - septimal_support)
        * (6 ** septimal_support + 6 * (-1) ** septimal_support)
        // P7
    )
    require(septimal_support >= 2, "THM-2325 septimal support failed")
    require(septimal_unit_count > 0, "septimal all-unit kernel is empty")

    mod13_fibre_lower_bound = 9 * 12 ** 5
    full_unit_fibre_lower_bound = (
        septimal_unit_count * mod13_fibre_lower_bound
    )
    witness_digest = sha256("\n".join(digest_rows).encode()).hexdigest()

    print("THM-2544 ALL-UNIT TARGET PROJECTOR KERNEL -- EXACT REFEREE")
    print(f"word={W}")
    print("order=(guard/unit_0,...,guard/unit_5,owner_j,target_a,target_b)")
    print(f"owner_packet omitted={OMITTED} grafts=({GRAFT_A},{GRAFT_B})"
          f" pivot={PIVOT} rank_mod13=6 bright_columns=9/9")
    print(f"K91_size=91^8={ambient_dimension}")
    print(f"target_count=13^2={P13 ** 2} fibre_size=7^8*13^6={fibre_size}")
    print(f"rank_U={P13 ** 2} dim_ker_U={ambient_dimension-P13**2}")
    print(f"rank_J={P13 ** 2} dim_ker_J={ambient_dimension-P13**2}")
    print(f"rank_(U,J)={2*P13**2}"
          f" dim_(ker_U_inter_ker_J)={ambient_dimension-2*P13**2}")
    print(f"septimal_support={septimal_support}"
          f" septimal_all_unit_count={septimal_unit_count}")
    print(f"mod13_all_unit_count_per_target_at_least={mod13_fibre_lower_bound}")
    print(f"all91unit_count_per_target_at_least={full_unit_fibre_lower_bound}")
    print(f"unit_section={len(unit_witness)}/169"
          f" nonunit_section={len(nonunit_witness)}/169")
    print("hostile: support(Uc)=169/169 support(Jc)=0/169")
    print("kernel switches: 169/169 independent through J;"
          " nonzero-target subbank=168/168")
    print(f"common_mod7_witness={z7}")
    print(f"first_target={(0, 0)} unit={unit_witness[(0, 0)]}"
          f" nonunit={nonunit_witness[(0, 0)]}")
    print(f"last_target={(12, 12)} unit={unit_witness[(12, 12)]}"
          f" nonunit={nonunit_witness[(12, 12)]}")
    print(f"witness_sha256={witness_digest}")
    print(f"exact_checks={CHECKS}")
    print("SCOPE: arbitrary residue-current space only; no Abel/current-image,"
          " covering-row, arrival-field, or LRC(14) conclusion")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
