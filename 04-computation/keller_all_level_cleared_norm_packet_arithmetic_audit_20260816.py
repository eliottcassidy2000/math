#!/usr/bin/env python3
"""Arithmetic audit for THM-3528's all-level cleared norm tower.

The geometric proof lives in THM-3528: regularity of the finite-algebra norm
on A[L^-1], the two complete divergent packet faces, and nonnegative finite
sheet valuation imply polynomiality after multiplying by L^e.  This companion
checks the complete recurrence arithmetic, invariant cone, Pell/Cassini laws,
named initial rows, and the exact finite-sheet defect bookkeeping.

It deliberately does not assert that the finite-sheet defect is zero after
R8, or that later raw cleared norms are image primes, irreducible, separable,
or new nonproperness components.
"""

from __future__ import annotations

from hashlib import sha256
import json


EXPECTED_SEMANTIC_SHA256 = (
    "a77811be1a53f2e0d0e0eeac3b4a4ecac358f79c8ed0b6ec5fa6f03f3bb0c826"
)
MATRIX = ((7, -2), (3, -2))
NAMED = ("L", "H", "J", "G", "R5", "R6", "R7", "R8")
EXPECTED_NAMED_ROWS = (
    (1, 0),
    (7, 3),
    (43, 15),
    (271, 99),
    (1699, 615),
    (10663, 3867),
    (66907, 24255),
    (419839, 152211),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def step(row: tuple[int, int]) -> tuple[int, int]:
    e, m = row
    return 7 * e - 2 * m, 3 * e - 2 * m


def quadratic(row: tuple[int, int]) -> int:
    e, m = row
    return 3 * e * e - 9 * e * m + 2 * m * m


def spinor(row: tuple[int, int]) -> tuple[int, int]:
    e, m = row
    return 6 * e - 9 * m, m


def main() -> None:
    require(MATRIX[0][0] * MATRIX[1][1] - MATRIX[0][1] * MATRIX[1][0] == -8,
            "matrix determinant")

    rows = [(1, 0)]
    for _index in range(15):
        rows.append(step(rows[-1]))
    rows = tuple(rows)
    require(rows[:len(NAMED)] == EXPECTED_NAMED_ROWS,
            ("named rows", rows[:len(NAMED)]))

    cone_record = []
    congruence_record = []
    pell_record = []
    cassini_record = []
    for index, (e, m) in enumerate(rows):
        require(e > 0 and 0 <= 2 * m <= e, ("packet cone", index, e, m))
        require(m % 3 == 0 and e % 3 == 1, ("packet congruence", index, e, m))
        # Symbolic cone identities used by the proof:
        #   m'=3e-2m >= 2e > 0,
        #   e'-2m'=e+2m >= 0.
        if index + 1 < len(rows):
            e_next, m_next = rows[index + 1]
            require(m_next == 3 * e - 2 * m and m_next >= 2 * e,
                    ("m cone transport", index))
            require(e_next - 2 * m_next == e + 2 * m,
                    ("cone identity", index))
        cone_record.append((index, e, m, e - 2 * m))
        congruence_record.append((index, e % 6, m % 6))

        x, y = spinor((e, m))
        require(x * x - 57 * y * y == 36 * ((-8) ** index),
                ("Pell norm", index, x, y))
        require(quadratic((e, m)) == 3 * ((-8) ** index),
                ("packet quadratic", index))
        pell_record.append((index, x, y, quadratic((e, m))))

        if index + 1 < len(rows):
            e_next, m_next = rows[index + 1]
            determinant = e * m_next - m * e_next
            require(determinant == 3 * ((-8) ** index),
                    ("Cassini", index, determinant))
            cassini_record.append((index, determinant))

    # Both scalar coordinates obey u_(n+2)=5u_(n+1)+8u_n.
    for index in range(len(rows) - 2):
        for coordinate in (0, 1):
            require(rows[index + 2][coordinate]
                    == 5 * rows[index + 1][coordinate]
                    + 8 * rows[index][coordinate],
                    ("order-two recurrence", index, coordinate))

    # The exact old-L bookkeeping.  For every formal nonnegative finite-sheet
    # valuation s, two divergent values -e/2 contribute -e in total, so the
    # norm and cleared valuations are (-e+s,s).  This is the theorem's defect
    # identity; s=0 is the strictly stronger unit/coprimality gate.
    defect_record = []
    for index, (e, _m) in enumerate(rows[:10]):
        samples = []
        for finite_order in (0, 1, 2, e, e + 1):
            norm_order = -e + finite_order
            cleared_order = e + norm_order
            require(cleared_order == finite_order and cleared_order >= 0,
                    ("clearing defect", index, finite_order))
            samples.append((finite_order, norm_order, cleared_order))
        defect_record.append((index, e, tuple(samples)))

    named_record = tuple(zip(NAMED, EXPECTED_NAMED_ROWS))
    next_rows = rows[len(NAMED):len(NAMED) + 3]
    record = (
        MATRIX,
        named_record,
        next_rows,
        tuple(cone_record),
        tuple(congruence_record),
        tuple(pell_record),
        tuple(cassini_record),
        tuple(defect_record),
        ("proved finite-sheet units through input R7/output R8", 7),
        ("later unit/coprimality/image/separability", "OPEN"),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== fixed Keller all-level cleared-norm packet arithmetic audit ==")
    print(f"matrix={MATRIX};det=-8;named_rows={named_record}")
    print(f"next_raw_packet_rows={next_rows}")
    print("cone_transport=(e>0,0<=2m<=e)->(m'>=2e,e'-2m'=e+2m>=0): PASS")
    print("congruences=(m divisible by 3,e=1 mod 3): PASS")
    print("recurrence=u_(n+2)=5u_(n+1)+8u_n for e and m: PASS")
    print("pell=3e_n^2-9e_nm_n+2m_n^2=3*(-8)^n through n=15: PASS")
    print("cassini=e_nm_(n+1)-m_ne_(n+1)=3*(-8)^n through n=14: PASS")
    print("valuation_identity=v_L(N(P_n))=-e_n+s_n;v_L(L^e_n*N(P_n))=s_n>=0: SYMBOLIC PASS")
    print("finite_unit_boundary=s_n=0 is proved through input R7/output R8 only;later s_n may be positive")
    print(f"semantic_sha256={semantic}")
    print("scope=arithmetic/valuation audit for raw cleared norms;no later L-coprimality,image,irreducibility,separability,new component,or general JC claim")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
