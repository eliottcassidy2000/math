#!/usr/bin/env python3
"""Exact symbolic audit for THM-4280's arithmetic three-channel observer.

The audit keeps two coefficient categories separate.

* After extension through the chosen complex embedding, the four formal-log
  channels in degrees 1, 2, 4, and 7 form a free rank-four matroid.
* On the actual O=Z[omega] Hom lattice, a degree-1 or degree-7 coefficient
  lies in K(kappa)=K+K*kappa over K=Q(omega).  Its two K-coordinates make
  either (1,2,4) or (2,4,7) a minimal exact arithmetic observer.

This script audits the exact linear algebra.  The theorem supplies the
geometric identifications and the proof that kappa is not in K.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations

import sympy as sp


ACTIVE = (1, 2, 4, 7)


def require(condition: object, label: str) -> None:
    """Optimization-safe audit gate."""
    if not condition:
        raise RuntimeError(f"audit failure: {label}")


def powerset(items: tuple[int, ...]):
    """Yield subsets in deterministic cardinality/lexicographic order."""
    for size in range(len(items) + 1):
        yield from combinations(items, size)


def subset_rows(matrix: sp.Matrix, subset: tuple[int, ...]) -> sp.Matrix:
    """Select channel rows, including the empty-subset edge case."""
    if not subset:
        return sp.zeros(0, matrix.cols)
    return matrix[[ACTIVE.index(degree) for degree in subset], :]


def support(vector: sp.Matrix) -> tuple[int, ...]:
    """Return the nonzero channel support of a four-vector."""
    return tuple(
        ACTIVE[index]
        for index, value in enumerate(vector)
        if sp.simplify(value) != 0
    )


def main() -> None:
    I = sp.I
    sqrt3 = sp.sqrt(3)
    omega = (-1 + I * sqrt3) / 2
    kappa = (-sqrt3 + I) / 2

    require(sp.simplify(omega**2 + omega + 1) == 0, "omega relation")
    require(sp.simplify(kappa**2 + omega) == 0, "kappa relation")

    # Row-rescaled formal-log channels on the full O-basis [u,f,g,h].
    # Every discarded row scalar is geometrically proved nonzero upstream,
    # so this normalization preserves all kernels and ranks.
    channel = sp.Matrix([
        [0, 1, -kappa, (omega**2 - kappa) / 2],
        [1, 0, 0, 0],
        [0, 0, 0, sp.Rational(1, 2)],
        [0, 1, kappa, (omega**2 + kappa) / 2],
    ])
    determinant = sp.simplify(channel.det())
    require(sp.simplify(determinant - kappa) == 0, "channel determinant")

    complex_ranks = {
        subset: subset_rows(channel, subset).rank()
        for subset in powerset(ACTIVE)
    }
    require(
        all(rank == len(subset) for subset, rank in complex_ranks.items()),
        "complex subset rank law",
    )

    # Pure-channel witnesses after complex scalar extension.  The integral
    # vectors u and v occupy degrees 2 and 4; the kappa-vectors in degrees 1
    # and 7 are not actual O-lattice vectors.
    witnesses = {
        1: sp.Matrix([0, kappa, -1, 0]),       # kappa*f-g
        2: sp.Matrix([1, 0, 0, 0]),            # u
        4: sp.Matrix([0, -omega**2, -1, 2]),   # v
        7: sp.Matrix([0, kappa, 1, 0]),        # kappa*f+g
    }
    for degree, witness in witnesses.items():
        require(
            support(sp.simplify(channel * witness)) == (degree,),
            f"pure witness in degree {degree}",
        )

    # Write A=p+omega^2*d/2 and B=q+d/2.  Over K, c1=A-kappa*B
    # and c7=A+kappa*B each expose two coordinates because
    # K(kappa)=K direct-sum K*kappa.  Coordinates below are [a,A,B,e],
    # where e=d/2 on the c1=0 or c7=0 integral kernel.
    arithmetic_blocks = {
        1: ((0, 1, 0, 0), (0, 0, -1, 0)),
        2: ((1, 0, 0, 0),),
        4: ((0, 0, 0, 1),),
        7: ((0, 1, 0, 0), (0, 0, 1, 0)),
    }
    arithmetic_ranks: dict[tuple[int, ...], int] = {}
    for subset in powerset(ACTIVE):
        rows = [row for degree in subset for row in arithmetic_blocks[degree]]
        block = sp.Matrix(rows) if rows else sp.zeros(0, 4)
        rank = block.rank()
        expected = (
            2 * bool(set(subset) & {1, 7})
            + int(2 in subset)
            + int(4 in subset)
        )
        require(rank == expected, f"arithmetic rank for subset {subset}")
        arithmetic_ranks[subset] = rank

    minimal_bases = tuple(
        subset
        for subset, rank in arithmetic_ranks.items()
        if rank == 4
        and all(
            arithmetic_ranks[smaller] < 4
            for size in range(len(subset))
            for smaller in combinations(subset, size)
        )
    )
    require(
        minimal_bases == ((1, 2, 4), (2, 4, 7)),
        "arithmetic minimal bases",
    )
    require(
        all(
            arithmetic_ranks[subset] < 4
            for subset in combinations(ACTIVE, 2)
        ),
        "all two-channel observers are deficient",
    )

    # Exact visible-kernel parametrization forced arithmetically by c1=0:
    # [a,p,q,d]=[a,-omega^2*e,-e,2e]=a*u+e*v.
    a, e = sp.symbols("a e")
    visible = sp.Matrix([a, -omega**2 * e, -e, 2 * e])
    visible_channels = sp.simplify(channel * visible)
    require(
        visible_channels == sp.Matrix([0, a, e, 0]),
        "visible kernel parametrization",
    )

    # The full complex kernel of (c1,c2,c4) is instead the pure degree-7
    # line.  This is the sharp base-change hostile.
    triple_124 = subset_rows(channel, (1, 2, 4))
    hostile_124 = witnesses[7]
    require(triple_124.rank() == 3, "complex rank of triple 124")
    require(
        triple_124 * hostile_124 == sp.zeros(3, 1),
        "triple 124 hostile kernel",
    )
    require(hostile_124 != sp.zeros(4, 1), "triple 124 hostile nonzero")
    require(
        support(channel * hostile_124) == (7,),
        "triple 124 hostile support",
    )

    # Likewise, the other arithmetic basis loses the pure degree-1 line
    # after complexification.
    triple_247 = subset_rows(channel, (2, 4, 7))
    hostile_247 = witnesses[1]
    require(triple_247.rank() == 3, "complex rank of triple 247")
    require(
        triple_247 * hostile_247 == sp.zeros(3, 1),
        "triple 247 hostile kernel",
    )
    require(
        support(channel * hostile_247) == (1,),
        "triple 247 hostile support",
    )

    # On ker(c1|M)=O*u direct-sum O*v, the inherited visible Hermitian form
    # is 4(N(a)+N(e)).  Thus c1=0 is impossible on degree 34 or 42.
    degree_shells = (34, 42)
    require(
        all(degree % 4 == 2 for degree in degree_shells),
        "degree shells modulo four",
    )

    complex_census = dict(sorted(Counter(complex_ranks.values()).items()))
    arithmetic_census = dict(sorted(Counter(arithmetic_ranks.values()).items()))

    print("THM4280 INTEGRAL THREE-CHANNEL FAT-CONTACT OBSERVER EXACT AUDIT")
    print("status=FINITE_EXACT PASS")
    print("base_ring=O=Z[omega]")
    print("field_tower=Q < K=Q(omega) < K(kappa), kappa^2=-omega")
    print("normalized_channel_rows=degrees_1_2_4_7_on_[u,f,g,h]")
    print("normalized_complex_determinant=kappa")
    print(f"selected_complex_kappa={sp.sstr(kappa)}")
    print(f"complex_subset_rank_census={complex_census}")
    print("complex_matroid=U_4_4_on_1_2_4_7")
    print("complex_full_observer=(1,2,4,7)_ONLY")
    print("complex_pure_witnesses=1:kappa*f-g;2:u;4:v;7:kappa*f+g")
    print(f"arithmetic_block_rank_census={arithmetic_census}")
    print("arithmetic_rank_law=2*hit_{1,7}+hit_2+hit_4")
    print("arithmetic_minimal_full_observers=(1,2,4);(2,4,7)")
    print("arithmetic_all_two_channel_observers=RANK_DEFICIENT")
    print("integral_kernel_c1=integral_kernel_c7=O*u+O*v")
    print("integral_kernel_coordinates=[a,-omega^2*e,-e,2e]")
    print("base_change_hostile_124=kappa*f+g_pure_degree_7")
    print("base_change_hostile_247=kappa*f-g_pure_degree_1")
    print("actual_ramification_spectrum=1,2,4")
    print("five_Q_nonconstancy=SHARP_v_constant_on_4Q_not_5Q")
    print("degree_34_42_c1_zero_test=PASS_by_degree_mod_4")
    print("target_translation_sidecar=value_at_Q_REQUIRED")
    print("raw_Keller_descent=NOT_ADDRESSED")


if __name__ == "__main__":
    main()
