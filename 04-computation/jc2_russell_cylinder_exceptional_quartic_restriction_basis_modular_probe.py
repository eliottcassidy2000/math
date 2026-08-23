"""Split-fibre restriction-basis sidecar for the THM-3683 quartic folds.

This is a FINITE-FIELD compiler-design probe, not a characteristic-zero
conductor theorem.  It checks every root in two completely split good primes.
"""

import hashlib
import importlib.util
from pathlib import Path
import subprocess
import types

from flint import nmod_mat


PREDECESSOR = Path(
    "04-computation/jc2_russell_cylinder_exceptional_quartic_modular_lift_thm3687.py"
)
EXPECTED_ROOTS = {137: (44, 82, 92, 134), 163: (32, 95, 109, 150)}
EXPECTED_PIVOT_HASH = "8fdc045de885e4a3f17b40cd93e99f36b5ddcf07012d801e35687ce74c5576cd"
EXPECTED_GAP_HASH = "0b7c480df6c67445d1acc0c016d81935f121d0aa0c0306c27508e7ab5d97daf4"
EXPECTED_GAPS = (
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35,
    37, 38, 40, 41, 43, 44, 45, 46, 47, 49, 50, 52, 53, 55, 56,
    58, 59, 61, 62, 64, 65, 67, 68, 70, 73, 74, 76, 77, 79, 80,
    82, 85, 86, 88, 91, 94, 95, 97, 98, 100, 103, 106, 109, 112,
    115, 116, 118, 121, 127, 130, 133, 136, 139, 148, 151, 157, 169,
)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def load_predecessor():
    if PREDECESSOR.is_file():
        spec = importlib.util.spec_from_file_location("thm3687_modular", PREDECESSOR)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    # Sparse worktrees may omit a tracked predecessor from the physical tree.
    source = subprocess.check_output(
        ["git", "show", f"HEAD:{PREDECESSOR.as_posix()}"], text=True
    )
    module = types.ModuleType("thm3687_modular")
    exec(compile(source, PREDECESSOR.as_posix(), "exec"), module.__dict__)
    return module


def integer_hash(values):
    return hashlib.sha256(",".join(map(str, values)).encode("ascii")).hexdigest()


probe = load_predecessor()
uniform_rows = []
for prime, expected_roots in EXPECTED_ROOTS.items():
    roots = tuple(probe.roots_mod_prime(prime))
    require(roots == expected_roots, f"split roots p={prime}")
    for root in roots:
        predecessor = probe.solve_j0(prime, root, 198)
        B, C, E, dB, dC, dE = predecessor["compiler"]
        monomials = probe.packet(
            (B, C, E), (dB, dC, dE), (30, 21, 18), 375, prime
        )
        restrictions = [item[1] for item in monomials]
        require(len(restrictions) == 1017, f"packet p={prime},r={root}")

        coefficient_matrix = nmod_mat(
            376,
            len(restrictions),
            [
                probe.coefficient(polynomial, degree)
                for degree in range(376)
                for polynomial in restrictions
            ],
            prime,
        )
        reduced, rank = coefficient_matrix.rref()
        require(rank == 287, f"restriction rank p={prime},r={root}")
        pivot_columns = [
            next(column for column in range(len(restrictions)) if reduced[row, column])
            for row in range(rank)
        ]

        # Reverse source degree before row reduction so pivot positions record
        # leading degrees of a canonical filtered basis rather than low jets.
        selected_reversed = nmod_mat(
            rank,
            376,
            [
                int(coefficient_matrix[375 - reverse_degree, column])
                for column in pivot_columns
                for reverse_degree in range(376)
            ],
            prime,
        )
        degree_reduced, degree_rank = selected_reversed.rref()
        require(degree_rank == rank, f"degree rank p={prime},r={root}")
        reverse_pivots = [
            next(column for column in range(376) if degree_reduced[row, column])
            for row in range(rank)
        ]
        basis_degrees = sorted(375 - value for value in reverse_pivots)
        gaps = tuple(degree for degree in range(376) if degree not in set(basis_degrees))
        conductor = min(
            degree
            for degree in range(376)
            if all(value in set(basis_degrees) for value in range(degree, 376))
        )

        # The coupled 395+174+395 selector uses precisely these 287 distinct
        # target restrictions; this is the bridge to a future exact basis solve.
        coupled = probe.structured_coupled_solution(
            prime, root, predecessor, 375
        )
        chosen_free = [
            coupled["free_columns"][index]
            for index in coupled["selected_direction_indices"]
        ]
        active_restrictions = sorted(
            set(
                column % len(restrictions)
                for column in (*coupled["pivot_columns"], *chosen_free)
            )
        )

        require(active_restrictions == pivot_columns, f"active basis p={prime},r={root}")
        require(gaps == EXPECTED_GAPS, f"gap list p={prime},r={root}")
        require(conductor == 170, f"candidate conductor p={prime},r={root}")
        require(integer_hash(pivot_columns) == EXPECTED_PIVOT_HASH, f"pivot hash p={prime},r={root}")
        require(integer_hash(gaps) == EXPECTED_GAP_HASH, f"gap hash p={prime},r={root}")
        uniform_rows.append((prime, root, rank, len(gaps), conductor))
        print(
            f"prime={prime};root={root};restriction_rank={rank};gaps={len(gaps)};"
            f"candidate_conductor={conductor};active_equals_pivots=True"
        )

require(len(uniform_rows) == 8, "eight split fibres")
print(f"pivot_hash={EXPECTED_PIVOT_HASH}")
print(f"gap_hash={EXPECTED_GAP_HASH}")
print("scope=FINITE_FIELD_ONLY;characteristic_zero_conductor=OPEN")
print("RESULT=PASS")
