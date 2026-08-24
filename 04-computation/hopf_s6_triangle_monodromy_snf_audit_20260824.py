#!/usr/bin/env python3
"""Exact finite audit of the integer linear algebra in Alpoge's S6 preprint.

This script checks the displayed matrices, their exterior powers, the
two-generator clutch presentation, and the finite conductor pushout conditional
on the manuscript's stated oriented opposite-side quotient.  It does not check
the analytic family, the toric/logarithmic fillings, realization of that
quotient, the global homology identifications, or the claim that the resulting
manifold is S^6.
"""

from __future__ import annotations

from itertools import combinations
from math import gcd

from sympy import Matrix, eye
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def exterior_power(matrix: Matrix, degree: int) -> Matrix:
    """Matrix of the induced action on the indicated exterior power."""
    n = matrix.rows
    if matrix.cols != n or not 0 <= degree <= n:
        raise ValueError("exterior_power expects a square matrix and 0 <= q <= n")
    if degree == 0:
        return Matrix([[1]])
    indices = list(combinations(range(n), degree))
    return Matrix(
        [
            [matrix.extract(target, source).det() for source in indices]
            for target in indices
        ]
    )


def smith_diagonal(matrix: Matrix) -> tuple[int, ...]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return tuple(abs(int(smith[i, i])) for i in range(min(smith.shape)))


def relation_matrix(
    a: int, b: int, ell_0: int, ell_1: int, ell_2: int
) -> tuple[Matrix, int]:
    """Abelian presentation in the basis (x,c), with det(R)=-p."""
    matrix = Matrix([[a, -ell_1], [b, ell_2 - b * ell_0]])
    p = a * b * ell_0 - b * ell_1 - a * ell_2
    assert matrix.det() == -p
    return matrix, p


def main() -> None:
    T1 = Matrix(
        [
            [1, 0, -6, 2],
            [0, -1, 1, 1],
            [0, -1, 0, 1],
            [0, 0, 0, 1],
        ]
    )
    T2 = Matrix(
        [
            [1, 6, 0, -3],
            [0, 0, -1, 1],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
        ]
    )
    I4 = eye(4)
    T0 = (T1 * T2).inv()
    N = T0 - I4
    A1 = T1.inv().T
    A2 = T2.inv().T
    M0 = T0.inv().T

    expected_T0 = Matrix(
        [[1, 0, 0, 1], [0, 1, -1, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    )
    expected_A1 = Matrix(
        [[1, 0, 0, 0], [6, 0, 1, 0], [-6, -1, -1, 0], [-2, 1, 0, 1]]
    )
    expected_A2 = Matrix(
        [[1, 0, 0, 0], [0, 0, -1, 0], [-6, 1, 0, 0], [3, 0, 1, 1]]
    )
    expected_M0 = Matrix(
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 1, 1, 0], [-1, 0, 0, 1]]
    )

    assert T0 == expected_T0
    assert A1 == expected_A1
    assert A2 == expected_A2
    assert M0 == expected_M0
    assert T1.det() == T2.det() == 1
    assert T1**3 == I4
    assert T2**4 == I4 and T2**2 != I4
    assert N**2 == Matrix.zeros(4)
    assert N.rank() == 2
    assert A1 * A2 * M0 == I4

    print("DISPLAYED MONODROMY MATRICES: EXACT")
    print("det(T1)=det(T2)=1; T1^3=I; T2^4=I; T2^2!=I")
    print("T0=(T1*T2)^-1=I+N; N^2=0; rank(N)=2")
    print("A1=(T1^-1)^t; A2=(T2^-1)^t; M0=(T0^-1)^t")
    print("A1*A2*M0=I")

    assert smith_diagonal(M0 - I4) == (1, 1, 0, 0)
    joint = (A1 - I4).row_join(A2 - I4)
    assert smith_diagonal(A1 - I4) == (1, 1, 0, 0)
    assert smith_diagonal(A2 - I4) == (1, 1, 0, 0)
    assert smith_diagonal(joint) == (1, 1, 1, 0)
    assert (M0 - I4).nullspace() == [Matrix([0, 0, 1, 0]), Matrix([0, 0, 0, 1])]

    print("\nLATTICE SMITH DATA: EXACT")
    print("SNF(M0-I)=(1,1,0,0); ker(M0-I)=<w_hat,delta_hat>")
    print("SNF(A1-I)=SNF(A2-I)=(1,1,0,0)")
    print("SNF([A1-I | A2-I])=(1,1,1,0): one free coinvariant coordinate")

    expected_invariant_ranks = (1, 2, 4, 2, 1)
    exterior_rows: list[str] = []
    for degree, expected_rank in enumerate(expected_invariant_ranks):
        action = exterior_power(M0, degree)
        response = action - eye(action.rows)
        diagonal = smith_diagonal(response)
        kernel_rank = action.rows - response.rank()
        assert kernel_rank == expected_rank
        assert all(entry in (0, 1) for entry in diagonal)
        assert diagonal.count(0) == expected_rank
        exterior_rows.append(
            f"q={degree}: size={action.rows}, rank(kernel)={kernel_rank}, "
            f"SNF={diagonal}"
        )

    print("\nEXTERIOR-POWER CUSP DATA: EXACT")
    for row in exterior_rows:
        print(row)

    wedge_T1 = exterior_power(T1, 2) - eye(6)
    wedge_T2 = exterior_power(T2, 2) - eye(6)
    expected_wedge_T1 = Matrix(
        [
            [-2, 1, 1, -6, 2, -8],
            [-1, -1, 1, -6, 2, -6],
            [0, 0, 0, 0, 0, -6],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, -2, 1],
            [0, 0, 0, 0, -1, -1],
        ]
    )
    expected_wedge_T2 = Matrix(
        [
            [-1, -1, 1, -6, 6, -3],
            [1, -1, 0, 0, 3, 0],
            [0, 0, 0, 0, 6, 0],
            [0, 0, 0, 0, -1, 0],
            [0, 0, 0, 0, -1, -1],
            [0, 0, 0, 0, 1, -1],
        ]
    )
    assert wedge_T1 == expected_wedge_T1
    assert wedge_T2 == expected_wedge_T2
    joint_wedge = wedge_T1.row_join(wedge_T2)
    wedge_coinvariant = Matrix([0, 0, 1, 6, 0, 0])
    assert smith_diagonal(joint_wedge) == (1, 1, 1, 1, 1, 0)
    assert joint_wedge.T * wedge_coinvariant == Matrix.zeros(12, 1)
    assert int(wedge_coinvariant.dot(Matrix([0, 0, 6, 1, 0, 0]))) == 12

    print("\nJOINT DEGREE-TWO COINVARIANTS: EXACT")
    print("SNF([wedge^2(T1)-I | wedge^2(T2)-I])=(1,1,1,1,1,0)")
    print("primitive coinvariant=(0,0,1,6,0,0); [uw]=6[gamma delta]")
    print("q=uw+6 gamma delta has coinvariant value 12")

    W2 = Matrix(
        [
            [0, -1, 0, 0, -1, 0, 0],
            [1, 1, 1, 0, 0, -1, 0],
            [0, 0, -1, 0, 0, 0, -1],
            [1, 0, 1, 1, -1, 0, 0],
            [0, 0, 0, -1, 0, -1, 0],
            [1, 1, 0, 1, 0, 0, -1],
        ]
    )
    assert smith_diagonal(W2) == (1, 1, 1, 1, 0, 0)
    assert W2.rank() == 4

    # Independent oriented-conductor presentation.  The first four rows are
    # H,E1,E2,E3 in H_2(dP6); the last three are the quotient double curves.
    # Columns are the six boundary curves around the anticanonical hexagon.
    conductor_pushout = Matrix(
        [
            [0, 1, 0, 1, 0, 1],
            [1, -1, 0, 0, 0, -1],
            [0, -1, 1, -1, 0, 0],
            [0, 0, 0, -1, 1, -1],
            [-1, 0, 0, -1, 0, 0],
            [0, -1, 0, 0, -1, 0],
            [0, 0, -1, 0, 0, -1],
        ]
    )
    kernel_1 = Matrix([-1, -1, 0, 1, 1, 0])
    kernel_2 = Matrix([1, 0, -1, -1, 0, 1])
    assert conductor_pushout * kernel_1 == Matrix.zeros(7, 1)
    assert conductor_pushout * kernel_2 == Matrix.zeros(7, 1)
    assert conductor_pushout.rank() == 4
    assert smith_diagonal(conductor_pushout) == (1, 1, 1, 1, 0, 0)

    # Hostile control: reverse the relative degree on only the second branch
    # of the first opposite-side pair.  Counts and Euler characteristic stay
    # fixed, but the integral response acquires 2-torsion.
    reversed_branch = conductor_pushout.copy()
    reversed_branch[4, 3] = 1
    assert reversed_branch.rank() == 5
    assert smith_diagonal(reversed_branch) == (1, 1, 1, 1, 2, 0)

    e_dp6 = 6
    e_hexagon = 6
    e_double_locus = 2
    e_W = e_dp6 + e_double_locus - e_hexagon
    assert e_W == 2

    # For the oriented quotient, the hexagon incidence loop maps trivially on
    # H_1(D).  The pushout LES and the saturated degree-two map then give the
    # free homology ranks in degrees 0,...,4.
    h1_boundary_rank = 1
    h1_double_locus_rank = 2
    h1_incidence_map_rank = 0
    h_W_ranks = (
        1,
        h1_double_locus_rank - h1_incidence_map_rank,
        (7 - conductor_pushout.rank())
        + (h1_boundary_rank - h1_incidence_map_rank),
        6 - conductor_pushout.rank(),
        1,
    )
    assert h_W_ranks == (1, 2, 4, 2, 1)
    assert sum((-1) ** q * rank for q, rank in enumerate(h_W_ranks)) == e_W

    alpha2 = Matrix(
        [
            [2, 1, 3, 0, 0, 0],
            [-4, -2, 0, 1, 0, 0],
            [-2, -2, -4, 0, 0, 0],
            [3, 3, 0, -1, 0, 0],
        ]
    )
    alpha2_annihilator = Matrix([4, 2, 3, 2])
    assert smith_diagonal(alpha2) == (1, 1, 1, 0)
    assert alpha2.T * alpha2_annihilator == Matrix.zeros(6, 1)

    print("\nDEGREE-TWO TOPOLOGY MATRICES: EXACT")
    print("W conductor/pushout map: rank=4, SNF=(1,1,1,1,0,0)")
    print("alpha2: rank=3, SNF=(1,1,1,0)")
    print("primitive left annihilator of alpha2=(4,2,3,2)")

    print("\nORIENTED CONDUCTOR PUSHOUT: FINITE-EXACT")
    print("oriented A: rank=4, SNF=(1,1,1,1,0,0), kernel rank=2")
    print("chi(dP6)+chi(D)-chi(hexagon)=6+2-6=2")
    print("with zero H1 incidence map: H_0..H_4 ranks=(1,2,4,2,1), all free")
    print("one reversed branch: rank=5, SNF=(1,1,1,1,2,0)")
    print("orientation hostile: same cell counts and Euler, but Z/2 appears")

    v1 = Matrix([1, 2, -4, 0])
    v2 = Matrix([-1, -3, 3, 0])
    assert A1 * v1 == v1
    assert A2 * v2 == v2
    assert (int(v1[0]), int(v2[0])) == (1, -1)
    print("\nTWIST VECTORS: EXACT")
    print("A1*v1=v1, A2*v2=v2, (ell1,ell2)=(1,-1)")

    chosen, chosen_p = relation_matrix(3, 4, 0, 1, -1)
    comparison, comparison_p = relation_matrix(3, 4, 0, 1, 1)
    assert chosen_p == -1 and smith_diagonal(chosen) == (1, 1)
    assert comparison_p == -7 and smith_diagonal(comparison) == (1, 7)

    print("\nCLUTCH PRESENTATIONS: EXACT")
    print(f"(ell0,ell1,ell2)=(0,1,-1): p={chosen_p}, SNF={smith_diagonal(chosen)}")
    print(
        f"(ell0,ell1,ell2)=(0,1,1): p={comparison_p}, "
        f"SNF={smith_diagonal(comparison)}"
    )

    checks = 0
    finite_groups = 0
    for a in range(2, 16):
        for b in range(2, 16):
            if gcd(a, b) != 1:
                continue
            for ell_0 in range(-3, 4):
                for ell_1 in range(-6, 7):
                    if gcd(ell_1, a) != 1:
                        continue
                    for ell_2 in range(-6, 7):
                        if gcd(ell_2, b) != 1:
                            continue
                        matrix, p = relation_matrix(a, b, ell_0, ell_1, ell_2)
                        checks += 1
                        assert gcd(p, a * b) == 1
                        first_divisor = gcd(
                            *(abs(int(entry)) for entry in matrix)
                        )
                        assert first_divisor == 1
                        if p:
                            finite_groups += 1
                            # For a full-rank 2 x 2 matrix the two Smith
                            # divisors multiply to |det|.  The first is 1.
                            assert abs(int(matrix.det())) == abs(p)

    print("\nCOPRIME (a,b) CLUTCH COMPILER: FINITE-EXACT CONTROL")
    print(f"admissible tuples checked={checks}; finite cokernels checked={finite_groups}")
    print("for every checked tuple: gcd(p,a*b)=1 and SNF(R)=(1,|p|) when p!=0")
    print("native response: Delta_ell0 p=ab, Delta_ell1 p=-b, Delta_ell2 p=-a")
    print("for (a,b)=(3,4), kernel move (ell0,ell1,ell2)+=(1,3,0) preserves p")

    print("\nSCOPE")
    print("VERIFIED: displayed integer matrices, exterior-power Smith data, twist determinant.")
    print("FINITE-EXACT GIVEN ORIENTED QUOTIENT: conductor homology and orientation hostile.")
    print("NOT VERIFIED: period-map analysis, fillings, global gluing, homology identifications, S6.")


if __name__ == "__main__":
    main()
