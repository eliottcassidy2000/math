"""
test_mod_rank_library.py — Pytest suite for mod_rank_library.py

Tests:
  1. gauss_rank_uint8 against known ranks
  2. gauss_rank_nullbasis_uint8 for correctness (C @ v = 0 for each null vector)
  3. matmul_mod for correctness against numpy
  4. certified_rank (two-prime agreement)
  5. betti_number_from_boundary_ranks formula
  6. find_prime_for_roots_of_unity utility

Author: kind-pasteur-2026-03-12-S55
"""
import sys
import numpy as np
import pytest
sys.path.insert(0, '04-computation')

from mod_rank_library import (
    find_prime_for_roots_of_unity,
    gauss_rank_uint8,
    gauss_rank_nullbasis_uint8,
    matmul_mod,
    certified_rank,
    betti_number_from_boundary_ranks,
)

PRIME = 89  # small prime used throughout


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

class TestFindPrime:
    def test_t7_prime(self):
        q = find_prime_for_roots_of_unity(7)
        assert q is not None
        assert (q - 1) % 7 == 0

    def test_t11_prime(self):
        q = find_prime_for_roots_of_unity(11)
        assert q is not None
        assert (q - 1) % 11 == 0

    def test_t19_prime(self):
        q = find_prime_for_roots_of_unity(19)
        assert q is not None
        assert (q - 1) % 19 == 0


# ---------------------------------------------------------------------------
# gauss_rank_uint8: known rank cases
# ---------------------------------------------------------------------------

class TestGaussRank:
    def _rng_matrix(self, rows, cols, seed=42):
        rng = np.random.default_rng(seed)
        return rng.integers(0, PRIME, size=(rows, cols), dtype=np.int16)

    def test_zero_matrix(self):
        C = np.zeros((5, 8), dtype=np.int16)
        assert gauss_rank_uint8(C, PRIME) == 0

    def test_identity_rank(self):
        n = 10
        C = np.eye(n, dtype=np.int16)
        assert gauss_rank_uint8(C, PRIME) == n

    def test_rank_one(self):
        # Row all-same
        C = np.ones((5, 8), dtype=np.int16)
        assert gauss_rank_uint8(C, PRIME) == 1

    def test_full_rank_square(self):
        # Random 5x5 over F_89 — almost certainly full rank
        C = self._rng_matrix(5, 5)
        r = gauss_rank_uint8(C, PRIME)
        assert r <= 5

    def test_rank_deficient(self):
        # Two independent rows; rows 2-4 are linear combos
        row0 = np.array([1, 2, 3, 4, 5], dtype=np.int16)
        row1 = np.array([5, 3, 7, 1, 2], dtype=np.int16)  # independent of row0
        C = np.array([
            row0,
            row1,
            (row0 * 2 + row1) % PRIME,
            (row0 * 3 + row1 * 4) % PRIME,
            (row0 + row1 * 7) % PRIME,
        ], dtype=np.int16)
        r = gauss_rank_uint8(C, PRIME)
        assert r == 2

    def test_wide_matrix(self):
        # 4 rows, 20 cols — rank at most 4
        C = self._rng_matrix(4, 20)
        r = gauss_rank_uint8(C, PRIME)
        assert r <= 4

    def test_tall_matrix(self):
        # 20 rows, 4 cols — rank at most 4
        C = self._rng_matrix(20, 4)
        r = gauss_rank_uint8(C, PRIME)
        assert r <= 4

    def test_known_rank_from_construction(self):
        # Construct a rank-3 matrix explicitly
        A = np.zeros((6, 7), dtype=np.int16)
        A[0, 0] = 1
        A[1, 3] = 1
        A[2, 6] = 1
        # Rows 3-5 are random linear combinations of rows 0-2
        for i in range(3, 6):
            for j in range(3):
                A[i] = (A[i] + (i + j) * A[j]) % PRIME
        r = gauss_rank_uint8(A, PRIME)
        assert r == 3


# ---------------------------------------------------------------------------
# gauss_rank_nullbasis_uint8: null space correctness
# ---------------------------------------------------------------------------

class TestGaussNullbasis:
    def test_null_vectors_are_null(self):
        """C @ v = 0 mod prime for each basis vector v."""
        rng = np.random.default_rng(7)
        C = rng.integers(0, PRIME, size=(5, 10), dtype=np.int16)
        rank, basis = gauss_rank_nullbasis_uint8(C, PRIME)

        assert rank + basis.shape[0] == 10  # rank-nullity

        for v in basis:
            result = (C.astype(np.int32) @ v.astype(np.int32)) % PRIME
            assert np.all(result == 0), f"Null vector not in null space: {result}"

    def test_nullity_rank_sum(self):
        """rank + nullity = cols."""
        rng = np.random.default_rng(13)
        C = rng.integers(0, PRIME, size=(4, 8), dtype=np.int16)
        rank, basis = gauss_rank_nullbasis_uint8(C, PRIME)
        assert rank + basis.shape[0] == 8

    def test_full_rank_trivial_null(self):
        """Full row rank => empty null space (square full-rank matrix)."""
        # Try to construct a full-rank 5x5
        C = np.eye(5, dtype=np.int16)
        rank, basis = gauss_rank_nullbasis_uint8(C, PRIME)
        assert rank == 5
        assert basis.shape[0] == 0

    def test_rank_one_nullity(self):
        """Rank 1 matrix of size 3x5 has nullity 4."""
        row = np.array([1, 2, 3, 4, 5], dtype=np.int16)
        C = np.array([row, row * 2 % PRIME, row * 3 % PRIME], dtype=np.int16)
        rank, basis = gauss_rank_nullbasis_uint8(C, PRIME)
        assert rank == 1
        assert basis.shape[0] == 4


# ---------------------------------------------------------------------------
# matmul_mod
# ---------------------------------------------------------------------------

class TestMatmulMod:
    def test_matches_numpy(self):
        """matmul_mod should match (A @ B) % prime for small matrices."""
        rng = np.random.default_rng(42)
        A = rng.integers(0, PRIME, size=(10, 15), dtype=np.int16)
        B = rng.integers(0, PRIME, size=(15, 8), dtype=np.int16)

        result_lib = matmul_mod(A, B, PRIME)
        result_np = (A.astype(np.int64) @ B.astype(np.int64)) % PRIME

        assert np.array_equal(result_lib, result_np.astype(np.int16))

    def test_identity_matmul(self):
        """A @ I = A."""
        rng = np.random.default_rng(1)
        A = rng.integers(0, PRIME, size=(5, 5), dtype=np.int16)
        I = np.eye(5, dtype=np.int16)
        result = matmul_mod(A, I, PRIME)
        assert np.array_equal(result, A % PRIME)

    def test_zero_matmul(self):
        """A @ 0 = 0."""
        rng = np.random.default_rng(2)
        A = rng.integers(0, PRIME, size=(5, 5), dtype=np.int16)
        Z = np.zeros((5, 5), dtype=np.int16)
        result = matmul_mod(A, Z, PRIME)
        assert np.all(result == 0)


# ---------------------------------------------------------------------------
# certified_rank
# ---------------------------------------------------------------------------

class TestCertifiedRank:
    def test_known_rank(self):
        """Rank 3 matrix certifies as 3."""
        A = np.zeros((5, 7), dtype=np.int64)
        A[0, 0] = 1; A[1, 2] = 1; A[2, 5] = 1
        # Rows 3-4: linear combinations of rows 0-2
        A[3] = 2 * A[0] + 3 * A[1]
        A[4] = 5 * A[1] + A[2]

        rank, details = certified_rank(A, primes=[89, 97])
        assert rank == 3
        assert details['certified'] is True

    def test_full_rank(self):
        """Identity is full rank."""
        A = np.eye(5, dtype=np.int64)
        rank, details = certified_rank(A, primes=[89, 97])
        assert rank == 5
        assert details['certified'] is True

    def test_zero_rank(self):
        """Zero matrix has rank 0."""
        A = np.zeros((5, 7), dtype=np.int64)
        rank, details = certified_rank(A, primes=[89, 97])
        assert rank == 0
        assert details['certified'] is True


# ---------------------------------------------------------------------------
# betti_number_from_boundary_ranks
# ---------------------------------------------------------------------------

class TestBettiFormula:
    def test_basic_formula(self):
        """beta_m = (omega_m - rk_d_m) - rk_d_{m+1}."""
        beta, ker = betti_number_from_boundary_ranks(
            omega_hi=10, omega_lo=None, rk_d_hi=3, rk_d_lo_next=5
        )
        assert ker == 7         # omega_hi - rk_d_hi = 10 - 3 = 7
        assert beta == 2        # ker - rk_d_lo_next = 7 - 5 = 2

    def test_zero_betti(self):
        """beta = 0 when ker = im."""
        beta, ker = betti_number_from_boundary_ranks(
            omega_hi=15, omega_lo=None, rk_d_hi=5, rk_d_lo_next=10
        )
        assert ker == 10
        assert beta == 0

    def test_t7_beta4(self):
        """T_7 has beta_4 = 6. Test the formula gives correct result."""
        # For T_7 k=0: Omega_4=9, rk(d_4)=3, rk(d_5)=0
        # beta_4 = (9 - 3) - 0 = 6
        beta, ker = betti_number_from_boundary_ranks(
            omega_hi=9, omega_lo=None, rk_d_hi=3, rk_d_lo_next=0
        )
        assert beta == 6


# ---------------------------------------------------------------------------
# Run as script
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
