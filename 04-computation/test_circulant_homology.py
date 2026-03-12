"""
test_circulant_homology.py — Pytest suite for circulant_homology.py

Tests:
  1. omega_dims() for T_3, T_7, T_11 against known values
  2. betti_numbers() (cached) for T_3, T_7, T_11
  3. betti_numbers() (computed) for T_3 and T_7
  4. euler_characteristic() consistency
  5. CirculantHomology general API
  6. path_count_sequence() for T_7
  7. Find-prime utility

Author: kind-pasteur-2026-03-12-S55
"""
import sys
import pytest
sys.path.insert(0, '04-computation')

from circulant_homology import (
    CirculantHomology, PaleyHomology,
    find_prime_for_roots, find_nth_root_of_unity, is_prime,
)


# ---------------------------------------------------------------------------
# Test utility functions
# ---------------------------------------------------------------------------

class TestUtilities:
    def test_is_prime(self):
        assert is_prime(2)
        assert is_prime(7)
        assert is_prime(11)
        assert is_prime(89)
        assert not is_prime(1)
        assert not is_prime(4)
        assert not is_prime(91)

    def test_find_prime_for_roots(self):
        # T_3: need prime q with 3 | (q-1)
        q3 = find_prime_for_roots(3)
        assert is_prime(q3) and (q3 - 1) % 3 == 0

        # T_7: prime should be 29 (28=4*7) or similar
        q7 = find_prime_for_roots(7)
        assert is_prime(q7) and (q7 - 1) % 7 == 0

        # T_11: should be 23 (22=2*11)
        q11 = find_prime_for_roots(11)
        assert is_prime(q11) and (q11 - 1) % 11 == 0

    def test_find_nth_root_of_unity(self):
        p = 11
        prime = find_prime_for_roots(p)
        omega = find_nth_root_of_unity(p, prime)
        # omega is a primitive p-th root of unity mod prime
        assert pow(omega, p, prime) == 1
        assert all(pow(omega, k, prime) != 1 for k in range(1, p))


# ---------------------------------------------------------------------------
# PaleyHomology: Omega dims (cached, fast)
# ---------------------------------------------------------------------------

class TestOmegaDimsCached:
    def test_t3_omega_cached(self):
        h = PaleyHomology(p=3)
        dims = h.omega_dims()
        assert dims == [1, 1, 0]

    def test_t7_omega_cached(self):
        h = PaleyHomology(p=7)
        dims = h.omega_dims()
        assert dims == [1, 3, 6, 9, 9, 6, 3]

    def test_t11_omega_cached(self):
        h = PaleyHomology(p=11)
        dims = h.omega_dims()
        assert dims == [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]

    def test_t7_chi(self):
        h = PaleyHomology(p=7)
        chi = h.euler_characteristic()
        assert chi == 1  # per eigenspace

    def test_t11_chi(self):
        h = PaleyHomology(p=11)
        chi = h.euler_characteristic()
        assert chi == 1


# ---------------------------------------------------------------------------
# PaleyHomology: Omega dims (computed, slower but correct)
# ---------------------------------------------------------------------------

class TestOmegaDimsComputed:
    def test_t3_omega_computed(self):
        h = PaleyHomology(p=3)
        dims = h.omega_dims(use_cache=False)
        assert dims == [1, 1, 0]

    def test_t7_omega_computed(self):
        h = PaleyHomology(p=7)
        dims = h.omega_dims(use_cache=False)
        assert dims == [1, 3, 6, 9, 9, 6, 3]

    def test_t11_omega_computed(self):
        h = PaleyHomology(p=11)
        dims = h.omega_dims(use_cache=False)
        assert dims == [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]

    def test_t7_palindrome(self):
        """T_7 Omega dims are palindromic in the interior degrees 1..n-1."""
        h = PaleyHomology(p=7)
        dims = h.omega_dims(use_cache=False)
        # Interior = dims[1:] (degrees 1..6)
        interior = dims[1:]
        assert interior == interior[::-1], f"T_7 inner dims not palindromic: {interior}"


# ---------------------------------------------------------------------------
# PaleyHomology: Betti numbers (cached)
# ---------------------------------------------------------------------------

class TestBettiCached:
    def test_t3_betti_cached(self):
        h = PaleyHomology(p=3)
        b = h.betti_numbers()
        assert b == [1, 1, 0]

    def test_t7_betti_cached(self):
        h = PaleyHomology(p=7)
        b = h.betti_numbers()
        assert b == [1, 0, 0, 0, 6, 0, 0]

    def test_t11_betti_cached(self):
        h = PaleyHomology(p=11)
        b = h.betti_numbers()
        assert b == [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0]

    def test_t7_betti_chi(self):
        """Total chi(T_7) = sum_m (-1)^m beta_m = p = 7."""
        h = PaleyHomology(p=7)
        b = h.betti_numbers()
        chi = sum((-1)**m * bm for m, bm in enumerate(b))
        assert chi == 7

    def test_t11_betti_chi(self):
        """Total chi(T_11) = 11."""
        h = PaleyHomology(p=11)
        b = h.betti_numbers()
        chi = sum((-1)**m * bm for m, bm in enumerate(b))
        assert chi == 11


# ---------------------------------------------------------------------------
# Betti (computed, no cache) — T_3 and T_7 only (fast)
# ---------------------------------------------------------------------------

class TestBettiComputed:
    def test_t3_betti_computed(self):
        h = PaleyHomology(p=3)
        b = h.betti_numbers(use_cache=False)
        assert b == [1, 1, 0]

    def test_t7_betti_computed(self):
        h = PaleyHomology(p=7)
        b = h.betti_numbers(use_cache=False)
        assert b == [1, 0, 0, 0, 6, 0, 0]

    def test_t7_betti_nonnegative(self):
        """All Betti numbers must be non-negative."""
        h = PaleyHomology(p=7)
        b = h.betti_numbers(use_cache=False)
        assert all(bm >= 0 for bm in b)

    def test_t3_betti_matches_cache(self):
        h = PaleyHomology(p=3)
        b_cached = h.betti_numbers(use_cache=True)
        b_computed = h.betti_numbers(use_cache=False)
        assert b_cached == b_computed

    def test_t7_betti_matches_cache(self):
        h = PaleyHomology(p=7)
        b_cached = h.betti_numbers(use_cache=True)
        b_computed = h.betti_numbers(use_cache=False)
        assert b_cached == b_computed


# ---------------------------------------------------------------------------
# CirculantHomology: general circulant (not Paley)
# ---------------------------------------------------------------------------

class TestCirculantGeneral:
    def test_n7_s124_omega(self):
        """T_7 = C_7^{1,2,4}, verify Omega dims match Paley."""
        h = CirculantHomology(n=7, S={1, 2, 4})
        dims = h.omega_dims()
        assert dims == [1, 3, 6, 9, 9, 6, 3]

    def test_path_count_sequence_t7(self):
        h = PaleyHomology(p=7)
        counts = h.path_count_sequence()
        # |A_0|=1, |A_1|=3, |A_6|=H(T_7)=? (all Ham paths divided by n)
        assert counts[0] == 1
        assert counts[1] == len(h.S)  # one diff seq per element of S
        assert counts[-1] > 0  # non-zero Hamiltonian paths

    def test_euler_characteristic_additive(self):
        """chi(full) = sum of chi per eigenspace = n * chi_per_eigenspace."""
        h = PaleyHomology(p=7)
        # Per-eigenspace chi = 1 for T_7
        # omega_dims() returns the k=0 eigenspace dims, chi=1
        chi = h.euler_characteristic()
        assert chi == 1

    def test_custom_circulant_n3(self):
        """3-cycle: S={1}."""
        h = CirculantHomology(n=3, S={1})
        dims = h.omega_dims()
        assert dims == [1, 1, 0]
        b = h.betti_numbers()
        assert b == [1, 1, 0]

    def test_custom_circulant_betti_nonneg(self):
        """Betti numbers should be non-negative for any circulant."""
        h = CirculantHomology(n=7, S={1, 2, 4})
        b = h.betti_numbers()
        assert all(bm >= 0 for bm in b)


# ---------------------------------------------------------------------------
# Run as script
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
