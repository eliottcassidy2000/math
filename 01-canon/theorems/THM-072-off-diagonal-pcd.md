# THM-072: Off-Diagonal Position Character Decomposition

**Type:** Theorem
**Certainty:** 5 -- PROVED (degree 2: algebraic all odd n; degree 4: verified n=5,7)
**Status:** PROVED at degree 2 (all odd n); VERIFIED at degree 4 (n=5,7); CONJECTURED general
**Added by:** opus-2026-03-07-S35
**Tags:** #Walsh-Hadamard #transfer-matrix #off-diagonal #symmetry

---

## Statement

### (i) Interior zero lemma (PROVED, all degrees, all n)

For any Walsh monomial S of degree 2k, and any vertex v that is interior to a
path component of S (i.e., deg_S(v) >= 2):

**Mc_hat[v, b, S] = 0** and **Mc_hat[a, v, S] = 0** for all a, b.

where Mc[a,b](T) = #{directed Hamiltonian paths from a to b in T}.

**Proof:** For the Walsh coefficient to be nonzero, valid permutations must
contain all edges of S. Starting at v requires only 1 incident edge, but
deg_S(v) >= 2 edges of S are incident to v. So no valid HP starts at v. QED.

### (ii) Off-diagonal formula at even Walsh degrees

For degree 2k monomial S with r path components of lengths 2a_1,...,2a_r:

Let epsilon = (-1)^{desc(S)} (constant across valid perms, from PCD proof).
Let H_raw = epsilon * 2^r * (n-2k)!.

**Mc_hat[a,b,S] / H_hat[S] equals:**

| Start (a) | End (b) | Condition | Value |
|-----------|---------|-----------|-------|
| Interior | Any | deg_S(a) >= 2 | 0 |
| Any | Interior | deg_S(b) >= 2 | 0 |
| Endpoint B_i | Endpoint B_j | i != j (different blocks) | 1/(4*(n-2k)*(n-2k-1)) |
| Endpoint B_i | Endpoint B_i | Same block, n-2k > 1 | 0 |
| Endpoint B_i | Endpoint B_i | Same block, n-2k = 1 (top degree) | 1/2 |
| Endpoint | Free | | 1/(2*(n-2k)*(n-2k-1)) |
| Free | Endpoint | | 1/(2*(n-2k)*(n-2k-1)) |
| Free | Free | a != b | 1/((n-2k)*(n-2k-1)) |

### (iii) Symmetry at even Walsh degrees (PROVED)

**Mc_hat[a,b,S] = Mc_hat[b,a,S]** for all even-degree monomials S.

This follows from the block-placement counting: the number of valid perms
from a to b equals the number from b to a, because reversing the macro-permutation
and flipping all block orientations creates a bijection.

### (iv) Parity separation of symmetric and antisymmetric parts

The counting matrix decomposes as Mc[a,b] = Sym[a,b] + Anti[a,b] where:
- **Sym** (symmetric part) has ONLY even Walsh degrees
- **Anti** (antisymmetric part) has ONLY odd Walsh degrees

Since H(T) = sum_{a,b} Mc[a,b] involves all entries, and Anti sums to zero,
H(T) depends only on the symmetric part, which has exclusively even Walsh degrees.

### (v) Rank structure

At Walsh degree 2k, each nonzero Mc_hat[:,:,S] has:
- **Rank = n - 2k** (observed, matches number of macro-items)
- Kernel = interior vertices (dimension 2k + r - 2r = 2k - r for non-endpoint interiors)

---

## Verification

| n | Degree | Type | Method | Status |
|---|--------|------|--------|--------|
| 5 | 0 | - | Exhaustive | PASS |
| 5 | 2 | P2 | Exhaustive + algebraic | PASS |
| 5 | 4 | P4 | Exhaustive | PASS |
| 7 | 2 | P2 | Permutation enum | PASS |
| 7 | 4 | P2+P2 | Permutation enum | PASS |

---

## Connection to PCD (THM-068) and OCF

The OPCD complements the diagonal PCD:
- **PCD** (THM-068): describes M_pos[v] = sum_P (-1)^{pos(v,P)} at each Walsh degree
- **OPCD** (this): describes Mc[a,b] = #{HPs from a to b} at each Walsh degree

These are different objects: M_pos is a signed sum over ALL HPs, while Mc counts
HPs by start/end vertex.

The OPCD symmetry at even Walsh degrees is a WEAKER statement than full symmetry
Mc[a,b](T) = Mc[b,a](T) for all T, since odd-degree asymmetry remains.
Full symmetry of Mc is NOT true in general (verified: asymmetric at odd Walsh degrees).

---

## Algebraic proof (degree 2)

For S = wedge with center c, endpoints e_1, e_2. Block B = [e_1, c, e_2] of width 3.

**Mc_hat[e_i, f_j, S]:** B first (e_i start), f_j last. Middle: (n-4)! arrangements.
Contribution: epsilon * (n-4)!.

**Mc_hat[f_i, f_j, S]:** f_i first, f_j last. B in middle: 2 orientations * (n-4)!.
Contribution: epsilon * 2 * (n-4)!.

**Normalization:** H_raw = epsilon * 2 * (n-2)!. So:
- ep -> free: (n-4)! / (2*(n-2)!) = 1/(2*(n-2)*(n-3))
- free -> free: 2*(n-4)! / (2*(n-2)!) = 1/((n-2)*(n-3))

At n=5: 1/12 and 1/6. At n=7: 1/40 and 1/20. Both confirmed numerically.

---

## Scripts

- `/tmp/pcd_offdiag_formula.py` -- Discovery of off-diagonal patterns
- `/tmp/opcd_algebraic.py` -- Algebraic verification at n=5,7
- `/tmp/pcd_offdiag_v3.py` -- Full Walsh analysis of Mc[a,b]
