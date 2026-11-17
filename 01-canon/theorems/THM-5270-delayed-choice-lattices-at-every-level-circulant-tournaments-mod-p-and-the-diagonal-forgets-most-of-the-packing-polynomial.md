---
id: THM-5270
title: "Packing-polynomial lattices and circulants. (A) At order n let P(q) = q^n + sum_{m in M} a_m q^m be the diagonal packing polynomial (THM-5269), M = {m = n (mod 2), 1 <= m <= n-2}, m_0 = min M. The integer kernel of the rows (c^m)_{m in M}, c = 1..k, has the Z-basis of band vectors (the coefficients of prod_{c=1}^k (X - c^2) on k+1 consecutive elements of M), so two tournaments of order n with equal P_1, ..., P_k differ in P_{k+1} by an element of L_n^{(k+1)} = (k+1)^{m_0} prod_{c=1}^k ((k+1)^2 - c^2) Z whenever n - 2 - m_0 >= 2k: 120 Z / 360 Z (level 3, n odd / even), 5040 Z / 20160 Z (level 4), 362880 Z / 1814400 Z (level 5); checked to n = 14. (B) A circulant tournament of prime order p has G_T(x, y) = y^p + ((p-1)/2) x (mod p) coefficientwise (the rotation fixes only the all-singleton packing and the (p-1)/2 constant-step Hamiltonian cycles), hence p | P_q(T) for every q, a_1 = -1 and a_m = 0 (mod p) for 1 < m < p; checked on every connection set of orders 3, 5, 7, 11, 13 (908 checks); the circulant counts are 3; 15; 175, 189; 92411, 93027, 93467, 95095; and six values from 3,669,497 to 3,711,175, with a_{p-2} and a_{p-4} constant across each order. (C) The diagonal P(q) forgets most of G: at orders 5, 6, 7, 8 there are 8, 24, 156, 1728 distinct diagonals against 8, 29, 203, 3007 distinct bivariate packing polynomials, and 0, 16, 196, 4503 classes whose diagonal is shared with a class of different cycle structure. (D) Typed: circulant Hadamard survivors {1, 4} (symmetry) and prime-count survivors {3} among vertex-transitive tournaments (symmetry: n | H and the floors) against the excluded counts {7, 21} (growth: floors and a finite kernel) are three finite exceptional sets of two kinds; the holes are the single-orbit counts 7·1 and 21·1 that orders 7 and 21 cannot have, as C_3's 3·1 is the only one they can."
status: >
  PROVED (klein-2026-11-16-circ): (A) by the triangular band basis of the Vandermonde kernel (lattices2.py builds the
  basis and compares with the formula for n <= 14; lattices.py found the level-3 values by an independent enumeration);
  (B) by Burnside on the rotation and Fermat's little theorem (circulant.py, every connection set, 0 failures). FINITE-EXACT:
  (C) (fibres2.py, packings enumerated directly on profiles8.npz). ANALOGY: (D), with the circulant Hadamard
  classification cited through the concurrent packet. Not independently audited. CITED: THM-5269, THM-002, THM-4861,
  THM-4747, THM-1370, and packing_hadamard_flatness_20261117 (the identity, the Hadamard defect in G, the orbit-phase
  reduction of which (B) is the circulant case, the (x - 2y) invisibility, the order-6 pair).
session: klein-2026-11-16-circ
date: 2026-11-16
source: 05-knowledge/results/circulant_lattices_fibres_20261116/LATTICES.md
scripts:
  - 05-knowledge/results/circulant_lattices_fibres_20261116/scripts/lattices2.py (+ receipts/lattices2.out)
  - 05-knowledge/results/circulant_lattices_fibres_20261116/scripts/circulant.py (+ receipts/circulant.out)
  - 05-knowledge/results/circulant_lattices_fibres_20261116/scripts/fibres2.py (+ receipts/fibres2.out)
related:
  - 01-canon/theorems/THM-5269-the-q-slot-cover-count-is-the-odd-cycle-packing-polynomial-at-2q-q-frozen-two-cover-fibres-delayed-choices-in-a-lattice-and-flexible-cores-under-the-one-port-alphabet.md
  - 01-canon/theorems/THM-002-ocf.md
  - 05-knowledge/results/circulant_lattices_fibres_20261116/HADAMARD.md
  - 05-knowledge/results/packing_hadamard_flatness_20261117/README.md
---
# THM-5270 — packing-polynomial lattices and circulant tournaments

**(A)** Write `m = m_0 + 2j`; row `c` of the constraint matrix is `c^{m_0}` times `(u_c^j)_j` with `u_c = c²`, so the kernel
is that of `k` Vandermonde rows in `u_1, …, u_k`. A band vector `v` (coefficients of `Π_c (X − c²)`) at start `j_0` satisfies
`Σ_i v_i u^{j_0+i} = u^{j_0} Π_c (u − c²) = 0` at every `u = c²`; the bands are triangular with leading entry 1, so they reduce
any integer kernel vector to one supported on the lowest `k` positions, which is zero by the square Vandermonde system.
The image of the band at `j_0` under `Δa ↦ Σ Δa_m (k+1)^m` is `(k+1)^{m_0+2j_0} Π_c ((k+1)² − c²)`; the gcd is the `j_0 = 0`
value. ∎ (`Π_c ((k+1)² − c²) = k!·(2k+1)!/(k+1)!`.)

**(B)** The rotation of `Z_p` permutes spanning packings; a packing containing a cycle of length `< p` lies in an orbit of
size `p` (no proper nonempty subset of `Z_p` is translation-invariant), so modulo `p` only the fixed packings count: the
all-singleton packing and the Hamiltonian cycles of constant step `d ∈ S`, one per element of the connection set. Hence
`G ≡ y^p + ((p−1)/2) x`, and `P_q = G(2q, q) ≡ q^p + (p−1)q ≡ 0` by Fermat. ∎

**(C), (D)** are LATTICES.md §2 and HADAMARD.md. Chan–Pak OPEN.
