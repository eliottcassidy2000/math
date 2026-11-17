# Pre-registration ledger — klein-2026-11-16-circ

Predictions in [prediction_circ.txt](prediction_circ.txt), two batches, each written before the runs it covers.

| # | prediction | outcome |
|---|---|---|
| H-1 | circulant congruence on every connection set of orders 3, 5, 7, 11, 13: `p | P_q`, `a_m ≡ 0` (`1 < m < p`), `a_1 ≡ −1` | **confirmed** (908 checks, 0 failures) |
| H-2 | `ΔP_3` lattice `120ℤ` (odd `n`), `360ℤ` (even `n`) for `7 ≤ n ≤ 14` | **confirmed** |
| H-3 | distinct diagonal polynomials at order 8 in [3000, 5000] | **refuted** (1728) |
| H-4 | one orbit of paths only for `C_3`; least circulant count at order 11 ≥ 5000; `P_11` in [10⁴, 10⁵] and divisible by 11 | **confirmed** (92,411; 95,095 = 11·8645) |
| H-5 | classes at order 8 whose diagonal is shared with a different `G` in [100, 688]; such classes at orders 5 and 6 | **refuted** (4503; none at order 5, 16 at order 6) |
| H-6 | `L_n^{(k+1)} = (k+1)^{m_0} Π_{c≤k} ((k+1)² − c²) ℤ` for `k = 2, 3, 4`, `n ≤ 14` | **confirmed** (band basis = formula; the small-box brute force found only multiples where its box was too small) |

Refuted clauses are recorded as refuted; none was revised after the fact. The parenthetical of H-3 ("= distinct packing
polynomials") was corrected in the second batch before any run on it: three evaluations fix the diagonal, not `G`.
