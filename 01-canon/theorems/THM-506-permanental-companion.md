---
id: THM-506
title: The permanental companion of the OCF — per(xI+A) is the unsigned twin of the char poly (the all-length vertex-graded face of Phi); (char,perm) determines H iff n<=7, breaking at n=8 via the D44<->D35 trade
status: PROVED (permanental coefficient theorem, det/per O/E split — proofs by permutation cycle decomposition / Ryser). VERIFIED computationally: coeff theorem + O/E split all samples n<=7; non-spectral first at n=6 (same wall as H), splits same 47 cospectral classes as H at n=7; (char,perm) DETERMINES H for n<=7 and does NOT for n=8,9 (within-class exact-Q-rank: rank[De_unsigned|DH] exceeds rank[De_unsigned] by 1 at n=8,9; break = D44<->D35). Engineering: (char,perm) strictly dominates H as a fingerprint (28->32 distinct vs 56 iso classes at n=6; H adds nothing beyond perm).
source: monad-explorer-2026-06-15-S6
depends_on:
  - THM-505   # OCF non-spectral defect: H = spectral skeleton + carriers; Sachs-basis form
  - THM-502   # closed-walk census; c6 non-spectral from onset; Witt transform
  - THM-499   # H = I(Omega,2); alpha_2 first non-spectral OCF ingredient at n=6
related:
  - HYP-2515
  - HYP-2514  # master cycle-packing polynomial Phi (spectrum & H as two specializations)
  - OPEN-Q-096
  - reflection: the-permanent-is-the-unsigned-twin-of-the-spectrum
  - reflection: the-master-cycle-packing-polynomial
---

# THM-506 — the permanental companion of the OCF

The master cycle-packing polynomial (HYP-2514, reflection `the-master-cycle-packing-polynomial`)
is `Phi(T;{y_k}) = sum_{L linear subdigraph} prod_{C in L} y_{|C|}`, with the spectrum the
**signed** vertex-graded face (`y_k=-x^{-k}`, Sachs) and `H=I(Omega,2)` the unsigned odd-only
fugacity-2 face. This theorem identifies the **unsigned, all-length, vertex-graded** face as a
second classical matrix function — the permanental polynomial — and uses it to bound exactly
how much of `H` is recoverable from the spectrum-plus-permanent.

Throughout, `A` is the 0/1 adjacency of a tournament `T` on `n` vertices; a *linear
subdigraph* `L` is a set of vertex-disjoint directed cycles (length `>= 3`, since tournaments
have no 2-cycles); `#cyc(L)` is the number of cycles; `|V(L)|` the number of covered vertices.

## Part 1 — the permanental coefficient theorem (PROVED)

> **per(xI + A) = sum_{m=0}^{n} e_m^{unsigned} x^{n-m},   e_m^{unsigned} = #{ L : |V(L)| = m }**,

the **unsigned, all-length** packing count covering exactly `m` vertices. Equivalently
`per(xI+A) = sum_L x^{n-|V(L)|}` — the unsigned twin of the Sachs coefficient theorem
`det(xI-A) = sum_L (-1)^{#cyc(L)} x^{n-|V(L)|}`, `e_m^{signed} = sum_{|V(L)|=m}(-1)^{#cyc}`.

**Proof.** `per(xI+A) = sum_{sigma in S_n} prod_i (xI+A)[i,sigma(i)]`. Decompose `sigma` into
cycles: a fixed point contributes `(xI+A)[i,i] = x`; a cycle of length `>= 2` contributes
`prod A[i,sigma(i)] in {0,1}`, equal to 1 iff that cycle is a directed cycle of `T`. In a
tournament there are no 2-cycles, so the surviving `sigma` are exactly the linear subdigraphs
`L` (non-fixed cycles = the cycles of `L`), each contributing `x^{#fixed} = x^{n-|V(L)|}`. ∎
(The same cycle decomposition gives `det(xI-A)`'s Sachs form, with `sgn(sigma)` producing the
`(-1)^{#cyc}`.) VERIFIED: Ryser `per(xI+A)` over `Z[x]` equals the independent packing
enumeration, all samples n=3..7 (`permanental_companion_monad.py` [1]).

## Part 2 — det/per is the cycle-parity SIGN channel; the O/E split (PROVED)

`det(xI-A)` and `per(xI+A)` are the **signed** and **unsigned** faces of `Phi`, differing
only by `(-1)^{#cyc} -> (+1)^{#cyc}`. Hence the pair determines, for every `m`,

> **O_m := (e_m^{unsigned} - e_m^{signed})/2 = #{ odd-cardinality packings on m vtx }**,
> **E_m := (e_m^{unsigned} + e_m^{signed})/2 = #{ even-cardinality packings on m vtx }.**

VERIFIED (all samples n<=7, `permanental_companion_monad.py` [2]). The permanent is the
*spectrum with the cycle-parity signs stripped*: `det` is a symmetric function of the
eigenvalues (spectral, polynomial-time), `per` is the non-spectral unsigned cycle-cover
count (Valiant: per of a 0/1 matrix is `#P`-hard — the ambient structural reason the
spectrum cannot see the unsigned data; not a hardness claim for the tournament-restricted
permanent).

## Part 3 — same wall, strictly finer fingerprint (VERIFIED)

The permanental polynomial is constant on cospectral classes for `n <= 5` and first splits a
cospectral class at **n = 6** (the disjoint-pair threshold `6 = 3+3`) — the **same wall** as
`H`, `c_6`, and every face of `Phi`. At `n = 7` it splits the **same 47** cospectral classes
`H` does (both 3/28 at n=6, 47/168 at n=7). But it is strictly **finer**: `(det, per)`
resolves the full carrier vector via `O_m, E_m` (e.g. `O_6 = c_6`, `E_6 = D_33`, recovering
`c_6` and `D_33` separately from the spectral `e_6^{signed}=D_33-c_6`), whereas
`H = I(Omega,2)` reads only the fugacity-2 functional `4c_6+2c_7+...`.

Exhaustive `n = 6` (all 32768 labeled tournaments; `permanental_companion_monad.py` [4]):
iso classes (A000568) = 56; distinct char poly = 28; distinct `(char, perm)` = **32**;
distinct `(char, perm, H)` = 32. The permanent recovers 4 distinctions the spectrum loses;
`H` adds nothing beyond the permanent. **The permanental polynomial dominates `H` as a
tournament fingerprint.**

## Part 4 — (char, perm) determines H iff n <= 7 (VERIFIED; break = D44 <-> D35)

> **The pair (characteristic polynomial, permanental polynomial) determines `H` for `n <= 7`
> and does NOT for `n >= 8`.**

Method: within a cospectral class the spectrum is fixed, so `H` and the perm coordinates
`e_m^{unsigned}` are both linear in the non-spectral carriers; `(char,perm) -> H` iff `H` is
affine in `e^{unsigned}`, i.e. `rank[ De^{unsigned} | DH ] = rank[ De^{unsigned} ]` over
within-class deltas (exact `Q`-rank, `perm_rank_test_monad.py`):

| n | rank[De_unsigned] | rank[De_unsigned | DH] | (char,perm)->H |
|---|---|---|---|
| 6 | 1 | 1 | determines |
| 7 | 2 | 2 | determines |
| 8 | 3 | 4 | does NOT |
| 9 | 4 | 5 | does NOT |

(thousands of within-class deltas per n; 47 classes split H at n=7, 657 at n=8, 1113 at n=9.)
The rank diagnosis at n=8: adjoining **either `D_44` or `D_35`** to the perm coordinates
restores the affine relation (interchangeable since `D_44+D_35 = E_8` is fixed). This is the
first **packing collision**: the disjoint odd pair `D_35` (which `H` needs — part of the
odd-pair level `alpha_2`) and the disjoint even-cycle pair `D_44` (which `H` does not) share
the vertex count 8 and the cardinality-parity (both 2 cycles), so the permanent's **vertex
grading** merges them while `H`'s **odd-length** grading must separate them. (At `n=9` the
analogous odd-cardinality collision `c_9 <-> T_333` adds a second deficit dimension.)

By-product: `rank[De^{unsigned}] = n - 5` for n=6..9 — the permanental polynomial's
within-cospectral-class non-spectral dimension is `n-5`, carried by the `n-5` coordinates
`e_6^{unsigned},...,e_n^{unsigned}`. This is the natural home of THM-505's *original*
"non-spectral dimension n-5" observation (the carrier-coordinate count `c_6,...,c_n`), as
opposed to `H`'s own functional dimension `floor(n/3)` (smaller; `H` reads level-sum
marginals).

## Files
- `04-computation/permanental_companion_monad.py` (+ `05-knowledge/results/permanental_companion_monad.out`)
  — Parts 1,2,3 + fingerprint completeness.
- `04-computation/perm_rank_test_monad.py` (+ `.out`, `perm_rank_test_n9.out`) — Part 4.

## Open
- The even-only face `I(Omega_even,.)` has no obvious single matrix function (possibly a
  Pfaffian on a derived graph). The permanental ROOTS as a non-spectral invariant. The
  general-n carrier deficit of `(char,perm)` past the first collision. (OPEN-Q-096; HYP-2515.)
