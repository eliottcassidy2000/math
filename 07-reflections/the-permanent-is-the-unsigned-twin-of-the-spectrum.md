# The permanent is the unsigned twin of the spectrum: det/per is the home of the non-spectral defect

*monad-explorer-2026-06-15-S6. Builds directly on the reflection
`the-master-cycle-packing-polynomial` (S5) and THM-505. That reflection identified the
master cycle-packing polynomial `Φ(T;{y_k}) = Σ_{packings} ∏ y_{|C|}` and showed the
spectrum (Sachs, signed `y_k=−x^{−k}`) and the OCF `H` (odd-only unsigned fugacity-2) are
two specializations of it. Its handoff (2) asked: "Is `I(Ω_all, 2)` a known graph
polynomial?" This reflection answers that — and finds the structural reason the spectrum
cannot see `H`.*

## The answer: the all-length unsigned face of Φ is the permanental polynomial

The characteristic polynomial is the **signed, vertex-graded** face of `Φ`. Its unsigned
twin — same packings, signs dropped — is a second classical matrix function:

> **det(xI − A) = Σ_L (−1)^{#cyc(L)} x^{n−|V(L)|}**   (SIGNED → spectrum / char. poly.)
> **per(xI + A) = Σ_L (+1)^{#cyc(L)} x^{n−|V(L)|}**   (UNSIGNED → permanental poly.)

where `L` ranges over linear subdigraphs (vertex-disjoint directed-cycle packings; length
`≥ 3` in a tournament). Reading off coefficients,

> `[x^{n−m}] det(xI−A) = e_m^{signed}   = Σ_{|V(L)|=m} (−1)^{#cyc}`   (Sachs / char poly),
> `[x^{n−m}] per(xI+A) = e_m^{unsigned} = Σ_{|V(L)|=m} (+1)^{#cyc}`   (**permanental** poly).

So the "all face" of `Φ`, graded by vertices covered, **is the permanental polynomial of
the tournament** — one of the most-studied non-spectral graph polynomials in the literature
(permanental polynomials, permanental roots). The proof is the permanent's own
cycle-decomposition: `per(xI+A) = Σ_σ ∏_i (xI+A)[i,σ(i)]`, where a permutation `σ`
contributes `x` per fixed point and `1` per non-trivial cycle iff that cycle is a directed
cycle of `T` (and `0` otherwise). Tournaments have no 2-cycles, so the surviving `σ` are
exactly the linear subdigraphs — matching `Φ` term-for-term. (Verified: Ryser's
`per(xI+A)` over `Z[x]` equals the independent packing enumeration, all samples n=3..7;
`permanental_companion_monad.py`.)

The determinant and the permanent of the **same** matrix (up to the sign on `A`) differ by
*exactly one thing*: the cycle-parity sign `(−1)^{#cyc} → (+1)^{#cyc}`. The reflection on
`Φ` named three channels through which information hides from the spectrum — (a) signs, (b)
the even/odd split, (c) length-vs-cycle-count grading. **The permanent strips channel (a)
and nothing else.** It is the *spectrum with the cycle-parity signs removed*.

## Why this is the home of the non-spectral defect

`det` is computable in polynomial time (eigenvalues / Gaussian elimination); `per` of a 0/1
matrix is `#P`-hard (Valiant 1979). The whole project slogan — "the spectrum is mean-field,
the OCF is correlation; the eigenvalues cannot see `H`" — acquires its natural ambient
explanation: the spectral invariants live on the **determinant** side (signed, tractable,
symmetric functions of eigenvalues), and the non-spectral carriers of `H` live on the
**permanent** side (unsigned cycle-cover counts). The non-spectrality of `H`, of `c₆`
(THM-502), of every face of `Φ` — these are shadows of the `det`/`per` gap. *(Honest scope:
the rigorous content is that `det` is a symmetric function of eigenvalues while `per` is the
non-spectral unsigned packing count; the `#P`-hardness is Valiant's for general 0/1
matrices and is invoked as the ambient structural reason, not claimed for the
tournament-restricted permanent.)*

## The same wall — n=6 — and the same split, but strictly finer

The permanental polynomial is non-spectral on the **same schedule as `H`**: it is constant
on cospectral classes for `n ≤ 5` and first splits a cospectral class at **n = 6** — the
disjoint-pair threshold `6 = 3+3`, exactly where `H`, `c₆`, and every face of `Φ` first turn
non-spectral. At `n = 7` it splits the **same 47** cospectral classes that `H` does
(`permanental_companion_monad.py` [3]: per and `H` both split 3/28 at n=6, 47/168 at n=7).

But `per` is strictly **finer** than `H` as a fingerprint. The pair `(det, per)` fixes, for
each vertex count `m`,

> `O_m := (e_m^{unsigned} − e_m^{signed})/2 = #{odd-cardinality packings on m vtx}`,
> `E_m := (e_m^{unsigned} + e_m^{signed})/2 = #{even-cardinality packings on m vtx}`,

i.e. it resolves the **full carrier vector** `(c₆, c₇, …)` — whereas `H = I(Ω,2)` reads only
the single fugacity-2 functional `4c₆ + 2c₇ + …`. Concretely at the first wall:
`O_6 = c₆` (the only odd-cardinality packing on 6 vtx is a single 6-cycle), `E_6 = D₃₃`
(the only even-cardinality one is a disjoint triangle pair), and combined with the spectral
`e_6^{signed} = D₃₃ − c₆` this recovers `c₆` **and** `D₃₃` separately. The spectrum gives
only their difference; the permanent supplies exactly the missing non-spectral carrier.

Exhaustive `n = 6` (all 32 768 labeled tournaments, 56 iso classes,
`permanental_companion_monad.py` [4]):

| fingerprint | # distinct values |
|---|---|
| iso classes (A000568) | 56 |
| char poly (det) | 28 |
| (char, perm) | **32** |
| (char, perm, H) | 32 |

The permanent recovers **4** distinctions the spectrum loses; `H` adds **nothing** beyond
the permanent (`32 → 32`). So `(char, perm)` **determines `H`** at `n = 6`, and more — the
permanental polynomial **dominates `H` as a tournament fingerprint**. (Engineering, domain
12: a single permanent `per(xI+A)`, `O(2ⁿ n)` by Ryser, is a strictly stronger fingerprint
than the Hamiltonian-path count `H`, and computes the whole non-spectral carrier vector at
once.)

## O_m, E_m vs H's carriers: where the two gradings agree and where they split

The pair `(det, per)` grades packings by **(vertex-count `m`, cardinality parity)**; `H`
grades them by **(odd lengths, exact cardinality, via `2^{#cyc}`)**. They coincide on the
carriers `H` actually needs, until two packings of different exact cardinality but the same
vertex count and the same cardinality-parity collide:

- A single odd cycle `c_k` is an odd-cardinality packing — `O_k = c_k` while `k` is the only
  vertex count reachable by one cycle, i.e. for `k ≤ 8`. At `m = 9` the triple `T₃₃₃` (three
  disjoint triangles, cardinality 3, also **odd**) shares `m = 9` with the 9-cycle `c₉`:
  `O_9 = c₉ + T₃₃₃`. The vertex grading **merges** them; `H` keeps them apart
  (weights `2·c₉` vs `8·T₃₃₃`).
- A disjoint **odd** pair (`D₃₃, D₃₅`, cardinality 2, **even**) shares its vertex count with
  an **even-cycle** pair: at `m = 8`, `E_8 = D₄₄ + D₃₅`. `H` needs `D₃₅` (it is part of the
  odd-pair level `α₂`); the vertex grading lumps it with `D₄₄`.

So `(det, per)` determines `H` while these collisions are absent, and the first collision
fixes where it must break.

**The verdict (within-cospectral-class exact-`ℚ`-rank test, `perm_rank_test_monad.py`).**
Within a cospectral class the spectrum is fixed, so `H` and the perm-poly coordinates
`e_m^{unsigned}` are both linear in the non-spectral carriers; `(char, perm)` determines
`H` iff `H` is affine in `e_m^{unsigned}`, i.e. iff `rank[Δe^{unsigned} | ΔH] =
rank[Δe^{unsigned}]`. Over thousands of within-class deltas per `n`:

| `n` | `rank[Δe^{unsigned}]` | `rank[Δe^{unsigned} | ΔH]` | `(char,perm) → H`? |
|---|---|---|---|
| 6 | 1 | 1 | **determines** |
| 7 | 2 | 2 | **determines** |
| 8 | 3 | **4** | **does NOT** |
| 9 | 4 | **5** | **does NOT** |

So **`(char, perm)` determines `H` exactly for `n ≤ 7`, and first fails at `n = 8`** — and
the rank-diagnosis pinpoints the cause: at `n = 8` adding either `D₄₄` *or* `D₃₅` to the
perm coordinates restores the affine relation (they are interchangeable because
`D₄₄ + D₃₅ = E_8` is fixed). This is precisely the predicted collision: the disjoint odd
pair `D₃₅` (which `H` needs, as part of `α₂`) and the disjoint even-cycle pair `D₄₄` (which
`H` does not) share the vertex count 8 and the cardinality-parity (both 2 cycles), so the
permanent's **vertex grading** cannot separate them while `H`'s **odd-length** grading
must. The break is the first place the two gradings of `Φ` genuinely diverge.

A clean by-product: `rank[Δe^{unsigned}] = n − 5` for `n = 6..9` — the permanental
polynomial's within-cospectral-class non-spectral dimension is `n − 5`, carried by the
`n − 5` coordinates `e_6^{unsigned}, …, e_n^{unsigned}`. This is the natural home of
THM-505's *original* "non-spectral dimension `n − 5`" observation: that count was always
the **carrier-coordinate** dimension (the simple-cycle counts `c_6,…,c_n`), which is exactly
what the perm poly resolves — not `H`'s own functional dimension `⌊n/3⌋`, which is smaller
because `H` reads only level-sum marginals.

## The lattice of faces, re-seen

`Φ` is the grand object; its faces are obtained by choosing **(sign, lengths, grading)**:

| face | sign | lengths | grading | object |
|---|---|---|---|---|
| char poly `det(xI−A)` | signed | all | vertices | **spectrum** |
| perm poly `per(xI+A)` | unsigned | all | vertices | **permanental polynomial** |
| `I(Ω,x)` | unsigned | odd | cycle-count (fugacity `x`) | independence poly of `Ω` |
| `H = I(Ω,2)` | unsigned | odd | fugacity 2 | **OCF / Rédei** |
| `I(Ω,−1) = −χ̃(Ind Ω)` | signed | odd | fugacity −1 | Euler char of packing complex |

The spectrum and the permanent are the two **vertex-graded** faces — the determinant and the
permanent of one matrix, differing by the cycle-parity sign. `H` and the Euler characteristic
are two **cycle-count-graded** faces of the odd sub-polynomial, differing by the fugacity
(`+2` vs `−1`). The project has been reading one tournament through five windows onto a
single packing polynomial; naming the windows `det` and `per` places the spectral/
non-spectral wall inside the oldest dichotomy in the theory of matrix functions.

## Handoffs

1. **The general-`n` carrier deficit of the perm poly** (the break is now located at
   `n = 8`, via `D₄₄ ↔ D₃₅`). Past the first collision, how many carriers must be adjoined
   to `(char, perm)` to recover `H` at each `n`? The collisions are the length-multisets
   that share a vertex count and cardinality-parity; the deficit grows as these multiply.
   Is the deficit `= (#packing-multisets on `≤ n` vtx) − (#perm coords) − (spectral)`?
2. **The even face.** `per(xI+A)` is the *all*-length unsigned face. Is there a matrix
   function for the **even**-only face `I(Ω_even, ·)`? (No obvious one — even-cycle packings
   are not a clean permanent restriction. Possibly a Pfaffian / half-determinant on a
   derived bipartite graph.)
3. **Permanental roots.** The roots of `per(xI+A)` ("permanental spectrum") are a genuine
   non-spectral invariant. Do cospectral-but-not-co-permanental tournaments have distinct
   permanental roots, and do those roots carry the `(c₆,c₇,…)` carriers more transparently
   than the coefficients?
4. **Fingerprint program (engineering, domain 12).** `(char poly, perm poly)` is a cheap,
   strictly-stronger-than-`H` fingerprint. How close to complete is it at `n = 7, 8`? (At
   `n = 6` it gives 32 of 56 iso classes — still far from complete, but it strictly dominates
   both the spectrum and `H`.) Pair it with one more carrier (`T₃₃₃`?) past the first
   collision.
