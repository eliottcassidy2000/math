# The walk function is the complement shift — why a tournament's walk counts are spectral

**Author:** monad-explorer-2026-06-15-S8
**Builds on:** THM-507 (this session), the S7 skew-determinant reflection
(`the-skew-determinant-is-the-signed-even-face-and-it-is-spectral-monad-s7`),
HYP-2517, THM-505/506 (det/per face-by-face program), and the project's recurring `−1/2`
(THM-055/059/080, central factorials).
**Status:** the closed form and all corollaries are PROVED (THM-507) and verified with
exact integer/rational arithmetic.

---

## The discovery in one line

S7 reduced "the whole A-affine pencil is spectral" to a single unproved fact: that the
**walk counts** `w_k = 1ᵀAᵏ1` of a tournament are spectral. They are, and the reason is a
clean shift identity:

> For a tournament, `A − J = −(Aᵀ + I)`. So the all-ones rank-1 perturbation `A ↦ A − J`
> — the operation that *generically* destroys spectrality of walk counts — is, for a
> tournament, a **transpose-and-shift** whose eigenvalues `{−1−λᵢ}` are *forced by the
> spectrum alone*. The walk-generating function collapses to a ratio of two characteristic
> polynomials:
>
> **`F(x) = 1ᵀ(xI−A)⁻¹1 = ∏ᵢ(x+1+λᵢ) / ∏ᵢ(x−λᵢ) − 1`,   i.e.   `1ᵀadj(xI−A)1 = (−1)ⁿcharA(−x−1) − charA(x)`.**

That `+1` shift inside the numerator is the **complement**. Everything below is a reading
of that single `+1`.

## Why this is not obvious — the cospectral walk obstruction

For a *general* graph, walk counts are **not** spectral. The smallest cospectral pair,
`C₄ ⊔ K₁` and `K_{1,4}`, both have spectrum `{2,0,0,0,−2}`, yet `w₂ = Σdegᵢ²` is `16`
versus `20`. Walk counts depend on the *main eigenvalues* and on how the all-ones vector
`1` is angled against the eigenspaces — data the spectrum does not carry. In the
"walk-matrix / main-eigenvalue" theory, the rank-1 update `A ↦ A − J` (or `A ↦ A + J`)
scrambles the spectrum in a way that *depends on those angles*.

The tournament miracle is that for `A + Aᵀ = J − I` (complement = converse), the angle
dependence never gets to act: `A − J = −(Aᵀ + I)` has eigenvalues pinned at `{−1−λᵢ}` no
matter what the eigenvectors do, because transposition preserves the spectrum and `+I` is a
scalar shift. The proof (THM-507) is two lines — matrix-determinant lemma plus this
identity — and it uses **nothing** about tournaments except `A + Aᵀ = J − I`.

This sharpens *which structural fact* makes tournaments special. It is not "0/1 entries"
and not "no 2-cycles" — it is precisely **complement = converse**.

## The complement fixed point

There is a general identity (Sherman–Morrison) for *any* `A' = J − I − A`:
`F_{A'}(x) = h(x)/(1 − h(x))`, with `h(x) = 1ᵀ((x+1)I + A)⁻¹1 = −F_A(−x−1)`.
For a tournament the complement equals the converse, `A' = Aᵀ`, and the walk function is
transpose-invariant (`F_{Aᵀ} = F_A`, it is a scalar). So the general complement relation
becomes a **self-reference**, and solving it gives the reciprocity

> **`(1 + F(x))(1 + F(−x−1)) = 1`.**

In words: *the tournament walk function is the fixed point of the graph-complement walk
map.* The reciprocity is its fixed-point equation. The closed form is the explicit fixed
point. (Verified exactly at rational points, `walk_counts_trace_formulas_monad_s8.py`.)

## The `−1/2` again

The involution `x ↦ −x−1` has fixed point `x = −1/2`. Recentre: `g(t) := 1 + F(t − 1/2)`.
Then the reciprocity becomes the clean odd symmetry

> **`g(t)·g(−t) = 1`,**   with   `g(t) = (−1)ⁿ q(−t)/q(t)`,   `q(t) := charA(t − 1/2)`.

So the natural variable for the walk function is `t = x + 1/2`, the characteristic
polynomial **recentred at `−1/2`**. This is the *same* `−1/2` that organises the
central-factorial / Eulerian structure of `W(r) = tr M(r)` (THM-055/059/080, via
`u = (r+1/2)(r−1/2) = r² − 1/4`). Two independent corners of the project — the transfer
matrix `M(r)` and the walk function `F(x)` — are both naturally centred at `−1/2`. The
common cause is the *same* relation `A + Aᵀ = J − I`: complementation sends `λ ↦ −1−λ`,
whose centre of symmetry is `−1/2`. **`−1/2` is the fixed point of complementation on the
spectral axis.**

## The walk counts wear their cycle counts on their sleeve

Expanding the closed form (using the tournament facts `tr A = 0` and `tr A² = 0`, no
2-cycles), the walk counts come out as

> `w₀ = C(n,1) = n`,  `w₁ = C(n,2)`,  `w₂ = C(n,3) + 2c₃`,  `w₃ = C(n,4) + (2n−3)c₃`, …

The leading term `C(n,k+1)` is exactly the walk count of the **transitive** tournament
(walks of length `k` = strictly increasing `(k+1)`-tuples). The corrections are the
**spectral cycle counts** `c₃ = tr(A³)/3, c₄, …`. So:

> **walk count = transitive baseline `C(n,k+1)` + spectral cycle-count correction.**

This is the open-walk analogue of the closed-walk census `tr(Aᵏ)` (THM-502): both the
open census `1ᵀAᵏ1` and the closed census `tr Aᵏ` live entirely on the spectral shore. The
contrast with `H = I(Ω,2)` is the whole point — `H` counts cycle *packings*, and *that* is
where non-spectrality enters (THM-505/506). Walks are spectral; packings are not.

## Where it sits in the det / per / spectral program

THM-506 drew the map: signed faces of the master packing polynomial `Φ` are determinants
(spectral, P); unsigned faces are permanents (non-spectral, #P). S7 showed every A-affine
*signing* lands on the determinant side, conjecturally. THM-507 makes "conjecturally"
into "provably": **no determinant of an A-affine matrix can see non-spectral content**,
because every such determinant reduces — by the matrix-determinant lemma — to `charA` times
the walk function, and the walk function is `charA` in disguise. The entire signed /
determinantal / linear-algebraic world is a closed (spectral) book. To reach `H` you must
cross to the permanent / #P side. The Valiant `det`∈P vs `per`∈#P wall is the
spectral/non-spectral wall, and THM-507 nails the easy shore shut.

## The striking refinement, explained

The score sequence — and every score power-moment `Σsᵢᵖ` for `p ≥ 3` — is *not* spectral
(it splits 14 cospectral classes at n=6, 85 at n=7). Yet the walk counts `w_k = 1ᵀAᵏ1`
*are*. The closed form says why: `w_k` is the *fully contracted* quantity `1ᵀAᵏ1`, which
the complement shift expresses through `charA`; the score moments `Σ(Aᵏ1)ᵢᵖ` are *diagonal*
data that retain exactly the eigenvector-angle information the contraction throws away. The
spectral/non-spectral line runs *between* `1ᵀAᵏ1` and `Σ(Aᵏ1)ᵢᵖ` — between contracting the
walk vector and reading its entries.

## Handoff / what to chase next

1. **The de-contraction wall — located (this session).** Probing how spectrality dies as
   one de-contracts the resolvent `R = (xI−A)⁻¹` (verified on all 3 genuine cospectral,
   different-`H` classes at n=6, `decontraction_wall_monad_s8.py`):

   | functional | spectral? |
   |---|---|
   | `1ᵀR1` (full contraction) | **yes** (THM-507) |
   | `tr R = charA'/charA` | **yes** (log-derivative) |
   | `M₂ = Σₐ(R1)ₐ²` (2nd moment of resolvent row sums) | **yes** |
   | `‖R‖_F² = tr(RRᵀ)` | **yes** |
   | `M₃ = Σₐ(R1)ₐ³` (3rd moment) | **NO** (splits all 3 classes) |

   The wall is between the 2nd and 3rd moment of the resolvent row-sum vector `r = R1` —
   the *exact resolvent-space echo* of the classical score wall "`Σsₐ²` spectral (Moon),
   `Σsₐᵖ` (`p ≥ 3`) not." And it is not a coincidence of `x`: expanding
   `(R1)ₐ = Σ_k(Aᵏ1)ₐ x^{−k−1}`, the leading nontrivial Taylor coefficient of `M₂` is
   `Σₐ sₐ² = (A1)ᵀ(A1)` — spectral by Moon — whereas `M₃` carries `Σₐ sₐ³`, the first
   non-spectral score moment. So `M₂ ⊃ Σsₐ²` and `M₃ ⊃ Σsₐ³`: the resolvent moments
   *contain* the score moments as their leading coefficients, and inherit the wall from
   them. (Note `tr(RRᵀ)` spectral is itself a small surprise — `RRᵀ` mixes `A` and `Aᵀ`,
   and stays spectral via the same complement mechanism; the two-sided walk counts
   `1ᵀ(Aᵀ)ⁱAʲ1` deserve their own closed form.) **Open:** is every two-sided walk count
   `1ᵀ(Aᵀ)ⁱAʲ1` spectral (a double complement-shift)? and where is the pointed
   `1ₐᵀ R 1_b` / `M[a,b]` (THM-080) relative to this wall? — that is the boundary between
   the spectral skeleton and the genuine eigenvector-angle data.
2. **Other complement-stable families.** The proof is purely "complement = converse." Which
   *other* matrices `A` with `A + Aᵀ = cJ + dI` inherit a spectral walk function? (c-tournaments
   already appear in `walk_gf_symmetry.py`.) This isolates the algebraic core away from
   tournaments entirely.
3. **The permanental roots** (still open from S6): the non-spectral content lives on the
   permanent side; are the *roots* of `per(xI+A)` a cleaner carrier of `(c₆,c₇,…)` than the
   coefficients? THM-507 guarantees the determinant side gives nothing new, so the entire
   fingerprint frontier is permanental.
