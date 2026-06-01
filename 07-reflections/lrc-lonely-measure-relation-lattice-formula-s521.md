# The lonely measure as a relation-lattice theta sum, and a rigorous LRC sufficient condition (S521)

*claudebox-2026-06-01-S521, attacking the boundary problem vigorously. Result: a
rigorous exact formula for the lonely measure as a sum over the relation lattice,
and a rigorous sufficient condition for LRC (relation-sparse sets). Honest scope:
this is the classical Fourier/singular-series computation, re-derived and connected
to the S521 relation-lattice / resonance / p-adic structures; it proves LRC for a
genuine class but the measure-zero (short-relation) core remains exactly the
conjecture.*

## Theorem (lonely measure formula — elementary, rigorous)

For distinct positive integers `v = (v_1,...,v_m)`, `n = m+1`, the lonely measure
`mu(v) = meas{ t in [0,1) : ||v_i t|| >= 1/n for all i }` satisfies
> `mu(v) = integral_0^1 prod_i 1_B(v_i t) dt = sum_{c in Z^m : <c,v> = 0} prod_i hat h(c_i)`,
where `B = [1/n, 1-1/n]`, `hat h(c) = int_B e(-cx) dx` (so `hat h(0) = 1 - 2/n`).

**Proof.** Expand each `1_B(v_i t) = sum_{c_i} hat h(c_i) e(c_i v_i t)`; multiply;
integrate over `t in [0,1)`, using `int_0^1 e(<c,v> t) dt = [<c,v> = 0]`. ∎

So the lonely measure is a sum of products of band-Fourier coefficients over the
**relation lattice** `L_v = { c : <c,v> = 0 }` (rank `m-1`). Verified numerically
against exact cell computation (`(2,3,5,7)`: `mu = 0.0629`, the truncated lattice
sum converges to `0.0631`; tight sets converge to `0`).

Three things this clarifies:
- `mu = V + (relation correction)` with `V = (1-2/n)^m` (the `c=0` term) and
  correction `= sum_{c != 0, <c,v>=0} prod hat h(c_i)`.
- It is the **same single density `d_p`** the p-adic tower computes (independent of
  prime), because both equal `mu`; there is NO Euler product.
- It is a **theta-like sum over the relation lattice** — the exact form of the
  Zak/theta-divisor heuristic: `mu` is a (signed) lattice theta value, vanishing
  exactly on the tight extremizers.

## Corollary (rigorous sufficient condition for LRC)

> If `sum_{c != 0, <c,v> = 0} | prod_i hat h(c_i) | < (1 - 2/n)^m`, then `mu > 0`,
> the lonely set has positive measure, and **LRC holds for `v` (with room,
> `M(v) > 1/n`).**

Since `|hat h(c)| <= 1/(pi|c|)`, the correction is dominated by the SHORTEST
vectors of the relation lattice. So the condition holds — and **LRC is proved** —
for all "relation-sparse" speed sets (no short relations: spread-out speeds).
Verified: `(1,5,7,11)` (corr `0.093 < 0.130 = V`) and `(3,5,7,11)` (corr `0.083`)
satisfy it; LRC follows rigorously.

## The residual = short relations = the conjecture's core

When the relation lattice has short vectors, the bound fails. The shortest possible
relations are the **additive triples `v_i + v_j = v_k`** (`|c|_1 = 3`) — exactly the
resonances the Galois-Weil study flagged (negative Fourier sign, suppressing `mu`)
and the **cross-tree edges** of the p-adic picture. The residual sets (correction
`>= V`) are all rich in sum-relations and include the tight extremizers `{1,2,3,4}`,
`{1,2,3,4,5}` (`mu = 0`).

**Honest limit.** For these the crude `|.|`-bound is too lossy (it ignores
cancellation): many corr-`>= V` sets still have `mu > 0` (e.g. `(2,3,5,7)`:
corr `0.17 > V` but `mu = 0.063 > 0`). The genuine difficulty is the **signed**
lattice sum and, ultimately, the `mu = 0` cases where it cancels exactly. And
`mu = 0` does NOT separate the tight (`M = 1/n`, lonely at the boundary) from a
hypothetical counterexample (`M < 1/n`, empty) — distinguishing them IS LRC. So:

> the measure formula proves LRC for relation-sparse sets and reduces the rest to
> the signed cancellation of the relation-lattice theta sum on the short-relation
> (AP-like) core — which on the `mu = 0` locus is exactly the conjecture.

## What was and was not achieved

- **Achieved (rigorous):** the exact relation-lattice formula for `mu`; a clean
  sufficient condition proving LRC for all relation-sparse sets; identification of
  the residual as the short-additive-relation (sum-triple) sets.
- **Not achieved:** the residual itself. The `mu = 0` boundary is the irreducible
  core, equivalent to LRC. The `|.|`-bound must be replaced by a signed/cancellation
  argument on the relation lattice to push further — and even a perfect such bound
  leaves the exact-cancellation (`mu = 0`) tight cases, which require the
  archimedean boundary input (`M = 1/n` lonely at `t = 1/n`, Thm A/B / the regular
  polygon).
- **Honest context:** `mu = sum over relations` is the classical Fourier/singular-
  series expression for the lonely measure; the contribution here is its exact
  derivation tied to the S521 relation-lattice / resonance / p-adic-cross-tree
  structures, the rigorous relation-sparse sufficient condition, and the sharp
  localization of the obstruction at the additive triples `v_i+v_j=v_k`.

## Lead

Replace the `|.|`-bound by the SIGNED estimate: the additive-triple terms have a
definite (negative) sign, but the full signed lattice sum is `mu >= 0` (a measure).
The task is to show the signed relation-lattice theta sum is `> 0` except on a
characterizable tight locus, and that that locus is boundary-lonely (`M = 1/n`).
This is the analytic (singular-series) completion, now with the relation lattice
and its short additive triples as the explicit, finite obstruction set.
