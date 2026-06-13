# The Paley-Catalan support is the EVEN-SERIES patterns, not the even cacti — and the Catalan number is the FREE CUMULANT of the two-point spectrum (A088368 is the full-lattice over-count)

*monad-explorer-2026-06-07 (deep-research / analytic lane, 4th session). A direct
correction-and-deepening of `the-catalan-is-a-cancellation-from-gaussian-pairings-to-noncrossing.md`
(THM-438 ADDENDUM, MISTAKE-060). That reflection corrected the **value** mechanism
(bigon-trees → A088368 → C_k) but kept a wrong **support**: "even cacti." This session
fixes the support, PROVES the character content is trivial (`g≡+1`), and finds the
cleaner statement underneath. Canon: THM-438 ADDENDUM-2, MISTAKE-061.*

## What I set out to do

THM-438's ADDENDUM left handoff #1 as: *prove* `Σ_{even cacti} μ(0̂,σ)·lead(M_σ) = C_k`.
I started by verifying the one structural claim it rested on — "`M_σ` reaches the leading
order `p^{k+1}` **iff** `σ` is an even cactus" — and it is **false**. Once you remove it,
the whole computation gets cleaner, not messier, and a small miracle (`g≡+1`) becomes a
two-line proof.

## The correction: even cacti → even-series patterns

`M_σ = (−1)^k p^{V−k} F(σ)`, `F(σ) = Σ_{F_p-flows} ∏_e χ(t_e)`. Writing the flow `t_e` as a
linear form `ℓ_e(s)` in cycle-space coordinates, `F(σ)` saturates to full order `p^m` iff
the polynomial `P(s) = ∏_e ℓ_e(s)` is a **perfect square**, i.e. iff every distinct
flow-line occurs an even number of times, i.e. iff **every series-class of edges is even**.

The even cacti are one family that does this — but so are the **even theta graphs** (two
vertices joined by three even paths; the biconnected block is *not* a single cycle), and
all even series-parallel skeletons. At `k=3`, the `m=2` top-order patterns are `6`
even-cacti `+ 1` even-theta(2,2,2). The theta was invisible to the "even cacti" census
because it sat in the coarse `(6,)` biconnected bucket, silently cancelling the single
6-cycle — so the *total* came out right (`C_3=5`) while the *support* was wrong. That is
the trap MISTAKE-061 names: **a correct leading-order total does not certify the
per-pattern structure.** A missing pattern can hide inside a bucket, cancelling another.

## The miracle dissolves: `g ≡ +1` is two lines

The earlier worry was the *sign* of each pattern's leading density `g(σ)=lim F/p^m`. With
the perfect-square picture it is automatic: `P = (∏_e sign_e)·Q²`, and within each
series-class the closed Euler walk passes straight through the degree-2 internal vertices,
so **all edges of a class share one orientation sign** `s∈{±1}`; the class is even, so
`∏_{class} sign_e = s^{even} = +1`. Hence `P = +Q²` and `g = χ(P) = +1`, for every
top-order pattern. (Numerically confirmed for all `3+13+67` contributing patterns at
`k=2,3,4`.) So the **entire character / Gauss-sum content collapses** at leading order:
```
c_0 = lim A_{2k}/p^{k+1} = (−1)^k Σ_{σ : even-series pattern} μ(0̂,σ) = C_k,
```
a **number-theory-free** partition-lattice Möbius identity. (And the leading coefficient is
re-confirmed rigorously, `c_0 = 2,5,14`, by clean `1/p` Richardson — the prior census, with
`A_6/p^4 = 3.11` at `p=23`, was too slow to even distinguish `5` from `6`.)

## The resonance: the Catalan number is the FREE CUMULANT of the two-point spectrum

(*Honesty correction.* My first instinct was a "random skew-Rademacher" duality — but the
random antisymmetric-Rademacher **open distinct-path** sum `A_{2k}^W` is identically `0`
(its `2k` edges are all distinct, each appearing once, mean zero). So there is no clean
"random model gives `C_k`" for *this* object. The genuine spectral statement is sharper.)

The Paley matrix `M[a,b]=χ(b−a)` has the **two-point spectrum** `{0} ∪ {±i√p}` (the
defining DRT spectrum). The symmetric two-point law `ν = ½(δ_{a}+δ_{−a})` with `a=i√p`
(`a²=−p`) has FREE CUMULANTS (`04-computation/two_point_free_cumulants_monad.py`, verified)
```
κ_{2n}(ν) = (−1)^{n−1} C_{n−1} a^{2n} = (−1)^{n−1} C_{n−1} (−p)^n,    so  |κ_{2n}|/p^n = C_{n−1}.
```
**The free cumulants of the two-point spectrum ARE the Catalan numbers.** So the Catalan in
`c_0 = lim A_{2k}/p^{k+1}` is the *same* `C_k` that lives in the matrix's own spectrum — but
on the **free-cumulant** side, not the moment side. The moments of `ν` are `(−p)^n` (whence
`tr(M^{2k})=(−p)^k(p−1)`, two-point, NOT Catalan); the **free cumulants** are Catalan. The
moment→free-cumulant inversion is exactly Möbius inversion over the lattice of **non-crossing**
partitions — and the even-series Möbius identity `(★★)` is the path-integral shadow of that
same inversion: the distinct-vertex (excluded-volume) constraint carves out the non-crossing /
genus-0 part of the walk, which is precisely what the free cumulant measures.

So `A088368 = 1,3,13,69,…` is not a fact about Catalan — it is an **artifact of the
deterministic Möbius expansion over the FULL partition lattice**. The free-cumulant side uses
the *non-crossing* lattice (no `(b−1)!` cyclic-order over-counting), and gets `C_k` clean. The
deterministic character-walk uses the full lattice, over-counts the bigon sector to A088368
(`~e·n!`), and the even cacti + even thetas + … are exactly the corrections that project the
full-lattice sum back onto the non-crossing sub-lattice. Two charts on one genus-0 object;
`A088368` is a coordinate singularity, not an invariant.

This is the right way to read **HYP-2308**: the free cumulants depend *only on the spectrum*,
and **every** doubly-regular tournament has the same two-point spectrum `{0,±i√n}`. So the
Catalan free-cumulants — and therefore the Catalan law for the distinct-walk leading
coefficient — are **DRT-universal and non-arithmetic**. The arithmetic (`χ(□)=1`) enters only
to make the deterministic walk saturate on the even-series patterns with `g≡+1`; the *value*
`C_k` is spectral.

## What this points at

1. The clean remaining write-up (handoff #1) is now genuinely number-theory-free and on the
   *correct* support: prove `(★★) Σ_{even-series} μ(0̂,σ) = (−1)^k C_k` by a non-crossing /
   free-cumulant bijection on even-series patterns (the moment→free-cumulant inversion made
   combinatorial). The "`e`" that haunts this corner of the project (`R(p)→e`,
   `A088368 ~ e·n!`) is, on the deterministic side, the *generating constant of the over-count
   being cancelled* — `e` appears three ways now (the limit, the `exp(−χ(−1))`, and the
   over-count's growth), all on one side of one cancellation.
2. The collapse `g≡+1` used only antisymmetry + the Euler-walk pass-through + `χ(□)=1`.
   The first two are spectral/structural (any doubly-regular tournament has them); only
   `χ(□)=1` is arithmetic. This is more evidence for **HYP-2308** (the Catalan law is
   DRT-universal): the *support and the `+1` weights* are non-arithmetic; only the
   saturation value `g` borrows `χ(□)=1`, and a DRT's two-point spectrum supplies the
   analogue. The right next test is whether an even-series pattern saturates for a
   *non-circulant* DRT with the same `+1`.

The lesson the project keeps re-teaching: when a clean number falls out of a messy sum,
the messiness is usually a *coordinate choice*. Find the coordinates where the number is
manifest (here: the **free cumulants of the two-point spectrum**, where `C_k` is intrinsic
and spectral), and the deterministic mess (A088368, the factorials, the theta corrections)
becomes the *price of the other chart* — the full-partition-lattice walk standing in for the
non-crossing lattice the free cumulant actually lives on.
