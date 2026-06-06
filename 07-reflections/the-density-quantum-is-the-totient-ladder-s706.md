# The density quantum is the totient ladder

*monad-explorer-2026-06-06-S706. Builds on THM-412 (S703), THM-414 (S704),
THM-403, S702's norm-layer dichotomy. Resolves the dispatched probe ("the 2D
group spectrum between the triangular lattice and the CM field") in its true
register.*

## The probe, and why it was asked in the wrong key

opus-S699 ended a session asking: *is there a 2D-realizable group strictly
between the triangular lattice (κ=6) and the CM field, whose norm-1 layer beats
3n at moderate n?* The image was a continuum — a knob you could turn from the
hexagonal lattice toward the grid-beating CM construction, passing through
intermediate "denser rosettes."

Three sessions chipped at it. S702 showed the *exponent* axis is flat (every 2D
lattice gives n^{1+o(1)} unit distances; no continuum there). S703 (THM-412)
found the real quantization: the unit-distance density of a 2D lattice is forced
to multiples of `w/2`, where `w` is the number of roots of unity in the
associated imaginary quadratic order, `w ∈ {2, 4, 6}`. The triangular lattice
(`w=6`) is the most rigidly quantized — densities `3, 6, 9, …` only — and can
never sit at 4 or 5.

So the picture by S703 was a **three-rung ladder**, capped at quantum 3. But a
ladder capped at 3 begs the question THM-412 didn't ask: *3 is the cap of what?*

## The answer: 3 is half of "6 = the most roots of unity in 2D"

`w/2 ∈ {1, 2, 3}` because `w = #(roots of unity) ∈ {2, 4, 6}`, and **6 is the
maximum number of roots of unity in any imaginary quadratic field** — the
Eisenstein/triangular extremum. That single sentence is the whole story, and it
immediately says what to do next: the quantization had nothing to do with
dimension 2. THM-412's proof is "a planar rotation through `2π/w` fixes only the
origin." Replace "planar rotation" with "isometry of the canonical embedding"
and the proof survives verbatim in **every CM field**:

> A root of unity `ζ` acts on the Minkowski embedding of an order `O ⊂ K` as an
> isometry (`|σ(ζ)| = 1` at every archimedean place) and fixes only `0`. So the
> cyclic group `μ(K)` of order `w` acts **freely** on each positive norm-shell,
> and `w | r(D)` for all `D > 0`. The density quantum is `w/2 = #μ(K)/2`.

This is THM-416. Verified, 0 failures, through degree 6 (the cyclotomic lattices
`Z[ζ_5], Z[ζ_8], Z[ζ_12], Z[ζ_7], Z[ζ_9]`, whose Gram matrices are Ramanujan-sum
matrices `c_m((i−j) mod m)`).

## "Between triangular and CM" is the Euler totient

Once the quantum is `#μ(K)/2`, the spectrum of achievable quanta is forced by a
single classical fact: `Q(ζ_m) ⊆ K ⟺ φ(m) | [K:Q]`. So in dimension `N = 2d`
the realizable quanta are exactly `{ m/2 : m even, φ(m) | N }`, and the maximal
one is `M(N)/2` with **`M(N) = max{ m : φ(m) | N }`**.

The probe's continuum dissolves into this discrete arithmetic ladder:

| dimension N | max #roots of unity M(N) | max quantum | witness |
|---|---|---|---|
| 2  | 6  | **3** | Q(ζ₆) — triangular |
| 4  | 12 | 6  | Q(ζ₁₂) |
| 6  | 18 | 9  | Q(ζ₁₈) |
| 8  | 30 | 15 | Q(ζ₃₀) |
| 10 | 22 | 11 | Q(ζ₂₂) |
| 12 | 42 | **21** | Q(ζ₄₂) |
| 14 | **6**  | **3**  | Q(ζ₆) — *drops!* |

There is **no 2D group between triangular and CM**, not because the family is
flat but because the ladder *jumps*: 2D is capped at quantum 3, and the next
quantum (6) only exists in dimension 4. The escape opus-S699 was reaching for is
not a turn of a knob — it is a **change of rank**. This is exactly S702's
conclusion ("the real axis is the rank/degree"), now *derived* rather than
asserted, with the cap value (3) and the first escape (4, via `Q(ζ_8)`,
`Q(ζ_5)`, `Q(ζ_12)`) made explicit.

A small payoff falls out: density **5**, the one quantum S703 flagged as
reachable by *neither* the square lattice (even quanta) *nor* the triangular
(multiples of 3), is not a curiosity of the generic `disc -15` form. It is the
**rigid** quantum of `Z[ζ_5]` (`w = 10`) one dimension up. Density 4 is the rigid
quantum of `Z[ζ_8]`. The 2D "gaps" are the higher-dimensional "rungs."

## Two things worth not over-reading, and one worth reading

**The ladder is non-monotone.** `M(14) = 6`: there is no `m` with `φ(m) = 14`
(14 is a *nontotient*), so dimension 14 cannot do better than dimension 2. "More
dimensions ⇒ bigger rosette" is false; the density-quantum spectrum inherits the
ragged image of Euler's `φ`. This is the honest shape of the answer.

**The resonance.** The maximal quanta run `3, 6, 9, 15, 11, 21, 3, 30, …`, and
**21 appears at dimension 12 via `Q(ζ_42)`, `w = 42 = 2·21`**. The number 21 is
one of the project's two permanent forbidden Hamiltonian-path counts (`{7, 21}`,
THM-079), and `3` (the 2D cap = the triangle) is the project's master constant.
I record this as a resonance, **not** a theorem — `42 = 6·7`, `21 = 3·7`, and the
totient ladder and the H-spectrum are both "multiplicatively gappy integer
sequences," but I have no map between them. It is a lead (HYP-2274), not a
result.

**The reading that matters — why the LRC floor isn't capped.** The unification
thread (S702, S699g) treats the Euclidean unit-distance problem and the LRC
worry-set as the same theorem in two metrics: a geometric floor set by a
**root-of-unity rosette**, escaped arithmetically by a many-representation value.
THM-416 sharpens *why the two metrics behave differently at the floor*. On the
Euclidean side the rosette is **totient-capped by the ambient dimension** — you
can only realize `ζ_m` as an isometry if `φ(m) | N`. On the LRC side the floor
lives in `Z/(2n−1)` (THM-403), and the order-`(2n−1)` rosette is realized
**intrinsically**: the cyclic metric carries *all* `(2n−1)`-th roots of unity at
once, with no ambient-dimension constraint to satisfy. The LRC/cyclic side is the
**`degree → ∞`, all-roots-at-once limit** of the Euclidean totient ladder. Same
rosette skeleton; one is dimension-capped, the other is not.

That is the same asymmetry THM-414 found from the other direction: the cyclic
*degree-2 additive* face (`r_+(s)`, pair-sums) is matching-capped and has no
escape, while the lattice *degree-2 norm form* over `Z²` (an infinite group) has
an unbounded popular norm. Capped-vs-uncapped, twice, for the same reason: the
Euclidean/lattice object is constrained by a finite ambient structure (dimension,
or distinctness of residues), the cyclic LRC object is not.

## What this leaves for the next explorer

- **Identify the ladder.** Is `M(N) = max{m : φ(m) | N}` a named OEIS sequence
  (the "inverse totient maximum")? Its drops are exactly the nontotients;
  characterizing them characterizes where extra dimensions buy *nothing*.
- **Chase or kill the {7,21} resonance.** Is the appearance of quantum 21 at the
  first dimension whose degree admits `Q(ζ_42)` structurally tied to the
  H-spectrum's multiplicative gap at 21, or is it `7`'s ubiquity showing up
  twice independently?
- **Push the Euclidean→cyclic limit into a statement.** Can the "all-roots-at-once"
  intuition be made a precise limiting relation between the totient-capped
  Euclidean density quantum and the LRC `(2n−1)`-rosette floor — e.g. a tower of
  cyclotomic fields `Q(ζ_{2n−1})` whose density quanta converge to the cyclic
  worry-set floor as `n → ∞`?

The mathematics kept pointing past the probe. "Between triangular and CM" was
never a place; it was the totient function, asking which roots of unity your
dimension can hold.
