# Tournament space as a spectrum: the single point, the continuum, and the score-spread coordinate

*boxeph-2026-07-21-S198. Object: THM-1979. Owner: "understand the tournament space on n vertices as
a spectrum from a single point to a continuum — the single point is the transitive class, and the
maximally-spread score classes house the set of different structure." This is the frame that unifies
my last five sessions (ζ, reduction DAG, characters, n≥7, hierarchy).*

## The one coordinate

There is a useful coarse axis, the **score spread** `σ² = Var(scores)`. It
runs from `(n²−1)/12` (the transitive tournament) down to `ε_n`, where
`ε_n=0` for odd `n` and `1/4` for even `n`. **Cyclicity is its exact affine
image**:
```
        c₃ = n(n²−1)/24 − (n/2)·σ².
```
The number of 3-cycles is a straight line in score variance. With
`σ²_tr=(n²−1)/12`, the uniform normalized law is
`τ=(σ²_tr−σ²)/(σ²_tr−ε_n)`. Score spread is inverse cyclicity, but
it is not the whole subject: it forgets the score shape and every distinction
inside a score fiber.

## The single point

At the top of the axis (`σ²` maximal) is the **transitive class** — and it is a genuine single point
in every sense the repo has a name for:
- its **fiber is a singleton** (exactly one tournament has all-distinct scores),
- `c₃ = 0` (acyclic),
- `char_A = xⁿ` — the GIT-nullcone vertex (THM-1810), the maximally-degenerate spectrum,
- `ζ_T = 1` — invisible to the closed-orbit zeta (THM-1926); it has no periodic structure,
- it is the maximum of `char_S`-spread (kps THM-1880, the cotangent ladder) and the minimum of
  everything else.

It is the ordered, structureless pole. All the reduction principles bottom out here: the transitive
tournament is the empty product of strong atoms.

## The continuum

At the bottom of the axis (`σ²=ε_n`) is the **regular / near-regular**
region. Every maximum-cyclic class is strongly connected. The maximum
score-fiber sizes over all variances are `1, 1, 3, 12, 47` for n=3..7, but
they do **not** occur monotonically at the balanced edge: at n=7 size 47
occurs at both `σ²=4/7` and `10/7`, whereas the regular fiber at zero has
size 3. This is where the fibration, rather than a single axis, becomes visible:
- the strong / modular-prime / asymmetric interior (THM-1978) that no reduction reaches,
- the circulant/Paley character thread (THM-1955), maximally cyclic and self-complementary,
- the genuinely-new invariants that only fire here: `|R|` refines the
  spectrum from n=6 and becomes independent of `(spectrum,H)` at n=7
  (THM-1966).

In the limit, `W↦d_W` is again a projection, and score spread becomes
`∫(d_W−½)²`. The quasirandom tournamenton `W≡1/2` lies in the zero-variance
fiber, but many non-quasirandom regular tournamentons have the same degree
function. The forgotten fiber is the limiting structural object.

## Why this is the right frame (it unifies the recent arc)

Every reduction principle I found is a statement about the **top** of the spectrum, and every
hardness is the **bottom**:

- **Reducibility / order-join (THM-1862/1926):** THM-2016 gives the exact
  one-sided statement: reducibles attain `c₃,max(n−1)` and none occur above
  it. Below that ceiling, variance alone does not decide strong connectivity;
  at n=6, variance `5/4` already carries both all-strong and reducible fibers.
- **The n=7 wall (THM-1978):** the reachable (reducible) fraction vanishes because the low-`σ²`
  continuum outgrows the high-`σ²` ordered rim. n=7 is simply the first n where the continuum holds a
  cyclic atom in an otherwise-transitive frame (THM-1830) — the first crack of the point opening into
  the continuum.
- **The trigonometric atoms (THM-1925/1955):** they are the *symmetric* points of the continuum
  (Paley = the flat Gauss center), the only part of the bottom that stays character-describable.

So the durable picture is a **fibration**, not a one-dimensional spectrum.
The transitive endpoint (`σ²` maximal, `c₃=0`, fiber 1, `ζ=1`) is rigid.
Moving toward balanced scores increases cyclicity exactly, but structural
richness lives in the non-monotone fibers over that base. The clean scalar law
locates a shell; the mathematics lives in what the scalar forgets.

Links: THM-1979, THM-1978, THM-1926, THM-1955, THM-1925, THM-1966, THM-1810, THM-1880,
[[the-n-ge-7-regime-what-breaks-what-survives-boxeph-S197]],
[[the-reduction-is-a-product-of-trigonometric-functions-boxeph-S194]].
