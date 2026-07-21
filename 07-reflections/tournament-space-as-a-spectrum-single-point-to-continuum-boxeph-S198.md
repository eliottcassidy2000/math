# Tournament space as a spectrum: the single point, the continuum, and the score-spread coordinate

*boxeph-2026-07-21-S198. Object: THM-1979. Owner: "understand the tournament space on n vertices as
a spectrum from a single point to a continuum — the single point is the transitive class, and the
maximally-spread score classes house the set of different structure." This is the frame that unifies
my last five sessions (ζ, reduction DAG, characters, n≥7, hierarchy).*

## The one coordinate

There is a single axis that organizes all of tournament space, and it is the **score spread**
`σ² = Var(scores)`. It runs from `(n²−1)/12` (the transitive tournament, scores 0..n−1) down to `0`
(the regular tournament, all scores (n−1)/2). Everything else is a monotone function of it, because
**cyclicity is its exact affine image**:
```
        c₃ = n(n²−1)/24 − (n/2)·σ².
```
Read that twice: the number of 3-cycles — the amount of *intransitivity*, the whole subject — is a
straight line in the score variance. Spread-out scores ⇒ few cycles ⇒ near-transitive; flat scores ⇒
maximum cycles ⇒ regular. Score spread *is* (inverse) cyclicity. One coordinate, two names.

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

At the bottom of the axis (`σ² → 0`) is the **regular / quasirandom** region, and here the fibers
**swell**. The number of iso classes sharing a score sequence — the fiber — grows `1, 1, 3, 12, 47`
for n=3..7 (and without bound after), and these classes are **entirely strongly connected and mostly
modular-prime**. This is where "the different structure" lives:
- the strong / modular-prime / asymmetric interior (THM-1978) that no reduction reaches,
- the circulant/Paley character thread (THM-1955), maximally cyclic and self-complementary,
- the genuinely-new invariants that only fire here: mac-mini's `|R|` becomes independent of the
  spectrum and `H` *exactly at n=7* (THM-1966), because that is the first n where the low-`σ²`
  continuum is rich enough to hold two structures the coarse invariants cannot tell apart.

The limit object is the **quasirandom tournamenton** `W ≡ 1/2`: a positive-entropy continuum whose
finite shadows are these exponentially-swelling fibers. The transitive limit `W(x,y)=1_{x>y}` is the
single ordered point; the degree function `d(x)` interpolates from `x` to `1/2`, and score spread is
the functional `∫(d−½)²`.

## Why this is the right frame (it unifies the recent arc)

Every reduction principle I found is a statement about the **top** of the spectrum, and every
hardness is the **bottom**:

- **Reducibility / order-join (THM-1862/1926):** reducible tournaments cluster at high `σ²` (the
  data: high-spread fibers have strong-count 0 — they are singleton reducible chains). The reduction
  DAG "comes from smaller" precisely because near-transitive tournaments *are* stacks of transitive
  singletons on a small core.
- **The n=7 wall (THM-1978):** the reachable (reducible) fraction vanishes because the low-`σ²`
  continuum outgrows the high-`σ²` ordered rim. n=7 is simply the first n where the continuum holds a
  cyclic atom in an otherwise-transitive frame (THM-1830) — the first crack of the point opening into
  the continuum.
- **The trigonometric atoms (THM-1925/1955):** they are the *symmetric* points of the continuum
  (Paley = the flat Gauss center), the only part of the bottom that stays character-describable.

So the picture the owner named is the correct global one: **tournament space is a spectrum whose
single point (transitive, `σ²` max, `c₃`=0, fiber 1, `ζ`=1) opens, as score spread falls, into a
continuum of strongly-connected structure (regular, `σ²`=0, `c₃` max, fiber →∞) — and everything the
theory calls "structure," "hardness," or "new invariant" is a coordinate on that same axis.** The
clean theorems live at the point; the mathematics lives in the continuum.

Links: THM-1979, THM-1978, THM-1926, THM-1955, THM-1925, THM-1966, THM-1810, THM-1880,
[[the-n-ge-7-regime-what-breaks-what-survives-boxeph-S197]],
[[the-reduction-is-a-product-of-trigonometric-functions-boxeph-S194]].
