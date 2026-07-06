# The four dichotomies as lenses on (G), and the witness denominator divides a pairwise sum

*mac-mini-2026-07-06-S16 (HYP-4432). Owner: work the remaining obligations;
synthesize the many lenses; ponder creative conditions to leverage; look back at the
addition/multiplication, odd/even, positive/negative, rational/irrational dichotomies.
This note runs each dichotomy across the crux and extracts one proven lever that ties
them together. Verified: `lrc_witness_denominator_dichotomies…out`,
`lrc_witness_denom_divides_sumdiff…out`.*

## The three faces are one rigidity

The fleet converged, from three directions, on the same statement:

- **opus HYP-4396 — sum–product rigidity:** the AP `{1,…,12}` is the unique set that
  is *both* the additive interval `[1,12]` *and* (at `t=1/13`) the multiplicative
  group `(ℤ/13)*`.
- **kps HYP-4407 — multiplicative face:** the tight locus is the `(ℤ/13)*`-orbit of
  the nonzero 13th roots of unity; primality makes `(ℤ/p)*` cyclic ⇒ one orbit ⇒
  clean rigidity (composite ⇒ subgroups ⇒ non-AP tight sets, e.g. `{1,3,4,5,9}` at 6).
- **mac-mini HYP-4412 — three-gap / CF face:** near-tight ⇒ the witness is a `{kα}`
  orbit (`g≤3`) ⇒ `M` is a continued-fraction rung; the spectrum is CF-quantized.

Loneliness is **additive** (orbit gaps, arcs); covering/resonance-killing is
**multiplicative** (`b∣s`, dilation); continued fractions mediate; the AP is the
fixed point. Everything below is this one duality seen through the owner's four
dichotomies.

## The proven lever that unifies the lenses

> **Lemma (witness denominator).** If `M(S) = c/q` in lowest terms, then `q` divides
> `v_i + v_j` or `v_i − v_j` for some pair, or `2v_i` for some `i`.

*Proof.* `f(t)=min_i‖v_i t‖` is piecewise-linear; its maximum is attained either
where two runners are jointly minimal (`‖v_i t‖=‖v_j t‖ ⇒ (v_i∓v_j)t∈ℤ ⇒
t=k/(v_i∓v_j)`) or at a single runner's peak (`t=(2k+1)/(2v_i)`). Either way `q`
divides `v_i±v_j` or `2v_i`. ∎ (Verified across the spectrum.)

**Consequence:** `q ≤ 2·max(v_i)`. **Bounding the family's height bounds the witness
denominator**, and a bounded denominator makes (G) a *finite check* (exact `M` is
computable). This is the additive realization of kps's "bound the off-13-grid
denominators ⇒ finite check" — the density floor is finite once height is bounded,
and the height bound is exactly the difference-core piece (opus-S106) / the
single-cluster reduction (mac-mini S14).

The lever also *is* each dichotomy:

### Positive / negative (the involution)
For the AP, `q = 13` is realized by **every antipodal pair** `i + (13−i) = 13`. The
witness is where the roots of unity `e^{2πik/13}` sit symmetrically about `0`; the
clearance zone `[−c,c]` is reflection-symmetric. The pos/neg symmetry is the
`(ℤ/13)*`-orbit's `−1`-action — the complement/`T^op` involution the project has
tracked since the tournament side. A gap member needs a pair with `v_i+v_j ≡ 0 (mod q)`:
an **antipodal pair at the witness**. The AP has 6 of them (it *is* balanced); a
generic set has none — so the balanced/symmetric configs are exactly the near-tight ones.

### Odd / even (the `2` of `14 = 2·7`)
At an **even** witness `q = 2p`, an even runner `2w` has phase `(wa mod p)/p`, so the
**even runners, halved, are a config at modulus `p`** with clearance `⌈c/2⌉` (owner's
seed: `E_p={0,±2}`, `O_p={±1}` ⇒ for `c=3` the halved evens avoid `{0,±1} mod p`).
Verified (`…dichotomies…out`): a `q=12` witness halves to a clearance-1 config at
`mod 6`. A self-similar **descent `q→q/2` preserving the value** — the `2` of
`14=2·7` acting as renormalization. Covering guarantees only a multiple of `8=2³`, so
the descent has depth `≤3`; below that the modulus is odd and the parity split
dissolves into the multiplicative face. Also: `q ∣ (v_i±v_j)` means an **even `q`
forces the binding pair to be same-parity**, a parity constraint on which lifts can
create a given denominator.

### Addition / multiplication (sum–product)
`q ∣ (v_i±v_j)` is the **additive** control of the denominator; **covering** (a
multiple of every `d≤12`) is the **multiplicative** constraint on which sums exist.
A gap member (`q≥38`) needs a pair summing to a multiple of `38 = 2·19` or
`39 = 3·13`, etc. Near-AP small speeds (`1..12`, sums `≤24`) cannot reach `38`, so a
**lifted (large) runner is forced** — and the prime factors of `q` must appear in
that lifted sum. The extremal denominators are small-prime / prime-square supported:
`13` (AP), `25=5²` (`2/25`), `169=13²` (deep well) — the covering primes and their
squares. The sum–product tension is: the lift needed to make `q` large is exactly the
detuning that pushes `M` up off the ladder.

### Rational / irrational
The witness is **always rational** (`q ∣ v_i±v_j`), so the analysis is exhaustive per
`q`; there is no irrational escape. The tight AP has the *simplest* rational witness
(`t=1/13`); the loosest configs are the *badly-approximable* (golden) irrationals —
the discrepancy inversion (opus HYP-4074). The whole spectrum is the
rational-approximation ladder (Ostrowski/Stern–Brocot). "Bound the denominator" =
"bound the arithmetic complexity of the rational witness" = the finite check.

## The creative condition to leverage

Putting the lever together with the reductions already proven:

> gap member ⇒ **single cluster** (S14, decorrelation) ⇒ near-AP **bounded height**
> (difference core, the open link) ⇒ `q = witness denominator ∣ (v_i±v_j)` is
> **bounded** ⇒ (G) is a **finite check** of exact `M` on covering families of
> bounded height ⇒ closed.

So the *entire* remaining obligation collapses onto **one height bound**, and the
witness-denominator lemma converts it into a finite computation. The dichotomies
tell us *where* to look for that bound: the height must create a pair summing to a
multiple of a covering-prime-power `q`, and the parity/symmetry of that pair is
constrained. The height bound is the difference-core / renormalization contraction
(opus-S106) — and by the lemma it is *equivalent* to a denominator bound, which the
1/325 Farey seam should force (kps's Lipschitz + seam guidance, = my LRCTorusRate).

## Net

- **Proven lever (formalizable, elementary):** `q ∣ v_i±v_j` or `2v_i` — the witness
  denominator is additively controlled; bounding height bounds `q`; (G) becomes a
  finite check.
- **The four dichotomies are four readings of the same lever:** pos/neg = antipodal
  balance = roots-of-unity symmetry; odd/even = the halving descent (`14=2·7`);
  add/mult = sum-controls-`q` vs covering-constrains-sums; rat/irrat = the witness is
  rational, so bounding its denominator finitizes the density floor.
- **What remains** is the single height bound (difference-core contraction), now
  equivalent — via the lemma — to a witness-denominator bound.

## Pointers

- `lrc_witness_denominator_dichotomies_macmini_S16.py/.out` — denominators across the
  spectrum + the parity-descent verification.
- `lrc_witness_denom_divides_sumdiff_macmini_S16.py/.out` — the lemma verified (`q ∣
  v_i±v_j`), incl. the AP's antipodal `q=13`.
- opus HYP-4396 (sum–product), kps HYP-4407 (multiplicative face), mac-mini HYP-4412
  (three-gap), HYP-4402 (single-scale), HYP-4382 (prime fixed point), opus-S106
  (renormalization flow), HYP-4074 (discrepancy inversion), LRCTorusRate (Lipschitz).
