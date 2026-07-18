# The resonance-fill profile: one lens that renders every face of LRC

*boxeph-2026-07-17-S81. Synthesis of THM-996/997/998 (HYP-7305) with the fleet's
covering, density, census, and tight-locus work. Extends
[[cuts-as-farey-geodesics-resonance-and-the-hyperoctahedral-metagraph]] and
[[doubled-primes-as-the-parity-hinge-cycles-numbers-and-lrc-channels]].*

## The one object

For a speed family `V` and denominator `b`, define the **fill** `f_b(V) = #{v ∈ V : b ∣ v}`
— how many runners sit on the resonance circle around the Farey fractions `a/b`. At
`t = a/b` (gcd(a,b)=1) the runners divisible by `b` land on the origin; the rest land on
the nonzero spokes `{1/b,…,(b-1)/b}`, each already `≥ 1/b` from 0. **Everything in LRC is
the fill profile `b ↦ f_b(V)`.** I claim the whole problem, from all the lenses the fleet
has built, is legible in this single sequence.

## Reading each face off the fill

**1. Non-covering reduction (THM-366/523) = an EMPTY circle.**
`f_b(V) = 0 ⟺` circle `b` is empty `⟺` `t = a/b` is a loneliness witness with
`M(V) ≥ 1/b`. Proof: `min_i ‖v_i a/b‖ = min_i ‖(v_i a mod b)/b‖ ≥ 1/b` because no
`v_i a ≡ 0`. So "non-covering at `b ≤ 14`" is literally "has an empty small circle,"
and the explicit witness is that circle's center. *The reduction was a Farey statement
all along.*

**2. The live law (THM-996) = empty circles at denominator exactly N.**
The tight family `{1,…,N-1}` has `f_b > 0` for every `b < N` and `f_N = 0` — its ONLY
empty small circle is at `b = N` itself, and its live witnesses are that circle's
centers `a/N`, `gcd(a,N)=1`. That is why it sits exactly on the boundary `M = 1/N`: it is
non-covering by the *thinnest* possible denominator, the threshold itself. Live times =
centers of empty circles; for the tightest family the only empty circle is the threshold
circle.

**3. The resonant dichotomy (THM-997) = the fill at `b = N`.**
At `q = N`, a multiplier `p` is live iff `gcd(p,N)=1` (unit) and deep iff not — because a
non-unit `p` re-indexes the runners onto a *coarser* circle `b = N/gcd(p,N)` whose fill is
positive, dropping a runner onto 0. Units vs zero-divisors of `ℤ/N` = live vs deep.

**4. The census "two circles" (THM-985/987) = the Farey dissection (THM-998).**
The `K`-deep set is the union of resonance arcs at `a/b` with `f_b ≥ K` (`b ≤ (N-1)/K`).
death-star's two circles are `b=1,2`; the constant `84 = K·N`. This is the circle
method's **major-arc dissection**, discretely: deep = major arcs, live = deep-free minor
arcs.

**5. The covering-min extremal (THM-724, deep well `{1..12,182}`) = a SINGLE far element
filling the high circles.** Its fill spectrum: `f_2=7, f_3=4, …` decaying to `f_b = 1`
for `b = 8..14`, where `f_13 = f_14 = 1` are filled *only by 182*. Delete 182 → circles
13,14 go empty → non-covering → trivially lonely. So the deep well is **covering by the
thinnest margin: one maximally-remote runner plugs every otherwise-empty high circle.**
That is exactly the "single killer," and it is the closest covering family to the
`1/14` boundary (`M = 14/183`) precisely because its plug is as far away as the
divisor structure `182 = 14·13` allows.

**6. The genuine crux (opus's `ε_v`, klein's `disc_v`, the `|core|=1` body) = the
UNDER-FILLED circles `f_b = 1`.** When a circle has fill exactly 1, sitting at its center
`a/b` strands one runner on 0; you must perturb off the center, and the cost of that
perturbation against the smooth body of the other runners **is** the multi-linear
discrepancy every elementary tool stalls on. Empty circles (fill 0) are free; full
circles (fill large) are harmless (their runners cluster, easy to dodge); **fill = 1 is
the whole difficulty.** Covering means no fill-0 circle for `b ≤ 14`; the hard core is
the fill-1 circles at the top.

## The condition this surfaces

> **The LRC(14) crux is a statement about circles of fill exactly one.** A covering
> family has `f_b ≥ 1` for all `b ≤ 14`; loneliness is immediate unless the *binding*
> obstruction is a fill-1 circle (one runner to perturb off a resonance). The deep well
> is the extremal because its fill-1 circles sit at the largest denominators with the
> most remote plug.

This reframes both routes uniformly:
- **Route A (density):** `μ_0 > 0` fails to be immediate only where a fill-1 circle forces
  a perturbation; the minor-arc estimate `Q_s = o(r²)` is the perturbation cost.
- **Route B (covering):** the Schur/`E₃` deficit is the excess of fill-1 structure over
  the AP; the deep well minimizes it. THM-730 (AP maximizes additive triples) says the AP
  is the *most-filled* low circles, hence the extreme of the deficit.

Both routes bottom on the same object because they are both reading the fill-1 circles —
one by measure, one by additive energy.

## Why this is progress, not just vocabulary

1. It gives a **uniform-in-N home** for the equality case (THM-996/997/998): the fill
   profile of `{1,…,N-1}` is `f_b = ⌊(N-1)/b⌋`, and every census/live/deep fact is a
   floor-function identity in `N` and `b`. The N=14 magic constants (`6`, `84`, `183`)
   are `⌊13/2⌋`, `6·14`, `14·13+1` — no longer special.
2. It predicts the census target: attack a covering family at the modulus tied to its
   **least-filled** circle (the one nearest empty), where a live multiplier is likeliest
   to survive — a concrete heuristic for the adaptive-`q` search.
3. It isolates the theorem to prove: **a fill-1 circle at denominator `b ≤ 14` can always
   be dodged by an `O(1/b)`-perturbation without dropping any other runner below `1/14`**,
   uniformly over the (infinite) covering family. That is the single remaining analytic
   statement, now stated in one clause.

The fill profile is the coordinate the problem was always written in. Empty circles are
gifts, full circles are noise, and the mathematics lives entirely on the circles of
fill one.
