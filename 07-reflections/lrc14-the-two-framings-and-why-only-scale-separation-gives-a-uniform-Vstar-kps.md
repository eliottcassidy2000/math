# The two framings of the apex ruler, and why only scale-separation gives a uniform V*

*kind-pasteur-2026-06-22-S34, THREAD 1. A reflection on rigorizing the NODE-1 three-gap
lemma (THM-565).*

## The coincidence that wasn't

The S33 closure (HYP-2853) stated the apex-ruler three-gap lemma and "verified the boundary
core {1,…,12,V}" with `meas(G)=319/560`, `arcCount=546`, `V*=958`. Three numbers. When I
recomputed them in the most literal framing — the small speeds {1,…,12} restricted by their
own safe set `G_P`, the apex as the lone fast clock — I got `meas(G)=6617/194040 ≈ 0.034`,
not `0.570`. A factor of 17 off. The arcCount was `12`, not `546`.

A discrepancy that large is never an arithmetic slip; it is two *different objects* wearing
the same name. Chasing it down turned out to be the whole content of the rigorization.

## Two ways to point the ruler

A covering tuple `S = P ∪ L`, apex `V = max S`. There are two ways to assign "who is the
slow clock and who are the teeth":

- **Framing B (all-as-teeth):** every non-apex speed becomes an offset tooth
  `e = V − u`, and `G = {x : maxgap{frac(e·x)} > 1/7}`. This gives `meas ≈ 0.57` for the
  boundary core — the S33 number. *But* the small-speed offsets `e = V − p ≈ V` are **huge**,
  so the teeth `frac(e·x)` oscillate at frequency `~V`, and the good set acquires
  `arcCount ~ V^{0.6}` arcs. As `V → ∞`, `V* = arcCount/meas → ∞`. **The threshold runs away.**

- **Framing A (scale-separated):** the small speeds `p ≤ 13` are handled by the *slow*
  condition `x ∈ G_P` (their phases barely move over one ruler period); only the genuinely
  clustered speeds `u ∈ L`, whose offsets `e = V − u` are **bounded**, become teeth. Then
  `G = G_P ∩ {maxgap{frac(e·x)} > 1/7}` is **literally independent of `V`** — `arcCount` and
  `meas(G)` are constants, and `V* = m/c` is a genuine fixed number (≈ 234 worst-case).

Both framings are *correct* descriptions of the same loneliness event. But only Framing A has
a uniform `V*`. The S33 lemma was stated in the language of Framing B (where `meas = 0.57`
looks comfortably large) but the bound `#good ≥ V·meas − arcCount` is only non-vacuous in
Framing A (where `meas = 0.034` looks alarmingly small but `arcCount = 12` is bounded).

The large `meas` of Framing B is a mirage: it is paid for, dollar for dollar, by an arcCount
that grows with `V`. The honest invariant is `meas/arcCount`, and that is maximized by
*refusing to treat a slow clock as a tooth*.

## Why scale-separation is forced, not chosen

This is the part that points beyond the particular lemma. The lonely-runner problem at a
covering tuple has a built-in **separation of time scales**: the apex `V` is the fast hand,
sweeping its phase `φ = frac(Vτ)` across `(1/14, 13/14)` within a single ruler period; the
slow hand is `x = j/V`, the period index. A speed `u` enters the analysis through *which hand
it tracks*:

- `u` close to `V` (the cluster) ⟹ `V − u` small ⟹ `u·τ` tracks the **fast** hand minus a
  small slow drift ⟹ it is a *tooth* in the fast phase.
- `u` small (the base) ⟹ `u·τ` tracks the **slow** hand essentially alone ⟹ it is a
  *constraint on `x`*, i.e. a `G_P` condition.

The drift between the two — `p·φ/V = O(maxP/V)` — is exactly the error term that the shrunk
safe set `G_P^δ`, `δ = maxP/(2V)`, absorbs. As `V → ∞` the two scales decouple perfectly and
`δ → 0`. **Framing A is not a clever choice; it is the eigenbasis of the slow-fast splitting.**
Framing B mixes the scales — it forces a slow clock (small `p`) to masquerade as a fast tooth
(offset `V − p`) — and pays with high-frequency oscillation it didn't need.

This is the same lesson as the triangle foundation's **Mode A vs Mode B** (CLAUDE.md): there
is a fast time scale (hypotenuse removal, `n → n−1`) and a slow one (both-legs removal,
`n → n−2`), and invariants are clean only when computed in the scale they live on. Here the
apex ruler *is* the fast scale and `G_P` *is* the slow scale; the three-gap lemma is clean iff
you keep them apart.

## The factorization that fell out

Once Framing A is fixed, the worst-case threshold has a clean closed shape:

> **V\* = (LRCArcComplexity arcCount) / (OPEN-Q-108 floor).**

The numerator is the *bounded* arc count of the cluster teeth (≤ 7·sumE, sorry-free in Lean).
The denominator is the *measure floor* `c = meas(G)`, which is minimized by the small parts
`P` that minimize `meas(G_P)` — exactly the quantity OPEN-Q-108 is about (`m_P = 14249/252252`).
So the finite-V gate (NODE 1) and the measure floor (NODE 3) are not two problems; they are
numerator and denominator of one ratio. The owner's intuition — "a single clean lemma knocks
over all three nodes" — is literally the statement that `V*` is `arcCount` over the floor: the
three-gap lemma *is* the division.

And the division is elementary. Every piece of the lemma except the floor `c > 0` is a finite,
exact, verified combinatorial fact (piecewise-linearity of maxgap; lattice-point counting;
the `1/7 ⟹ 1/14` reach midpoint). The hardness was never in turning `ρ*` into a finite-V
witness — that step is a ruler and a floor function. The hardness is, and remains, the floor.

## The moral

When two computations of "the same thing" disagree by a factor of 17, do not reconcile them by
splitting the difference. One of them is measuring an object whose defining structure (here:
which speeds are fast) was left implicit. The disagreement is the structure announcing itself.
The right framing is the one in which the invariant you actually need (`meas/arcCount`, the
non-vacuous bound) is *stable* — and stability under the natural deformation (`V → ∞`) is the
tell. Scale-separation wins because it is the only framing in which the slow-fast change of
variables is an isometry rather than a shear.
