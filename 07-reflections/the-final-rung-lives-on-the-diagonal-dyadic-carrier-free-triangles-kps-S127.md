# The final rung lives on the diagonal — a dyadic carrier over free triangles

*kind-pasteur-2026-07-09-S127. Owner: "work the endgame, think quantitative stability and how it connects
to any possible past concepts in the repo." This note is the synthesis. It reframes the one remaining
LRC(14) obligation — the quantitative Freiman-stability rung — through the concepts the corpus already
holds, and lands on a single sentence: the obstruction is 2-adic (the diagonal of E₃), and the odd-cyclic
(triangle / Schur) content is free. That sentence unifies six independent threads.*

---

## Where quantitative stability actually stands

The endgame reduces (grand assembly + my `lrc14_from_B5` / `from_liveness`) to one obligation: a covering,
compressed, distinct, difference-primitive, no-detune, **not-near-AP** family of 13 speeds is lonely. Via
THM-681 this is the exact-load dichotomy: `W₀ ≤ 0.08` ⟹ live by the classical `B ≥ 2n−3` bound; `W₀ >
0.08` ⟹ the Freiman ladder ⟹ near-AP (excluded) or GAP (dispatched). The **final rung** is the middle.

THM-682 has since collapsed that rung almost entirely. Every core family has `B ≥ 32` (the restricted-sumset
ladder `diam ≤ B−11`, exhaustive through `B ≤ 31`, with the `B = 3k−7 = 32` twin-AP escape dispatch-owned).
And its part (d) is the pivot of this whole note:

> **THM-682(d): the only support-2 global exact relations are doublings `v_b = 2 v_a`.**  Schur triples
> are nearly weightless (line weight `0.0027`); `W₀ > 0.08` needs `≥ 3` doublings, hence an even-rich,
> 2-adically coherent family.

So the remaining rung is not an interval of exotic sets — it is the **doubling-rich corner**, and the
dispatch that owns it is 2-adic (LEM-019 descent, LEM-021 depth-4, the `g = 2` detune).

## The redirect I had half-right: E₃ is the honest invariant (opus-S182)

opus-S182 is the most important thing the sweep surfaced. Loneliness is **dilation-invariant but not
translation-invariant**, and among additive functionals only `E₃` (the Schur count `#{a+b=c}`) shares that
exact symmetry group — `E₂` (additive energy `a+b=c+d`) and the classical burden `|A +̂ A|` are
translation-invariant and *cannot* separate the tight AP from its loose translate (the sum-free `{1,3,…,25}`
has maximal `E₂` but `E₃ = 0` and is loose). So the honest stability invariant is `E₃`, the **anchored**
quantity — exactly the axis my peel ladder lives on, and *not* the burden axis the classical Freiman
`3k−4` machinery would use. The burden ladder (opus/mac-mini) is a correct but wrong-symmetry proxy that
works because the residual class is small enough to exhaust; the symmetry-honest object is `E₃`.

## The bridge to the project's origin: a Schur triple is a triangle (THM-446)

Here is the connection that makes "everything is the triangle" literal for the endgame. THM-446 records
the **additive-relation ↔ cycle-length ladder**: on the additive Cayley graph, a Schur 3-term relation
`a+b=c` is a **triangle — the smallest odd cycle**; a Sidon 4-term `a+b=c+d` is a 4-cycle; additive
energy is the `C₄`-count; and the **power-of-two relation `d ↦ 2d` is the dyadic hard core** (Erdős-64).

So `E₃` — the invariant opus-S182 says is the honest one — is a **triangle count on the speed set**, and
the project's founding subject is the **Odd-Cycle Collection Formula**: counting odd cycles in tournaments.
The triangle is the OCF's atom (`Φ₆(n) = n²−n+1 = 2T_{n−1}+1`, "the triangle doubled and shifted"). The
LRC endgame did not wander away from tournaments; it walked back into the smallest odd cycle.

## My contribution this session: the diagonal split (Lean, sorry-free)

Put THM-682(d) and THM-446 together in the `E₃` ledger. `E₃ S = #{(a,b) ∈ S² : a+b ∈ S}` splits by its
**diagonal**:

  `E₃ S = doublingCount S + (off-diagonal Schur count)`,   `doublingCount S = #{a ∈ S : 2a ∈ S}`

(`LRCSchurPeel.schurCount_eq_doubling_add_offDiag`, sorry-free; verified `E₃ = D₂ + 2T` exhaustively). The
two pieces are exactly the two rungs of THM-446's ladder:

| `E₃` piece | additive relation | cycle | endgame role |
|---|---|---|---|
| **diagonal** `(a,a)`, `2a∈S` (`doublingCount`) | doubling `d ↦ 2d` | dyadic hard core | **the `W₀`-carrier** (THM-682(d)) |
| **off-diagonal** `(a,b)`, `a+b∈S` | Schur triple `a+b=c` | **triangle** (odd cycle) | nearly weightless (`0.0027`) |

The computation confirms the diagonal *is* the 2-adic axis: doubling-rich ⟺ even-rich (mean even speeds
climbs `5.75 → 8.63` as `D₂: 0 → 6`), and the doubling graph is a **forest of geometric-ratio-2 chains**
(`D₂ = |S| − #chains`). A long doubling chain `a, 2a, 4a, …` is a rank-1 GAP *in the multiplicative group* —
which is exactly what klein's multiplicative character frame (HYP-5835) sees at prime rulers.

## The one sentence: the obstruction is 2-adic, the triangles are free

Six threads, one statement.

- **THM-682(d):** only doublings carry `W₀`.
- **My split:** doublings = the diagonal of `E₃`; Schur triples = the off-diagonal.
- **THM-446:** off-diagonal = triangles = odd cycles; diagonal = the dyadic hard core.
- **LEM-019/020/021 (mac-mini):** the doubling corner is dispatched 2-adically; the first live parity
  layer is 2-adic depth 4; all-odd sets are `½`-lonely.
- **klein-S222:** the single irreducible obstruction is a *signed bit* — a 2-adic parity phenomenon, not a
  triangle phenomenon (the absolute triangle/Schur tail is free to diverge because it cancels).
- **opus-S182:** the honest invariant is `E₃`, and its symmetry is dilation — under which the diagonal
  (doublings, `2a`) is rigid and the off-diagonal (translates of the sum relation) is not.

**The LRC(14) final-rung obstruction lives on the diagonal of `E₃` — it is 2-adic / dyadic — and the
odd-cyclic (triangle, Schur-triple) content is free.** That is why every road "merges back into the one
calibrated wall" at the 2-adic dispatches (LEM-021's `16 < 28 < 32` knife-edge), why the odd-cycle
machinery never sees a lone obstruction, and why the absolute triangle mass is allowed to diverge. The
triangle — the object the whole project is built from — turns out to be the *free* part of the last rung;
the load is carried by its dyadic shadow.

## Two proof-method echoes (the tournament move is already here)

1. **Rank bounds the local count.** My peel ladder's Rung B, `repCount S (max) ≤ |S| − 1`, is the same
   move as the tournament score / cut-space hierarchy (the vertical leg of the staircase `δ_{n−2}`): a
   rank/level bound caps a local incidence count. mac-mini's LEM-018 ω-map (hole-accounting injection) and
   Rédei-parity counting (LEM-020) are the same signature move — existence by exact counting.
2. **The doubling tower has a tournament twin.** THM-455's skew-Sylvester doubling tower (transitive
   subtournaments, the `+1` law, the Mersenne 7–31 window) is the tournament-side image of the speed set's
   doubling forest — the same 2-adic climb, Mode B of the triangle foundation (Cayley–Dickson descent).

## What this buys the endgame, concretely

- The `E₃`-axis contribution now hands the carrier to the right lane: **the additive content that matters
  for the final rung is `doublingCount` (the diagonal), which is the 2-adic corner mac-mini/monad own.**
  The off-diagonal (my peel ladder's bulk, opus's burden) is the *free* part — good to have formalized as
  the extremal anchor, but not where the last bit hides.
- It says the remaining work is **not** a general Freiman/BSG/`3k−4` import (Tier-4: none exists in canon,
  and opus-S182 says it is the wrong symmetry anyway). It is the finite 2-adic corner: the `B = 33`
  boundary (`native_decide`), the flagged LEM-016 sliver, and the doubling-chain dispatch — all dyadic.
- Honest caveat (MISTAKE-115): do **not** re-frame the corner as a "bounded-defect" Freiman exclusion —
  the governing parameter is the order/2-adic depth, not a defect count. The diagonal split respects this:
  it counts doublings (2-adic events), not deficits.

*Files: `LRCSchurPeel.schurCount_eq_doubling_add_offDiag` (sorry-free), `lrc14_e3_diagonal_split_kps_S127.py`
/`.out`. Builds on THM-681/682 (monad), THM-446 (Sidon/cycle ladder), opus-S182 (E₃ symmetry-match),
LEM-019/021 (mac-mini 2-adic), klein-S222 (the signed bit). See
[[the-two-axes-share-a-threshold-e3-peel-ladder-kps-S126]],
[[the-interval-is-the-shared-extremizer-schur-triples-and-lrc-loneliness-kps-S113]],
[[the-final-rung-is-signed-not-absolute-kps-S124]].*
