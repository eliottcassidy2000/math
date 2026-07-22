# Where GMC(2) reaches LRC(14) — the CT-functional — and where it stops — the volume ceiling

*boxeph-2026-07-21-S211. Owner: ponder creatively, with multiple pulls into past/incoming threads, how the
GMC(2) proof can be leveraged in combination with other things toward an LRC(14) proof. Builds on codex
THM-2047 (phase-height carrier; corrects my S209 — MISTAKE-223), mac-mini S157 (the obstruction), boxeph
S205 (Frobenius is rank-1), S206 (Wall A = joint-order AP-extremality), S208 (Vandermonde localization),
S210 (minimax = saddle), THM-1840 (shared seed), THM-1820 (relation-lattice pairing), THM-730 (AP maximizes
triples). Seven pillars verified in `04-computation/ct_functional_bridge_gmc_to_lrc_boxeph_S211.py`.*

## The honest starting point

A *direct* GMC(2)⇒LRC(14) implication is dead: the functionals differ (mac-mini S157 — factorial-monotone
vs. sinc-oscillating), and Frobenius amplification is rank-1 while LRC(14) is rank-12 (S205). Only the rank-1
**seed** transfers. So the question is what GMC(2)'s *techniques*, in **combination** with the LRC-native
machinery, actually buy — and, just as importantly, where they hit a ceiling. Both answers turned out to be
sharp.

## Where it reaches: the constant-term functional is the shared object

GMC's Gaussian moment factors through the **polar bridge `E = L ∘ CT`** (THM-1645, proved/Wick-verified):
an *angular* constant-term part `CT` (the DvdK/Duistermaat–van der Kallen integral, `CT_u[Λ_s^m]`) and a
*radial* Laplace part `L` (`L(s^k)=k!`). Two facts from that theorem matter here. First, the angular half is
**DvdK-closed** (THM-1630 = Duistermaat–van der Kallen 1998), so *the entire open part of GMC(2) is the
radial Laplace descent* — "GMC(2) is obstructed by **Laplace determinacy, not tori**" (THM-1645). Second,
and decisively (deathstar-S77): **LRC's covering moments live in GMC's already-*proved* CT half; GMC's open
radial half is *perpendicular* to circle-covering.** So the transferable part of GMC(2) is exactly its
*proved* part, and the angular `CT` is precisely what LRC needs, because (verified, P1):

> **additive `m`-energy of a speed set `S` = `CT[Pᵐ P̄ᵐ]` = `‖P‖_{2m}^{2m}`**, where `P(z)=Σ_{v∈S}z^v`.

So the additive energy — the object Wall A is about — **is a CT-moment**, and `CT` is literally the angular
half of the GMC functional. (opus-S153 already had the general form: the order-`j` energy `E_j` is a
support-`2j` relation-lattice / constant-term count — I re-derive the identity cleanly and read it through
the polar bridge.) This makes the GMC↔LRC kinship concrete rather than analogical:

- **The shared seed (THM-1840) is a CT statement** (verified P2): `[u⁰](au^p+bu^{−q})^m` has a single
  uncancellable term at the coprime return `m₀=(p+q)/gcd` — a constant-term extraction, functional-agnostic
  (`E`, Laplacian, sinc alike).
- **Wall A is a CT-moment extremality.** S206 recast LRC(14)'s crux as "the AP jointly maximizes additive
  energy at *all* orders." Verified (P3), extending THM-730 (the `m=2`/triple case) to `m=3,4`: over `k`-
  subsets, **every** additive-`m`-energy maximizer is an AP, for `m=2,3,4`. A disproof would have to *peel*
  higher-order energy from triple energy; none does.
- **The AP is the cyclotomic extremal** (verified P4): `P_AP=(z^k−1)/(z−1)=∏_{d|k,d>1}Φ_d`. Two distinct
  cyclotomic structures meet at `n=14` and I should keep them separate (per the mining pull): the *argmax
  value* `t*=14/183=n/Φ₆(n)` is **Eisenstein** (`Φ₆=x²−x+1`, the `ζ₆`-line / `ℤ[ζ₆]`), while the *speed-set
  apex hardness* localizes to **`Φ₇`** — `x¹⁴−1=Φ₁Φ₂Φ₇Φ₁₄`, and the-recursion-modes note pins LRC(14)'s
  hardness on the single apex `Φ₇` factor (`14=2·7`). The AP's cyclotomic factorization is the common
  scaffold; the value is `Φ₆`, the apex is `Φ₇`.

So GMC(2)'s proven machinery — the seed, and the Vandermonde/localization technique of S208 — lives on the
CT/angular side, and that side is genuinely shared. This is the reach.

## Where it stops: the volume ceiling (codex THM-2047, correcting my S209)

Pulling the incoming thread was decisive. **codex's THM-2047** proves the exact LRC carrier is the *signed
phase-height cell complex* `vt±δ∈ℤ` (not the standard toric arrangement I proposed in S209 — MISTAKE-223,
which I accept and adopt), and it establishes the crucial fact that guillotines the whole CT/volume program:

> **Safe volume is only a strict-exit certificate.** At the tight threshold the good set `G_{M}(S)` can be
> **nonempty but measure-zero**; Haar volume reads `0` and *misses the tight witness*. The correct iff
> criterion is `χ(G_δ(S)) > 0` (the number of components), which *sees* the isolated tight phase.

The CT-moment / additive energy / `|G_δ|` **is a volume functional** — precisely the quantity that goes
blind here. Verified head-on (P7): for `S={1,2,3}`, as `δ↑¼` the volume `|G_δ|→0`, hitting exactly `0` at
the threshold `δ=¼=M(S)` — yet the runner **is** lonely there, at the two isolated phases `t=¼,¾`, so
`χ(G_{¼})=2>0`. Volume says "covered"; topology says "lonely." **LRC(14)'s crux — Wall A, the tight AP —
lives exactly in this measure-zero boundary that the CT/volume functional cannot resolve.** codex also
showed the S209 block-localization has the wrong local dimension after restriction to the LRC orbit; the
exact first-order object is his boundary-layer law `|G_{M−ε}| = ε·Σ_{maximizers}(1/s₊ + 1/(−s₋))`, whose
linear vanishing in `ε` *is* the volume→0 collapse P7 measures.

This is the *measure-theoretic* face of a fact the repo already knew *analytically*: additive energy is
**necessary but not sufficient** for loneliness (opus-S181) — the translated AP `{1,3,…,25}` has maximal
energy yet is loose, because loneliness `L` is **dilation-invariant but not translation-invariant** while
energy is translation-invariant; and THM-730's "open mile" is exactly that the Schur/energy deficit does
not track `L` (the witness `{1,…,11,13,84}` has the largest `E₃`-deficit yet the smallest `L`). Volume-
blindness at the saddle (P7) and energy-insufficiency (opus-S181) are one obstruction seen two ways: the
CT-moment is the wrong *final* invariant, even though it is the right *shared* one. And there is a reason
the correct invariant is topological and finite: LRC's danger count `X∈{0,…,13}` has a **bounded alphabet**,
so detection depth is finite (klein-S389, Bonferroni terminates at depth 5) — which is exactly why a finite
topological count, `χ(G_δ)=#components`, is the natural closing invariant where the unbounded-depth volume
expansion is not.

## The bridge across the ceiling: the minimax is a saddle

Why does volume degenerate exactly at the crux, and what is the right object there? Because
**`M(S)=max_t min_v ‖vt‖` is a minimax — a saddle-point value** (from S210; this optimization fact is
untouched by codex's MISTAKE-224, which only retracts S210's stronger "torus needs saddles" equivalence —
I keep only the minimax and the Morse-index bookkeeping). At the tight witness the piecewise-linear `f_S`
attains its max at a **critical/saddle phase** where an active rising wall (`+v_i`) meets an active falling
wall (`−v_j`) — codex's pair-sum law `q | v_i+v_j` (THM-2047 §2). That is a *measure-zero critical point*,
so volume cannot see it and `χ` can. The **volume→χ gap is a saddle degeneration**, and the correct
invariant there is `χ(G_δ)=#components` (codex THM-2047 §4).

## The picture: two complementary halves, and the residual named precisely

```
   GMC(2) / CT-moment  ──►  VOLUME of the good set          (strict-exit certificate; measure-positive)
        (angular half E=L∘CT; seed THM-1840; AP=cyclotomic extremal, verified)
                                    │  minimax = SADDLE (S210)  ← the bridge
   codex THM-2047 / phase-height ──►  χ(G_δ)  (the tight, measure-zero boundary; the iff criterion)
```

GMC(2)'s CT leverage and codex's phase-height `χ` are **complementary**: volume handles the strict
(measure-positive) loneliness, `χ` handles the tight (measure-zero) boundary, and the minimax/saddle
structure is the hinge between them. Neither closes LRC(14) alone; together they **partition** it, and the
residual is pinned exactly: *detect the measure-zero tight witness — force `χ(G_{1/14}(C))>0` for every
13-speed covering core* — which is Wall A (HYP-7310), now with a precise topological form (`χ`, not volume)
and a precise reason the CT/volume half cannot reach it (measure-blindness at the saddle).

## What the session actually establishes

- A concrete, verified GMC↔LRC bridge: **additive energy = CT-moment = the angular half of the GMC
  functional**, with the seed as a shared CT statement and Wall A as a CT-moment extremality (AP maximizes
  at orders 2–4, cyclotomic/Eisenstein).
- A precise **ceiling**: that leverage is a *volume* functional and (codex THM-2047) volume is measure-blind
  to the tight witness — verified on `(1,2,3)` (`vol→0`, `χ=2`).
- The **bridge across it**: `M(S)` is a minimax/saddle (S210), so the residual is a `χ`/saddle-degeneration
  problem — codex's phase-height complex is the right carrier, not my S209 toric complement.

Honest scope: no new implication toward a proof. The value is the **decomposition** — GMC(2) buys the
volume/strict half and the shared seed; the LRC(14) crux is the complementary `χ`/tight-boundary half — and
the explicit identification of *why* (measure-blindness at the minimax saddle) and *what carries the rest*
(codex THM-2047's phase-height `χ`). It folds S205/S206/S208/S210 and codex's THM-2047 into one map, and
adopts codex's correction to my S209.

Links: HYP-8840, THM-2047 (codex), THM-1645 (polar bridge), THM-1630 (DvdK), THM-1840, THM-1820, THM-730,
HYP-7310 (Wall A); deathstar-S77 (CT is GMC's proved half, radial half ⊥ LRC), opus-S153 (E_j = support-2j
CT counts), opus-S181 (energy necessary-not-sufficient), klein-S389 (finite detection depth),
[[what-an-lrc14-disproof-must-be-and-why-fibonacci-is-the-foil-boxeph-S206]],
[[antisymmetry-is-the-hinge-tori-odd-functions-saddles-and-tournaments-boxeph-S210]].
