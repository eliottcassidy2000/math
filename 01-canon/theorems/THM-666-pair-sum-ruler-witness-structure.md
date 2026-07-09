# THM-666 — The pair-sum ruler theorem: the witness always lives at t = p/(v_i+v_j), and LRC(14) is a band statement on pair-sum moduli

**Status:** PROVED (parts 1–4, elementary, proofs below) + VERIFIED (300 random sets: pair-sum max
= full-breakpoint max, 0 mismatches) + the realization corollary's uniform-slack refinement is
HYP-5720 (OPEN, exact tests this session).
**Source:** mac-mini-2026-07-09-S65.
**Precedents (honest):** THM-420 (opus-S700) observed *empirically* that the shell-free residual's
witnesses sit at pair-sum resolutions; mac-mini-S64 and klein-S206 *used* the local-max enumeration
as an exact-computation device without proof. This file supplies the proof and the structural
consequences (witness-ruler boundedness, the event-realization criterion, the Schur-triple kill
rule).

**Convention:** `S = {v_1 < … < v_k} ⊂ Z_{>0}` primitive, `m(t) = min_i ‖v_i t‖`,
`M(S) = max_{t∈[0,1)} m(t)`. LRC(n) for `k = n−1` runners ⟺ `M(S) ≥ 1/n`. Here `n = 14`, `k = 13`.

---

## Part 1 — Local maxima only at pair-sum rationals

> **Theorem.** `m` is continuous piecewise linear with every slope in `{±v_i}` (never 0, so no
> plateaus), and every local maximum point `t*` of `m` satisfies
> `t* = p/(v_i + v_j)` for some `1 ≤ i ≤ j ≤ k`, `p ∈ Z`
> (the case `i = j` being the single-runner peaks `t* = p/(2v_i)`, `p` odd).

*Proof.* Each `‖v_i t‖ = min(frac(v_i t), 1 − frac(v_i t))` is piecewise linear with slopes `±v_i`;
a finite min of such is piecewise linear with slopes in the union, and no slope is 0.
At a local max `t*` of `m`, the left derivative is `> 0` and the right derivative is `< 0`.
So some active constraint `i` (with `‖v_i t*‖ = m(t*)`) is **rising** on the left:
`frac(v_i t*) = m(t*) ∈ (0, 1/2]`, and some active `j` is **falling** on the right:
`frac(v_j t*) = 1 − m(t*)`.
- If `i = j`: `frac = m = 1 − m` forces `m(t*) = 1/2`, i.e. `v_i t* ∈ 1/2 + Z`: a peak
  `t* = p/(2v_i)`, `p` odd. ∎(case)
- If `i ≠ j`: `frac(v_i t*) + frac(v_j t*) = 1`, hence `(v_i + v_j)·t* ∈ Z`. ∎

(Difference-crossings `frac(v_i t) = frac(v_j t)` — both rising or both falling — never give a
local max: the min keeps a one-signed slope through them. Kinks at `‖v_i t‖ = 0` are minima.)

## Part 2 — Witness-ruler boundedness

> **Corollary.** `M(S)` is attained at `t* = p/q` with `q = v_i + v_j ≤ 2·Vmax` (some `i ≤ j`);
> `M(S) ∈ Q` with denominator dividing some pair sum; and
> `M(S) = max_{i≤j} max_{0<p<q} min_l dist(v_l·p mod q, {0, q})/q, q = v_i+v_j` —
> a finite INTEGER computation (grid-free; MISTAKE-130-proof; native_decide-shaped).

This resolves klein-S207's observed mystery ("the witness is on a DIFFERENT ruler: `M = 3/13` at
`τ* = 11/39` and `39 ∤ 91`"): lawful, because `11/39 = 22/78` and `78 = 24 + 54` is a pair sum of
that cluster; the active pair at the max is `{24, 54}` with residues `≡ ±9 (mod 39)`. The
Vmax-ruler can never host the witness (klein-S207: ruler points are never lonely) — **the pair-sum
rulers always do.** Note also: mac-mini-S64's "no bounded-q reduction" (adversarial min-q ≥ 37,
unbounded over the family) is compatible — the bound `q ≤ 2Vmax` is *relative* to the cluster, not
absolute.

## Part 3 — The event-realization criterion (the finite form of "good gap ⟹ lonely")

Fix `i < j`, `q = v_i + v_j`. Call `t` with `q·t ∈ Z` an **(i,j)-event**. At an event, either
`frac(v_i t) = frac(v_j t) = 0` (a *zero event*, `p ∈ (q/g)Z`, `g = gcd(v_i, v_j)`; these are
spaced `≥ 2` events apart since `g ≤ v_i < q/2`), or `frac(v_i t) + frac(v_j t) = 1`, which forces
**`‖v_i t‖ = ‖v_j t‖ = w_{ij}(t)/2`** where `w_{ij}(t) := ‖v_i t‖ + ‖v_j t‖`
(the *pair-safety function*; in co-offset language: the observer phase sits exactly at the
midpoint of the arc between teeth `i` and `j` through it, of width `w_{ij}`).

> **Event-Realization Lemma.** If there is an interval `I` with `|I| ≥ 2/(v_i+v_j)` on which
> (a) `w_{ij}(t) ≥ 1/7` and (b) `‖v_l t‖ ≥ 1/14` for every `l ∉ {i,j}`,
> then `I` contains a nonzero event `t*` with `m(t*) ≥ 1/14` — a lonely time.

*Proof.* `|I| ≥ 2/q` contains two consecutive events, at least one nonzero. There
`‖v_i t*‖ = ‖v_j t*‖ = w_{ij}(t*)/2 ≥ 1/14` by (a); the rest by (b). ∎

**This is the finite, per-pair replacement for the equidistribution ρ_K → ρ* node**: no drift, no
Vmax-ruler, no measure theory — the fast phase is *pinned to the midpoint by the event*, and the
only question left is whether some pair keeps its safety function above `1/7` for two event
spacings while the other eleven runners stay `1/14`-safe. The known partial discharges are special
cases: klein-S205's drift embed (`Vmax > 1.41·spread`) constructs such an `I` inside one
Vmax-ruler cell; kps-S106's scale separation constructs it inside the slow-safe window; kps-S112's
continuum bridge is the `q → ∞` limit where events are dense. The tight AP `{1..13}` is the
degenerate boundary instance: at `t* = 1/14` the pair `(1,13)`, `q = 14`, has `w = 1/7` EXACTLY at
its event — a tangency, `|I| = 0` (the pinch/lemniscate node of opus-S177, in closed form).

## Part 4 — The Schur-triple kill rule (why E3 governs tightness)

> **Proposition.** If `(v_i + v_j) | v_l` for some runner `l` (in particular if `v_i + v_j = v_l`
> is a **Schur triple** in `S`), then EVERY (i,j)-event `t` has `‖v_l t‖ = 0`: the pair-sum ruler
> `q = v_i+v_j` is **dead** (no witness can ever sit on it).

*Proof.* `v_l·(p/q) ∈ Z` for all `p` when `q | v_l`. ∎

So each Schur triple `a + b = c` in `S` removes the `(a,b)` ruler from the witness supply — the
exact pair-sum analog of klein-S207's "the observer kills its own ruler" (`e_0 = 0` there; here
runner `l` sits at the origin at every `(i,j)`-event). **This is the mechanism behind
opus-S182's finding that the resonance sum / tightness is governed by E3 (Schur triples), not
additive energy E2:** the AP `{1..13}` maximizes E3 = 78 among nonzero 13-sets ⟹ a maximal number
of dead rulers ⟹ the witness supply is squeezed to the single tangency event `1/14` on the ruler
`14 = 1+13` (note `14 ∤ v_l` for all `l` — the one ruler the AP does NOT kill, precisely because
the AP is non-covering at `q = 14`!). Dilation-invariance matches (Schur triples are
dilation-invariant, not translation-invariant), and the sum-free separator (`{1,3,…,25}`, E3 = 0:
no dead rulers, loose) is explained.

**The two sieves are one sieve.** The non-covering dispatch (THM-369/THM-523: `q ∈ {2..14}`
dividing no speed ⟹ `t = 1/q` lonely) and the pair-sum witness structure are the same statement at
two scales: a modulus `q` hosts a witness iff its residue configuration `{v_l mod q}` admits a
multiplier `p` with all residues in the middle band `[q/14, 13q/14]`. Small `q ≤ 14` with no zero
residue: `p = 1` works trivially. General witness: some pair-sum `q = v_i + v_j` works — and by
Part 1 this is EXACT:

> **LRC(14) ⟺ for every primitive covering 13-set `S`, there exist `i ≤ j` and `p` such that
> every `v_l·p mod (v_i+v_j)` lies in `[⌈(v_i+v_j)/14⌉, ⌊13(v_i+v_j)/14⌋]`.**

A dead ruler (`q | v_l`) has a zero residue and can never host; covering says exactly that the
*constant* rulers `q ∈ {2..14}` are all dead-at-`p=1`-in-the-weak-sense (some zero residue). The
conjecture is that the pair sums always supply a live ruler.

## HYP-5720 — the realization slack (OPEN, tested this session)

Define `σ(S)` = (length of the lonely component `{m ≥ 1/14}` at the witness) × (active pair sum).
`σ(AP) = 0` (tangency). **Conjecture: `σ(S) ≥ c₀ > 0` uniformly over primitive covering 13-sets**
(equivalently: on covering clusters the witness is never a tangency — the strict-cushion
counterpart of klein-S206's strict good-period margin ≥ 1.2353 and mac-mini-S64's covering
`min M = 1/12 > 1/14`). Verification: `04-computation/lrc14_pair_sum_ruler_macmini_S65.py`
(+ `.out` in 05-knowledge/results/).

## Consequences for the Lean endgame

1. Any finite dispatch (bounded window, `Vmax ≤` bound, census legs) should enumerate pair-sum
   events, not grids: exact integers, `native_decide`-shaped, and COMPLETE (the max provably lives
   there) — upgrading the certified-rational direction monad-explorer-S1 recommended (their (d)).
2. `hembed`'s open corner can be restated pair-locally (Event-Realization Lemma): the hypothesis
   to discharge becomes "some pair keeps `w_{ij} ≥ 1/7` for `2/(v_i+v_j)` while others stay safe",
   which is a *finite intersection of sawtooth conditions* — no equidistribution quantifier.
3. The Schur-triple kill rule gives the a-priori bridge target a combinatorial form: covering
   13-sets cannot kill all their pair-sum rulers (E3 budget vs. 91 rulers) — a candidate for the
   sum-free/Schur extremal literature (opus-S182's two-step target, now with a mechanism).

**Depends on:** none (elementary). **Related:** THM-420, THM-369/THM-523, klein-S205/S206/S207,
kps-S106/S112, opus-S177/S182, monad-explorer THM-665, LEM-012/013, MISTAKE-130.
