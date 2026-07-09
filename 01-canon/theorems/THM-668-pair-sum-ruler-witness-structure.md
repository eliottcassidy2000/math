# THM-668 — The pair-sum ruler theorem: the witness always lives at t = p/(v_i+v_j), and LRC(14) is a band statement on pair-sum moduli

**Status:** PROVED (parts 1–4, elementary, proofs below) + VERIFIED (300/300 random sets: pair-sum
max = full-breakpoint max, exact; klein-S207's 91-cluster witness reproduced and explained; AP
tangency σ = 0 exact). The uniform-slack refinement HYP-5720 is OPEN (data below — honest scope
note included: it reformulates the covering cushion rather than independently delivering it).
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

## HYP-5720 — the realization slack (OPEN; first data this session, with an honest scope note)

Define `σ(S)` = (length of the lonely component `{m ≥ 1/14}` at the witness) × (active pair sum).
`σ(AP) = 0` (tangency). Conjectured: `σ(S) ≥ c₀ > 0` uniformly over primitive covering 13-sets.
**Data (exact):** 966 exhaustive covering sets `[1,18]`: min `σ = 0.049` (> 0, at the
`min M = 1/12` set); adversarial hill-climbs: min `σ = 0.124` (cap 60), **`0.0343` (cap 120)**,
`0.091` (cap 200). The cap-120 dip below the exhaustive minimum means σ is NOT yet evidently
bounded away from 0 — OPEN, trend unclear.
**Honest scope note (self-caught):** since the witness IS an event (Part 1), the lonely component
at the witness always contains its own event; so `σ > 0 ⟺ M(S) > 1/14` strictly, and HYP-5720 is a
scale-adjusted reformulation of the covering cushion (klein-S206 margin / mac-mini-S64
`min M = 1/12`), NOT an independent realization mechanism. The Event-Realization Lemma's
long-interval hypothesis (`|I| ≥ 2/q`) is sufficient but empirically not how most witnesses
realize (e.g. the 91-cluster witness component spans only 0.44 event spacings). The counting
route to the a-priori bridge must therefore use the sharper fact that events sit AT the local
maxima (attracted, not independent) — see Part 1, and the ruler-supply question in Part 4.
Verification: `04-computation/lrc14_pair_sum_ruler_macmini_S65.py` (+ `.out`).

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

---

## ADDENDUM (mac-mini-S65 cont., HYP-5730) — the live-ruler certificates and the two-domain structure

Attacking the Schur-budget statement head-on. Four elementary PROVED certificates for "ruler `q`
is live", then the census. Notation: `m := ⌈q/14⌉ − 1` (danger radius), `D := {d : |d| ≤ m}`,
`B_l := {p ∈ Z/q : v_l·p mod q ∈ D}` (runner `l`'s bad multipliers; `0 ∈ B_l` always).

> **C0 (window).** The ruler `q = Vmin + Vmax` is live at `p = 1` iff `Vmax ≤ 13·Vmin`.
> *Proof:* `v_l/q ∈ [Vmin/(Vmin+Vmax), Vmax/(Vmin+Vmax)] ⊆ [1/14, 13/14] ⟺ r ≤ 13`. ∎
> (This is kps-S28's `spread13_lonely` witness recognized as a pair-sum event — one more instance
> of Part 1's universality.)

> **C1 (gcd-exact ledger).** With `g_l = gcd(v_l mod q, q)`:
> `|B_l| = g_l·(2⌊m/g_l⌋ + 1)`, and `B_l = B_k` whenever `v_l ≡ ±v_k (mod q)` (D symmetric) —
> for `q > Vmax` the merges are exactly the `r(q) ≤ 6` pair representations of `q`. Since every
> `B` contains 0: `|bad ∩ [1, q−1]| ≤ Σ_classes (|B| − 1)`. **If that sum `< q − 1`, `q` is
> live.** Rule-of-thumb form at `g ≡ 1`: fires iff `2·r(q) + (cheapening) ≥ 13`.

> **C2 (divisor descent).** If `k | q`, `k > 14`, and some `s ∈ [1, k−1]` has all
> `v_l·s mod k ∈ [⌈k/14⌉, k − ⌈k/14⌉]`, then `p = (q/k)·s` is banded mod `q`.
> *Proof:* `v_l·p mod q = (q/k)·(v_l·s mod k)`, and `(q/k)·⌈k/14⌉ ≥ q/14`. ∎
> (For `14 < k ≤ 28` the mod-`k` condition is exactly "avoid `{0, ±1}`". Recursion on the divisor
> lattice; the base is THM-420-Lemma-B-style counting.)

> **C3 (six-pair prime).** If `q` is prime, `Vmax < q`, `q ≡ t (mod 14)` with `t ≥ 3`, and `q`
> has SIX pair-sum representations, then `q` is live with `≥ t − 1` multipliers.
> *Proof:* 6 merges leave `13 − 6 = 7` distinct `B`-classes, each `|B| = 2m+1` (units);
> `|bad ∩ [1, q−1]| ≤ 7·2m = 14(⌈q/14⌉−1) = q − t`. ∎ (The `13 = 2·6+1` pairing wall is beaten
> exactly by the pair-merge; this is the sharpest pure-counting case.)

### Census (all exact; scripts `lrc14_live_ruler_certificates_…` / `…_certificate_stress_…` / `…_blocking_configs_…`)

- **Soundness:** 0 unsound firings over 500 random sets × all rulers.
- **Covering `[1,18]` (966 exhaustive): C1 alone certifies 100%.** (C2: 12.7%, C3: 8.7%.)
- **Random covering cap 60 (600): C2 alone certifies 100%;** C1 79.5%; 0 residuals.
- **Structured adversaries all certified:** monad-S2's detuned harmonics (the family that defeated
  the φ-interval composition) fall to C2 at `k = 23, 25`; covering near-AP blocks fall to C1 at
  `q ∈ [17, 26]`; the S65 min-supply set fires C1 at `q = 23`; worst7Struct@91: 34 certified rulers.
- **Defeaters of pure counting exist, and they are SMALL-SCALE:** hill-climbs found covering sets
  with 0 counting-certified rulers, e.g. `{1,2,3,5,6,8,…,14,23}` (live at `q=21, p=5` via the
  residue collision `23 ≡ 2`, which the union bound cannot see) and `{9,10,14,…,29}` (live at
  `q ∈ [32,36]`, `p = 1` — C0-shaped, all speeds in band). All defeater live-rulers sit at
  `q ≤ 36`.
- **The genuinely open annulus is certificate-saturated:** in the sliver `r > 13`, `Vmin ≥ 18`
  (all pair sums > 36 — counting-only territory): 250/250 random certified; targeted hill-climbs
  (5 restarts × 200 steps) could not push below **38 certified rulers** — the adversary cannot
  even approach a defeater at scale.
- **Blocking census** (13-subsets of `(Z/q)\{0}`, exact, up to dilation): blocked fraction falls
  `100% (q=15) → 80% (17) → 5.7% (23) → 7.1% (26)`; lex-first blocked classes are all
  near-intervals (longest-AP 12–13). Small rulers are classification territory, not counting
  territory.

### The two-domain structure (the finding)

**The live-ruler theorem factors:**
`[q ≤ Q₀ ≈ 36: a BOUNDED, ABSOLUTE domain — exact blocking-classification / exhaustive sweep]`
`+ [q > Q₀: counting certificates C0 ∨ C1 ∨ C2 (∨ C3), empirically total and adversarially
robust]`.
Every observed defeat of counting happens below Q₀; every small-Vmax set is decidable by finite
check (and `[1,18]` is already exhaustively covered, klein-S206/mac-mini-S64). The remaining
PROOF obligations, explicitly:
1. **Small domain:** exhaust primitive covering 13-sets with `Vmax ≤ ~30` (C-code scale, one run —
   the analog of klein's 966 at 18) — turns the domain into a machine-checked lemma.
2. **Large domain:** prove "covering + all pair sums > Q₀ ⟹ some C1/C2 certificate fires" — the
   candidate mechanism is C2 divisor abundance (91 pair sums, each even sum has `k = q/2 > 14`;
   blocked fraction per `k` is 5–7% and falling) + C1 representation-richness. This is a finite
   residue-geometry statement; the sliver data (min 38 certified) suggests enormous slack.
3. **Lean:** all four certificates are one-page elementary lemmas over `Z/q` — native_decide-free,
   grid-free.

**Honest scope:** the certificates prove liveness per set/family; the infinite compressed stratum
still needs obligation 2 as a theorem (or the fleet's analytic Kronecker route) — the two
approaches meet exactly at the R1 box, and the certificates decimate it from the arithmetic side.

**Depends on:** none (elementary). **Related:** THM-420, THM-369/THM-523, klein-S205/S206/S207,
kps-S106/S112, opus-S177/S182, opus-S183/LEM (Schur extremal), monad-explorer THM-665/THM-666,
boxeph-S1 LEM (P-separated composition), death-star-S1 (pure-cluster corner), LEM-012/013,
MISTAKE-130.
