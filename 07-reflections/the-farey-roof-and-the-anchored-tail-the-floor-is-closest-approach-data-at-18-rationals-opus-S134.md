---
source: opus-2026-07-07-S134
status: synthesis + new mechanism (Farey roof PROVED-shape; F7 finitization PROVED; A″-F6 candidate VERIFIED adversarially, OPEN)
tags:
  - lonely-runner
  - LRC14
  - density-floor
  - three-gap
  - farey
  - anchored-tail
  - observer-lens
  - history-audit
---

# The Farey roof, and the anchored tail: the density floor is closest-approach data at 18 rationals

**Owner mandate this session:** *understand the past LRC(14) work deeply — how the proof state
changed as corrections surfaced, what is valid, what we've been missing; uncover forgotten
factoids that were true but not load-bearing, extend and connect them; aim at how/why, not just
if.* This note is the synthesis. Everything computational is in THM-637 and the four S134 scripts.

## 1. What the history actually says (validity audit, compressed)

Tracing every era (sieve prehistory → p0 skeleton → q-witness → covering-min/Φ₆ saga → safe-band
map → Route 2 rise & kill → Route 1 present) yields one clean conclusion: **the survivors are all
sup-side reductions plus one analytic object.** GREEN and load-bearing: the q-witness sieve
(oracle-S18), the coarse reduction (kps-S53/THM-636), primitive-saturated equivalence (kps-S56),
far-peel; verified-not-formal: the conjugate witness (klein-S152/opus-S131). The single analytic
residual is stated two ways — the MOAT (sup: single-scale non-AP ⟹ M ≥ 1/13) and the FLOOR
(average: μ_{1/7} ≥ …) — and mac-mini-S40 proved these need *different tools* (signed/location vs
moments). The recurring failure modes: wrong-object errors (MISTAKE-116/117), bounded-range and
dilation artifacts (kps-S56), weak-adversary sampling (MISTAKE-100/101/102, and again TODAY:
klein-S153's "AP minimizes E[maxgap], descent converges to AP" — 60 random families + local
descent — is refuted exactly by death-star/boxeph/monad's structured parity-interlaced families),
and tail-vs-mean conflation (kps-S58's "(1) AP-minimality of E[maxgap]"). **monad-S1's audit is
the current sharpest correction:** the mean route's honest bar is T* = 1/7 + (6/7)m_P ≈ 0.19127,
the record adversary sits at 0.19699 (margin +0.0057 and shrinking), and the mean serves only the
k=13/P=∅ leg. **The load-bearing open lemma is (A′): per-k tail minimality of μ_{1/7} at the AP.**
So this session attacked the *tail*, not the mean.

## 2. The Farey roof — the mechanism behind every exact constant (what we'd been missing, part 1)

klein-S153 noticed numerically that the AP's origin gap equals its max gap a.e. That observation
is a **theorem with a three-line proof** (THM-637a): removing 0 from the with-0 config merges the
two flanking gaps `‖q_L x‖ + ‖q_R x‖`, which is the maximal three-distance value. Consequently on
every Farey-k cell `(p/q, p′/q′)`:

> **maxgap(AP_k, x) = q(x−p/q) + q′(p′/q′−x)** — the "Farey roof", linear on each cell,
> node value `1/q` at each Farey point p/q.

Everything the fleet computed piecemeal falls out of this one function:
- `E[maxgap(AP_k)] = Σ 1/(qq′²)` (consecutive Farey pairs) — reproduces 93/440, 43/140, … and
  gives all k; `E[min_i frac(ix)] = half` (kps-S58's reflection identity is the roof's mean).
- `μ_θ(AP_k)` = the roof's superlevel measure — reproduces 477/1078 and the whole canon table.
- The good set of the AP is **exactly** the q ≤ 6 Farey neighborhoods: the AP's floor is *pure
  resonance-window measure, zero bulk*. The roof nodes `1/q` are the Ostrowski/Kravitz rung
  heights — **the M-ladder and the μ-floor are two readings of one piecewise-linear function.**
- `q = 7` windows are exactly marginal (node = 1/7): **the apex prime 7 is invisible to the
  density floor.** The composite-14 hardness (klein-S151) lives entirely on the sup side. This is
  a structural reason the floor is tractable while the moat is not — beyond sup-vs-average.
- Small numerology with content: pigeonhole covers k ≤ 7 (1/k ≥ 1/7); for k ≥ 8 the iid mean
  H_k/k clears 2/(k+1) precisely because ln k > 2 for k ≥ 8 — the crossover e² ≈ 7.39 sits exactly
  at the pigeonhole boundary. The threshold 1/7 = 2·(1/14) is ALSO the iid size-biased mean gap
  2/(k+1) at k = 13 exactly — the break-even that killed every single-statistic bound (kps-S58).

## 3. The anchored tail — finitizing (A′) (what we'd been missing, part 2)

The failed reductions all tried to catch the max gap with ONE statistic (origin gap, length-biased
mean, fixed φ). The repair is not a cleverer single statistic but a **finite family of anchors**:

> **F₇ identity (exact, every E):** the Farey-7 points have max spacing 1/7, so any gap > 1/7
> contains one; hence `{max_{a∈F₇} gap∋a > 1/7} = {maxgap > 1/7}` **pointwise**, and
> **μ_{1/7}(E) is identically an 18-anchor statistic** — the joint law of the config's closest
> approaches to 18 fixed rationals. (A′) is *exactly* an inhomogeneous-approximation statement
> at rational targets, governed by residue profiles mod small q and the diameter.

> **(A″-F₆) candidate (tight, stronger, structured):** with the 12 Farey-6 anchors,
> `t_{F₆}(E) ≥ μ_{1/7}(AP_k)` for every affine-normalized k-set, equality iff AP. Survived:
> normalized corpus, dedicated descent (k = 13, 12, 10 — all converge back to the AP), exhaustive
> 1-swap and sampled 2-swap tight-locus scans. Zero violations; zero (A′) violations anywhere.

Two mechanisms make F₆ the principled set:
- **Window-exactness:** near x = p/q ± δ (q ≤ 6) the config clusters at the q-grid and each
  cluster spreads *one-sidedly* by e·δ, so every inter-cluster gap contains its flanking grid
  point j/q. Inside windows the anchored max IS the max gap. At the AP (pure window), t_{F₆} = μ
  exactly — the bound is tight precisely at the conjectured minimizer, the signature of a correct
  extremal formulation.
- **Escape confinement:** F₆'s only spacings exceeding 1/7 are (0, 1/6) and (5/6, 1); an
  unanchored >1/7 gap must hide there with length in (1/7, 1/6). Empirical bulk capture ≈ 0.99–1.0.

And one hard-won lesson repeated on a new axis: **all naive-anchor failures were affine images of
the AP** ({2..14} = AP+1, odds = 2AP+1, {3+7i} = 7AP+3). μ is affine-invariant; anchored tails are
not; so the lemma must be stated over affine-normalized representatives. This is kps-S56's
dilation-artifact lesson recurring on the anchor side — the invariance group of the *statistic*
must be quotiented out of the *domain* before any extremality claim is meaningful.

## 4. Forgotten factoids, now wired in (the connective tissue)

- **The observer lens** (memory: it cracked the covering-min): anchor 0 IS the observer;
  `gap∋0 = m⁺ + m⁻` is the orbit's two-sided closest approach to the observer, and the AP
  saturates its max gap AT the observer (THM-637a). The whole anchored frame is the observer lens
  multiplied over 18 rational vantage points.
- **THM-591 (inhomogeneous-AP linear law, forgotten Cluster-A fact):** exact control of the AP's
  shifted-target loneliness — precisely the AP-side data the anchored joint law needs at each
  a = p/q. It was proved for the *shifted problem* nobody was using; the anchored tail is where
  it plugs in.
- **THM-580 (2-adic odd/even exact recursion) + THM-590 (apex-7 cyclotomic gap):** the anchor ½
  splits E by parity (odd elements see ½ as their own observer under x ↦ x+½·(odd)), so the
  anchored laws inherit the 2-adic descent; and the F₆/F₇ boundary (q=7 marginal) is the apex-7
  again — the floor stops exactly where the cyclotomic hard core begins.
- **oracle-S18's residue sieve:** saturation (covering residues mod q) was a sup-side tool; in the
  anchored frame it reappears on the *average* side as the thing that sets window masses — the
  window at p/q is fat iff E's residues mod q are few/clustered; covering all residues mod q
  shrinks it. The additive↔multiplicative tension (mac-mini-S15) becomes, anchor by anchor, a
  finite tradeoff: window mass (correlation) vs bulk capture (decorrelation).
- **boxeph-S1's co-offset audit:** speed-config vs offsets-with-0 differ by a rotation by
  (base)·x — the same translation-covariance that forced affine normalization above. One
  phenomenon, three sightings (kps-S56 dilation, boxeph co-offset, S134 anchors).

## 5. Where this leaves the program (honest)

- (A′) per-k remains open and unbreached — now with an exact finite reformulation (F₇) and a
  tight extremal candidate (A″-F₆). The natural proof shape: window-exactness (rigorous,
  lemma-sized) + per-window mass comparison vs the AP (residue profile + diameter, sieve-flavored)
  + bulk capture via the two escape arcs (second-moment friendly — it's an average, mac-mini-S40
  licenses moments here). None of that is done; each piece is now finite-dimensional.
- The mean sidecar (E[maxgap] vs T*) is monad's razor; my independent descent reconverged to
  their record family {2,…,22 even, 11, 13} = 0.19699 and found trisection variants *worse*
  (0.204–0.216) — mechanism: mod 3 has two nonzero residues to fill but only one spare slot
  budget-neutral way; mod 2 has one. Bisection is the uniquely cheap interlacing. If the sidecar
  dies below T*, nothing above is lost.
- Corrections filed: klein-S153's "AP minimizes E[maxgap]" (INDEX note; weak-adversary artifact),
  kps-S58's part-(1) framing (already flagged by S133; INDEX note added).
- NOT touched: the moat (sup side) — mac-mini-S40's division stands; anchors are an average-side
  tool. The Part-A finite-Vmax bridge (o(Vmax) arc bound, empirically ~S^0.45) remains the other
  open link and is untouched by everything above.
