# New strategies for closing hp0cap and hpartA — the Farey-level-14 program

**mac-mini-2026-07-01-S96 (HYP-3852).** The two named analytic sorries carrying LRC(14)
(`hp0cap`: `p0(E) ≤ cap_k`; `hpartA`: `witnessG2 > 0 ⟹ Mreach ≥ 1/14`) re-attacked with
the S93–S95 machinery (THM-592 radius-derivative, THM-593 unit-residue, THM-594 pair law +
Mirsky–Newman, THM-596 bands). Everything below is pure mathematics; verification in
`lrc_farey_level14_decode_macmini_S96.py` (exact rationals throughout).

---

## 0. The organizing discovery: one Farey level, four appearances

The number `n = 14` appears as the SAME Farey-structural level in four previously separate
places, and the two open obligations sit inside this one structure:

1. **Pair overlaps branch at `p + q = 14`** (THM-594(B)): `|D_p ∩ D_q| = 2r/q` for
   `p+q ≤ 14`, the wrapped form `(q+2p−14)/(7pq)` beyond.
2. **The collapse constant is the level-14 Mordell–Tornheim slice** (NEW, exact at
   n = 6, 8, 14, 20):
   `c_AP(n) = Σ_{a∈(ℤ/n)ˣ} 1/(a(n−a)) = 2n · Σ_{p+q=n, p<q, gcd=1} 1/(p·q·(p+q))`.
   klein-S87's mean-loneliness layer-cake `E = Σ_Farey 1/(qq′(q+q′))` is the FULL MT sum;
   the cusp slope is its level-n slice. (Each profile merge-kink at `d/(v+w)` has jump
   `(v+w)/(vw)`; jump × radius = `d/(vw)` — the layer-cake literally integrates the Farey
   grid, term by term.)
3. **The final-window kink bands are the numerator strata of the window's Farey/Stern–
   Brocot tree** (THM-596 refined): reduced `d′/q*` in `(1/15, 1/14)` with `q* ∈ (14d′, 15d′)`:
   `d′=2: {29}` (the MEDIANT `2/29 = 1/15 ⊕ 1/14`), `d′=3: {43,44}` (its two children),
   `d′=4: {57,59}` (58 excluded by reducedness), `d′=5: {71,72,73,74}`, `d′=6: {85,89}`.
4. **The 23520/392 decode** (the S95 correction constants): `392 = 2·14²`
   (`27/392 = 1/14 − 1/392`); `Λ_GW − Λ_AP` has slope exactly `29/60 = 2(1/5 + 1/24)`
   on `[1/15, 2/29]` — the GW outlier pair `(5, 24)` and its mirror; `5 + 24 = 29 = 2n+1`,
   and `1/23520 = (29/60)(2/29 − 27/392)`, `23520 = 60·392`. Bonus universal fact: the GW
   outlier is `2(n−2)`, so its mediant-merge partner is the runner **5 for every n**
   (`2(n−2) + 5 = 2n+1`), with jump `2(1/5 + 1/(2n−4)) → 2/5`.

Reading: the window `(1/(n+1), 1/n)` is a Farey parent pair; its mediant `2/(2n+1)` is
simultaneously GW's identity onset, the first exposed band, and the band-emptying anchor.
The ladder should CLIMB THE STERN–BROCOT TREE: each rung anchored at a mediant, each
rung's defect confined to strictly deeper tree levels with `q*` growing and exposure
requiring 13 simultaneous residue clearances mod `q*` — a convergent defect series.

---

## 1. Strategy for hp0cap: exact multi-overlap arithmetic replaces Vitali

**State:** THM-534 proves `p0 ≤ L_y` (moment-LP dual). Open half: "consecutive maximizes
`L_y`", split (S28/S29) into [bounded-spread finite check] + [far elements drop `L_y` by
≥ 0.044 — currently a Vitali/marginal-uniformity estimate].

**New route.** `L_y = q0 + q6 + q3/10` is a linear functional of the sector miss-
distribution, and every `q_j` is an inclusion–exclusion over sector-hit events whose
d-fold intersection measures are **exactly computable in closed form** by the multi-
overlap generalization of THM-594(A/B): the Fourier coefficients of each `1_{D_v}`-type
indicator are `sin`-kernels on divisor frequencies, so d-fold overlap integrals are
finite cosine-series evaluations = **Bernoulli-polynomial expressions** (the classical
`Σ_s cos(sθ)/s^d = (2π)^d B_d(θ/2π)`-family; the pair case d = 2 is THM-594(B), already
proved two-branch). Consequences:

- (i) **The far-element drop becomes exact arithmetic.** For `E = B ∪ {w}`, every term of
  `L_y(E)` involving `w` is a finite sum of overlap deficits with explicit branch-2
  closed forms; the "≥ 0.044" decorrelation is then a rational inequality per base shape
  — no measure theory, no Vitali. (Probe: consec_9 `L_y = 2111/4410 ≤ cap_9 = 1979/4004`
  with margin 0.0156; far swaps drop `L_y` by ~0.086, with visible arithmetic
  oscillation in `w` — exactly what an exact per-pair formula handles and a uniform
  estimate fights.)
- (ii) **The bounded-spread census is exact and finite.** THM-529 caps the spread at 30;
  the miss-distribution of every bounded shape is a rational vector via the same closed
  forms; "consec maximizes L_y" on the compact set is a finite exact comparison.
- (iii) **Proof shape:** hp0cap = [finite exact census (ii)] + [per-shape branch-2
  deficit bound (i)] + THM-534. The only lemma to write is the d-fold overlap closed
  form (d ≤ 7 suffices — 7 sectors), which is the Bernoulli/MT layer of §0.2. Expected
  side effect: the `cap_k` rationals (2243/5880, 1979/4004, 55/91, 66/91) should emerge
  AS Bernoulli/MT-slice combinations — worth checking; if they do, cap and slope are
  finally the same species (opus §7.4's ζ(2)↔ζ(−1) note made quantitative).

## 2. Strategy for hpartA: windowed Mirsky–Newman

**State:** `hpartA` = THM-527 Part A: positive two-scale witness density (`G2 > 0`, the
small-part good-period measure) implies actual reach `≥ 1/14`. The danger: the large
cluster's arcs could cover the whole positive-measure window.

**New route.** The `G2` window is a finite union of rational-endpoint intervals (the
S93 breakpoint grid gives its exact structure). Killing the reach requires the large
cluster to cover **every component** of the window at radius 1/14 — a LOCAL exact-cover
event. THM-594(C) says finite distinct-speed systems never tile globally; what hpartA
needs is the **windowed quantitative version**:

> **Target lemma (windowed MN).** For any interval `I` and finite speed set `F` whose
> danger arcs at radius `r` have total mass `A_I` inside `I`: the uncovered measure
> `|I ∩ {C_F = 0}| ≥ |I| − A_I + (windowed Parseval defect)`, where the defect isolates
> a divisor-minimal frequency against a Beurling–Selberg majorant of `1_I` (band-limited
> window ⟹ the frequency-`w` coefficient survives with an explicit `O(deg⁻¹)` loss).

Ingredients all exist in-repo: the BS majorant machinery (HYP-3765's signed
Beurling–Selberg smoothing), the divisor-minimal coefficient (THM-594(A/C)), the
critical-mass bookkeeping (THM-594(E)). The point: Part A stops being an
equidistribution estimate and becomes **rigidity + window Fourier bookkeeping** — the
same one theorem (MN) read through a window. Where the cluster is divisor-chained (the
only near-tiling direction, THM-594(C) corollary), the 2-adic descent (THM-580) has
already quotiented it away — the fixed-locus alignment of the S94 synthesis, now doing
work inside a named sorry.

**Fallback within the same frame:** if the BS loss fights, anchor the window at the
Stern–Brocot mediants of its components (§0.3): each component's first exposed band is
its mediant; the peel + ladder (THM-592(v), kps/opus union floor) then transports the
G2 mass to `1/14` with the defect confined to depth-≥2 tree levels — the
"anchor-vs-band tradeoff" of THM-596(v) run recursively.

## 3. Tangents of tangents (recorded for future sessions)

- **The SB-renormalization ladder:** rung k anchored at the depth-k mediant; exposure at
  depth k needs 13 residues clear mod a `q*` that grows Fibonacci-fast along the tree;
  the defect series converges geometrically; the deep-well CF signature `[0; 13, 14]` is
  a PATH in this tree — "danger zone = short tree path" would unify HYP-3792 with
  THM-596.
- **Universal 5:** the GW mediant partner is 5 for all n. Is there a runner that plays
  this role for every sporadic tight family (the q=8 relift beater's partner at its
  mediant `2/17`... check)? If the sporadics' mediant partners are bounded, the
  tight-locus finiteness (fattening) gets a new finite handle.
- **jump×radius = `d/(vw)`:** the layer-cake E as a literal sum over profile kinks —
  differentiating klein's `E = Σ 1/(qq′(q+q′))` BY LEVEL gives the collapse slice (§0.2);
  differentiating BY TREE DEPTH should give the defect series of the SB ladder. One
  generating function, three derivatives.
- **`13 = (29−3)/2`... resisted numerology:** kept only what has structure; the decode
  items in §0.4 all have mechanism attached. (Note: the owner's "23530" read as 23520;
  no structure found at 23530 = 2·5·13·181 beyond 181 = the deep-well resonance
  denominator — flagged, not pursued.)

-> THM-592/593/594/596, THM-534, THM-527, HYP-3834/3950 (peel), HYP-3765 (BS),
HYP-3852, OPEN-Q-108.
