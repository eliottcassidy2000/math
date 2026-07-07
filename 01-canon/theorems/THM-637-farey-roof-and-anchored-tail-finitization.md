# THM-637 — The Farey roof (AP max-gap mechanism) and the anchored-tail finitization of μ_θ

**Source:** opus-2026-07-07-S134
**Status:** parts (a),(b),(c) PROVED (elementary, from the three-distance theorem; machine-verified k ≤ 13 to 1e-15);
part (d) PROVED (3 lines); part (e) VERIFIED-computational (candidate lemma, OPEN).
**Scripts:** `lrc_farey_roof_apmaxgap_opus_S134.py`, `lrc_anchored_tail_A2_opus_S134.py`,
`lrc_farey6_anchored_tail_opus_S134.py`, `lrc_F6_tightlocus_swapscan_opus_S134.py` (+ .out in 05-knowledge/results).

Throughout: `config(E, x) = {frac(e·x) : e ∈ E}` on the circle, `maxgap` = largest circular gap,
`gap∋a` = the circular gap containing the point `a`, `μ_θ(E) = Leb{x : maxgap > θ}`.

---

## (a) Origin-gap saturation (sharpens klein-S153's a.e. numeric to an identity)

**For the AP `E = {1,…,k}` and every irrational `x`: `gap∋0 = maxgap`.**

*Proof.* In the with-0 config `{frac(i·x) : 0 ≤ i ≤ k}`, the two gaps flanking 0 have lengths
`δ₁ = ‖q_L x‖` and `δ₂ = ‖q_R x‖`, where `q_L, q_R ≤ k` are the one-sided best-approximation
denominators of `x` at height `k`. By the three-distance theorem every gap of the with-0 config
has length in `{δ₁, δ₂, δ₁+δ₂}`. Deleting the point 0 merges its two flanking gaps into one gap
of length `δ₁+δ₂` — the maximal possible value — and leaves all other gaps unchanged. ∎

Consequence: `maxgap(AP_k, x) = m⁺(x) + m⁻(x)` where `m⁺ = min_{i≤k} frac(ix)`,
`m⁻ = min_{i≤k}(1−frac(ix))` — the **two-sided closest approach of the orbit to the observer**.
(The observer-lens object; kps-S58's identity `E[gap∋0] = 2E[min]` is its mean.)

## (b) The Farey roof

**On each Farey-order-`k` cell `(p/q, p′/q′)`:
`maxgap(AP_k, x) = q·(x − p/q) + q′·(p′/q′ − x)` — linear, interpolating `1/q → 1/q′`.**

*Proof.* On the cell, `m⁺(x) = q(x − p/q)` and `m⁻(x) = q′(p′/q′ − x)` (the classical first-return
facts: the closest approach from above/below at height `k` is realized by the flanking Farey
denominators); add and use (a). ∎ (Machine-verified to 1e-15, k = 3..13, 4000 random x each.)

The roof's node values are `1/q` at each Farey point — **the roof endpoint heights ARE the
Ostrowski/Kravitz rung heights**; the M-spectrum ladder and the μ-floor are two readings of one
piecewise-linear function.

## (c) Exact closed forms (all k)

- **`E[maxgap(AP_k)] = Σ 1/(q·q′²) = (1/2) Σ (q+q′)/(q²q′²)`** over consecutive Farey-k
  denominator pairs. Verified against the canon table: `43/140, 47/168, 19/72, 151/630,
  796/3465, 93/440` (k = 8..13) — all match; new values k = 1..7: `1, 3/4, 7/12, 1/2, 5/12,
  23/60, 1/3`.
- **`E[min_{i≤k} frac(ix)] = (1/2)·E[maxgap(AP_k)]`** (reflection symmetry); e.g. `93/880` at k=13.
- **`μ_θ(AP_k)` = the roof's superlevel measure**, a per-cell linear-tail sum. Reproduces the
  canon exactly: `691/735, 247/294, 38/49, 1381/2205, 13823/24255, 477/1078` (θ = 1/7, k = 8..13),
  and gives e.g. `μ_{2/7}(AP_13) = 829/4620` exactly.
- Mechanism readout: for the AP, `{maxgap > 1/7}` is **exactly the union of the Farey
  neighborhoods of rationals with `q ≤ 6`** (roof node `1/q > 1/7` ⟺ `q ≤ 6`). The AP's density
  floor is pure resonance-window measure; its bulk contribution is zero. Note `q = 7` windows are
  exactly marginal (roof node `= 1/7`, never counted): **the apex prime 7 is invisible to the
  density floor** — its hardness lives entirely on the sup/moat side.

## (d) The anchored-tail finitization (exact, every E)

Let `F_Q = {p/q ∈ [0,1) : q ≤ Q}`. The maximal spacing of `F₇` is `1/7`, so any circular gap of
length `> 1/7` contains a point of `F₇`. Hence for **every** finite integer set `E`:

> **`{x : max_{a∈F₇} gap∋a > 1/7} = {x : maxgap > 1/7}` pointwise (up to the measure-zero
> boundary), so `μ_{1/7}(E) = P_x(max_{a∈F₇} gap∋a > 1/7)` — an 18-anchor statistic.**

(Verified: zero disagreement on 160k-point grids for AP, parity-interlaced, and spread families.)
The load-bearing lemma (A′) `μ_{1/7}(E) ≥ μ_{1/7}(AP_k)` is therefore **exactly** a statement
about the joint law of the config's closest approaches to 18 fixed rationals — inhomogeneous
approximation data at rational targets, governed by E's residue profile mod lcm(2..7)=420 (mod 60
for the F₆ part) and its diameter. THM-591 (inhomogeneous-AP linear law) supplies exact AP-side
control of precisely such shifted targets.

## (e) OPEN candidate (A″-F₆), verified adversarially — tight at the AP

`F₆` (12 anchors) has max spacing `1/6 > 1/7`; escapes are only gaps of length in `(1/7, 1/6)`
hiding inside `(0, 1/6) ∪ (5/6, 1)`. Since for the AP the good set is pure q≤6-window and inside
a window every inter-cluster gap contains its base Farey point (**window-exactness**: for
`x = p/q ± δ` the q clusters spread one-sidedly by `e·δ`, so each inter-cluster gap contains the
flanking `j/q`), we get `t_{F₆}(AP_k) = μ_{1/7}(AP_k)` **exactly**.

> **(A″-F₆):** for every affine-normalized k-set `E` (min = 1, gcd of differences = 1, k = 8..13):
> `t_{F₆}(E) := P_x(max_{a∈F₆} gap∋a > 1/7) ≥ μ_{1/7}(AP_k)`, equality iff `E = AP_k`.
> Since `t_{F₆} ≤ μ_{1/7}` pointwise and `μ` is affine-invariant, (A″-F₆) ⟹ (A′).

Affine normalization is NECESSARY: every naive-anchor failure found was an affine image of the AP
(`{2..14} = AP+1`, odds `= 2AP+1`, `{3+7i} = 7AP+3`) — μ is affine-invariant, anchored tails are
not. (The anchor-side analogue of the kps-S56 dilation-artifact lesson.)

**Evidence:** normalized corpus (14 families incl. parity-interlaced record, GW, spread, primes,
deep well, missing-residue adversaries) — worst is the AP itself; dedicated normalized descent at
k = 13, 12, 10 converges to the AP at the bar; exhaustive 1-swap scan (351 families, k=13; 324,
k=12) and 2-swap samples (468/396): **zero below the bar**. Also zero μ-violations anywhere —
(A′) itself remains unbreached.

**Honest status:** (A″-F₆) is *logically stronger* than (A′); its value is structural — it
replaces the full order statistic by 12 closest-approach statistics with sieve-friendly rational
targets, and its window/bulk decomposition (window mass ↔ residue coverage & diameter;
bulk capture ≈ 1 with escapes confined to two explicit arcs) is a finite-dimensional version of
the additive↔multiplicative tradeoff (mac-mini-S15). Not a proof of (A′); a reformulation (d)
plus a tested extremal candidate (e).
