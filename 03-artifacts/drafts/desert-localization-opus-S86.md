# Desert Localization: the exact mathematics of the cluster sliver (HYP-4146)

**opus-2026-07-05-S86.** The one open l ≥ 7 case (same-scale clusters, c ≥ 7, fee-dead
by the six-top ceiling, descent-dead by ratio) reduced to: one provable theorem + one
finite rational-incidence table. Exact-math-first per the owner directive; numerics in
`cluster_deserts_opus_S86.out`.

## 1. The phase decomposition

Cluster C = {v_1 < ... < v_c}, scale W = v_1, ratio v_c/v_1 < 26/11, band ρ = 1/13.
On a short interval around t₀, each phase splits as

    v_j·t = [v_1·t] + [d_j·t₀] + d_j·(t − t₀),      d_j = v_j − v_1 ≤ (15/11)W,

a FAST common rotation (rate W), FROZEN offsets {d_j t₀ mod 1}, and a SLOW drift
(rate ≤ 1.36W but applied to |t − t₀| ≤ desert length). The union of the c combs covers
a neighborhood of t₀ iff for every position θ of the common rotation some offset lands
in the band: **iff the offset set {0, d_2 t₀, ..., d_c t₀} is a 2/13-net of the circle.**

## 2. The resonance census (exact combinatorics — proved)

If the offsets form a 2/13-net and t₀ ≈ p/q, the offsets sit near the q-grid restricted
to the residues {d_j mod q}, whose max circular gap g must satisfy g/q ≤ 2/13:

* q ≤ 6: max gap ≥ 1/q ≥ 1/6 > 2/13 — **immune**: no desert can center at low-q rationals;
* q = 7: g = 1 forced — needs ALL 7 residues hit by {0, d_j mod 7}: c ≥ 7 and the
  difference multiset covering ℤ/7;
* 8 ≤ q ≤ 12: g = 1 forced (2q/13 < 2) — needs all q residues: **c ≥ q**;
* q = 13: g ≤ 2 gives slack 2/13 − 2/13 = **ZERO** — marginal nets have no room to
  persist, and g = 1 needs all 13 residues (c ≥ 13 > 12): **q = 13 is slack-dead** for
  every admissible cluster;
* q ≥ 14: slack > 0 requires g < 2q/13 with g ≥ ⌈q/c⌉: for c = 7 only q ≡ 0 (mod 7)
  (slack exactly 1/91, same as q = 7); measured: even these did NOT materialize beyond
  q = 7 itself (the common-phase interplay adds constraints) — the per-cluster resonance
  set is SPARSER than the census bound.

**Sharpened consequence (the incidence correction):** a naive all-q danger table covers
the circle (0/179 base points clash-free with worst-case radii) — but PER CLUSTER the
danger zone is tiny: a c = 7 cluster resonates only at q ∈ 7ℤ (measured: q = 7 alone,
two deserts, total measure 0.028); q ∈ [8,12] needs c ≥ q; q = 13 never. The correct
statement: danger(C) = the localization theorem's finite desert list for C, of total
measure ~ 0.03-0.16 (c = 7..12), against ≥ 50 spread base margin points (merge grid,
denominators ≤ 24). Generic configurations trivially admit a clean base point; the
residual is the THIN CO-INCIDENCE VARIETY (base's few margin points all inside C's few
deserts simultaneously) — finitely parametrized (resonant residue pattern × scale tuning
× base shape), the right target for a box sweep or the S85 exact machinery.

Verified: consecutive c=7 (differences {1..6}: covers ℤ/7) deserts EXACTLY at 1/7 and
6/7, length 0.0138 each, tail total 0.028; random c=7: **zero** components > 5/W;
all-differences-≡ 0 (mod 7): gap 1.0, no desert (concentration = alignment = harmless —
the dangerous degeneracy is residue-COVERING, not residue-concentration).

## 3. The localization theorem (the one lemma to prove)

**Theorem (desert localization — conjectured form, numerics 100%):** every covered
component of length > K/W (K ≈ 5) is contained in an interval around a q-resonant
rational p/q (census above), of half-length at most

    slack(q, C) / d_rms-ish,   slack = 2/13 − g(q, C)/q,

i.e. the desert persists only while the offset drift has not eaten the net's slack.
Proof route: net-drift rigidity — over a desert of length ℓ every offset moves d_j·ℓ;
a 2/13-net of the circle by c ≤ 12 points cannot survive total displacement beyond its
slack; kps's pair walk (S8) and triple walk (S9) are exactly the c = 2, 3 cases of this
rigidity (alternation/exposure arguments), and the c-walk cascade they offered is its
general engine. The band-parametric constants transpose (their 2/25 arithmetic → 1/13:
pair balance 11/6, seam to descent (11/6, 26/11) closed by fee-blocks below).

## 4. What survives for the assembly

1. **Generic clusters** (no netable q): all components ≤ K/W — the window works at ANY
   base point once W ≥ K·468/7 ≈ 67K ≈ 350 (K = 5); below that, bounded box.
2. **Resonant clusters**: deserts at KNOWN rationals p/q, q ∈ [7,13], with EXPLICIT
   lengths. The base window sits at citation/merge-grid points with denominators ≤ 24.
   The clash question is a FIXED FINITE INCIDENCE TABLE: (base point, resonance point)
   pairs with |t_base − p/q| < desert half-length + window. Computable once, exactly,
   for all of: ~50 base points × ~46 resonance rationals. Where no clash: done. Where
   clash: the cluster is heavily structured (residue-covering mod q AND scale-tuned):
   near-consecutive — and the S85-concurrent exact machinery (M = a/(2a+r−1) saturation,
   the full-consecutive family is loose OUTRIGHT) plus residue-pinning eat exactly these.
3. **c ≥ q ≥ 8 resonances need c ≥ 8**: at l ≥ 7 with ≤ 5 unlifted, a c ≥ 8 cluster
   leaves ≤ 4 runners outside — the family is nearly ALL one cluster: the S85 exact
   forms + near-AP perturbation own this corner.

## 5. Honest status

* Proved here: the phase decomposition; the resonance census combinatorics; immunity of
  q ≤ 6; the two-sided numerics (localization 100%, lengths, tails).
* To prove: net-drift rigidity (the c-walk cascade at 1/13) — REAL mathematics, kps's
  cascade is the tool; and the incidence table computation (mechanical, exact).
* Then the sliver is CLOSED: generic by short-deserts, resonant by incidence/structure,
  and the l ≥ 7 unbounded leg assembles: composition descent (S84) + fee-blocks (c ≤ 6)
  + desert localization (c ≥ 7) + bounded boxes (mac-mini sweeps).
