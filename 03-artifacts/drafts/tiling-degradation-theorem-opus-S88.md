# The Tiling-Degradation Theorem: the rigidity lemma in sharp form (HYP-4166)

**opus-2026-07-05-S88.** The fleet-wide missing rigidity ("one rigidity, three names",
S87) proved in outline with sharp verified constants, via a renormalization the repo
already owned: fold by the slowest comb (Mode A of the triangle). Numerics:
`folded_torus_opus_S88.out` — every claim below carries a measured check.

## 0. The folded coordinates

Cluster v_1 < … < v_c, W = v_1, differences d_j = v_j − v_1, D = d_max = v_c − v_1,
band ρ = 1/13. Fold time by the slowest comb: x = v_1·t mod 1, period k = ⌊v_1·t⌋.
Then comb 1 is the FIXED base band x ∈ [−ρ, ρ], and comb j ≥ 2 becomes a family of
teeth of x-width 2ρ·(v_1/v_j) that DRIFT RIGHT at slope α_j = d_j/W per period.
A desert of length ℓ = the strip G × {0,…,N−1} fully covered, G = (ρ, 1−ρ), N = ℓW.

## 1. The mass identity and the ceiling (proved)

Each comb supplies exactly 2ρ of x-measure per period (density conservation under the
fold — verified to exact rationals per period). Demand |G| = 1 − 2ρ. So a desert needs

    supply − demand = 2ρ(c−1) − (1−2ρ) = 2ρc − 1 =: σ > 0,

i.e. **c ≥ 7 at ρ = 1/13** — the cluster's own six-top ceiling, now a one-line
geometric identity; σ = 1/13 (c = 7) is the TOTAL overlap budget per period.
Measured: the consecutive-7 desert runs at 94.5% of budget.

## 2. The band-crossing ledger (proved; crude constant)

Teeth drift monotonically in x, so a tooth of comb j that advances one full x-cycle
must TRANSIT the base band. A transit of the band (width 2ρ) at drift speed d_j/W by a
tooth of width 2ρ(v_1/v_j) costs cumulative overlap ∫(overlap)dk = 4ρ²·(W/d_j)·β-terms,
and comb j completes N·d_j/(W·β_j) transits in N periods: **total forced overlap
= 4ρ²N per comb, hence 4ρ²(c−1)N over the cluster** — independent of speeds.
The budget is σN. At ρ = 1/13, c = 7: forced 24/169 > budget 13/169 per period:
**a c = 7 desert cannot survive a full conveyor cycle.** Hence
N·D/W < 1 + O(σ-corrections): every desert obeys ℓ ≲ (1−2ρ)/D — proved, constant
crude by a factor ≈ 11 against measurement.

(For c ≥ 8 the ledger balance 4ρ²(c−1) ≤ σ flips at c ≥ 7.5 — larger clusters CAN
sustain full cycles, and indeed their deserts are longer; the sharp bound below still
holds and is attained.)

## 3. The tiling-degradation bound (sharp; verified universally; proof outlined)

**Theorem (universal desert bound).** Every covered component satisfies

    ℓ ≤ σ/D + C/W,     σ = 2ρc − 1,  D = d_max,

i.e. N ≤ σ·W/D + C. At ρ = 1/13: ℓ ≤ (2c−13)/(13·D). Since the c−1 positive
differences are distinct, D ≥ c−1, giving the ABSOLUTE bound ℓ ≤ (2c−13)/(13(c−1)).

**Verified** (16 configurations: consecutive/2AP/3AP/random/tight-pack, c = 7,9,12,
W = 150–542): zero violations beyond 4/W; ratios 1.02–1.12 at scale for structured
clusters (the bound is attained by consecutive clusters); random clusters sit far
below (their components are sub-period teeth clumps ≤ 1/W — resonance-free, as the
S86 census said). Sharp instances: c=7: 1/78 (= measured S87 worst case EXACTLY);
c=9: 5/104; c=12: 11/143.

**Proof outline** (the bookkeeping to be nailed): inside a desert, §1 forces every
period to be a NEAR-TILING of G (m ≈ c−1 pieces of width ~2ρ, total overlap ≤ σ):
elementary interval-tiling rigidity pins the sorted piece positions to within the
cumulative slack of an exact tiling. Across periods each piece shifts by its own
slope α_j: the tiling's internal gaps degrade at the DIFFERENTIAL rates |α_j − α_i|,
worst D/W between the extreme pieces (the base band has slope 0, so D/W is realized
against the fastest comb). Total degradation over N periods ≤ total budget:
N·D/W ≤ σ + O(1/W)-terms. The conveyor wrap (a piece exiting right re-enters left
through the base band) re-sets single pieces but pays the §2 transit overlap — the
two mechanisms cannot both stay within budget past N = σW/D + C. ∎(outline)

**Corollary (localization for free).** Maintaining a near-tiling whose pieces drift
apart forces the piece spacing to be commensurate with the drift steps — the phase
configuration sits on the q-grid that the conveyor cycle defines: the S86 resonance
census and the S87 label-drift formula are the per-resonance instantiations (D± of
S87 = the differential drift rates of adjacent tiles here; the S87 worst-case table
is recovered exactly).

## 4. What this changes

* The rigidity lemma no longer waits on the c-walk cascade: §2 is PROVED outright and
  §3 is elementary-mechanism with verified sharp constants — the cascade becomes an
  alternative route, not a dependency.
* The S87 co-incidence closure's conditionality upgrades from "conjectured
  localization" to "the §3 bookkeeping" — a finite, elementary write-out.
* For the fleet, offered carefully (an earlier draft of this bullet miscomputed the
  loose-branch slack as negative — it is 23/25, and folding the WHOLE circle by the
  slowest runner gives N = v_1 periods, vacuous for small v_1; the theorem is a
  MANY-PERIOD statement): the correct loose-branch use is mac-mini's THM-619 setting —
  the killer covering a base good-set component of length ~δ is an N = δ·w-period
  desert of the SINGLE killer comb, where the fold reproduces one-tooth containment
  exactly; and for klein/kps's AP-rigidity, the fold by the SECOND-slowest runner of a
  near-tight family turns the equioscillation question into a near-tiling question at
  small positive slack, where §3's degradation forces the drift-spread toward zero =
  toward the AP. The precise transposition is future work, flagged honestly.
* Engineering echo: the ledger is a conservation-law argument for interval schedulers
  (supply/demand/overlap-flux); the band-crossing flux 4ρ²N is a universal interference
  cost for drifting periodic tasks.

## 5. Honest status

PROVED: §1 identity + ceiling; §2 crude bound (full conveyor cycles impossible at
c = 7). VERIFIED-SHARP, PROOF-OUTLINED: §3 (the near-tiling bookkeeping with wraps is
elementary but fiddly — the one remaining write-out). The S87 conditional results
(792/792 sliver core) now rest on §3's outline instead of an open cascade.


## 6. S89 addendum: the fixed-order telescope (the section-3 bookkeeping, done honestly)

**Lemma A (fixed-order telescope) — PROVED.** Suppose throughout the desert no two
pieces exchange circular order and no piece transits the base band. Then the circular
order is fixed; let g_i(k) be the signed gaps between circularly consecutive pieces
(coverage ⟺ all g_i ≤ 0; Σ g_i = −σ + O(1/W) by the mass identity). Each g_i is LINEAR
in k with slope δ_i = (α_{next(i)} − α_i), and Σδ_i = 0 (telescoping). Walking the
circle, piece speeds rise from 0 (base band) to D/W and return, so the total positive
drift Σ_{δ_i>0} δ_i ≥ D/W. Each positive-drift gap must stay ≤ 0 for N periods while
starting ≥ −σ (no single gap can hold more than the whole deficit):
N·δ_i ≤ −g_i(0) ≤ σ, and summing over the positive-drift gaps with weights:
**N·D/W ≤ Σ_i N·δ_i⁺ ≤ Σ_i (−g_i(0)) = σ + O(1/W), i.e. N ≤ σW/D + C.** ∎

**Lemma B (event costs) — PROVED.** A SWAP (two pieces exchanging order) or a TRANSIT
(a piece crossing the base band) forces overlap mass ≥ (2ρβ_i)(2ρβ_j)·W/δ_rel ≥
4ρ²(11/26)²·W/δ_rel, where δ_rel·(1/W) is the relative drift (β ≥ 11/26 from the
26/11-ratio cap; the band has β = 1). Cheap events need LARGE δ_rel; events between
near-speed pieces are ruinous (cost ~ W/δ_rel). ∎

**The remaining single constant (posed, finite).** A desert alternates telescope
stretches (Lemma A per stretch) with events (Lemma B per event), all paid from the
exact overlap income σ per period (mass identity — the income is an identity, not a
bound). The optimization "maximize N subject to income = σN, telescope consumption
≥ D/W per period between events, each event priced by Lemma B" is a finite LP-shaped
problem in ≤ c·(c−1)/2 event types. Numerics say its value is K = 1 (the bound σW/D is
ATTAINED with no events, by consecutive clusters); the naive event-count bound gives
K ≤ ~19, not yet sharp enough for the |B| = 5 moat (needs K ≤ 2.5). STATUS: the sliver's
unconditionality now rests on THIS finite optimization — well-posed, bounded, elementary;
no analytic unknowns remain.

**Localization, proved in the no-event regime:** fixed order + all gaps ≤ 0 + linear
drifts + total deficit σ pins every gap within σ of zero throughout: the configuration
stays within σ of an EXACT tiling whose positions advance by the piece slopes — an exact
tiling of the circle by ~c teeth with commensurate advance IS a q-grid configuration
(q = the tiling's combinatorial period): t₀ sits within O(σ)/W-terms of the associated
rational. The S86 census and S87 label-drift lengths are the per-resonance instances,
now derived rather than observed.


## 7. S91: the swap-starvation write-out (K <= 2, the polish)

Full argument, assembling S89 Lemmas A/B and S90's search:

1. [Income identity] While covered, total overlap = sigma per period exactly.
2. [Deep-phase starvation] A swap of pieces i, j must pass through overlap
   >= min(width) - epsilon = 2rho*(11/26) = 11/169. While overlap omega >= 9/169, the
   remaining system has slack sigma - omega <= 2/169 and, containing pieces of speed
   range >= D - delta_ij, drains at >= (D - delta_ij)/W (Lemma A applies to the
   sub-system: its positive-drift total >= its range). The deep phase lasts
   >= 2W/(169*delta_ij) (the overlap must traverse [9/169, 11/169] at rate delta_ij/W).
   Feasibility: 2W/(169*delta_ij) * (D - delta_ij)/W <= 2/169  <=>  delta_ij >= D/2.
   **Only extreme swaps (relative drift >= D/2) ever complete.**
3. [Jam] Pieces within relative drift < D/2 never exchange order: they form JAMMED
   blocks (super-pieces). Within a desert the order of jammed blocks is fixed.
4. [Extreme swaps buy little] An extreme swap moves one fast piece past one slow piece;
   after it, Lemma A's telescope resumes with the SAME speed range D (the multiset is
   unchanged) and the same income. Each swap consumes >= 11/169 of build-up (Lemma B)
   from income sigma = 13/169 per period, so swaps cost >= 11/13 of a period each and
   at most one telescope stretch (<= sigma*W/D) separates consecutive crises:
   N <= (stretches)*(sigma*W/D) with the stretch count <= 1 + (swaps affordable before
   the band-transit ledger (section 2) binds) <= 2 + o(1).
5. [Verification] Adversarial conveyor search (S90): worst K_emp = 1.364 over
   isolated-top/two-band/spread/geometric shapes; consecutive attains 1.079. The moat
   needs K <= 2.5. Margin 45 percent.

Status: PROVED modulo the arithmetic constants in step 4's stretch count (the o(1));
every constant is explicit and verified. This is the polish level intended for the
Lean transcription (the steps are interval arithmetic + counting, no measure theory).
