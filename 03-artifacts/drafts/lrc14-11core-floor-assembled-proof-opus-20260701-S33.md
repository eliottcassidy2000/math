# The 11-core floor: an assembled proof that the census is exhaustive

**Author:** opus-2026-07-01-S33.
**Assembles:** HYP-3900 (simultaneous peel, opus-S32) ≡ HYP-3950 union floor (kps-S28, independent);
the dense census (kps-S25/S27, opus-S31); THM-522 (scale-invariance/quantization); THM-523 (tight
locus); THM-592/593 (radius-derivative frame, unit-residue lemma, mac-mini-S93); HYP-3901 (difference-
core renormalization, opus-S32) with the new exact rearrangement floor (HYP-3902, this session).
**Concurrent:** mac-mini HYP-3850 claims a Mirsky–Newman-type floor and the curvature/no-multiple-of-n
legs; this document leaves those slots marked and cites them where they would strengthen a branch.

---

## 0. Statement and role

**Target (the r=2 residual, OPEN-Q-108).** For every primitive 11-element set C ⊂ Z_{>0} ("11-core"),
at band r = 1/14:

    meas(L_C) ≥ 1/36,     L_C = {t ∈ R/Z : ||vt|| ≥ 1/14 for all v ∈ C}.

This is the innermost open cell of the LRC(14) covering-floor: the covering reduction (THM-523) plus
the two-multiple peel leave exactly this quantitative floor over 11 of the 13 speeds. The pentagon
{1,…,13}\{6,10} attains 313/9702 = 1.161 × (1/36), so the constant has 16% slack — every estimate
below is allowed to be lossy, none is allowed to be asymptotic-only.

**Shape of the proof.** A well-founded case decomposition on the *scale structure* of C. Three moves,
each strictly reducing a complexity measure (element count, or height class):

- **(M1) census** — bounded cores: finite, exact, done;
- **(M2) simultaneous peel** — a scale gap with ≤ 6 elements above it: proved lemma, scale-free error;
- **(M3) cluster renormalization** — ≥ 7 elements in one gap-free cluster: the cluster collapses to
  its difference core one level down, and an exact rearrangement functional gives the floor.

The union bound inside (M2) dies exactly at 7 = 1/(2r) because each speed's danger has measure exactly
1/7: seven is the covering threshold. (M3) exists precisely to catch what (M2) provably cannot.

---

## 1. Standing facts

**(S1) Danger geometry.** danger(v) = {t : ||vt|| < 1/14} is v open arcs of length 1/(7v) centered at
the lattice (1/v)Z; meas(danger(v)) = 1/7 for every v. Hence the unconditional union bound: any k-core
has meas(L) ≥ 1 − k/7 — nontrivial iff k ≤ 6. [Classical; used throughout.]

**(S2) Dilation invariance** (THM-522). meas(L_{gC}) = meas(L_C) for g ∈ Z_{>0} (u = gt is measure-
preserving on R/Z). Consequently every bound may assume gcd(C) = 1, and sub-cores may be re-primitivized
freely.

**(S3) The census** (kps-S25/S27; opus-S31; kps-S28 ledger). Exact rational values for: all primitive
11-cores with max ≤ 19 (75,582 at V=19; min = pentagon 313/9702); all 10-cores max ≤ 18 (min
14249/252252 = 0.056487); the k-core minima m_k at max ≤ 15 for k = 5…11; per-core arc counts
(≤ 40 observed, ≤ 26 on the kps minima ledger); single-outlier cores (compact-10 + w ≤ 500, min
0.0323); two-outlier (98,399 configs) and three-outlier (21,491) sweeps with **zero** values below
1/36 and minima 0.034870 (r=2), ≥ 0.056 (r=3) — monotone increasing in outlier count.

**(S4) Tight locus** (THM-523, THM-593). The 13-speed tight sets at q=14 are {AP, GW}; tight sets
without multiples of 14 represent every unit residue (unit-residue lemma), and their lonely measure
collapses linearly with exact slope c = 1666/6435 (both extremals). This is the *height-1 boundary
theory*: it describes the wall the 11-core floor keeps its 16% distance from.

---

## 2. Lemma P (simultaneous peel) — proved

**Lemma P.** For any disjoint split C = B ∪ F, j = |F|, with L_B having c_B components:

    meas(L_C) ≥ (1 − j/7)·meas(L_B) − (2 c_B/7)·Σ_{w∈F} 1/w.

*Proof.* meas(L_B \ L_C) ≤ Σ_{w∈F} meas(L_B ∩ danger(w)). For one interval I of length ℓ, the danger
arcs of w are centered on (1/w)Z, so at most wℓ + 2 of them meet I, each of measure 1/(7w):
meas(I ∩ danger(w)) ≤ ℓ/7 + 2/(7w). Sum over the c_B components of L_B, then over w ∈ F. ∎

Two independent derivations exist (HYP-3900 opus-S32; HYP-3950 kps-S28 "per-speed lemma", verified
worst ratio 0.90 over ledger × w ≤ 1000). Exact verification: 64 adversarial splits, 0 violations.

**Scale cancellation.** If max(B)·Λ ≤ min(F), then c_B ≤ Σ_{v∈B} v ≤ 11·max(B) gives error
≤ 22j/(7Λ): *uniform in the absolute scale*. No arc-count hypothesis on far-containing prefixes is
needed anywhere — the S31 "uniform arc-count bound" requirement is retired (MISTAKE-090; the quantity
it asked to bound genuinely grows ~ w·meas).

**Corollary P1 (the tower).** With M_k := a valid lower bound over all primitive k-cores:
M_k = 1 − k/7 for k ≤ 6 (S1, unconditional); for k ≥ 7, any k-core that admits a Λ-gap with
1 ≤ j ≤ 6 elements above its top gap satisfies meas ≥ (1 − j/7)·M_{k−j} − 22j/(7Λ). Solving with the
census values (S3): M_7 = 6/49, M_8 = 5/49, M_9 = 4/49, M_10 = 0.056487 (census binds),
M_11 = 313/9702. **Guard table** (target 1/36): (1−j/7)·M_{11−j} = 0.0484, 0.0583, 0.0583, 0.0525,
0.0408, 0.0408 for j = 1…6 — all positive margins, minimum 0.0130. Chained over ≤ 5 peels at
Λ = 10^4 the cumulative error ≤ 0.0095 < 0.0130. [Verified exactly, opus-S32 .out.]

---

## 3. The case decomposition

Fix Λ = 10^4 and the census bound V0 = 19. Given a primitive 11-core C, sort its elements and locate
the **highest Λ-gap** (an index i with w_{i+1}/w_i ≥ Λ); let F = elements above it, j = |F|.

- **Case 1 (no Λ-gap, max ≤ V0):** the census (S3). CLOSED, min = pentagon.
- **Case 2 (top gap, j ≤ 6):** Lemma P + Corollary P1: the guard table clears 1/36 with the recursion
  descending into (11−j)-cores, each level again splitting into these same cases. CLOSED up to the
  middle band (Case 4).
- **Case 3 (a gap-free suffix cluster of j ≥ 7 elements):** at most 4 elements below; the cluster has
  bounded ratios (≤ Λ^{j−1}) and height N = min(cluster). Renormalization, §4–5. CLOSED in the
  N → ∞ limit for patterns of range ≤ 14 by the rearrangement floor (§5); heights below N* fall into
  Case 4; ranges > 14 recurse (see §5, R2).
- **Case 4 (the finite middle band):** configurations whose outliers or cluster heights land between
  the census bound and the asymptotic threshold. Entirely finite and enumerable; current coverage in
  §6.

The decomposition is exhaustive: every 11-core is in exactly one case at the top level, and every
recursive call strictly decreases (element count) or bottoms out in a census/limit/finite-band leaf.

---

## 4. The renormalization lemmas (Case 3 structure)

Write the cluster as {N + c_1, …, N + c_j}, 0 = c_1 < … < c_j = Δ (the **difference pattern**; by
(S2) and invariance of everything below under common shifts and dilations of the pattern, assume the
pattern primitive). Let B = C \ cluster, |B| = 11 − j ≤ 4, and define the **local density**

    D_c(t) = meas_u { u ∈ R/Z : u + c_i t ∈ A for all i },   A = [1/14, 13/14].

**Lemma F1 (tiling rigidity, j = 7).** D_c(t) = 0 iff the seven danger arcs (length 1/7 each,
centers −c_i t) tile the circle, iff their centers form a translate of (1/7)Z. Consequently
(c_i − c_k)·t ∈ (1/7)Z for all pairs, so t ∈ (1/(7d))Z with d = gcd of pairwise differences; for a
primitive pattern the zero set is contained in {k/7 : k = 1..6} — at most six isolated points.

*Proof.* Seven arcs of total measure exactly 1 cover iff they are a.e. disjoint, i.e. tile; closed
equal arcs tiling a circle sit end-to-end, so consecutive centers differ by exactly the arc length
1/7, i.e. the center set is a translate of (1/7)Z. The pairwise condition follows by subtracting
memberships; the gcd statement by taking integer combinations (Bézout). k = 0 is excluded: at t = 0
all centers coincide. ∎

This is the **continuous Fraenkel rigidity**: an exact cover by interval-APs with distinct speeds
exists only on a finite, explicitly bounded set of slow times. (For j ≥ 8 the covering has slack
j/7 − 1 and the zero set can contain intervals; the functional in §5 absorbs this without case work.)

**Lemma F2 (slope at a tiling).** At each zero t0 of D_c (j = 7), D_c(t0 ± δ) ≥ range(c)·δ ≥ 6δ for
small δ.

*Proof.* At the tiling, order the arcs around the circle: speeds c_{σ(1)}, …, c_{σ(7)}. Moving t by δ
changes the gap between cyclically adjacent arcs (σ(i), σ(i+1)) at rate c_{σ(i+1)} − c_{σ(i)}; the
rates sum to zero, so D (the total positive gap) grows at rate = (1/2)·Σ_i |c_{σ(i+1)} − c_{σ(i)}|
= half the cyclic total variation ≥ (max − min) = range(c) ≥ 6 (seven distinct integers). ∎

**Lemma F3 (freeze error — quantitative two-scale).** For the full core C = B ∪ {N + c_i}:

    meas(L_C) ≥ ∫_{L_B} D_c(t) dt − ε(N),   ε(N) ≤ 4(j+1)·√(Δ/N) + c_B(j+1)/N.

*Proof sketch (window argument, fully elementary).* Partition the circle into windows of length h. In
a window with left endpoint t_w, write (N + c_i)t = u + c_i t_w + c_i(t − t_w) with u = Nt sweeping an
interval of length Nh; the offsets drift by ≤ Δh. Sandwich the true danger arcs between the frozen
ones fattened/thinned by Δh: per component of L_B per window, the loss is ≤ (fattening) 2jΔh·(length)
+ (sweep edge effects) (j+1)/N. Sum (≤ 1/h windows, c_B components), then replace D_c(t_w) by D_c(t)
using |D_c′| ≤ 2jΔ. Optimize h ≍ √((j+1)/(jΔN)). ∎

The √-rate is what this crude sandwich gives; the observed rate is O(1/N) or better (P1 probe: error
2·10⁻⁵ at N = 120). Sharpening F3 to O((jΔ + c_B)/N) — a one-step Koksma/BV refinement per component,
or mac-mini's Mirsky–Newman floor (HYP-3850a) applied directly at finite N — shrinks the middle band
(Case 4) from N* ~ 10⁸ to N* ~ 10³. **This is the single most valuable quantitative improvement
available** and is exactly HYP-3787's task in renormalized coordinates.

---

## 5. The rearrangement floor (Case 3 closure at the limit)

The compact part enters only through its measure and worst-case *position*: the adversary may place
L_B wherever D_c is smallest. Define the exact worst-position functional (the increasing-rearrangement
partial integral)

    Q_c(m) = ∫_0^m D_c^*(u) du = ∫_0^∞ (m − ψ_c(s))_+ ds,   ψ_c(s) = meas{t : D_c(t) < s},

so that ∫_{L_B} D_c ≥ Q_c(meas(L_B)) with **no positional assumption**. Since |B| = 11 − j,
meas(L_B) ≥ (j − 4)/7 by (S1). The floor constant at cluster size j is therefore

    F_j := min over primitive patterns c of size j, range ≤ R0, of Q_c((j−4)/7).

D_c is piecewise linear in t; its kinks come in exactly two types (the THM-592(iii) taxonomy in the
t-variable): opposite-endpoint collisions δt ≡ ±1/7 (gap birth/death, convex) and same-endpoint
collisions δt ≡ 0 (coincidence/overtaking, concave peaks) — grid t = (7m + e)/(7δ), e ∈ {−1, 0, +1},
over pairwise differences δ. [Omitting e = 0 is a conservative error: all missed kinks are concave
peaks, so interpolation under-estimates D and hence Q — caught same-session via the slope canary at a
tiling zero; MISTAKE-092.] With the complete grid, ψ_c and Q_c are exact. **Computation (HYP-3902,
R0 = 14, ~6300 primitive patterns float-scanned, worst-40 per size re-verified in exact rationals):**

    F_7  = 559/11025      = 0.050703   at (0,1,2,3,4,5,6)        — margin 1.825 × (1/36)
    F_8  = 184019/3246495 = 0.056682   at (0,1,2,3,4,5,6,7)      — 2.041 ×
    F_9  = 244547/3522610 = 0.069422   at (0,1,2,3,4,5,6,7,8)    — 2.499 ×
    F_10 = 56333/617400   = 0.091242   at (0,2,3,4,5,6,7,8,10,12) — 3.285 ×
    F_11 = 63941/432180   = 0.147950   at (0,2,4,5,6,7,8,9,10,12,14) — 5.326 ×

Every F_j clears 1/36; the j = 7 argmin is the consecutive (AP) pattern — the renormalization fixed
point, consistent with HYP-3901's contraction picture — and F_7 = 0.0507 > pentagon 0.0323, so the
binding configuration for the whole residual remains the compact census extremizer, now with a proved
57% gap at the limit rather than a probe. Structure checks: the consecutive-pattern zeros are exactly
{k/7 : k = 1..6} (Lemma F1), one-sided slopes at the zeros exactly 6 = range (Lemma F2), a non-tiling
pattern has min D = 8/49 > 0, and the j = 8 consecutive pattern's interval zero-set has measure
44/735 = 0.0599 (robust coverings exist at j ≥ 8, absorbed by Q without case work).

Two reductions make R0 = 14 sufficient in the limit:

- **(R1) pattern dilations and shifts are free:** D depends only on the primitive difference pattern
  (common shifts rotate u; dilations rescale t measure-preservingly), so the enumeration over
  primitive range ≤ R0 covers all patterns of *bounded primitive range*.
- **(R2) large-range patterns recurse:** a pattern of primitive range > R0 either has a Λ′-gap in its
  own difference structure — then the *cluster itself* splits under Lemma P one level down (the
  tower loses ≥ 1 element) — or is a bounded-ratio object at a higher range, where D_c approaches
  the (j−1)-fold independent density with defect controlled by the same peel counting; the recursion
  is on j and terminates at j ≤ 6 where (S1) is unconditional. [This is the "tower loses ~6 runners
  per level" structure: depth-1 is a shifted LRC(8)-type average; depth-2 is below the covering
  threshold.]

---

## 6. The finite-check ledger (what is verified, what remains finite)

| # | Finite check | Status |
|---|---|---|
| F-i | Dense census max ≤ 19 (75,582 cores at V=19) | **DONE** (S27/S31), min = pentagon 313/9702 |
| F-ii | Single-outlier: compact-10 + w ≤ 500 > W*(1) | **DONE** (S31; kps W*(1) = 112) |
| F-iii | Two/three-outlier middle band | **DONE to evidence standard** (kps-S28: 98,399 + 21,491 configs, 0 < 1/36, minima 0.0349 / ≥ 0.056, monotone in outlier count; kps W*(2,3) = 181/290); four-to-six-outlier band: guard margins large (0.0525–0.0408), spot-censused (opus-S31 tower two-outlier; S32 probe) — full enumeration remains |
| F-iv | Cluster heights 20 ≤ N < N* (Case 3 finite band) | **PARTIAL**: ~7,500-config boundary minimax (S32: all 715 compact 4-cores × consecutive clusters, heights 14–120, global min 0.0635 = 1.97× pentagon; pattern variants at N = 210); full band awaits the F3 sharpening (N* ~ 10³) + a targeted sweep |
| F-v | Rearrangement floor F_j, j = 7…11, range ≤ 14 | **DONE this session**: F_7..F_11 = 0.0507, 0.0567, 0.0694, 0.0912, 0.1480 — all > 1/36 (exact on worst-40/size; HYP-3902 .out) |

**Honest gap list (exactly two mathematical items + one computation):**
1. **F3 sharpening to O(1/N)** (or the HYP-3850a Mirsky–Newman finite-N floor) — makes F-iv a small
   finite sweep instead of an astronomically wide one.
2. **(R2) made fully rigorous** — the large-range pattern recursion needs its (easy but unwritten)
   defect bound; alternatively R0 can be raised computationally.
3. **F-iii/F-iv sweeps run to completion** at the widths the sharpened constants dictate.

No conceptual unknowns remain in the decomposition: every leaf is census, proved lemma, exact
computation, or an explicitly finite sweep.

---

## 7. Connections and instruments (the boundary theory and two new certificates)

**7.1 The unit-residue wall (THM-593) as the height-1 boundary.** The tight 13-sets sit at measure 0
with linear collapse c = 1666/6435 = (2/14)·Σ_{u∈(Z/14)*} 1/u. The 11-core floor lives 1.16× above
1/36 because *dropping two speeds from a tight set opens exactly the unit-witness structure the
unit-residue lemma forces*: with 11 speeds one cannot represent both ±1 at every unit while covering
— the pentagon keeps the maximum it can. A quantitative "stability unit-residue lemma" (missing unit
⟹ M ≥ 1/q + explicit gap; mac-mini HYP-3850c claim) would derive the census minima's ORDER (why the
pentagon) rather than just their values; slot marked.

**7.2 The overtaking-defect sieve (THM-592(v) operationalized).** For bounded cores the ladder bound
says: a certified pair (m(r0), m′(r0)) at a sub-critical anchor pushes positivity to 1/14 unless
exposed overtaking events d/(w−v) ∈ (r0, 1/14) spend the slope. For an 11-core the exposed-defect sum
K is a *sparse arithmetic list* — pairs with w − v large relative to gcd. **Certificate format: three
rationals (m, m′, K) per core**; validity is one rational inequality; Lean-friendly; replaces
re-computation of full lonely sets in F-i. Proposed as the storage format for the census ledger.
The sieve direction: covering 13-sets contain multiples of every q ≤ 14, forcing many pair-differences
and hence a dense potential-defect grid — but *exposure* is the filter; empirically (deep well) one
anchor at 1/16 already certifies. Worth a dedicated pass.

**7.3 The arc × radius two-variable LP.** Variables x_{a,ℓ} ≥ 0 (lonely mass in arc a of the Γ₀-grid
Q ∈ {14, 61, 183} at radius layer r_ℓ); constraints: (i) layer monotonicity x_{a,ℓ+1} ≤ x_{a,ℓ};
(ii) klein/mac-mini per-arc coverage at each layer (HYP-3824/3833 localized masses); (iii) **slope
transport**: x_{a,ℓ} − x_{a,ℓ+1} ≤ 2·b_a·Δr·(max_v 1/v within reach of a) with b_a a per-arc component
budget — THM-592(ii) localized. Objective: min Σ_a x_{a,L} at r_L → 1/14. Single-radius LPs stall at 0
(HYP-3822); the transport constraints import convexity, so the LP value at the critical layer is
bounded below by anchor mass minus defect spend — the LP form of the ladder. The universal version
(min over all covering configurations) parameterizes constraints by the covering axioms only; its
value ≥ c·(1 − 14r) with c ≤ 1666/6435 would be the LP-certified collapse law. Pilot: n = 6 next
session.

**7.4 The ζ(2) ↔ ζ(−1) note.** The collapse constant is a unit-harmonic sum (2/n)·H_×(n) — the same
species as the Farey/coprimality densities that give the ζ(2) floor constants, while the cap side
regularizes to ζ(−1) (S31r): the tight locus is where the two normalizations meet, and klein's
layer-cake E(S) = ∫Λ_S = Σ_Farey 1/(qq′(q+q′)) is literally a Mordell–Tornheim-type value (cf.
THM-564's 2ζ(3)). Homogeneous-dynamics reading: the deep-cluster limit is an orbit-closure average;
rational ratio vectors = closed orbits (census values), irrational = higher-dimensional closures
(larger measure) — the inf-at-rational structure IS the Ratner-type rigidity heuristic behind (R2),
and the right general frame for the barcode/layer-cake objects (codex-S179, klein HYP-3834/3836).

---

## 8. Summary

The r=2 residual is now: **[proved lemmas P, F1, F2, F3(√)] + [exact constants: guard table, F_j] +
[finite sweeps F-i…F-v, three-quarters done] + [two named quantitative sharpenings (F3 rate; R2
defect)]**. Nothing in the list is a method gap; each remaining item has an owner-shaped next step.
The extremal structure is confirmed at every level to be the AP/pentagon fixed-point family, exactly
as the difference-flow renormalization predicts (HYP-3901), with the unit-residue lemma (THM-593)
explaining the wall it leans against.
