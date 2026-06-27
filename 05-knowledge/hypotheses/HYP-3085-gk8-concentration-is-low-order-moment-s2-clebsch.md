---
id: HYP-3085
title: The gK8/Delsarte concentration bound (CRUX 1 of the covering route) is a LOW-ORDER moment-LP statement led by the pairwise S2 term, carried by the Clebsch "support-six" biplane (= cut-space Cayley graph of K5); and this BOUND — not the AP/GW census — is exactly what the proof needs
status: VERIFIED localization (moment-driver decomposition + S2 extremality, this session) + SYNTHESIS reduction (literal design-Hodge route REFUTED, replaced by reflection-Perron route). Not a proof of CRUX 1.
source: mac-mini-2026-06-27-S60
note: renamed from HYP-3084 (ceded to kind-pasteur-S31af, the level-7 sieve / margin refutation) to avoid namespace collision.
related:
  - HYP-3084   # kps level-7 sieve + the dilation/margin correction (integrated in "the sharp-bound correction" below)
  - HYP-2829   # gK8 binding is single-far (the r=0/r=1/r>=2 split) — this S2-localization sharpens it
  - HYP-2823   # variance reframe (Var(N)=S1+2S2-S1^2) — S2 is the covariance layer
  - HYP-2809   # the dichotomy / census route — NOT needed (off critical path)
  - HYP-2892   # Clebsch / Bruhat / unital exact finite carriers ("support-six")
  - HYP-2890   # support-six residual leak
  - THM-534    # coverage sieve = moment-LP (the working frame)
  - THM-563    # apex-periodicity (single-far supremum is a finite window-max)
  - OPEN-Q-108 # the covering-moment / bounded-core positivity crux
reflections:
  - the-proof-lives-at-the-exponential-variation-covering-bound-not-census   # S59 redirect
  - the-cut-side-is-classical-clebsch-and-the-permutohedron                  # Clebsch = cut-space K5 (S40)
  - multiplication-is-repeated-addition-the-lrc-hyperoperation-tower-s548    # the tower
---

# HYP-3085 — The covering-moment bound is a low-order (S2-led) moment-LP on the Clebsch carrier

## The reconciliation first (what the proof actually needs)
The S59 redirect (the proof = the covering bound; the AP/GW census = the easy case's extremal,
OFF the critical path) is sometimes misread as "the p0/cap work is irrelevant." It is NOT. The
covering route's load-bearing open node is the **bounded-core positivity** — O2 / GAP-G2 /
OPEN-Q-108 / "covering-moment discharge": every bounded covering core has a positive-margin safe
set, i.e. `p0(E) <= cap_k` with margin. That is a **BOUND**, and it is exactly what the gK8 /
Delsarte machinery delivers (`10*p0(E) <= L_yK8(E)`, proved sorry-free; `max_E L_yK8(E) <= 10*cap_k`,
CRUX 1). What is NOT needed is the **CLASSIFICATION** of the tight locus (the census AP/GW =
HYP-2809 dichotomy). **Proof needs the bound, not the classification.** So HYP-2829 (gK8 binding is
single-far, comfortable margin) sits directly on the critical path; the census does not.

## The sharp-bound correction (integrating kps HYP-3084 — my S59 "margin" framing was WRONG)
kps-S31af refuted the tempting reading that the covering bound is "strictly weaker / has free margin":
`2·{1..13} = {2,…,26}` **is covering** (contains `14`) and **tight** (`M = 1/14`, by scale-invariance
`M(cS)=M(S)`), so the covering bound `M ≥ 1/14` is **achieved** — no slack. The resolution sharpens,
rather than weakens, this file: the **dilation `×2` that moves the AP into the covering case is exactly
the H2 multiplicative face** (`14 = 2·7`; multiply by the 2-part and `7 ∈ AP` becomes the `14`). So the
census tight configs, *in dilated form*, ARE the tight covering sets — the equality locus the bound runs
up against. This is consistent with the gK8 route precisely because **`cap_k` is the SHARP extremal value**
(`p0 = cap_k` is attained at the tight configs): a sharp bound is the right tool for a tight inequality.
What stays true from S59: the **classification** (which configs are tight) is not logically required —
only the sharp inequality `p0 ≤ cap_k`, achieved with equality at the dilated census. So: *bound (sharp),
not classification*; and the census is the equality locus of that sharp bound, reached through the H2
dilation — not an irrelevant easy-case artifact.

## The new localization (VERIFIED this session)
The Delsarte duals, in factorial-moment form (`S_r = E[C(N,r)] = sum_t C(t,r) q_t`, the r-th binomial
moment of the inner-sector miss-count `N = #{j in 1..6 : sector j empty}`):
```
L_yK8  = 10 S0 - 10 S1 + 10 S2 -  9 S3 + 6 S4     (k=8)
L_yK9  = 18 S0 - 13 S1 +  8 S2 -  3 S3            (k=9,10)
L_yK11 =  6 S0 -  3 S1 +    S2                    (k=11,12,13)
```
**All binding duals are supported on `S0..S4` — a LOW-ORDER moment functional.** And the
concentration GAP (the binding consec value minus a wide value) is driven by the **pairwise** term
`+S2` (`04-computation/lrc_gk8_moment_decomposition_macmini_S60.py`):

| k | dual driver of the consec−wide gap | per-term ΔS·c (S0..S4) vs single-far(21) | S2 extremal? |
|---|---|---|---|
| 9 | **S2** (+8) | `[0, -0.39, +2.17, -1.18, 0]` | consec maximizes S2 (0/42 beat it) |
| 10 | **S2** (+8) | `[0, -0.22, +2.07, -1.09, 0]` | consec maximizes S2 (0/42) |
| 11 | **S2** (+1) | `[0, -0.08, +0.25, 0, 0]` | — |
| 8 | S2/S3/S4 balance | `[0, +0.49, +2.75, -3.88, +1.85]` | near (1/42 beat it) |

Reading: `S2 = sum over the 15 = C(6,2) inner-sector PAIRS of meas{both empty}` — the pairwise
co-emptiness (covariance) layer of the miss-count. At `k=9,10,11` the binding is a clean **S2
maximization** (consec = maximal positive pairwise sector-correlation = the three-gap/additive
extremal). At `k=8` (the tightest cap, the perennial "finite exception") the cubic/quartic
corrections `-9 S3 + 6 S4` also bite — `S2` up, `S3` down, `S4` up — so k=8 is a low-order *balance*,
not a pure S2 max. Either way: **the bound lives at moment order <= 4, with S2 leading.**

## The Clebsch carrier (the H4 / cut-space face)
The pairwise object is exactly the project's **Clebsch "support-six" carrier**:
- `15 = C(6,2)` inner-sector pairs `= 2^4 - 1` = the nonzero vectors of `(Z/2)^4`.
- `(Z/2)^4` = the **cut space of K5** = the **Clebsch graph** (folded 5-cube, SRG(16,5,0,2)),
  VERIFIED as the cut-space Cayley graph of K5 (reflection `the-cut-side-is-classical-clebsch-and-the-permutohedron`, S40).
- The Clebsch closed-neighborhoods form the **2-(16,6,2) biplane** (block size 6 = "support-six"),
  the pair-balanced design carrier in HYP-2890/HYP-2892.

So CRUX 1 (the covering-moment bound) is a statement about the **pairwise covariance of 6 sector
indicators**, and its natural certificate is the **design-Hodge decomposition of that 15-dim
pairwise space on the Clebsch cut-space geometry**. This is the H4 (2-adic / cut-side) face of the
hyperoperation tower (S548), and the binding config (consec) is the H1 (additive / three-gap)
extremal sitting on it — a cross-level weld, the shape S548 predicts for LRC.

## The proposed reduction (the improved argument for CRUX 1)
```
max_E L_y(E) <= scale_k * cap_k         (CRUX 1, the covering-moment bound)
  <==  (i)  it is a moment functional of order <= 4 (verified: duals on S0..S4);
       (ii) the leading binding term +S2 is the pairwise co-emptiness, maximized at consec
            (verified k=9,10; near at k=8) — an S2-extremal certified by the Clebsch/biplane
            design-Hodge decomposition (pairwise covariance has the 4I+2J / SRG eigenstructure);
       (iii) the bounded (r=0) finite check + single-far (r=1) THM-563 periodicity + r>=2-lower
            of HYP-2829 supply the wide tail;
       (iv) the k=8 row additionally bounds the -9 S3 + 6 S4 correction (the only non-pairwise
            obligation; the "finite exception" row).
```
This replaces the razor-thin p0 dichotomy (margin 0.13) and the full census by a **low-order
moment-LP with a pairwise (Clebsch design-Hodge) extremal certificate** — comfortable margins
(HYP-2829: 0.90–1.44) and a finite, structured obligation.

## What is solid vs to-verify
- SOLID (this session): moment form supported on S0..S4; S2 leads the binding gap (k=9,10,11);
  consec maximizes the pairwise S2 (k=9,10); the 15-pair = (Z/2)^4 = Clebsch cut-space count.
- SOLID (project): Clebsch = cut-space Cayley graph of K5 (S40); Clebsch closed-nbhd = 2-(16,6,2)
  biplane = the support-six carrier (HYP-2892).
- TO VERIFY (the next computation): that the 6-sector pairwise covariance matrix literally carries
  the biplane / SRG eigenstructure (`N^T N = 4I + 2J`), i.e. the design-Hodge decomposition
  diagonalizes the S2 form and certifies its consec extremum. This is the concrete step that would
  turn the synthesis into a certificate.

## VERIFIED STRUCTURE (S60, `lrc_gk8_pairwise_covariance_structure_macmini_S60.py`) — corrects the naive design ID
Forming the exact `7x7` matrix `M[i][j] = meas{sectors i,j both empty}` at `consec_k` (sector 0 pinned
empty-measure 0 by `e=0`), the two clean guesses are BOTH **refuted**:
- M is **NOT circulant** on the inner 6 sectors (no `Z/7` Gauss-sum/Fano structure at consec — the
  base-path AP breaks the cyclic symmetry; `QR/NQR` test does not even apply).
- M's spectrum is **NOT** a clean design `{a, b(mult m)}` (so not a literal biplane `4I+2J` Gram).

What IS true and is the real certificate target:
- **Reflection symmetry (Z/2):** the diagonal is palindromic, `[0, .272, .244, .239, .244, .272, .335]`
  (k=8) — invariant under the involution `s ↦ 6−s` on inner sectors. This is the complement/reality
  reflection `x ↦ −x` (HYP-2657), the **2-adic** symmetry — confirming the H4/2-adic FACE while
  refuting the literal Clebsch-design realization. (The Clebsch graph is the right H4 *avatar* — the
  cut-space of K5 — but the covariance is its own reflection-symmetric matrix, not the design Gram.)
- **Dominant Perron mode:** one eigenvalue (`0.867 / 0.747 / 0.655` at k=8/9/10) is well-separated and
  carries ~55% of the trace; it is the positively-correlated **concentration mode** (≈ all-ones on the
  empty sectors). This single mode is *why* `S2` is large at consec and *why* the bound is a
  concentration statement.

**Revised route to CRUX 1:** not "diagonalize by a fixed design," but **bound the top (Perron)
eigenvalue of the reflection-symmetric pairwise co-emptiness matrix M.** The `Z/2` reflection
block-diagonalizes `M` into a symmetric (3-dim) and antisymmetric block — the half-tiling halving — and
the Perron mode lives in the symmetric block. So `S2 = (1/2)(1^T M 1 − tr M)` is controlled by a
**3x3 symmetric-block eigenvalue bound** (exact rationals), the genuinely tractable certificate.

## Next
1. Block-diagonalize `M` by the `s↦6−s` reflection; extract the `3x3` symmetric block; bound its Perron
   eigenvalue (hence `S2`) over admissible bounded configs — the finite, low-dim certificate for CRUX 1.
2. Combine with the `−9S3+6S4` k=8 correction (the only non-pairwise obligation) and the HYP-2829 tail
   (bounded finite + single-far THM-563 + r≥2-lower).
3. Feed to Lean: extend `LRCGk8SingleFar.lean` with the moment-order-≤4 reduction and the `3x3`
   reflection-block Perron bound (exact-rational eigenvalues).
NOTE (honest): the literal "Clebsch `4I+2J` certifies S2" claim is REFUTED at consec; the surviving,
sharper claim is the reflection-symmetric Perron bound. The H1↔H4 cross-level reading (a consec/additive
extremal on a 2-adic/reflection carrier) stands.
