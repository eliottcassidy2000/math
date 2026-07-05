# The gap violator wears handcuffs at every prime

**kind-pasteur-2026-07-05-S3 (HYP-4105)**

## The setting

The last open predicate of LRC(14) is `TightLooseDichotomy`: every compressed
covering primitive family's base is dilated-AP-shaped or carries margin 2/25.
Its truth is equivalent to the **spectral gap**: no primitive 12-set has
`M(W) ∈ (1/13, 2/25)`, plus the exact-tight rigidity `M = 1/13 ⟹ AP`.
This session built the formal interrogation kit for a hypothetical gap violator.

## What ¬loose forces — the filter package (all now kernel-pure Lean)

Failure of the loose branch means: at EVERY time t some runner is within 2/25 of
an integer. Evaluate at t = a/q and clear denominators (`not_loose_eval`):
some runner has `|w·a − q·m| < (2/25)·q`. Integer conclusions:

- **q ≤ 12** (`(2/25)q < 1`): `q ∣ w·a` for some runner — for every direction `a`.
- **13 ≤ q ≤ 25** (`(2/25)q ≤ 2`): `w·a ≡ 0, ±1 (mod q)` for some runner —
  for every direction `a`.

The second family is the striking one. At a PRIME `p ∈ {13, 17, 19, 23}` with no
p-multiple in W, running over all directions `a` and inverting
(`pair_pinning_13`-style): **W's residues must hit every unit ±pair class mod p.**
Count the handcuffs:

| p  | pair classes to cover | elements available |
|----|----------------------|--------------------|
| 13 | 6                    | 12                 |
| 17 | 8                    | 12                 |
| 19 | 9                    | 12                 |
| 23 | **11**               | **12**             |

At p = 23 a gap violator without a 23-multiple must scatter its twelve residues
over ALL eleven ±pair classes — at most one doubled class. Simultaneously it
must cover every modulus ≤ 12 in every direction, be spread (`w_max > 11.5·w_min`,
the S2 gate), and be primitive. The dilated AP `c·{1..12}` passes every filter
with zero slack (its ±residues mod 23 are `±c·{1..12}` = all eleven pairs
exactly once — the unique perfect scatter, up to the one doubling). The filters
say: *a gap violator must impersonate the AP at every prime simultaneously.*

## The height side — the merge exclusion

THM-592 (mac-mini): M(W) is attained where the last lonely component dies — a
merge event at radius `d/(v+w)`, `v, w ∈ W`. Values in the gap `(1/13, 2/25)`
of the form `d/(v+w)` require:

- `d = 1`: `v+w ∈ (12.5, 13)` — **empty**;
- `d = 2`: `v+w ∈ (25, 26)` — **empty**;
- `d = 3`: `v+w ∈ (37.5, 39)` — first candidate `3/38`.

So any gap violator's binding merge is **at depth d ≥ 3 with v + w ≥ 38, hence
w_max ≥ 19**. The gap census with all filters active (P4, [1,24]) returns zero
violators; the low-height regime is closed by exhaustion, and every sweep the
fleet has run agrees.

## Why this transcends the session

1. **Necessary conditions are cheap; the geometry is one theorem away.** Every
   filter took ~10 lines of Lean (evaluation + omega). What remains open is not
   more congruences — it is the covering-efficiency statement that a 12-runner
   arc system covering [0,1] at radius < 2/25 wastes nothing, which only the AP
   achieves. The filters shrink the suspect pool; they do not convict. The
   conviction lives on the THM-592 merge geometry, where the fleet's tower
   (klein S134) and CRT sweeps (mac-mini S51/S52) are already digging.

2. **The certificate language is now closed both ways.** `rational_point_margin`
   (S2) turns integer conditions into margins; `cert_of_margin` (S3) turns
   margins into integer conditions with modulus `s ≤ B/(2(β'−β)) + 1`. At the
   verified rigidity slack (second value 14/169 vs margin 2/25) this is
   `s ≈ 176·B`. Loose is no longer an analytic property — it is a bounded
   integer search, formally.

3. **The asymmetry of effort inverted.** A month ago the loose side was the
   "analytic mystery" and the tight side looked combinatorial. Now the loose
   side is decidable-with-certificates and the tight side (rigidity) is the
   analytic residue. Formalization did not follow the mathematics here; it
   reorganized which mathematics is hard.

Related: [[the-loose-witness-is-one-decidable-inequality-kps-S2]],
[[the-tight-locus-rigidity-of-lrc13-kps-S1]].
