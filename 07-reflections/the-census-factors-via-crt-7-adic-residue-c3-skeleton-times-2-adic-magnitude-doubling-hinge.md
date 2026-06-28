# The census factors via CRT 14 = 2·7: the 7-adic residue layer (C₃ unit-skeleton, the rigid seed) × the 2-adic magnitude layer (the one doubling hinge, the blind data) — the direction is to attack the 2-adic layer

*mac-mini-2026-06-28-S84. The owner: even more creative reframes, test ideas, be inspired by concurrent work,
push/pull toward all-rigorous. This session FUSES three concurrent threads that independently converged this day:
my C₃ (S83), kps S256's census split (rigid unit-skeleton + one Jacobsthal hinge), and codex's rank-3 unit
nullspace (the unit index is blind to the height/magnitude data). They are ONE structure: the census FACTORS via
CRT `14 = 2·7`. Builds on [[the-hidden-c3-the-witness-space-is-one-galois-orbit-leverage-the-proved-binding-pair-plus-c3]],
[[the-census-splits-rigid-unit-skeleton-and-one-doubling]] (kps), [[lrc14-unit-equioscillation-nullspace-codex-20260628]] (codex).*

## Three threads, one structure (the convergence of 2026-06-28)
- **My S83 (C₃):** the 3 binding-pair witnesses `{1,13},{3,11},{5,9}` are one C₃-Galois orbit; C₃ = (ℤ/7)*/{±1} =
  Gal(ℚ(cos 2π/7)). The proof should be C₃-equivariant from the proved single binding pair (HYP-2909).
- **kps S256 (census split):** binding runners = the units `(ℤ/14)*` = a RIGID residue skeleton; covering runners
  = evens + apex-7 = the MAGNITUDE layer; the entire tight-locus freedom is ONE Jacobsthal doubling `12→24`.
- **codex (rank-3 nullspace):** the six unit gradients have RANK 3 (= the C₃ pairs), and are BLIND to the
  non-unit height/magnitude directions (AP, GW, decoy `2→16` share the unit projection but differ in blind data).

**These are the same fact seen three ways, and CRT `14 = 2·7` is why.** Writing each runner `s` as `(s mod 2,
s mod 7)`:

## The factorization
```
  BINDING = units (ℤ/14)* = {1,3,5,9,11,13} = (odd, nonzero mod 7)   ← 7-ADIC RESIDUE LAYER (the C₃ skeleton)
  COVERING = {2,4,6,8,10,12} ∪ {7}          = (even) or (7-ramified) ← 2-ADIC MAGNITUDE LAYER (the blind data)
```
- **The 7-adic residue layer = the C₃ skeleton.** The units' `mod 7` parts are the 6 nonzero residues, paired by
  `±1` into the 3 binding pairs = the C₃-orbit (S83). This is the rank-3 unit index (codex) = the rigid skeleton
  (kps). It is handled by the C₃-Galois + the proved binding pair (HYP-2909): one orbit-point, spread by C₃.
- **The 2-adic magnitude layer = the doubling hinge = the blind data.** The covering runners and the lone flex
  `12 → 24 = 2·12` (a DOUBLING — contains the old 12-kill lattice and adds interleaved 24-kills, while staying
  14-free). Crucially `12 → 26`
  (residue-preserving, `26≡12`) is NOT tight (M jumps to 1/12): so the magnitude freedom is the **doubling**
  (a 2-adic-flavored ×2 move), not a residue shift. This is exactly codex's "blind height data" — invisible to
  the 7-adic unit index — and kps's Jacobsthal hinge.
- **The apex 7 is the RAMIFIED prime, in the covering layer.** `7 = (1 mod 2, 0 mod 7)`: odd (would be a unit)
  but 7-ramified, so it covers rather than binds. kps's aphorism: *the prime that names the field (7) is the prime
  that does the covering.* In CRT terms, 7 is where the 7-adic factor degenerates (`0 mod 7`), so it falls into
  the magnitude layer despite being odd.

## Why this is the right reframe: two layers, two tools, one CRT
The census `M(S)=1/14 ⟹ S ∈ {AP, GW}·(dilations)` factors into two independent sub-rigidities, one per prime of 14:

| layer | prime | object | freedom | tool | status |
|---|---|---|---|---|---|
| **residue** | 7 | C₃ unit-skeleton (rank-3) | rigid (0 flex) | C₃-Galois + HYP-2909 | **the SEED — nearly done** |
| **magnitude** | 2 | covering runners + doubling | 1 flex (`12→24`) | deformation / 3-gap rigidity | **the BLIND residual — hard** |

This is the descent `14 → 7 → 2` (kps's HYP-3087 CRT bridge / the c=7 and c=2 lifts) made into a clean
factorization of the *rigidity* itself. The 7-part is the residue layer (the C₃, the de Moivre, the units —
the part I, kps, and codex all independently mapped). The 2-part is the magnitude layer (the doubling, the
Jacobsthal, the blind data — the part that is genuinely hard).

## The direction (the creative leap on where to head)
> **Stop spending effort on the 7-adic residue layer — it is the SEED and is essentially in hand** (the C₃ is
> rank-3, the units are forced as residues by HYP-2909 + C₃, the skeleton is rigid). **Concentrate the remaining
> hard rigor on the 2-adic MAGNITUDE layer:** prove the cover's deformation space is exactly 1-dimensional (the
> single `12→24` doubling hinge). This is codex's "blind data" and kps's "deformation space of the cover" — the
> magnitude-level rigidity the residue field ℚ(√−7)/ℚ(cos2π/7) cannot see.
The natural tool for the 2-adic layer is **infinitesimal rigidity** (bar-joint style, kps's pointer): the cover by
13 arithmetic progressions has a rigidity matrix whose nullspace (the flex) must be shown 1-dimensional. The
C₃-symmetry block-diagonalizes the residue part out (it is rigid), leaving the 2-adic magnitude flex as the
single open eigen-direction.

## Codex extension: the C₆ subfield guardrail
The executable HYP-3311 audit, now read as a sidecar complement to HYP-3310, adds one important correction to the
slogan.  In `Gal(ℚ(ζ₇)/ℚ)=C₆=C₂×C₃`, the C₃
quotient is the real-cubic binding-pair orbit:

```
(1,13) → (3,11) → (5,9) → (1,13)
```

But the quadratic `ℚ(√−7)` character is transverse to that quotient: every binding pair contains one quadratic
residue and one quadratic nonresidue mod 7.  So the proof packet has three legal coordinates, not one:

```
C₃ unit-pair skeleton
ℚ(√−7) / χ_7 transverse sidecar
2-adic height/flex ledger on 2U ∪ {7}
```

This is the controlled-forgetting rule.  A quotient may forget the quadratic sidecar or the height ledger only
after proving the LRC predicate is constant on that fiber, reconstructing the lost coordinate, dual-annihilating
it, descending through a HYP-3300 Morse chamber, or routing it to named residual debt.

Incoming HYP-3301 now gives the clean theorem-language for this guardrail: the first coordinate lost by such a
quotient should be recorded as a first-obstruction cocycle or as a Farey-cusp boundary-transfer kernel.  In this
packet, the three possible first losses are exactly the C3 unit-pair orbit, the transverse `chi_7`/`Q(sqrt(-7))`
sidecar, and the 2-adic height/flex ledger.  The sidecar audit is therefore the small exact chart that HYP-3301's
sheaf-exactness theorem needs to overlap with HYP-3310's broader C6 residue-magnitude synthesis.

Incoming HYP-3400 says the same thing in conservation language: these three coordinates are shadow-charge
reservoirs.  Any scalar shadow of the census must preserve them, transfer them to a dual observable, or emit
named charge debt before it can be used as proof currency.

## Honest status
- **CONFIRMED (structural):** the census factors via CRT `14=2·7` into the 7-adic residue layer (C₃ unit-skeleton,
  rank-3, rigid) and the 2-adic magnitude layer (the one doubling hinge `12→24`, the blind data). The apex 7 is
  the ramified prime in the covering layer.
- **SYNTHESIS:** my C₃ (S83) = kps's binding skeleton = codex's rank-3 unit index = the 7-adic factor. kps's
  Jacobsthal hinge = codex's blind height data = the 2-adic factor. One CRT.
- **DIRECTION (the leap):** the 7-layer is the seed (nearly done via C₃ + HYP-2909); the remaining hard rigor is
  the 2-adic magnitude layer — prove the deformation space is 1-dimensional (the doubling hinge). NOT a proof;
  LRC(14) open, but the census is now factored, and the hard residual is localized to one prime (2) and one
  object (the cover's deformation flex).

Related: HYP-3311 (this), HYP-3400 (shadow-charge conservation), HYP-3310 (C6 residue-magnitude frame), HYP-3301 (sheaf exactness / cusp transfer), HYP-3265 (contact-graph case split), HYP-3257 (codex, unit nullspace), HYP-3258 (kps,
census split), HYP-3259 (my S83 C₃, renamed from 3257 — ceded to codex), HYP-3300 (observability/Morse chamber
route), HYP-2909 (proved binding pair), HYP-3087 (kps CRT bridge), LTI-361, LTT-261, OPEN-Q-108.
