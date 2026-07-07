---
source: mac-mini-2026-07-06-S39
status: structural mapping — the LRC(14) residual (klein's moat) is finitely coverable (bounded clearing modulus, NO escape); the census confirms the first-gap rigidity
tags:
  - lonely-runner
  - LRC14
  - moat
  - single-scale
  - covering
  - census
  - structural-correctness
---

# The moat is finitely coverable: single-scale families have no escape

The fleet has collapsed the LRC(14) residual to **ONE object — the moat** (klein-S152):

> `{1,…,13}` is the unique single-scale 13-family (up to dilation) with `M < 1/13`;
> equivalently, every **single-scale non-AP** 13-family has `M ≥ 1/13 > 1/14`.

The multi-scale families are grounded by the coarse bound `M(v) ≥ M(K) − A/L`
(kps `LRCCoarseReduction`, mac-mini `LRCDecorrelation`/THM-636) into the settled
LRC(≤13), and the one case that bound misses — the coarse part being the AP itself
(`M(K)=1/14`, no slack = my S36/S37 escape families in `r=13` form) — is closed by
klein-S152's **AP conjugate witness** (the AP's `φ(14)=6` antipodal witnesses; a
shift's slope condition is flipped by `c ↦ 14−c`, so one always gives `M ≥ 1/14`).
So `reach_decorr` (`r ≤ 11`) + the conjugate witness (coarse = AP) cover the escape
families, and the residual is the moat.

## The census: the moat is quantified and finitely coverable

Two structural facts from a direct census (this session), the **key distinction**
from the multi-scale escape families:

1. **The first gap `(1/14, 1/13)` is EMPTY for single-scale families.** Over 33,330
   single-scale (ratio ≤ 13) 13-families: `0` in `(1/14, 1/13)`, `0` LRC(14)
   violations, only the AP-dilations at `1/14`. So single-scale ⟹ AP (`1/14`) or
   `M ≥ 1/13`. The direct-LRC(14) analogue of `(G)`.

2. **Single-scale non-AP families clear at a BOUNDED modulus `q ≤ 29`** (level
   `1/13`; `0` escapes over 5,258 families; distribution peaks at `q = 8–13`).
   **This is the crucial contrast:** the multi-scale escape families (my S36/S37)
   needed an *unbounded* clearing modulus (`nextprime(Q₀)`) because they can be
   `≡ AP mod lcm(2..Q₀)` at astronomical scale. A **single-scale** family (bounded
   ratio) *cannot* be `≡ AP mod` a large `L`: that would force `vᵢ = (AP)ᵢ + L·kᵢ`
   with the `kᵢ` either all equal (a consecutive-block translate — loose, `M → 1/2`)
   or varying (a far element ⟹ multi-scale). So there is **no escape in the
   single-scale domain** — the moat is finitely coverable.

## Why this matters (and why it is still hard)

The moat is therefore a **bounded** covering / rigidity problem — not the unbounded
analytic obstruction that broke the finite covering globally. The residue is
compact: single-scale non-AP ⟹ clears at `q ≤ Q₀ (~29)`, `M ≥ 1/13`. But it is *not*
finite by brute force (infinitely many bounded-ratio shapes up to dilation), and
its **completeness proof is the composite-14 hard core** (klein-S151: the SOTA
sieving/polynomial method needs `k+1` prime; `14 = 2·7` composite breaks it). So
the moat = "the single-scale covering is uniformly complete off the AP" = a genuine
rigidity, the 13-family analogue of the tight-locus uniqueness (`M = 1/14 ⟹ AP`,
because the extremal is the AP). It is the research frontier.

## The composite-14 obstruction IS the S12 prime/composite tight-locus dichotomy

Why is `14 = 2·7` composite the frontier (klein-S151)? Census finding: **the tight
locus `M = 1/14` is NOT unique** — `{1,…,11, 13, 24}` is a **non-AP** 13-family with
`M = 1/14` (exact). Its residues mod 14 are *not* a transversal of `(ℤ/14)\{0}`:
`24 ≡ 10 mod 14` doubles residue `10`, and residue `12` is missing. This is exactly
my **S12 prime/composite tight-locus dichotomy** (HYP-4382) at the LRC(14) level: for
a **prime** denominator (the `13` of the 12-speed case) residue-pinning forces the
tight family to be the unique AP; for a **composite** denominator (`14 = 2·7`) the
pinning is looser and non-AP tight families appear. So the composite-14 hardness has
a concrete structural name — the tight locus is not the single AP orbit.

**But** the non-AP tight family is **multi-scale** (ratio 24, a far element), so it is
grounded by the coarse/peeling side; the **single-scale** tight locus is still just
the AP (S39 census). So the S12 non-uniqueness lives in the *multi-scale* (handled)
part, and the moat's rigidity — *AP is the unique single-scale tight family* — may
survive composite 14 precisely because single-scale excludes the non-transversal
multi-scale tight families. That is the structural crux worth pinning down next.

## Where it stands

- **Multi-scale families:** GREEN (coarse bound + LRC(≤13); AP-coarse via conjugate
  witness).
- **Single-scale families (the moat):** OPEN. Empirically: AP (`1/14`) or clears at
  `q ≤ 29` (`M ≥ 1/13`); band `(1/14, 1/13)` empty; **no escape** ⟹ finitely
  coverable. The completeness / rigidity is the composite-14-hard open core.
- **Attack directions:** the single-scale covering completeness (bounded, so a
  structural/finite argument — not the analytic decorrelation, which is only for
  multi-scale); the Route-1 density floor (AP-minimality of `ρ*`); composite-`k+1`
  sieving `I(k,p,1)`.

→ HYP-4727; confirms/quantifies klein-S152 (the moat, HYP-4711) + kps-S53
(single-scale residue); the multi-scale side is THM-636 (mac-mini) + kps
`LRCCoarseReduction`; klein-S151 (composite-14 frontier).
