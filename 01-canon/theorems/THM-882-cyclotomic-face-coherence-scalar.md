---
id: THM-882
title: THE CYCLOTOMIC FACE AND THE COHERENCE SCALAR — (Φ) the tight locus of LRC(14) is the root set of the 14th cyclotomic polynomial: Φ₁₄(z) = z⁶−z⁵+z⁴−z³+z²−z+1 = Φ₇(−z) exactly, so tight times = 1/2 + {±1,±2,±3}/7 — THE HALF-TURN TIMES THE HEPTAGON STAR: the 14 = 2·7 factorization is geometric (the "2" is z ↦ −z, the "7" is the F₇ deck), identifying the μ₆ spine (THM-873) with the heptagon/token machinery (THM-773) as ONE object; (C) THE COHERENCE SCALAR answering S181's open call: saw(S) = Σ_{pairs} [ρ(a,b) − 4/169] is the first repo scalar that, like L itself, is DILATION-invariant but NOT translation-invariant — it separates S181's matched-energy twins (AP {1..13} vs 2·{1..13}−1, both E = 1469: saw = +0.602 vs +0.116, 5×; maxtree 0.668 vs 0.416) because additive energy counts sum-coincidences (translation-blind) while the sawtooth evaluates WHERE the sums sit mod 13; the GAP case confirms Freiman dimension remains the second axis (saw-sum alone does not separate GAP from AP) — S181's "resonance-lattice geometry = dimension + coherence" now has coordinates: (Freiman dim, saw-sum), with the Y*-spectrum as the resonance-depth profile
status: Φ PROVED (classical identity Φ₂ₙ = Φₙ(−z) for odd n, verified symbolically; the tight-time translation p/14 − 1/2 = (p−7)/14 with |p−7| ∈ {2,4,6}); C VERIFIED on S181's exact twin battery (4 sets, exact ℚ: M, E, saw, maxtree, Y*-spectrum) with the mechanism proved (translation-variance of the sawtooth)
source: opus-2026-07-16-S321 (owner: more coordinatizations; merge past unrelated threads)
depends_on:
  - THM-873 (the tight locus this gives the polynomial face to)
  - THM-863/864 (the sawtooth pair law the scalar is built from)
related: [S181 (additive-energy correction — the open call answered), S179/S175 (the energy thread), THM-773 (the F₇ tokens now identified with the tight locus's heptagon component), kps a-dictionary-of-loop-maps (OPUC/Blaschke — the tight locus as a Blaschke zero-set is the named next coordinate)]
verification: 05-knowledge/results/resonance_coordinates_s181_merge_opus_S321.out
---

# THM-882 — the cyclotomic face and the coherence scalar

## (Φ) The cyclotomic face

THM-873's tight locus {p/14 : gcd(p,14) = 1}, viewed on the unit circle
z = e^{2πit}, is exactly the root set of **Φ₁₄(z) = z⁶−z⁵+z⁴−z³+z²−z+1** —
degree φ(14) = 6, the μ₆ spine's defining polynomial. The classical identity
Φ₂ₙ(z) = Φₙ(−z) (n odd) gives **Φ₁₄(z) = Φ₇(−z)**: every tight time is a
half-turn away from a 7th root of unity,

> tight times = 1/2 + {±1, ±2, ±3}/7.

The 14 = 2·7 factorization is thereby GEOMETRIC: the factor 2 contributes
the half-turn z ↦ −z; the factor 7 contributes the heptagon. The deck's
F₇ token machinery (THM-773: coverage ⟺ X⁷ − X | Π(X − k_a)) and the μ₆
tight-locus spine (THM-873) are the same object in two coordinate systems —
the deck sees the heptagon component, the maximin sees its half-turn.
S52's Eisenstein lever (14 a primitive 6th root mod Φ₆(14) = 183) is the
third face: Φ₆, Φ₇, Φ₁₄ triangulate the problem's cyclotomic core.

## (C) The coherence scalar (S181's call, answered)

S181 proved additive energy E(S) is necessary-not-sufficient for tightness
and called for "resonance-lattice geometry (dimension + coherence)" — with
no coherence coordinate available. The sawtooth pair law (codex-S14 +
THM-863) supplies it:

> **saw(S) = Σ_{pairs} [ρ(a,b) − 4/169]** — dilation-invariant (ρ is
> scale-free), NOT translation-invariant (the sawtooth evaluates a+b mod 13
> POSITIONS, not coincidence counts).

Exact battery (S181's twin sets):

| set | E | M (exact) | saw | maxtree |
|---|---|---|---|---|
| AP {1..13} | 1469 | 1/14 (TIGHT) | **+0.602** | 0.668 |
| 2·{1..13}−1 | **1469** | 1/2 (loose) | **+0.116** | 0.416 |
| GAP {1..7}∪{10..16} | 1546 | 2/17 (loose) | +0.619 | 0.715 |
| near-AP {21..32} | 1156 | 21/53 (loose) | +0.043 | 0.304 |

The matched-energy twins (E = 1469) are separated 5× by saw — the exact
failure mode of E (translation-blindness: E(S+c) = E(S)) is the exact
strength of saw. The GAP row shows saw alone does not replace Freiman
dimension: **the S181 geometry now has coordinates (Freiman dim, saw), with
the Y*-spectrum as the per-pair resonance-depth profile.** (Note: S181's
L-column was good-set measure; the M-column here is the exact maximin —
verdicts agree, normalizations differ.)

## Named next

- The OPUC/Blaschke coordinate (kps's loop-map dictionary): the tight locus
  as the zero set of a degree-6 Blaschke product; Verblunsky coefficients of
  the safe-set measure as wall coordinates — the next coordinatization to
  build.
- saw(S) along the well atlas (mediant vs deep wells) and the signed
  selectors: is the coherence scalar the (Z/14)*-character's absolute value?
- The Φ₇-transport: rewrite the deck's token criterion at the tight locus
  via z ↦ −z and see what the half-turn does to token walks (THM-779).
