---
id: HYP-6830
title: Scale uniformization of Claim B — the two-regime splice (sheets above c=43, bounded band inflation below), the fragmentation⟺divisibility complementarity, and the discrete-LRC recursion probes
status: OPEN — THM-761 proves the large-scale regime; the splice composition, the complementarity claim, and the r>=7 tiling residue are the open content
source: opus-2026-07-14-S299
depends_on:
  - THM-755   # capped envelope, v* = r_P / (pi |G'_P|)
  - THM-760   # r=1 coprime sheet dodge
  - THM-761   # multi-exception sheet covering bound (this session)
  - HYP-6780  # v*(cP) = c v*(P): the scale covariance that killed raw-height bands
related: [THM-756, THM-757, THM-758, HYP-6785, HYP-6815, HYP-6820]
---

# HYP-6830 — scale uniformization of the ≥4-far covering endgame

## The claim to prove (the splice)

For every covering 13-speed family V and every scale c ≥ 2, write V = cP ⊔ W
(canonical: P = {v/c : c|v}, r = |W|). Define c*(V) = the largest c with |P| ≥ 7
(equivalently r ≤ 6), or 1 if none. Then:

1. **Large scale (PROVED, THM-761):** if c*(V) ≥ 43 (or the exact per-(r,c)
   criterion fires, or the gcd budget Σ g_a(⌊c/(7g_a)⌋+1) ≤ c−1 holds), V is closed:
   M(V) ≥ 1/14, by a free witness sheet. No enumeration, no witness search.
2. **Small scale (OPEN — composition):** if c*(V) ≤ 42, then every dilated
   sub-structure inflates the capped-envelope band edge by at most the SAME bounded
   factor (v*(cP) = c·v*(P) ≤ 42·v*(P), HYP-6780 used positively), so the ≥4-far
   band protocol (THM-755/756) runs over a genuinely bounded normalized domain.
   TO PROVE: the composition statement — every ≥4-far covering family with
   c*(V) ≤ 42 either safe-peels (THM-753), fires the capped envelope, or lies in an
   explicitly bounded normalized band family. This is the corrected, quantifier-honest
   replacement for THM-758 Claim B's "finite check" and for the REFUTED q≤25 finish
   (codex-S3, family 26·{1..12} ∪ {339}, first witness q=27).
3. **The complementarity claim (OPEN — the load-bearing sentence):** the good-set
   fragmentation r_P of a 12-core is large ONLY in the presence of divisibility
   structure. Quantitatively: there is an explicit function B(c) such that a core
   with c*(core) ≤ 42 has r_P ≤ B(42) — i.e., component count is controlled by the
   maximal dilated-part scale, not by raw height. (The mechanism: G'_{cQ} = the
   c-fold shrunken copies of G'_Q, so components multiply exactly under dilation;
   the claim is that NOTHING ELSE multiplies them. If true, regimes 1+2 cover
   everything and the covering endgame is uniformly finite.)
   FALSIFIER SHAPE: a scale-free core (no c ≥ 2 dividing ≥ 7 elements) with r_P
   unbounded along a family. Adversarial search must include near-dilates
   (MISTAKE-140 discipline) and multi-gcd hybrids.

## The residual inside regime 1 (named by THM-761)

- **r ≥ 7** (core ≤ 6 elements): union bound structurally cannot fire (7 arcs of
  length 1/7 can tile). The wall is realized (c=7 instance, all sheets bad) but the
  realized family is still lonely (M = 1/7) by a non-sheet route. Routes to close:
  (a) sweep t0 over the core's whole safe set, not one optimum (the origin-cap move
  on the deck); (b) signed/spectral refinement: rotated arcs of length exactly 1/7
  tile Z_c only under rigid divisibility (Newman-type cyclic tiling conditions —
  vanishing of arc Fourier coefficients on subgroup frequencies); classify tiling
  residue profiles as the AP/GW-analogue corners of the deck. (c) capped envelope
  applied to the exceptions within a sheet.
- **Deep gcd entanglement** (Σ g_a over budget): recursive descent c → c/g;
  bookkeeping to write; termination is clear (c strictly decreases).

## The recursion (the structural content, from the S299 reflection)

The sheet residual at scale c IS an inhomogeneous discrete lonely-runner instance on
Z_c: runners = exceptional residues w_a mod c, offsets = w_a·t0/c, radius 1/14,
lonely time = free sheet. The tight case (r = 7, arcs tile) is the 7-clock partition
(THM-754) one level down: tight instances are tilings at every level. The underlying
object is a pointed circle with burned arcs, self-similar under scale descent; the
free sheet is the next level's basepoint (the observer-lens principle survives
descent).

## Probes filed (testable, not claimed)

- **FI cubic probe** (arXiv:2605.29035 structural analogy): on the 8,260-family band
  bank, test whether an exact third-moment (cubic) functional of the good-set
  indicator decides the 19 direct-L closures that the quadratic disc certificate
  (THM-731/732) misses — if yes, the band protocol becomes one uniform inequality
  (strengthen-then-deduce, the Frank–Ivanisvili shape: bulk by low moments, sharp
  constant one moment deeper).
- **Rédei deck-parity probe:** does the sheet deck carry a parity statement (odd
  number of free sheets under a suitable weighting) refining bare existence? Aim at
  equality/boundary structure only (guardrail C18).

## Verification and tooling

- THM-761 battery: 04-computation/lrc14_multi_exception_sheet_bound_opus_S299.py
  (+ .out): 50,964 exact counting instances, exact failure sets, end-to-end 13-speed
  closures including the q25-refutation family, r=7 wall realization.
- Terminal certificate in the library: `sheet_certificate(speeds, c)` in
  lrc14_certificates.py (self-test 15/15) — speed-only, O(r), exact.
- Regime-2 composition and the complementarity claim are the natural next sessions:
  (i) prove r_P ≤ B(c*) or find the falsifier; (ii) write the bounded normalized
  enumeration for c* ≤ 42; (iii) the r ≥ 7 deck classification.
