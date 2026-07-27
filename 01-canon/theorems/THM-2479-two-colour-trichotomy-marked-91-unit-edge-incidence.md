---
id: THM-2479
title: "Two-colour trichotomy: 91-unit c3-edge incidence with a word-marked endpoint on every positive middle-owner word stratum"
status: >
  PROVED, with one explicitly flagged extraction: the load-bearing
  Lemma A (unconditional coherent fibre) is THM-2323 Section 7.2's own
  construction (32)-(38), which two independent audits this session
  verified to be free of hypothesis (31) gcd(a,91)=1 -- (31) enters
  only at the final unit conclusion (39)/(40).  Owner (mac-mini)
  countersign of that extraction is requested before downstream use as
  a proved dependency.  VERIFIED-EXACT companion: 444,000-pair exact
  residue census + six fibre configurations in Z[zeta_364], zero
  failures.  NOT claimed: a two-word unit edge when 7|a; cubic/address
  participation; target-plane gain; terminal-component phase; any
  scalar-row decrement.  The ledger stays 165; LRC(14) remains OPEN.
source: >
  mac-mini-2026-07-25 (MSG-2153 composition; THM-2326 quantifier repair
  MSG-2152); audited (two independent adversarial lenses) and canonized
  by opus-2026-07-26.
depends_on:
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
  - THM-2326-vertexwise-septimally-primitive-c3-degree
related:
  - THM-2293 (91-unit edge existence, unlocated)
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
  - THM-2315 (target-gain corolla)
  - HYP-9031-d5-h1-dictionary-lrc-word-current-vs-jc-flux (prediction 1)
script: 04-computation/msg2153_trichotomy_hostile_probe_opus.py
output: 05-knowledge/results/msg2153_trichotomy_hostile_probe_opus.out
---

# THM-2479 -- the two-colour trichotomy

## Statement

Let a strict canonical row among the `150` carry a positive-return
middle-owner word stratum: owner `j=2`, positive canonical word `Q` with
`measure(E_Q)>0`, `c_2=13^b u_2`, `c_3=13^c u_3`, `c>b`. Put
`g=gcd(c_2,c_3)`, `a=c_2/g`, `d'=c_3/g`; then `13|d'` and `13` does not
divide `a` (THM-2323 (30), using `c>b`). Fix any nonzero root character
`kappa` mod 13. Then the exact-grade-`b`, root-`kappa` spectrum of the
bare middle source `1_{E_2}` contains a `c_3` edge

```text
X - Y = M c_3,     gcd(M, 91) = 1,     0 < |M| <= a(1092 S^2 - 1) + 14S - 1,
```

with both endpoints of exact 13-valuation `b` and root character
`kappa`, and **at least one endpoint simultaneously a bare atom and a
literal positive-word atom of `1_{E_Q}`** (a THM-2323 common bare/word
atom). Precisely:

- **(i) `7` does not divide `a`.** Then `gcd(a,91)=1` (13-nondivision is
  automatic) and THM-2323 Section 7.2 as published gives the edge with
  BOTH endpoints word-marked and `|M| <= a(1092S^2-1)`.
- **(ii) `7 | a`.** Lemma A gives common bare/word atoms
  `A=c_2 q_{z_1}`, `B=c_2 q_{z_2}` with `B-A=m_2 c_3`, `m_2 = a beta`,
  `gcd(beta,91)=1`; hence `7|m_2` and `13` does not divide `m_2`.
  Lemma B applies at `A` and yields a bare canonical atom
  `C = A + m_1 c_3` with `1 <= m_1 <= 14S-1`, `7` not dividing `m_1`,
  same grade and root character. Trichotomy:
  - if `13` does not divide `m_1`: the edge `(A,C)` has
    `gcd(m_1,91)=1`, endpoint `A` word-marked;
  - if `13 | m_1`: then `C != B` and the edge `(B,C)` has multiplier
    `m_1 - m_2`, nonzero mod 7 (`= m_1`) and nonzero mod 13 (`= -m_2`),
    hence `gcd(m_1-m_2, 91)=1`, endpoint `B` word-marked.

  In case (ii) the third vertex `C` is bare-only and is NOT a
  `c_2`-multiple (THM-2323 (40): a `c_2`-multiple partner would force
  `a | m_1`, contradicting `7`-freeness of `m_1`): it is exactly the
  affine-coset escape that (40) proves necessary. **No word coefficient
  at `C` is used or asserted** (the care point of THM-2323 S7.1 and
  MSG-2153). Edges are undirected; signs are immaterial.

## Lemma A (unconditional coherent fibre; extracted from THM-2323 S7.2)

On every stratum as above and every prescribed nonzero root character,
the modulus-`N=91d'` fibre construction (32)-(38) of THM-2323 S7.2 --
`K_z = K_0 + d'z` primitivity, the primitive same-gauge landing of the
common bare/word atoms `q_z = K_z + N h_z`, and the factorization (38)
-- holds for EVERY `a`: it supplies adjacent common bare/word atoms
`A = c_2 q_{z_1}`, `B = c_2 q_{z_2}` of exact grade `b` and root
`kappa` with `B - A = a beta c_3`, `gcd(beta, 91) = 1`,
`0 < |beta| <= 1092 S^2 - 1`. Hypothesis (31) `gcd(a,91)=1` is used in
S7.2 only to convert `m = a beta` into a `91`-unit at (39); it is NOT
used in the construction. (Verified independently by two audit lenses
this session, reading (32)-(38) directly: `H_Q, H_E` are supported in
`D_1` by S6, so the plain one-arc same-gauge theorem applies at
`N = 91d'` for every `a`; `K_z` primitivity and the `K_7 x K_13` /
`K_6 x K_13` adjacency use only `d'`. Countersign requested.)

## Lemma B (THM-2326 in its MSG-2152 repaired form)

For a marked atom `A` with `(1_{E_2})_hat(A) != 0` and
`nu_13(A) = b < c = nu_13(c_3)`, with `E_2` disjoint from `D_{c_3}`
(strict-row exclusivity), THM-2326 yields `1 <= m_1 <= 14S-1` with `7`
not dividing `m_1` and `(1_{E_2})_hat(A + m_1 c_3) != 0`; the marked-
atom hypothesis guarantees `A + m_1 c_3 != 0` and preserves exact grade
and root character. `A = c_2 q_z` satisfies the hypotheses
(`13` does not divide `q_z`). MSG-2152 should be folded into the
THM-2326 canon file by its owner; this file cites the repaired form.

## Hostile companion (exact)

`msg2153_trichotomy_hostile_probe_opus.py`: (PART A) all `444,000`
exact residue pairs `(m_1, m_2)` with `7 not| m_1`, `7|m_2`,
`13 not| m_2`: zero cases with `m_1 = m_2`, zero cases where neither
`(A,C)` nor `(B,C)` is a `91`-unit edge; all `53,760` grade/root
retention rows pass. (PART B) six fibre configurations
`(a,d') in {(7,13),(7,143),(14,13),(7,2431),(1,91),(2,91)}`: `K_z`
distinct mod `N`, retained fibre full (`78` or `91` per parity of
`7|d'`), never edgeless; `Z[zeta_364]` arithmetic exact. MSG-2152's
sharp signed-step boundary hostile does NOT extend to carry a
`7|m_2, 13 not| m_2` pair against the composition (checked in-probe).

## What this closes and what remains

This discharges the "factor 13 / incidence composition" residual named
in THM-2323's connection contract (42) and THM-2326's Section 4: on
every positive middle-owner word stratum, gain selection now has an
unconditional 91-unit *marked* edge at the residue-label level. The
single remaining open frontier of the bispectrum thread is the
**address/gain participation** (THM-2321 (51), THM-2315 (30)) plus
terminal-component phase -- i.e. exactly the `H^2`/incidence half in
HYP-9031's language (prediction 1 of that file, existence half, is
hereby fulfilled modulo the Lemma A countersign). No scalar row is
removed; LRC(14) remains open.
