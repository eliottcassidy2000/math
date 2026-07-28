---
id: THM-2742
title: "Full two-target present sheet and deepest-source semantic current"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the canonical
  typed row, restoring the second lawful target coordinate to the source-one
  present packet changes the THM-2720 deepest-source wall into exactly 936
  positive sections.  Every one of those sections also meets the genuine
  E3 -> D^6 -> Q_(3,{1,2}) semantic handoff, and every nonzero character of
  the restored C13 target coordinate survives.  The aggregated target
  current is target-coordinate inversion-even.  This is a target-active
  semantic mass current, not a physical-deck/target identification, paired
  left relation current, row exclusion, or LRC(14) conclusion.
source: root/full-two-target-present-sheet-2026-07-28
audit: full-target-sheet-hostile-audit-2026-07-28 (independent periodic-antiderivative reconstruction, witness audit, character and replay checks)
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2407-owner-or-source-deletion-target-current-dichotomy
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2720-unshifted-deepest-source-present-packet-global-disjointness
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2712-semantic-following-congruence-lock-and-address-coboundary-descent
script: 04-computation/lrc14_two_target_present_semantic_attachment_probe_20260728.py
output: 05-knowledge/results/lrc14_two_target_present_semantic_attachment_probe_20260728.out
script_sha256: 062b352f4db12a5f01822b293cdbb10629632dacc5fa27b406d8dd321e550709
output_sha256: d7159343e91b593d4be670cd7a53b89b5423ea077f78fc038eb7766c43939c03
hash_basis: LF-normalized bytes
---

# THM-2742 -- the full target sheet repairs the deepest-source present wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Work on the canonical typed row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5).                  (1)
```

The one-target present grammar of THM-2720 is the zero slice of a lawful
two-target object.  Restoring its missing coordinate repairs the uniform
deepest-source disjointness, reaches the actual terminal semantic fork, and
produces every primitive character of that target coordinate.

## 1. Lawful two-target sheet

Write `d_v(theta)` for the ordinary width-`1/14` danger comb of speed `v`
centred at `theta`, and put `g_v(theta)=1-d_v(theta)`.  The width-`1/7`
guard-safe set is denoted by `C_H`; it is not an ordinary safe comb.
The lawful target differences from THM-2365/2407 are

```text
eta    = e_c2-e_q1,
lambda = e_c3-e_q2.                                    (2)
```

For `ell in Z/7` and `(s,t) in F_13^2`, define the full source-one present
section

```text
F_(ell,s,t)
 = d_c1(ell/7)
   * 1_C_H * product_(i=3)^5 g_qi(0)
   * g_q1(-s/13) * g_c2(s/13)
   * g_q2(-t/13) * g_c3(t/13).                         (3)
```

Thus `(s,t)` acts with opposite signs on the two graft pairs in `(2)`.
Equation `(3)` is the inherited source-deleted `U_(s,t)` packet with the
owner probe `d_c1(ell/7)` reinserted.  It neither deletes nor relabels any
further factor.  The old present packet is exactly `F_(ell,s,0)`.

Let

```text
A0 = C_H intersect product_(i=1)^5 g_qi(0),

E_j = A0 intersect d_cj(0)
             intersect product_(i!=j) g_ci(0).          (4)
```

The one-target slice has the complete compatibility matrix

```text
E1 intersect F_(ell,s,0) != empty  iff ell=0,
E2 intersect F_(ell,s,0) != empty  iff ell!=0 and s!=0,
E3 intersect F_(ell,s,0)  = empty  for all ell,s.       (5)
```

Its positive-section counts are `13,72,0`.  The last equality is precisely
THM-2720: `E3` contains `d_c3(0)` while the zero slice contains `g_c3(0)`.

On the full sheet the sharp law is instead

```text
E3 intersect F_(ell,s,t) != empty
        iff ell!=0 and t!=0,                            (6)
```

with no condition on `s`.  Exactly

```text
6*13*12=936                                             (7)
```

of the `7*13^2=1183` sections are therefore positive.  On the exact integer
grid used by the companion their unnormalized masses range from

```text
242609986080  at (ell,s,t)=(1,9,1)

to

1362887826300 at (ell,s,t)=(3,0,2).                    (8)
```

The zero slice supplies the inherited hostile control, while every label in
`(6)` is a changed-coordinate positive control.

## 2. Attachment to the terminal semantic fork

Source compatibility alone would not recover the blocker-word handoff.  Put

```text
Q=Q_(3,{1,2})
 =A0 intersect d_c1(0) intersect d_c2(0) intersect g_c3(0). (9)
```

For every label satisfying `(6)`, the exact periodic prefix-correlation
identity gives

```text
H_(ell,s,t)
 = measure(E3 intersect F_(ell,s,t) intersect D^(-6)Q)>0. (10)
```

All `936` values are positive, and their exact extrema are

```text
114819491/12545122758259
 <= H_(ell,s,t) <=
672887730/12545122758259,                               (11)
```

again attained at `(1,9,1)` and `(3,0,2)`, respectively.  A strict rational
witness for the first lexicographic positive section is

```text
(ell,s,t)=(1,0,1),
x=15991693680925/100360982066072,
D^6 x=3179229/20792408.                                (12)
```

Direct factor evaluation places `x` in `E3 intersect F_(1,0,1)` and
`D^6x` in `Q_(3,{1,2})`; its source word is exactly `(3,{1,2})`.  Hence the
new mass attaches to the genuine `E3 -> D^6 -> Q_(3,{1,2})` semantic cospan,
not merely to an isolated source point.

## 3. Every primitive target character survives

Fix `ell!=0` and `s`.  The rational thirteen-vector

```text
t |-> H_(ell,s,t)                                      (13)
```

is zero at `t=0` and strictly positive at all twelve nonzero values.  If one
primitive Fourier coefficient vanished, its rational coefficient polynomial
of degree at most twelve would be divisible by

```text
Phi_13(X)=1+X+...+X^12.
```

All thirteen entries would then be equal, contradicting `(13)`.  Galois
conjugacy therefore proves

```text
sum_(t in F_13) H_(ell,s,t) zeta_13^(tau*t) != 0
        for every tau in F_13^*.                       (14)
```

There are `6*13*12=936` nonzero marked coefficients in `(14)`.  Summing over
all `(ell,s)` preserves the zero-at-zero/positive-away-from-zero profile, so
all twelve primitive deepest-target characters survive after aggregation as
well.  The aggregate profile satisfies the exact symmetry

```text
A_t=A_(-t).                                             (15)
```

Its coefficients are consequently real and agree in `tau,-tau` pairs.  The
current is target-active but target-coordinate inversion-even; it does not by
itself orient the two conjugate target directions or supply a geometric
reflection.

## 4. Type boundary and surviving debt

The `t`-character `tau` in `(14)` is the endpoint-dipole character
`lambda.(u-v)` for `lambda=e_c3-e_q2`; at `m_deep=1`, its full target residue
is `q_b=1+tau`.  It is not the root-six physical deck character selected by
the THM-2712 atom atlas, and it is not the left relation residue required by
THM-2334.  This theorem supplies no comparison map between those coordinates.
The separately VERIFIED, noncanonical calculation
`07-reflections/lrc14-two-index-carrier-is-spectral-not-a-packet-label-20260728.md`
finds all-depth physical/target transversality; in particular, a common
numeral modulo `13` is not an identification.

Thus this theorem supplies

```text
deepest semantic source and terminal word
 + full lawful two-target present sheet
 -> positive semantic mass table
 -> all primitive deepest-target characters.           (16)
```

It does not supply a common affine half/`C_221` cycle, a private-root or
primitive-unit endpoint current, a fixed-`X` physical/target intertwiner, or
the paired left relation residue.  No all-row statement, scalar-row
exclusion, or LRC(14) conclusion follows.

## 5. Exact reproduction and audit

Run

```bash
python3 04-computation/lrc14_two_target_present_semantic_attachment_probe_20260728.py
python3 -O 04-computation/lrc14_two_target_present_semantic_attachment_probe_20260728.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_two_target_present_semantic_attachment_probe_20260728.out
```

with the LF-normalized SHA-256 values declared in the frontmatter.  The
companion pins the audited canonical interval carrier, rebuilds the three
sources and terminal fork, scans all `1183` sections, applies the exact
`D^6` prefix identity, constructs witness `(12)`, and reduces the target
coefficients in the power basis of `Q(zeta_13)`.  It contains no truth-bearing
Python `assert`.

An independent audit rederived the signs in `(3)` from `(2)`, rebuilt all
cells by a separate periodic-antiderivative implementation, and reproduced
the `936/1183` laws and exact extrema.  It checked `(12)` factor by factor,
replayed all `936/936` marked and `12/12` aggregated primitive-character
certificates, matched the carrier and artifact hashes, and verified normal,
optimized, and stored transcript equality.

QED.
