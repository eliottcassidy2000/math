---
id: THM-3669
title: "LRC typed-control all-packet three-twist defects"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the canon-recorded THM-2334 typed non-cover row, every one of the 120
  owner-pivot packet choices has a nonzero support-minimal three-twist target
  defect at the trivial center.  The packets collapse to thirty ordered graft
  charts and only thirteen distinct present-set geometries; both exact
  cyclotomic embeddings certify every defect.  One chart additionally has a
  nonzero mass defect, an exclusively owned boundary jump, and a nonzero first
  marked-product Fourier jet.  This is a positive control on a typed non-cover
  row, not a covering-row transfer or an LRC(14) proof.
source: kps-s193 / THM-3666 physical positive-control continuation, 2026-08-21
audit: >
  PASS -- Bohr reproduced normal and optimized transcripts and the semantic
  digest, then independently checked all thirty charts, their fourfold
  omitted-label multiplicity, the thirteen present geometries, signs, deep
  phases, primality/order certificates, canonical residuals, mass valuations,
  exclusive boundary, and first marked-product jet.  The compact 8269-byte
  packet ledger reproduced SHA-256
  9521c4654c9fddb43094e1194b36ef75f6118b8b216918c3f80787270c6db55a.
  No correction was required; the non-cover scope remains explicit.
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-3665-lrc-support-minimal-three-twist-target-detector
  - THM-3666-lrc-owner-pivot-dual-pair-swap-twist-basis
related:
  - THM-2327-two-colour-marked-unit-c3-triangle
  - THM-2350-owner-pivot-dual-dipole-normal-form
script: 04-computation/lrc_positive_control_all_pair_swap_defects_thm3669.py
output: 05-knowledge/results/lrc_positive_control_all_pair_swap_defects_thm3669.out
script_sha256: 757d7c8521a1756ec15ee5785f37fac89c20b1c9cfc0a21003b2f822b8ef29f4
output_sha256: 1b951cdd9863c1bd65f5e80263f2e99fab000fabed2ef35b19e1254c7527bd43
semantic_sha256: a33d438898b3c714557bc652add55c82e7d648ded6e9416a5ddb0a674282e27f
hash_basis: raw bytes
---

# THM-3669 -- every typed-control owner packet has a local target defect

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem evaluates THM-3665's support-minimal detector in the
physical pair-swap coordinates of THM-3666.  Its purpose is to separate
detector availability from the still-open covering-row transfer.

## 1. The fixed typed control

Use the canon-recorded primitive row and word

```text
w=(1,14,27,40,53,66,13,2197,742586),
owner=c1, targets a=c2 and b=c3,
word={a}, R=169,
(X,m,Y)=(13,1,742599).                               (1)
```

This row has strict target profile `(1,3,5)` and satisfies every local typing
condition used by THM-2334.  It is explicitly **not** a covering row.

Let an owner-pivot packet choose an omitted unit `u_0` and ordered distinct
graft units `(k_a,k_b)` among the remaining five.  THM-3666 identifies the
target-twist basis as

```text
alpha=e_a-e_(k_a),
beta =e_b-e_(k_b).                                  (2)
```

Write `H(x,y)` for the full marked current, including the deepest-target
phase, at twist `x alpha+y beta`.  The local detector at the trivial center is

```text
D_(u0,ka,kb)=H(0,0)+H(-1,0)-2H(0,-1).              (3)
```

By THM-3665, any nonzero value in (3) proves nonzero unrestricted mod-13
target mass for that packet.

## 2. Thirty charts exhaust all 120 packets

The omitted label does not occur in either vector in (2).  Hence fixed
`(k_a,k_b)` gives the same three twists for each of the four legal omitted
labels.  There are therefore

```text
6*5=30 ordered graft charts,
30*4=120 owner-pivot packets.                        (4)
```

The three-term census needs only the untwisted present set, six negative
`alpha` shifts and six negative `beta` shifts: thirteen distinct interval
geometries in total.  For each geometry the companion independently refines
the exact half-open Boolean intervals, intersects the 169 word preimages, and
computes the two endpoint factors in the exact embeddings

```text
Z[zeta_NN] -> F_(p_j),
NN=50334435734703120,

(p1,zeta_NN)=(352341050142921841,435817657216),
(p2,zeta_NN)=(956354278959359281,153943385426666320). (5)
```

The result is

```text
D_(u0,ka,kb) != 0 in both embeddings
for all 120 legal packet choices.                    (6)
```

One nonzero finite-field image already proves that the cyclotomic integer is
nonzero; the second is an independent exact control.  The compact 120-row
ledger has 8269 bytes and digest

```text
9521c4654c9fddb43094e1194b36ef75f6118b8b216918c3f80787270c6db55a. (7)
```

An independent blob-pinned implementation produced the same digest.

## 3. Canonical chart and local boundary anatomy

For `(u_0,k_a,k_b)=(0,1,2)`, the three exact current factors `gamma`, before
multiplication by their common nonzero scalar, are

```text
gamma(0)      =(310354333794505177,928156429768775202),
gamma(-alpha) =(173156527261826497,288184598848954178),
gamma(-beta)  =(317831260334985463,423756868609874985). (8)
```

Thus (3) has images

```text
(200189390529282589,368827291397979410).            (9)
```

None of the thirteen marked products is empty.  In the canonical chart the
three overlap numerators and component counts are

```text
overlaps  =(60084076348296,54135630512964,61887542465528),
components=(188056,169431,186674).                  (10)
```

Their signed mass defect is

```text
-9555378069796/NN
 =-183757270573/967969917975060 !=0,                (11)
```

with raw 13-adic valuations `(2,1,1)` and defect valuation `1`.

The marked products themselves give two finer local witnesses.  On the
`NN` endpoint grid, the first endpoint (in the deterministic term/endpoint
ordering) owned by exactly one of the three jump maps is

```text
t=7791347595131100/NN=1609245/10396204.             (12)
```

It is a left endpoint of the untwisted product, with raw and three-term jump
both `+1`; neither shifted product jumps there.  “Owned by exactly one” in
(12) is local to that endpoint and does not assert that no other exclusively
owned endpoint exists.  The first marked-product Fourier jet is already
nonzero at frequency one, with images

```text
(212565399985231344,197131116389696829).            (13)
```

The diagnostics (11)--(13) concern the marked product.  Only the full
cyclotomic residual (9) includes the separate bare-`Y` endpoint factor and
deep phase; no inference from a diagnostic alone is used to prove (6).

## 4. What the positive control changes

THM-2334 previously proved full-bank nonconstancy on this row.  The genuinely
new information is physical sparsity and uniformity:

```text
every owner packet works,
the same three-site detector works at the trivial center,
thirteen interval geometries replace 169 generic twists.           (14)
```

Consequently detector choice and packet choice are not the obstruction on
the typed control.  The live LRC obligation is a transfer theorem forcing one
such defect, boundary jump, valuation mismatch, or equivalent nonconstancy on
a hypothetical covering row while retaining the all-`91`-unit projector,
visible height and terminal phase.

No ancestry-digit identification is made, and (6) is only for the unrestricted
mod-13 target aggregate.  LRC(14) remains open.

## 5. Exact companion

Reproduce with

```bash
python3 -B 04-computation/lrc_positive_control_all_pair_swap_defects_thm3669.py
python3 -B -O 04-computation/lrc_positive_control_all_pair_swap_defects_thm3669.py
```

The assertion-free companion source-pins both independent THM-2334/3666
engines, checks the thirteen geometries, all thirty charts and all 120 packet
rows, reproduces the independent compact ledger byte contract, and verifies
the mass, jump and first-jet controls exactly.  **QED.**
