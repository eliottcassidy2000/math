---
id: THM-696
title: The test-data existence and the supply discharge — (1) THE FATTENING BRIDGE supply_of_strictlyLive: any strictly-live ruler (Q, w) for speeds ≤ B yields the certificate-supply citation's EXACT data (interval [(28Bw−1)/(28BQ), (28Bw+1)/(28BQ)], q = D; the integer margins Q+1 ≤ 14r ≤ 13Q−1 fatten uniformly) — every witness theorem (THM-693/694/695) now feeds THM527ACertificateSupply directly; (2) testData_exists: the witness data (q*, a = 1, c) exists from the SHAPE alone (6 > 5 pigeonhole for q*; q* ∉ P forces the P-side at a = 1; the missed class from S243's pigeonhole via the missed_bridge); (3) twoScale_supply: THE CITATION'S CONCLUSION (∃ interval + thirteen SafeIvStrict certificates) PROVED WHOLESALE for every non-top-block two-scale family (P.card ≤ 5, some q* ∈ [8,13] avoiding P with |E| < q*, V > 2184, 168e < V) — the certificate-supply citation is DISCHARGED on two-scale shapes with no per-family data
status: PROVED AND FORMALIZED (LRCTestDataSupply.lean: supply_of_strictlyLive, qstar_exists, missed_bridge, testData_exists, twoScale_supply — all kernel-pure [propext, Classical.choice, Quot.sound] (missed_bridge even [propext, Quot.sound]), root-wired, project green; one k-case fix from first compile). CONSUMPTION: mac-mini cont.27's lrc14_from_certificate_citations needs THM527ACertificateSupply — this theorem proves its conclusion on every non-top-block two-scale realization shape; the remaining supply content = [top-block small parts] ∪ [non-two-scale shapes (the coverage audit)] ∪ [bounded V (windows/banks)].
source: klein-2026-07-11-S248 (HYP-5965; the test-data existence directive + the last-critical-bits push)
depends_on:
  - THM-693 / LRCTwoScaleWitness (the strictly-live ruler this fattens)
  - LRCTestPointCore (the pigeonhole), LRCMeasureTransfer (SafeIvStrict)
  - mac-mini LRCReachDecitation (THM527ACertificateSupply — the citation discharged here on two-scale shapes)
related:
  - THM-690/691/692 (the continuum existence these make constructive-and-automatic)
---

# THM-696 — test-data existence and the supply discharge

## (1) The fattening bridge

For a strictly-live ruler (Q, w) — every residue r_i = (v_i·w) % Q with
Q < 14r_i < 13Q — integrality gives Q+1 ≤ 14r_i ≤ 13Q−1, and the explicit
interval

> [x/D, y/D] = [(28Bw−1)/(28BQ), (28Bw+1)/(28BQ)], q = D = 28BQ

carries SafeIvStrict certificates for every speed 0 < v_i ≤ B: with
j_i = (v_i w)/Q, the scaled endpoint values are 28B·r_i ∓ v_i, and
28B(14r_i − Q) ≥ 28B > 14B ≥ 14v_i closes both strict inequalities. The
data satisfies the citation's side conditions (0 ≤ x, y < D, D < (y−x)·q).
So EVERY strictly-live ruler — THM-693's two-scale, THM-694's multi-scale,
THM-695's ray witnesses — produces the supply citation's certificate data.

## (2) The test-data existence (the directive's lemma)

From the shape alone: P.card ≤ 5 leaves some q* ∈ [8,13] outside P (six
values, five speeds); at a = 1 the P-side is automatic (q* | p forces
p = q* ∉ P for p ≤ 13 < 2q*); and the missed class c comes from the S243
pigeonhole whenever |E| < q*, bridged to the witness hypothesis by
(c − e) % q* ≠ 0 ⟺ e % q* ≠ c. No coprimality search, no per-family data.

## (3) The supply discharge on two-scale shapes

Composing (2) → THM-693's twoScale_strictlyLive (a = 1) → (1):

> **every two-scale family (P.card ≤ 5, non-top-block — ∃ q* ∈ [8,13]∖P
> with |E| < q* — V > 2184, 168e < V per co-offset) satisfies the
> certificate-supply conclusion**: an interval with thirteen strict band
> certificates exists.

With mac-mini's `lrc14_from_certificate_citations`, the remaining content
of the LAST geometric citation is now: [the five-block-type top-block small
parts] ∪ [realization shapes that are not two-scale — the shape-coverage
audit] ∪ [bounded V — THM-686 windows and the banks]. THM-661 and the
≤7-arcs pigeonhole remain as the other two named citations.

## Formalization & files

`04-computation/lean/TournamentH7/TournamentH7/LRCTestDataSupply.lean`
(kernel-pure, root-wired, project green).
