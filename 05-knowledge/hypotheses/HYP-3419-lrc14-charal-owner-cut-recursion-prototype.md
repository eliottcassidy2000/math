---
id: HYP-3419
title: LRC14 charal owner-cut recursion prototype
status: SYNTHESIS / finite cut-recursion API prototype; auxiliary to HYP-3415 critical-path map
source: codex-2026-06-28 second pass over incoming HYP-3410, rebased after HYP-3415/S258
tangent: T1380
technique: LTI-380
tournament_technique: LTT-280
script: 04-computation/lrc14_charal_owner_cut_recursion_prototype_codex_20260628.py
result: 05-knowledge/results/lrc14_charal_owner_cut_recursion_prototype_codex_20260628.out
reflection: 07-reflections/lrc14-charal-owner-cut-recursion-prototype-codex-20260628.md
related:
  - HYP-3418
  - HYP-3417
  - HYP-3416
  - HYP-3415
  - HYP-3410
  - HYP-3409
  - HYP-3408
  - HYP-3407
  - HYP-3406
  - HYP-3405
  - HYP-3404
  - HYP-3402
  - HYP-3401
  - HYP-3311
  - HYP-3310
  - HYP-3301
  - HYP-3265
  - HYP-3124
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3419: LRC14 Charal Owner-Cut Recursion Prototype

## Claim

The HYP-3410 Bring/Schwarz-Christoffel/BDH/Menger/charal atlas points to a
more concrete finite lemma target:

```text
charal quotient -> mixed theorem-exit fiber -> owner cut sidecar
-> binary cut recursion -> terminal theorem exit or named debt
```

The first two HYP-3406/HYP-3410 mixed fibers are separated by one owner label.
The newer `10->20` frontier is not: it needs a size-`3` owner cut and an exact
depth-`3` decision tree.  So the next theorem should not be a one-label owner
theorem.  It should be a bounded owner-cut recursion theorem.

After the concurrent HYP-3415/S258 critical-path map, this prototype is not a
competing completion route for LRC14.  The main proof spine is still:

```text
q-witness for non-covering sets
-> LRC<=13 induction on Q
-> one uniform decorrelation-floor inequality |SPEC| < product
```

The role of HYP-3419 is auxiliary but concrete: it can classify finite
exception packets, certify owner cut-code purity, and name the residual
accessory/height/off-grid debts that must not be hidden inside the final
decorrelation floor.

After rebasing over HYP-3416 and HYP-3417, the local division of labor is:
HYP-3416 supplies the abstract recursive quotient ladder; HYP-3417 supplies
the signed owner-current certificate language; HYP-3419 supplies the concrete
charal decision-tree and cut-purity harness.  In particular, HYP-3419 should
feed HYP-3417-style certificates when a finite cut is found, and should record
the exact quotient failure when no bounded tree appears.

S259/HYP-3418 adds one sharper warning: the covering floor appears 2-adic,
with even speeds as the binding obstruction.  In this light, the `10->20`
frontier tree ending in the even-cover label `2:g2`, and HYP-3417's frontier
current `{2:g2,11:g1,13:g1}`, should be treated as a live even-cover
sidecar.  The charal owner-cut recursion should separate "binding labels" from
"2-adic/even-cover labels" before importing any apex-7, Galois, or scalar
floor story.

## Executable Readout

Script:

```text
04-computation/lrc14_charal_owner_cut_recursion_prototype_codex_20260628.py
```

Stored output:

```text
05-knowledge/results/lrc14_charal_owner_cut_recursion_prototype_codex_20260628.out
```

The exact owner-cut recursion over the three HYP-3410 fibers is:

```text
height_leak_12_family:
  min_owner_cut_size = 1
  core = ('5:g1',)
  optimal tree:
    test 5:g1
      yes -> positive-Haar-open
      no  -> unit-petal-named

persistent_owner_leak_26_40_54_family:
  min_owner_cut_size = 1
  core = ('1:g1',)
  optimal tree:
    test 1:g1
      yes -> unit-petal-named
      no  -> positive-Haar-open

height_persistent_owner_leak_10_20_drop_add_family:
  min_owner_cut_size = 3
  min_owner_cut_count = 5
  core = empty
  top labels across minimum cuts:
    13:g1 appears in 4/5 cuts
    11:g1 appears in 3/5 cuts
    2:g2 appears in 3/5 cuts
  one optimal tree:
    test 13:g1
      yes -> positive-Haar-open
      no  -> test 11:g1
        yes -> positive-Haar-open
        no  -> test 2:g2
          yes -> positive-Haar-open
          no  -> unit-petal-named
```

This is the strongest new signal.  The frontier does not merely say "owner
support matters"; it says the owner proof must allow a bounded multi-label cut
with no single mandatory label in the intersection of all minimum cuts.

## Charal Purity

On the represented HYP-3410 fibers, several charal levels already split the
visible rows.  The useful distinction is not that turn words fail on this
small sample.  The useful distinction is proof legality:

- a turn word is a Schwarz-Christoffel-style shadow;
- owner labels are the accessory parameters that can certify a cut;
- a full charal signature is safe for finite classification but too close to a
  row label unless the proof names the sidecars it needs.

The prototype therefore favors a small owner-cut sidecar over full row-level
charal retention.  Known exact cuts are `1,1,3`; that is small enough to become
a finite lemma target and large enough to rule out the overly optimistic
one-label theorem.

## Guardrail Integration

The requested outside ideas are used as gates, not as imported proof engines:

```text
Bring radical:
  five-exit branch alphabet; not a quintic formula chase.

Schwarz-Christoffel:
  contact turn word plus accessory owner debt; turns alone are demoted.

Barban-Davenport-Halberstam:
  finite owner-channel variance; ranks candidate labels but does not replace
  the exact cut certificate.

Menger:
  bounded endpoint-owner separator; this is the current theorem-shaped core.

Krasner:
  owner/contact stability gate for same-residue or p-adic-near lifts.

Sophie Germain:
  quartic height/flex debt can split into two quadratic channels, but only
  after the recursion names a live height sidecar.

Hermite-Lindemann-Weierstrass:
  no-scalar-shadow firewall; transcendental shadows cannot certify finite
  rational packets.

Ramanujan-Soldner:
  zero-level normalization hygiene for branch potentials, not a mixed-fiber
  separator.

Meissel-Mertens:
  tail-entropy label after finite exceptional owner cuts are named.
```

## Candidate Lemma

For every expanded-bank mixed charal fiber after residue/height sidecars:

1. either a bounded endpoint-owner cut makes `theorem_exit` a function;
2. or the fiber admits a dual owner-current / Farkas certificate;
3. or a stable charal ladder preserves the positive-Haar-open exit;
4. or it terminates at AP/GW, strict-open mass, q-witness, H7/state lift,
   Schwarz-Christoffel accessory debt, exact-period/BDH exception,
   height-factor debt, or a newly named finite residual.

The first finite target is to prove, or refute by first counterexample, a
bounded owner-cut theorem on the next expanded HYP-2963 bank.  The known bound
on represented leaks is `3`.

## Tournament Analysis

Vertices are proof modules and guardrail gates, not runners, arcs, or named
constants.

Pairwise observable:

```text
finite exactness + cut power + recursion - scalar risk
```

Fingerprint:

```text
vertices = 8
score_hist = {28:1, 31:1, 49:1, 57:1, 60:1, 70:1, 79:1, 80:1}
directed_3cycles = 0
hamiltonian_path =
  M00 -> M01 -> M02 -> M03 -> M04 -> M05 -> M07 -> M06
```

Priority path:

```text
bounded owner-cut theorem
-> charal decision-tree API
-> finite BDH label variance
-> Krasner owner-stability gate
-> Schwarz-Christoffel accessory reconstruction
-> Sophie-Germain height-factor channel
-> Soldner-HLW-Mertens scalar firewall
-> Bring branch alphabet
```

The order is deliberate: Menger/charal finite exactness beats branch
symbolism, and finite BDH is a label-ranking tool rather than a replacement
for the cut certificate.

## Assumption Challenge

Alternate vertices considered:

```text
runners, owner labels, cut sets, charal turns, BDH labels, Bring branches,
Schwarz-Christoffel prevertices, algebraic constants, proof modules.
```

The chosen vertices are proof modules and guardrail gates.  The preserved LRC
predicate is theorem-exit purity after controlled forgetting.  Raw data that is
destroyed by a quotient is allowed back only when a later mixed fiber demands
it as a named sidecar.
