# HYP-2445 - Diagonal Frobenius obstruction is a support gate

**Status:** OPEN synthesis; exact numerology confirmed by a small atlas.
**Source:** codex-2026-06-12.
**External source:** Benjamin Church, arXiv:2508.14876,
https://arxiv.org/abs/2508.14876.
**Companions:** HYP-2425, HYP-2430, HYP-2441, HYP-2443, HYP-2444, OPEN-Q-069.
**Computation:** `04-computation/product_quotient_support_gate_atlas_codex.py`;
stored output `05-knowledge/results/product_quotient_support_gate_atlas_codex.out`.

## Statement

Church's product-quotient obstruction should be imported into the repo as a
support-gate template:

```text
scalar property passes  +  retained side channel descends  =>  obstruction
```

In the paper, Shioda supersingularity is the scalar property.  It is not
enough: the examples are simply-connected and Shioda supersingular at
infinitely many good primes, yet not unirational.  The proof-bearing channel is
instead the existence, on every asymmetric partial Frobenius twist, of diagonal
symmetric differential forms.  These force rational and elliptic curves into
finite exceptional types or into a partial-Frobenius descent where a projection
degree drops.

That is the same proof shape as the recent repo support gates:

```text
Type II [72,36,16]:
  scalar Gleason enumerator passes;
  support/design/matroid realization is the obstruction.

LRC14:
  scalar "q blocked" ledgers are too coarse;
  Q27 classes, 13-clock, divisor fibers, and Bprime targets carry the proof.

Order-5 [72,36,16] branch:
  projected Type II [16,8] type is too coarse;
  marked tetrad support kills d16+ and forces an e8+e8 heptad split.
```

## Exact Numerology Worth Keeping

The atlas records an unexpectedly tight bridge around the number `1092`:

```text
|PSL2(F_13)| = 1092 = 84*(14-1) = 13*84 = 2^2*3*7*13.
```

This is simultaneously:

- the Hurwitz-bound order for a genus-14 curve in Church's construction;
- the LRC14 single-stranger dominance cutoff `13*84` from HYP-2444;
- `12*C(14,2) = 14*C(13,2) = 3*C(14,3)`.

The three subgroups used by the paper produce the sharper indices:

```text
D6, A4: index 91 = C(14,2) = the LRC14 Q27 rescue q=91.
D7:     index 78 = C(13,2) = the [72,36,16] minimum-design lambda_5.
```

These are search beacons, not proofs.  In particular, `78` should not be read
as an automorphism claim for a putative length-72 code.  It is a hint that the
minimum-word design ledger and the genus-14 Hurwitz subgroup ledger may share a
useful incidence arithmetic.

The first explicit supersingular prime norms in the paper also echo the LRC14
residue language:

```text
1091, 2339, 6551, 7643 are all -1 mod 13.
2339 and 6551 lie in the LRC14 missing +/-10 class mod 27.
```

Again, this is only a beacon.  The atlas keeps it separate from the actual
obstruction mechanism.

## Proposed Transfer Lemma

Try to prove an abstract support-gate lemma with the following parts:

1. A scalar quotient `Q` that is necessary or naturally visible but lossy.
2. A retained side channel `S` that survives the relevant twists, fibers,
   markings, or deletion/carry operations.
3. A descent or finite-exception statement: any bad object either lands in a
   named exceptional type or descends to a smaller resource.

For Church's surfaces, `Q` is Shioda supersingularity, `S` is the diagonal
form package on all partial Frobenius twists, and descent lowers projection
degree.

For LRC14, the proposed analogue is:

```text
Q  = blocked denominator / no plain witness,
S  = Q27 resource vector
     (shell-27 unit class, 13-clock, divisor fiber, Bprime target),
D  = each blocked shell spends an independent resource or opens an exceptional
     Bprime/owner-private deletion route.
```

For `[72,36,16]`, the proposed analogue is:

```text
Q  = the healthy Type II scalar enumerator,
S  = minimum-word support as a 5-(72,16,78) design / binary matroid,
D  = automorphism or local-neighborhood gates either force a support split
     or create forbidden low weight.
```

## Tournament Analysis

The computation uses proof routes as vertices, not curves, runners, or
codewords.  The pairwise observable is:

```text
(retained_channel, theorem_status, repo_bridge, computability, -risk)
```

The tie Hamiltonian path is declaration order.  The stored run is transitive:

```text
score histogram: {0:1, 1:1, ..., 7:1}
directed 3-cycles: 0
SCCs: singleton
Hamiltonian paths: 1
leader: diagonal_forms_all_partial_frobenii
```

Assumption challenge: alternate vertex sets considered were curves, quotient
surfaces, diagonal forms, Frobenius twists, subgroup indices, supersingular
prime residues, LRC runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, Pisano classes, code supports, matroid circuits, and
proof obligations.  The chosen quotient preserves proof-channel strength and
computability.  It destroys the actual surface geometry, explicit Magma group
data, code support realizability, and multi-stranger LRC interactions.

## HYP-2469 Addendum

HYP-2469 upgrades this bridge using the newer LRC14 finite atlases.  The `1092`,
`91`, and `78` coincidences are now explicitly treated as beacons rather than
proof.  The proof-bearing transfer is the Church grammar:

```text
scalar quotient + retained side channel on every twist/fiber
+ finite exceptions or strict descent.
```

For LRC14, HYP-2469 identifies the retained channel as Q27 obligations plus
`13`-clock, deleted-core, shell-27, divisor-fiber, owner/Bprime, and support-load
data.  HYP-2463, HYP-2464, and HYP-2465 now certify the one-stranger, hard-stack,
two-stranger residual, and near-core regimes.  The live problem is reduced to
two portals: below-nine-core descent and outside-window normalization.

## Next Moves

1. Use the `D7` index `78` as a design-ledger search seed for the
   `5-(72,16,78)` minimum layer.
2. Use the `D6/A4` index `91` as the geometric analogue of the LRC14
   `q=91` fibered rescue.
3. Try to model HYP-2444's Q27 resource descent after Church's partial
   Frobenius descent: each obstruction step should lower a resource or land in
   a finite exceptional type.
4. Build a small incidence toy: compare the subgroup coset geometries of
   `D7`, `D6`, and `A4` inside `PSL2(F_13)` against the support parameters
   of the `[72,36,16]` minimum design, explicitly checking what is preserved
   and what is destroyed.
