---
id: THM-2900
title: "Flag-conditioned rank-selective partition closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.  Combining
  THM-2895's forced H4 pair with THM-2897's odd rank selection replaces the
  child B2+q1 cap by B2+q3.  A fresh independent census closes all 784
  four-root child residuals, eliminating the recursive H2 layer.
source: root-2026-07-29
depends_on:
  - THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
related:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2898-unique-max-gate-five-parity-matching-closure
verification:
  - 04-computation/lrc14_j6_child_rank3_pair_partition_closure_thm2900.py
  - 05-knowledge/results/lrc14_j6_child_rank3_pair_partition_closure_thm2900.out
---

# THM-2900 -- flag-conditioned rank-selective partition closure

**PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.**

## 1. Composition lemma

Let a five-cover problem on a carrier `C` satisfy the THM-2895 singleton
cap, so every hypothetical cover contains an unordered pair from its finite
core `H_4`.  For every such pair `L`, form the literal residual

```text
R_L=C minus union_(w in L) D_w
```

with the actual excluded-prefix sidecar.  Let `q_3(R_L)` be the third
largest allowed singleton coverage of `R_L`, and let `B_2(R_L)` be any
global allowed pair-union cap.

If every nonempty residual satisfies

```text
q_3(R_L)+B_2(R_L)<|R_L|,                                  (1)
```

then no five-cover of `C` exists.

Indeed, THM-2895 forces a pair `L`.  Any three labels covering the residual
have a member of singleton coverage at most `q_3(R_L)`; their other two
labels cover at most `B_2(R_L)`.  THM-2897's rank-selective odd partition
bound makes their total coverage strictly less than `|R_L|`, a
contradiction.

The same composition applies at every odd residual size:

```text
(2m+1)-cover after a forced flag:
q_(2m+1)(R)+m B_2(R)<|R|.                                 (2)
```

It is a sufficient certificate, not an equivalence.  The literal residual
and excluded labels remain essential.

There is also a parent-carrier form.  Suppose a `k`-cover of `C` is forced
to contain a flag `L`, let `P` be the prior excluded set, and let
`F^C_j(P union L)` be any THM-2897 partition cap computed on `C` after
forbidding the flag labels.  Then

```text
U_C(L)+F^C_(k-|L|)(P union L)<|C|                          (3)
```

also excludes the cover.  Indeed, write the hypothetical cover as the
disjoint union of `L` and its remaining labels, use subadditivity once, and
apply the partition cap to the remainder.  At the three-slot leaf, `(3)`
becomes

```text
U_C(L)+q^C_3(P union L)+B^C_2(P union L)<|C|.              (4)
```

This inherited form can reuse a parent cap bank and avoids reconstructing
the child interval set.  It may be weaker or stronger than a fresh child
calculation; neither direction is automatic.

The child certificate `(1)` and the parent five-slot certificate
`q_5(C)+2B_2(C)<|C|` are incomparable.  Subtraction shrinks pair caps but
can remove two low-ranked labels and make the child's third remaining
singleton larger than the parent's fifth rank.  A proof pipeline should try
the parent certificate first because it is amortized once, then use `(1)`
only on the forced flags that survive.

## 2. Four-root compression

On the `25` marked nonscalar branches of THM-2895, the H4 cores produce
exactly `784` pair residuals.  The independent census finds

```text
q_3(R)+B_2(R)<|R|
```

strictly on all `784`.  The smallest margin is

```text
1720879997/753665913000
```

at

```text
E=(2,8,9,10,11,13,14), rank=1, apex=19, H4 pair=(37,125).
```

Thus the rank-selective child cap closes the four roots after one H4
descent.  None of THM-2895's five recursive H2 rows is needed.

## 3. Verification and scope

```text
04-computation/lrc14_j6_child_rank3_pair_partition_closure_thm2900.py
SHA-256 21a3cc54ab91d289128bede2c09c5639f787e9891dee4ab0ce7e215f033ba7e7

05-knowledge/results/lrc14_j6_child_rank3_pair_partition_closure_thm2900.out
SHA-256 7f8c0de2d02fa5f9564b3f89ea61fb45240b2ecc8d54dad7d8b2eb24ebb32d74
```

The verifier imports the independently reconstructed rational interval
engine from THM-2895's audit, not the primary LRC implementation.  Ordinary
and optimized replays are byte-identical, the source contains no Python
`assert`, and the complete child ledger digest is

```text
60ad47208695c5aa22ec79ca428eccfdd5e214cd7cf4a20a85629aa5023a6fa2.
```

The abstract composition lemma is general.  The finite-exact application
covers only the four roots already proved by THM-2895; it compresses their
proof and the future pipeline but does not close another root, the
seven-body rung, or LRC(14).
