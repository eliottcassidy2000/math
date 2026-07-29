---
id: THM-2900
title: "Flag-conditioned rank-selective partition closure"
status: >
  RESERVED / UNPROVED PROOF CANDIDATE UNDER LOCKED-ARTIFACT AUDIT.  Combining
  THM-2895's forced H4 pair with THM-2897's odd rank selection replaces the
  child B2+q1 cap by B2+q3.  A fresh independent census closes all 784
  four-root child residuals, but promotion awaits a hash-pinned replay.
source: root-2026-07-29
related:
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2898-unique-max-gate-six-seed-parity-closure
verification:
  - 04-computation/lrc14_j6_child_rank3_pair_partition_closure_thm2900.py
  - 05-knowledge/results/lrc14_j6_child_rank3_pair_partition_closure_thm2900.out
---

# THM-2900 -- flag-conditioned rank-selective partition closure

**RESERVED / UNPROVED PROOF CANDIDATE UNDER LOCKED-ARTIFACT AUDIT.**

## 1. Candidate composition lemma

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

The child certificate `(1)` and the parent five-slot certificate
`q_5(C)+2B_2(C)<|C|` are incomparable.  Subtraction shrinks pair caps but
can remove two low-ranked labels and make the child's third remaining
singleton larger than the parent's fifth rank.  A proof pipeline should try
the parent certificate first because it is amortized once, then use `(1)`
only on the forced flags that survive.

## 2. Candidate four-root compression

On the `25` marked nonscalar branches of THM-2895, the H4 cores produce
exactly `784` pair residuals.  The candidate independent census finds

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

Thus the rank-selective child cap would close the four roots after one H4
descent.  None of THM-2895's five recursive H2 rows would be needed.

No proved result is claimed by this reservation until the ordinary and
optimized replays agree, the complete ledger is locked, and the artifact
hashes are recorded.  Even after promotion, the exact application would
cover only four roots, not the seven-body rung or LRC(14).
