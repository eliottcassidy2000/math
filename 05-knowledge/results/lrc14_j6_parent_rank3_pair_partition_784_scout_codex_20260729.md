# Inherited-parent versus literal-child rank-three cap

**Status: FINITE-EXACT SCOPED NEGATIVE CONTROL.  NO ALL-ROOT OR LRC(14)
CLOSURE.**

## Question

THM-2900 closes a forced H4 pair `L` by reconstructing the literal child

```text
R_L=C minus union_(w in L) D_w
```

and proving

```text
q_3(R_L)+B_2(R_L)<|R_L|.
```

Its abstract parent-carrier alternative is

```text
U_C(L)+q_3^C(P union L)+B_2^C(P union L)<|C|.             (1)
```

Formula `(1)` is attractive at all-root scale because one carrier and its
cap bank can be reused across every forced pair.  This scout determines
whether that amortization can replace literal child reconstruction.

## Exact verdict

On the complete `784`-pair control universe from the four THM-2895 roots,
the inherited-parent certificate closes only

```text
204/784
```

rows.  The literal-child THM-2900 certificate closes `784/784`.  Parent
closures by root are

```text
(2,8,9,10,11,13,14):       48/225
(1,3,9,10,11,12,14):       46/118
(2,5,9,11,12,13,14):       83/351
(2,3,4,5,6,7,8):           27/90.
```

The worst parent margin is

```text
-135564641/3577774200
```

at root `(2,5,9,11,12,13,14)`, rank one, apex `16`, flag `(20,25)`.
The closest failure is already strictly negative:

```text
-17489/1873980108
```

at the same root, rank three, apex `17`, flag `(20,38)`.

Thus the parent form is a useful cheap prefilter, but it cannot be used as
a wholesale replacement for literal residuals: it loses precisely the
overlap with the forced flag that makes child singleton and pair caps
smaller.  On this control, `580` obligations still require the child
carrier or another sidecar such as matching overlap credit.

## Verification and scope

The verifier imports the independent rational-interval engine used to
audit THM-2895, reconstructs all 25 scalar-hard marked branches and 784
forced pairs, globally tail-seals every parent pair cap, and checks the
known child theorem as a positive control.  It pays `7,310` finite parent
pair-union evaluations.  Ordinary and optimized outputs are byte-identical.

```text
04-computation/lrc14_j6_parent_rank3_pair_partition_784_scout_codex_20260729.py
SHA-256 007d5da5bd5336fbc33216bc93e6ad81eaa0d42a94dc8b31b1300eb9b2a29813

05-knowledge/results/lrc14_j6_parent_rank3_pair_partition_784_scout_codex_20260729.out
SHA-256 0daa702f97ab1971f917cc53dbabc9606c8444fc064c775d0c05c14cb056ac7c

complete ledger
SHA-256 952b436fe47f9227de7791c7ee788bbb2b9e00b6214f8170c7d449178c51c298
```

The computation is deliberately scoped to the four-root control.  It
proves neither that `(1)` closes exactly 204 rows on another universe nor
that any of the remaining 3,423 seven-body roots closes.
