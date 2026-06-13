---
source: codex-2026-06-03-S608
status: finite quotient certificate + proof-program
tags: [LRC, n14, Res27, pinch, coimage, Yoneda, THM-401, HYP-2161, HYP-2164]
---

# n=14 `Res_27` Pinch Certificate

The useful new fact from S608 is small enough to say plainly:

```text
In the least-positive C=27 quotient, after unit-shell coverage plus the
canonical D/N blocker gates, pair-sum pinches leave only AP, V*, and 2*AP
at exact floor.  There is no subfloor quotient row.
```

The computation is exact and integer-only.  It enumerates `340,928` deduped
rows that hit all nine unit shells mod `27`, keeps the `27,733` rows that pass
the canonical D/U/N ledger, and checks pair-sum candidates `t=a/(v_i+v_j)`.
The output is:

```text
27,730 strict rows
3 floor rows
0 below rows
```

The floor rows are AP, `V*=(1,2,3,4,5,6,7,8,9,10,11,13,24)`, and the
nonprimitive doubled AP.  Exact maximin confirms all three are truly at
`1/14`.

## Why This Matters

HYP-2161 reframed the problem as:

```text
coimage + Yoneda; C=2n-1 shell probes are the finite observers of cancellation.
```

S608 makes that concrete at `n=14`.  The first `C=27` coimage layer does not
hide an army of new primitive tight rows.  It has exactly the expected two:
AP and the known prime-3 shell swap `V*`.

The incoming HYP-2162/THM-407 result now explains the shell geometry behind
this: the `13` shell face folds to `3` prime-3 strata under `G=<2,-1>`.  The
incoming HYP-2163 proof-pipeline result also separates the no-multiple clock
witness from the remaining multiple/Cprime branch.  S608 is the complementary
census rather than a rival reduction.  It checks the full least-positive
quotient and confirms that no extra primitive floor row appears inside those
strata.

So the open problem becomes cleaner.  We no longer need to wonder whether the
base `Res_27` quotient itself has another primitive floor family.  The burden
moves to lift/CRT conservativity: can a lifted row with the same shell data
change the D/N ledger or endpoint-owner behavior enough to create a genuinely
new floor or subfloor branch?

## Caveat

This is not a proof of LRC at `n=14`.  The script applies D and N gates to the
least-positive representatives in `{1,...,26}`.  Divisibility and `14`-clock
blocking are not determined by shell data alone for arbitrary lifts.  That
caveat is the point: the result identifies exactly what the next theorem must
control.

## Assumption Challenge

The tempting vertex set would be runners or raw residues.  S608 uses
proof-obligation types instead:

```text
missed nonunit shells,
doubled shells,
primitive/nonprimitive status.
```

This preserves the quotient predicate "does the canonical row survive and does
pair-sum prove it?" while forgetting lift sizes and phase order.  Tournament
Analysis on the `148` resulting types is transitive with singleton SCCs and no
directed 3-cycles, which says there is no cyclic residual left at this quotient
level.  The cyclic obstruction, if any, must live one layer higher in the
lift/CRT data.
