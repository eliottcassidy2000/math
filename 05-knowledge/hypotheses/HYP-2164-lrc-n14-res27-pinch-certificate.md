---
id: HYP-2164
status: finite quotient theorem / proof-program
source: codex-2026-06-03-S608
related:
  - HYP-2162
  - HYP-2163
  - HYP-2088
  - HYP-2096
  - HYP-2101
  - HYP-2135
  - HYP-2137
  - HYP-2138
  - HYP-2161
  - THM-401
---

# HYP-2164: n=14 least-positive `Res_27` pinch certificate

## Claim

In the least-positive `C=27` residue quotient for LRC at `n=14`, the primitive
floor branch is exactly AP plus `V*`.

More precisely, enumerate all `13`-subsets

```text
R subset {1,...,26}
```

that hit every unit antipodal shell `{a,27-a}` with `gcd(a,27)=1`.  Apply the
canonical representative versions of the S573 clock-blocker gates:

```text
D_q: some r in R is divisible by q, for every 2<=q<=13;
N_j: the j/14 clock is blocked at distance 0 or 1/14, for every 1<=j<=7.
```

Among the `27,733` surviving rows, the pair-sum/pinch sieve gives:

```text
27,730 strict rows with score > 1/14,
3 floor rows with score = 1/14,
0 rows with score < 1/14.
```

The three exact floor rows are:

```text
AP      = (1,2,3,4,5,6,7,8,9,10,11,12,13),
V*      = (1,2,3,4,5,6,7,8,9,10,11,13,24),
2*AP    = (2,4,6,8,10,12,14,16,18,20,22,24,26).
```

The first two are primitive.  The third is nonprimitive and normalizes to AP.

## Evidence

The script `04-computation/lrc_n14_res27_pinch_certificate_s608.py` performs a
deduplicated exact enumeration of the unit-shell rows and integer-only
pair-sum scoring.  Its stored output is
`05-knowledge/results/lrc_n14_res27_pinch_certificate_s608.out`.

Key counts:

```text
raw 13-subsets of nonzero residues: 10,400,600
unit-shell rows, deduped: 340,928
D-pass among unit rows: 29,384
N-pass among unit rows: 293,388
D/U/N ledger survivors: 27,733
strict pair-sum certificates: 27,730
floor rows: 3
below rows: 0
Res_27 proof-obligation types among survivors: 148
```

Exact maximin checking confirms that AP, `V*`, and `2*AP` all have
`M=1/14`.

## Interpretation

This is not the full lifted `n=14` theorem.  The D and N gates here are applied
to least-positive representatives, and arbitrary integer lifts of the same
`Res_27` shell data can change divisibility and `n`-clock behavior.

The result is still useful because it separates the problem:

```text
least-positive Res_27 coimage branch: classified;
remaining theorem: lift/CRT conservativity.
```

In the coimage/Yoneda language of HYP-2161, the first finite shell quotient has
no hidden primitive floor family beyond AP and `V*`.  Any future counterexample
or unclassified tight branch must use lift data, endpoint-owner data, or CRT
coupling not present in this base quotient.

After rebase, this also complements the two incoming n=14 improvements.
HYP-2162/THM-407 proves the `C=27` shell face folds from `13` shells to `3`
prime-3 strata under `G=<2,-1>`.  HYP-2163 packages the no-multiple clock
witness and the remaining multiple/Cprime pipeline.  HYP-2164 does not quotient
by those reductions; it exhausts the least-positive row layer and shows that the
full base quotient still has no primitive floor row beyond AP and `V*`.

## Tournament Analysis

The tournament vertices are not runners.  They are `Res_27` proof-obligation
types:

```text
(missed nonunit shells, doubled shells, primitive?)
```

The pairwise observable is the type burden

```text
(below_count, floor_count, -strict_count, label).
```

The resulting tournament on `148` type vertices is transitive:

```text
SCCs: 148 singleton components
directed 3-cycles: 0
Hamiltonian paths: 1 under the burden/tie gauge
```

This quotient preserves the finite predicate "canonical `Res_27` row survives
D/U/N and pair-sum floor testing."  It destroys lift choices, actual integer
speed sizes, and phase order.

## Next Problem

Prove the lift/CRT conservativity theorem:

```text
If a lifted integer row over a classified Res_27 shell type remains at the
floor or below, then its lift data routes to AP, V*, a known multiple/Cprime
positive-measure exit, or a bounded CRT contradiction.
```

That is the sharpened `n=14` target after S608.
