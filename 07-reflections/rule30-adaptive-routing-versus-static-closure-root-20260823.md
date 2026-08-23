# Rule 30 adaptive routing versus static closure

**Status (2026-08-23): PROVED ambient static-closure theorem plus FINITE-EXACT
adaptive audit.**  The bounded routed-observer computation is now paired with
an all-depth theorem: the least literal uniform two-stage point bank contains
the whole odd half.  This still proves no bounded first-return gap, physical
signalizer recurrence, finite graph, center-column statement, balance
statement, or prediction lower bound.

## Inheritance and question

THM-3511 gives the exact first-return section formula.  The incoming bounded
observer signal says that a depth-`D+B` portrait either recovers a first gap
`d<=B` and successor portrait or reports honest overflow.  The unresolved
operation question was whether the next selected bank could be compiled from
a small fixed extension of the current bank.

The canonical hostile is the pair `CA,CCAC`: both have the same depth-four
selected bank `(7,10,1)`, but their gaps are `3,10` and their successor banks
differ.  Thus a selected bank is not closed under first return.

## Exact adaptive compiler

Let `pi` be the source permutation at depth `H=D+B`, and let the selected
rays be `r=0,1,2`.  Once the zero chain reveals `d<=B`, the inherited
section identity is

```text
next(r) = (pi(pi(r 2^d)) >> d) mod 2^D.                         (1)
```

At `D=B=4`, the current bank plus the three routed chains

```text
C_0=(0,pi(0),pi^2(0)),
C_1=(2^d,pi(2^d),pi^2(2^d)),
C_2=(2^(d+1),pi(2^(d+1)),pi^2(2^(d+1)))                        (2)
```

therefore recovers the exact next bank.  The zero chain supplies the gap and
new ray zero; the other two chains supply rays one and two.  After the current
bank is known this is at most six additional point evaluations, with fewer
when routed addresses collide.

Exhausting active words of lengths one through nine gives `14,762` distinct
depth-eight portraits: `13,853` exact branches and `909` honest overflows.
Equation (1) agrees with direct renormalization on every exact branch.  Two
explicit collision pairs show that deleting either `C_1` or `C_2` loses the
corresponding next output even after retaining the current bank, gap, zero
chain, and other off-ray chain.  This is minimality only in that declared
whole-chain deletion family, not among arbitrary symbolic compressions.

## Why static compilation explodes

Before the gap is known, the possible first-stage inputs are

```text
X_B={0,2,4,...,2^(B+1)}.                                       (3)
```

A literal nonadaptive two-stage point-query bank must retain `pi` on `X_B`
and on every possible routed address `pi(x)`.  The wreath recursion contains
the self-replication certificate

```text
A|_0=A,                    (B^(-1)AB)|_0=B.
```

Hence the root stabilizer surjects by sections onto `<A,B>`, so this group is
level-transitive.  On each finite level the positive monoid has the same
orbits, and active words send every even ray onto exactly the entire odd half.
Therefore, for every `D>=2,B>=1`, the literal static closure is

```text
S_B=X_B union {all odd depth-(D+B) rays},
|S_B|=2^(D+B-1)+B+2.                                           (4)
```

The earlier finite sizes `19,36,69,134,263,520` are instances of this theorem.
At the central test
`B=4,H=8`, static closure consumes `134` of `256` coordinates.  The routed
observer avoids this cost because it asks the second query only after its
state-dependent address is known.

This is the useful connection to the operation-generated Rule 30 tariff:
closure is cheap in an adaptive query category but expensive in a frozen
spatial-coordinate category.  The preserved predicate is exact selected-ray
first return; the destroyed information is the rest of the portrait; the
indispensable sidecar is the routed intermediate address `pi(x)`.

## Boundary and next pull

The all-depth proof ranges over ambient active words, not the physical
signalizer path, and its minimality is only for literal point banks.  It does
not show that the physical path realizes either hostile pair or that two
chains are optimal against a richer Mealy-state representation.  The first
live question remains whether `pi(x)` itself has a bounded recursive
carry/section description along the physical path.  Such a recursion could
retain adaptive closure without storing half a portrait.

Exact artifacts:

```text
04-computation/rule30_adaptive_routed_query_closure_20260823.py
05-knowledge/results/rule30_adaptive_routed_query_closure_20260823.out
script SHA256  ba00dbf59505ccb1373de0fde4a456e7d839ace78def29ffb27642a1f6b8b849
output SHA256  bfbfef19cd3194e2941798e42fbffed2956fb7d24bcc6471e4f27af5a5641b1f
semantic SHA256 0608bc5298aa7855dabb4784dd56df3c2df6cba06a7bfe0e2a60b754f1a3ca63
```

Normal and optimized runs byte-match the frozen LF output.

The universal addendum and its import-free hostile audit are canonical at

```text
04-computation/rule30_universal_static_closure_thm3511.py
05-knowledge/results/rule30_universal_static_closure_thm3511.out
04-computation/rule30_universal_static_closure_independent_audit_thm3511.py
05-knowledge/results/rule30_universal_static_closure_independent_audit_thm3511.out
```
