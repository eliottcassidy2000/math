---
id: THM-1131
title: The PG(2,13) parameter identity is real, but the proposed one-twelfth loneliness gap is false
status: PROVED CORRECTION.  The primitive Covering family 2*{1,...,12} union {13} has exact M=1/13, by an elementary matching upper and lower bound.  A second primitive Covering family has exact M=3/37, strictly between 1/13 and 1/12.  Therefore neither the claimed empty interval [14/183,1/12) nor the target non-AP implies M>=1/12 can be used for INV.  The identities 183=13^2+13+1 and 14=13+1 survive as a parameter coincidence; no predicate-preserving map to PG(2,13), a Singer difference set, or the tournament metagraph has been supplied.
source: codex-2026-07-18-S67 (audit of boxeph-S110 / HYP-7635)
depends_on:
  - THM-668    # exact pair-sum peak enumeration
  - THM-724    # single-killer covering-min scope
  - THM-737    # prior exact compressed-near-dilate value
related: [HYP-7635, HYP-7604, MISTAKE-141]
script: 04-computation/lrc14_pg213_gap_referee_codex_S67.py
output: 05-knowledge/results/lrc14_pg213_gap_referee_codex_S67.out
---

# THM-1131 — the proposed `1/12` PG(2,13) gap is refuted

The numerical identities

```text
183 = 13^2 + 13 + 1,       14 = 13 + 1
```

are exact.  They identify the deep-well denominator and the LRC band parameter
with the point count and line size of `PG(2,13)`.  What they do **not** prove is
that the LRC family is a projective plane, a Singer difference set, or a vertex
of the tournament metagraph.  More decisively, the proposed quantitative gap is
already contradicted by canonical exact families in this repository.

## An elementary endpoint counterexample

Let

```text
V = {2,4,6,...,24} union {13}.
```

It has thirteen distinct positive speeds and is primitive.  It is `Covering`:
for every `2<=q<=12`, the speed `2q` is present; `13` carries modulus 13; and
`14=2*7` carries modulus 14.

The twelve-speed block `2*{1,...,12}` has loneliness value `1/13`, so adding
the speed 13 gives `M(V)<=1/13`.  At `t=1/26`,

```text
||2i*t|| = ||i/13|| >= 1/13      (1<=i<=12),
||13*t|| = 1/2.
```

Hence `M(V)>=1/13`, and therefore

```text
M(V)=1/13 < 1/12.                                      (1)
```

This row is not the deep-well AP family.  Equation (1) alone refutes both
“every non-AP covering row has `M>=1/12`” and the assertion that
`[14/183,1/12)` contains only the deep-well value.

## A counterexample in the open interval

The exact pair-sum peak replay also checks

```text
W = {2,3,5,8,9,11,12,13,14,15,17,20,23},
M(W)=3/37,
```

where `W` is again primitive and `Covering`.  Since

```text
1/13 < 3/37 < 1/12,                                    (2)
```

even loose non-AP rows populate the alleged gap.  The exact value in (2) was
already present in the bounded divisor-complete census; the new referee
recomputes it by testing every numerator at every THM-668 pair-sum denominator.

## What survives, and what INV actually asks for

The `PG(2,13)` parameter coincidence remains a potentially suggestive analogy.
But the S110 computation evaluated only five hand-picked cores, several of
whose completed rows were not `Covering`; it did not enumerate non-AP covering
families.  It also did not construct a map preserving `Lonely 13`, dominance,
or the `Covering` predicate between runner families, projective-plane objects,
Singer difference sets, and tournament classes.

The kernel capstone asks only for

```text
not exists t, Lonely 13 v t
  => some speed is at least 13 times every other speed.       (INV)
```

That statement is not equivalent to the false stronger target
`non-AP => M>=1/12`.  A viable stability reformulation must respect the
`1/13` boundary families and must specify which normal form and predicate its
quotient preserves.

## Tournament Analysis

Ordering the deep well, the compressed near-dilate, and `W` by their exact
`M` values gives a transitive three-vertex tournament: score histogram
`0,1,2`, no directed cycle, singleton SCCs, and one Hamiltonian path.  This
preserves only scalar `M` order.  It destroys divisor carriers, exact residues,
maximizer times, and the dominance conclusion, so it cannot transport
transitive-class isolation to INV.  Divisor obligations or exact maximizer
cells are more faithful candidate vertices than runners or analogy classes.

The correction does not refute INV or LRC(14).  It removes one false proposed
strengthening and leaves the genuine `M<1/13` dominance problem open.

The dependency-free referee uses exact `Fraction` arithmetic and explicit
optimization-stable checks.  Normal and `python -O` runs match the frozen
output byte-for-byte.  SHA-256:

```text
source  24a10aec74c3b5f8e238967311690df92d20bf18cf19c0b86863357417509569
output  2ea698e874fc1ec6960ea6a1c9bd2e0631522af1e3fadcbe53eec326e2939b47
```

∎
