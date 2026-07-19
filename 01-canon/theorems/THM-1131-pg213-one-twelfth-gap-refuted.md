---
id: THM-1131
title: The PG(2,13) parameter identity is real, but the proposed one-twelfth loneliness gap is false
status: PROVED CORRECTION.  Exact primitive Covering rows occur at M=1/13, at M=3/37 strictly between 1/13 and 1/12, and at M=28/365 strictly below 1/13.  Moreover two primitive Covering integer lifts with the same residue subset modulo 183 have different global maxima.  Thus the proposed one-twelfth gap, global one-thirteenth floor, and PG quotient transport are false.  The parameter identities survive only as an analogy.
source: codex-2026-07-18-S67 (audit of boxeph-S110 / HYP-7635)
depends_on:
  - THM-668-pair-sum-ruler-witness-structure  # exact pair-sum peak enumeration (not the colliding detuned-harmonic file)
  - THM-724    # single-killer covering-min scope
  - THM-737    # prior exact compressed-near-dilate value
related: [HYP-7635, HYP-7604, MISTAKE-141, MISTAKE-166]
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

## The strict sub-`1/13` branch is nonempty

The second deep-well rung is

```text
U = {1,...,12,364},             M(U)=28/365<1/13.       (3)
```

It is primitive and `Covering`, and its twelve-speed core is the same AP as the
first rung.  Thus a claim that every primitive Covering family has
`M>=1/13` is also false.  What survives is the narrower structural question:
after the already-lonely and normalization branches are typed correctly, does
strict `M<1/13` force the AP/dominance mechanism?  Equation (3) is compatible
with that question; it refutes the purported empty gap, not the AP-core route.

## The `Z/183` quotient does not preserve `M`

Let

```text
B  = {1,2,...,12,182},
B' = {2,3,...,12,182,184}.
```

Both rows are primitive and `Covering`, and their residue subsets modulo 183
are identical because `184=1 (mod 183)`.  Nevertheless exact pair-sum replay
gives

```text
M(B)=14/183,                    M(B')=13/93.             (4)
```

Therefore even the complete unlabelled residue support in `Z/183` does not
preserve global loneliness.  It records the 183-grid phases but forgets the
integer lift, all other rational charts, and hence the optimization over every
real time.  A PG/Singer quotient cannot support an inverse theorem until that
forgotten fibre is controlled.

## What survives, and what INV actually asks for

The `PG(2,13)` parameter coincidence remains a potentially suggestive analogy.
But the S110 computation evaluated only five hand-picked cores, several of
whose completed rows were not `Covering`; it did not enumerate non-AP covering
families.  It also did not construct a map preserving `Lonely 13`, dominance,
or the `Covering` predicate between runner families, projective-plane objects,
Singer difference sets, and tournament classes.

The surviving research target in this branch must be stated after its actual
domain is fixed—for example primitive rows in the strict `M<1/13` structural
branch with Covering rederived after gcd reduction.  The first formal capstone
asked universally for no-`Lonely 13` to imply dominance and was refuted by
`{1,...,13}`.  Adding Covering produced the historical premise

```text
Covering(2..14) v and not exists t, Lonely 13 v t
  => some speed is at least 13 times every other speed.       (INVcov)
```

over positive rows.  Although `{1,...,13}` misses modulus 14, THM-1158 shows
that its dilation `2*{1,...,13}` covers every modulus through 14 and still has
maximum `1/14` and no dominant speed.  Thus `INVcov` itself is refuted.  The
Lean surface also exposes the exact outer no-`Lonely 14` residual explicitly:

```text
Covering(2..14) v and not exists t, Lonely 14 v t
  => some speed is at least 13 times every other speed.       (ResidualINV)
```

The conditional implication from `INVcov` remains valid but has a false
premise.  With the AP bridge already assumed by the formal route,
`ResidualINV` is logically equivalent to the
LRC(14) target, rather than a reduction to the known twelve-runner theorem.
Neither formal target is equivalent to the false stronger assertion
`non-AP => M>=1/12`.  A viable stability reformulation must respect the `1/13`
boundary, the AP tower below it, and the integer-lift fibre.

## Tournament Analysis

Ordering the five referee rows by their exact `M` values gives a transitive
five-vertex tournament: score histogram `0,1,2,3,4`, no directed cycle,
singleton SCCs, and one Hamiltonian path.  This
preserves only scalar `M` order.  It destroys divisor carriers, exact residues,
maximizer times, and the dominance conclusion, so it cannot transport
transitive-class isolation to `INVcov`.  Divisor obligations or exact maximizer
cells are more faithful candidate vertices than runners or analogy classes.

The suggested Singer/tournament identification also fails at the parameter
level.  A Singer `(183,14,1)` set has 14 elements; a tournament connection set
on 183 vertices must be a skew half-set of size 91, and a doubly regular
tournament there has parameters `(183,91,45)`.  A regular tournament on 14
vertices is impossible.  The stationary-augmented deep-well residue set
`A={-1,0,...,12}` has `|A-A|=27` and additive energy 1834, while a Singer set
has `|D-D|=183` and energy 378.  These are contrasting additive objects, not
the two poles of one proved LRC-preserving tournament map.

The correction does not refute a genuinely primitive inverse problem or
LRC(14).  It refutes the literal covering-`M<1/13` dominance proposition and
shows that gcd/divisor fibres are part of the object.  The no-`Lonely 14`
qualifier types the exact `ResidualINV` target (MISTAKE-170).

The dependency-free referee uses exact `Fraction` arithmetic and explicit
optimization-stable checks.  Normal and `python -O` runs match the frozen
output byte-for-byte.  SHA-256:

```text
source  aa364b067c8ce25100ae24d0fa4e728ec5d56dbf59dc8da092314a7f08bf2dd6
output  ae8fd76c6449f7100829dcede6b799e404a3b6d9a73574ef405fc374e5efcd41
```

∎
