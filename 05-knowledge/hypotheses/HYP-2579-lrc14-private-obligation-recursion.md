# HYP-2579 - LRC14 covering residuals should reduce to private-obligation recursion

**Status:** OPEN proof program  
**Source:** codex-2026-06-17-S5  
**Computation:** `04-computation/lrc14_covering_residual_classifier_codex.py`  
**Stored output:** `05-knowledge/results/lrc14_covering_residual_classifier_codex.out`

## Claim

After THM-523, THM-524, and THM-526, the remaining LRC(14) gap-side proof
should be organized around **private q-obligations** of parked runners.

Let `S` be a primitive q-covering 13-set, so every `q in {2,...,14}` divides at
least one speed.  Let `w=14m` be a parked runner and write `A=S\{w}`.  THM-526
certifies `S` when the widest level-`1/14` safe arc of `A` has width
`W(A)>1/(98m)`.  The residual should have this shape:

1. `w` is small enough that arc-width does not certify it;
2. deleting `w` uncovers at least one private q-obligation, meaning some
   `q<=14` divides `w` and no other speed in `S`;
3. the exact THM-524 optimum is controlled by a binding crossing whose flank
   and crossing index are determined by those private obligations.

Equivalently, if no parked runner has private q-obligation, then the row should
either be discharged by arc-width, or admit a recursive deletion/descent to a
smaller parked-runner obligation packet rather than form a genuinely new hard
case.

## Evidence

The classifier sampled `103` primitive q-covering 13-sets from principal
single-drop towers, small two-drop/one-replacement cores plus a parked runner,
and seeded random covering rows.

Exact results:

| statistic | value |
|---|---:|
| sampled q-covering primitive rows | 103 |
| THM-526 arc-width certified rows | 94 |
| arc-width residual rows | 9 |
| LRC breaks `M<1/14` | 0 |
| arc residuals with no parked private obligation | 0 |

The closest row in the principal single-drop q-covering tower is not the older
unit-residue-complete `{1,...,11,13,84}` example, but

`S={1,2,3,4,5,6,7,8,9,10,11,12,182}`,

with

`M(S)=14/183`, slack `14/183 - 1/14 = 13/2562`,

attained at the sum-binding pair `(1,182)` with `D=183`, `j=14`, and
`14j-D=13`.  This row is q-covering, but it is missing unit residue `13` modulo
14.  Thus the older `7/89` champion remains the natural champion for the
stronger unit-residue-complete principal family, while the q-covering floor has
a different visible pressure point.

The principal table from the run:

| dropped `q` | best inserted `w` | best `M` | missing unit residues |
|---:|---:|---:|---|
| 12 | 84 | `7/89` | none |
| 13 | 182 | `14/183` | `(13,)` |
| 11 | 154 | `14/157` | `(11,)` |
| 9 | 126 | `14/131` | `(9,)` |

This separates two notions that had been blending together:

- **q-covering:** the THM-523 necessary condition for a counterexample;
- **unit-residue completeness:** a stronger mod-14 section condition that is
  useful but not equivalent.

## Tournament Analysis

This hypothesis challenges the default runner-vertex model.

- Vertices: sampled q-covering rows.
- Pairwise observable: exact hardness by `M(S)-1/14`, with lower slack harder.
- Switch/gauge: orient the harder row toward the easier row.
- Tie Hamiltonian path: lexicographic order on the speed tuple.
- Comparison gauge: private/blocking-height proxy
  `fragility=sum_q 1/count_q`, number of private parked obligations, arc-width
  margin, then unit-residue coverage.

Fingerprints on the `96` hardest sampled rows:

- exact hardness tournament score histogram is `0,1,...,95`;
- directed 3-cycles: `0`;
- pressure proxy agreement with exact hardness: `2963/4560 = 0.650`;
- edge flips versus pressure proxy: `1597/4560 = 0.350`.

Readout: exact hardness is scalar/transitive once `M` is known, so the useful
information is the flip rate.  Blocking pressure helps but is not the invariant.
The missing invariant is more local: which parked runner is uniquely paying
which q-obligation, and which binding flank that private obligation forces.

## Assumption Challenge

Considered vertex sets:

- runners,
- q-obligations,
- fixed sections,
- cover arcs,
- endpoint events,
- binding pairs,
- safe components,
- proof obligations.

Chosen quotient: parked-runner private q-obligations plus binding crossings.  It
preserves the THM-523 q-covering predicate, the THM-526 large-runner discharge,
and the THM-524 exact LRC predicate at the residual crossing.  It destroys most
global timing away from the crossing and loses unit-residue section data unless
that data is explicitly recorded as a side channel.

The challenged assumption is that "more divisibility pressure" is the right
hardness order.  The tournament flip rate shows that global blocking height is
too coarse.  Private obligation at the parked runner is the sharper pressure
coordinate.

## Proof Target

Prove a private-obligation crossing lemma:

> Let `S=A union {14m}` be a q-covering LRC14 row not certified by THM-526.  If
> `14m` has a private q-obligation, then the THM-524 binding crossing forced by
> that private obligation has crossing index `j` satisfying `j >= D/14`, where
> `D` is the binding denominator.

In principal single-drop families this reduces to the already observed
closed-form/monotone-collapse law.  The next step is to prove the analogous
flank law for mixed cores, or show that mixed cores recurse until they enter
the principal/private packet.

Cross-links: THM-523, THM-524, THM-526, HYP-2577, HYP-2578, HYP-2571, OPEN-Q-108,
T843, T844.
