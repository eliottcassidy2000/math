# Rule 30 left-permutive trace torsor and the seed boundary

**Research synthesis, 2026-08-15.**  This is the proof-and-loss ledger for the
proved and independently audited
[THM-3456](../01-canon/theorems/THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary.md)
package and its exact companion.  The current official page was checked on 2026-08-15: it still
lists all three prizes and accepts submissions until a satisfactory solution
is achieved, so the problems are treated here as open.

External anchors are Stephen Wolfram's 2019
[Rule 30 prize announcement](https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/),
the [official prize page](https://www.rule30prize.org/), and Chernikov's
[SOP2=SOP3, arXiv:2608.13291v1](https://arxiv.org/abs/2608.13291).

## Outcome

Every binary left-permutive radius-one cellular automaton has the same free
one-column trace space: after fixing the initial right half-line, varying the
initial left half-line produces every infinite binary temporal trace exactly
once.  At finite depth this is a triangular bijection

```text
(x_-1,...,x_-n)  <-->  (X_1(0),...,X_n(0)).
```

This exact universality is an obstruction, not a randomness theorem.  Rule 30
and Rule 60 have identical free trace space, while Rule 60's distinguished
single-seed center is the constant word `111...` and Rule 30's is not.  The
free trace quotient has discarded the selected initial ancestry.

Two adjacent temporal columns repair a different loss on a finite spatial
cylinder: Rule 30's sideways inverse reconstructs the whole spacetime, so the
adjacent-pair trace period equals the global orbit period.  One column does
not.  Width seven is the first nonconstant finite-cylinder alias after an
exhaustive width-three-through-seven audit.

Finally, a width-`n` single-seed cylinder agrees with the infinite single-seed
center only for the guaranteed prefix `0<=t<n`.  Every finite cylinder is
eventually periodic; that periodic tail is created after the periodic copies
of the seed can enter the center light cone.  It is not an infinite-seed tail.

## Inheritance pass

| role | inherited object |
|---|---|
| closest mechanism | Rule 30's classical sideways/left-permutive law |
| canonical hostile | Rule 60: the same trace torsor but constant single-seed center |
| corrected near miss | finite-cylinder period data treated as evidence about an infinite single-seed tail |
| least-used sidecar | the chosen initial left wing and the spatial boundary/lift |
| method cards | `Separate observer type, recurrence class, and finite head`; `Controlled forgetting and unlabeled quotients require a sidecar` |

There was no existing Rule 30 theorem, computation, or dedicated reflection in
the repository at startup.  During this session the independent sibling
[Rule 30 depth-observer no-go](lrc14-rule30-depth-observer-typed-phase-no-go-codex-20260815.md)
arrived.  It proves that sampling the seed trace at the `165` LRC valuation
depths collapses to four fibres and that an additive binary observer cannot
carry `F_13^2` phase.  That result and the trace torsor have the same loss
diagnosis—binary observation forgets ancestry/typed phase—but no LRC map is
created by that analogy.  The nearby historical `random031` files concern an
LRC row, not elementary cellular automaton Rule 30.

## 1. Exact algebra and the triangular trace torsor

Write the binary local rule as `f(l,c,r)`.  Rule 30 has

```text
f_30(l,c,r) = l XOR (c OR r).                         (1)
```

A binary rule is left-permutive when

```text
f(0,c,r) != f(1,c,r)       for every (c,r).
```

Equivalently, for a unique Boolean function `g`,

```text
f(l,c,r) = l XOR g(c,r).                              (2)
```

There are exactly sixteen such elementary rules:

```text
15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240.
```

Let `X_t(i)` be the spacetime value at time `t` and site `i`, with
`X_0(i)=x_i`.  Causality and (2) give, for every `t>=1`,

```text
X_t(0) = x_-t XOR H_t(x_(-t+1),...,x_t)               (3)
```

for a Boolean function `H_t` depending on the rule.  The proof is induction.
For `t=1`, this is (2).  For the next step,

```text
X_(t+1)(0) = X_t(-1) XOR g(X_t(0),X_t(1)).
```

The translated induction hypothesis puts `x_(-(t+1))` with coefficient one
in `X_t(-1)`, while both other terms depend only on initial cells with index at
least `-t`.  This proves (3).

Fix `x_0,...,x_n`.  Equation (3) is triangular in the order

```text
x_-1, x_-2, ..., x_-n.
```

Given a desired trace `y_1,...,y_n`, recover `x_-t` uniquely from `y_t` after
the previous `t-1` left cells have been recovered.  Thus

```text
Phi_n : {0,1}^n -> {0,1}^n,
Phi_n(x_-1,...,x_-n)=(X_1(0),...,X_n(0))              (4)
```

is a bijection for every `n`, every fixed right block, and every binary
left-permutive rule.  The maps commute with truncation.  Passing to the inverse
limit makes `Phi` a homeomorphism from the free left half-line to the one-sided
full binary trace shift.  It is a triangular homeomorphism, not a claim that
the spatial CA is conjugate to the temporal shift.

This is naturally a torsor statement: fixing a desired trace selects exactly
one left wing relative to the marked right half-line.  Forgetting the wing
makes every trace merely “possible” and destroys the distinguished seed.

## 2. Rule 30 versus Rule 60: the ancestry hostile

Rule 60 is also left-permutive:

```text
f_60(l,c,r)=l XOR c.                                  (5)
```

Start from `x_0=1` and `x_i=0` for `i!=0`.  Every site strictly left of zero
remains zero under (5), so

```text
X_(t+1)(0)=0 XOR X_t(0)=1.
```

Its center trace is `111...`.  Rule 30 begins

```text
11011100110001011001001110101110...
```

Nevertheless, (4) says that after freeing the left wing, both rules realize
every finite word once and every infinite trace once.  Therefore none of the
following can distinguish the selected Rule 30 seed from the Rule 60 hostile:

- the set of attainable trace prefixes;
- the number of left wings realizing each prefix;
- the abstract prefix tree or its branching statistics; or
- “balance across all free inputs.”

The first lost coordinate is the actual left wing, not a subtle statistic of
the trace language.  The repair is to retain the initial-state selector and
compatibility across the chosen orbit.

## 3. The SOP2 tree is real but generic

For a finite word `sigma`, fix the right half-line and define the trace
cylinder

```text
C_sigma = {left wings whose center trace begins with sigma}.
```

The triangular homeomorphism proves:

```text
C_sigma is nonempty;
C_tau subset C_sigma when sigma is a prefix of tau;
C_sigma intersect C_tau is empty when sigma,tau are incomparable.
```

Thus the cylinders form a literal full binary prefix tree as a set system.
If one explicitly takes a two-sorted first-order structure containing left
wings, finite binary words, and one uniform incidence relation

```text
R(x,sigma)  iff  x lies in C_sigma,
```

then `R` has SOP2: a branch is jointly consistent and incomparable nodes are
pairwise inconsistent.  Chernikov's Definition 2.1 gives precisely this
pattern (using `omega^(<omega)`, equivalently the binary tree), and his Theorem
1.1 proves `SOP2 => SOP3` for complete first-order theories.

This is a typed bridge only for the stipulated incidence structure.  The bare
cellular automaton language has not been shown to uniformly define arbitrary
finite trace-prefix parameters, and no canonical SOP3 formula for a natural
Rule 30 structure is produced here.  More importantly for the prizes, the
same cylinder tree occurs for Rule 60 and all fourteen other left-permutive
elementary rules.  The tree records input freedom, not selected-seed
randomness, density, or computational irreducibility.

## 4. Two adjacent traces reconstruct a Rule 30 cylinder

Equation (1) can be solved sideways:

```text
X_t(i-1) = X_(t+1)(i) XOR (X_t(i) OR X_t(i+1)).       (6)
```

Knowing the complete temporal columns at `i` and `i+1` therefore determines
the column at `i-1`; iterating (6) determines every column to the left.  On a
spatial cylinder `Z/nZ`, after at most `n-1` iterations this is the whole
spacetime.

Consequently, on a temporally periodic cylinder orbit, the minimal period of
any adjacent-pair trace equals the minimal period of the full orbit.  One
divisibility is immediate because the trace is an observation of the orbit.
For the reverse divisibility, a period of the pair propagates through (6) to
every column, so it is a period of the full state.

One column has neither property.  With words written in labelled order
`x_0 x_1 ... x_6`, the three width-seven states

```text
A = 0100000,
B = 0111000,
C = 0000111
```

lie on period-four Rule 30 orbits and all have the same labelled-site-zero
trace

```text
(0110)^infinity.                                      (7)
```

Their full cycles are

```text
A: 0100000 -> 1110000 -> 1001001 -> 0111111 -> A,
B: 0111000 -> 1100100 -> 1011111 -> 0010000 -> B,
C: 0000111 -> 1001100 -> 1111011 -> 0000010 -> C.
```

The separate width-seven state `0100110` has full period four but center-trace
period two, `(01)^infinity`; its adjacent-pair trace again has period four.
Exhausting every periodic state finds no nonconstant one-column alias for
widths `3,4,5,6` and ten alias classes at width seven.  Widths one and two are
excluded from the minimality claim because cyclic neighbor identifications
degenerate the radius-one neighborhood.

This is the sharp finite-cylinder loss ledger:

```text
source:       full labelled periodic spacetime
target:       one labelled temporal column
map:          restriction to site zero
preserved:    that column's exact infinite temporal word
destroyed:    neighboring column and global orbit phase/state
sidecar:      either adjacent temporal column
hostile:      (7), plus 0100110 for period loss.
```

## 5. A finite cylinder gives a prefix, not an infinite tail

The width-`n` cylinder with one black cell at residue zero lifts to the
infinite periodic initial configuration

```text
x_i = 1 iff i is in nZ.
```

Compare it with the actual infinite single seed, whose only black cell is at
zero.  A radius-one center value at time `t` depends only on initial sites in
`[-t,t]`.  For `0<=t<n`, that interval contains no nonzero multiple of `n`
other than zero.  Hence the two center traces agree exactly for

```text
0 <= t < n.                                           (8)
```

The bound cannot be uniformly extended to `t=n`: the width-three cylinder
first differs from the infinite seed exactly at `t=3`.  Width eight is a useful
cancellation control—its first difference is delayed to `t=10`, showing that
boundary arrival need not immediately change the observed bit.

Every width-`n` cylinder map acts on only `2^n` states, so every orbit is
eventually periodic and its transient plus period is at most `2^n`.  This
finite-state fact begins precisely where (8) no longer protects the selected
infinite seed.  Increasing `n` extends the certified prefix, but the eventual
cycles live in different periodic-lift systems and do not form an inherited
infinite tail.

The exact single-seed audit through width 24 happens to find center-trace
period equal to full orbit period for every width `4..24`, matching the
finite-size phenomenon described in the prize announcement.  It is
`FINITE-EXACT`, not an all-width theorem and not evidence that can cross (8).

## 6. Keep the three prize problems separate

The official prize page asks about the infinite center sequence generated from
one nonzero cell.

### Problem 1: absence of eventual periodicity

The adjacent-column theorem cannot be applied because the prize supplies one
column, and the width-seven alias proves that one column is not a faithful
finite-cylinder observable.  Finite cylinders are necessarily eventually
periodic and agree with the infinite seed only through (8).  No result
excluding eventual periodicity follows.

### Problem 2: limiting frequency one half

The full free trace language contains every binary word with exact uniform
finite multiplicity one, but that balance is across initial conditions.  The
prize asks for a time average along one selected initial condition.  Rule 60's
constant seed trace is the decisive countercontrol.  No density or even
density-existence statement follows.

### Problem 3: an `Omega(n)` computation lower bound

The triangular inverse is an existence and reconstruction statement when a
desired trace is supplied and the left wing is free.  It is not an algorithm
for the Rule 30 seed bit and supplies no lower bound in a specified machine
model.  Conversely, finite-cylinder lookup or cycle acceleration changes the
initial system after the `n`-step prefix boundary.  No complexity lower bound
or sublinear algorithm follows.

## Concept board and connection ledger

| object | representation | predicate | operation | lost coordinate | cheapest test |
|---|---|---|---|---|---|
| Rule 30 local law | Boolean polynomial/formula | left permutivity | solve sideways | right neighbor if only one trace is kept | width-seven alias |
| free trace torsor | triangular inverse system | every prefix once | forget left wing | seed ancestry | Rule 60 |
| trace cylinders | binary prefix tree | branch consistency/incomparability | first-order coding | natural definability | explicit two-sorted `R` |
| finite cylinder | functional graph on `2^n` states | eventual period | periodic spatial quotient | remote-seed boundary | width three at `t=3` |
| selected seed | one infinite orbit | the three prize properties | take prefixes/limits | none may be discarded | preserve all-zero wings |

The exact connection contracts are:

| source -> target | map | preserved | destroyed | sidecar | decisive hostile |
|---|---|---|---|---|---|
| left wings -> traces | triangular `Phi` | full prefix and inverse-limit topology | distinguished wing after existential projection | initial ancestry | Rule 30/60 |
| trace cylinders -> SOP tree | incidence coding | nested/incomparable set pattern | natural-language definability | explicit word sort and `R` | same tree for Rule 60 |
| two columns -> cylinder | sideways recursion (6) | full spacetime and period | nothing on a labelled finite ring | both adjacent traces | delete either column |
| periodic lift -> infinite seed | center restriction before boundary arrival | times `0..n-1` | all later ancestry | spatial lift/boundary label | width three |

## Exact companion and controls

The stdlib-only companion is

```text
python3 04-computation/left_permutive_trace_torsor_rule30_thm3456.py
python3 -O 04-computation/left_permutive_trace_torsor_rule30_thm3456.py
```

It freezes:

- the truth tables and complete list of sixteen binary left-permutive ECA;
- every fixed right block and every left wing through depth six for all
  sixteen rules;
- every fixed right block and left wing through depth eight for Rules 30 and
  60;
- `502,400` forward trace evaluations and `5,568` constructive inversions;
- every periodic state at widths `3..7`, including alias minimality;
- adjacent-pair/full-period equality on that periodic-state universe;
- the cylinder/infinite prefix comparison for widths `1..64`; and
- exact single-seed orbit data through width `24`.

Normal, optimized, and stored transcripts must be byte-identical before the
package is promoted.  Hashes are pinned in the matching result file after the
final replay.

## Honest frontier

The useful theorem is a no-go with a repair: free trace languages and abstract
trace-cylinder trees are maximally rich for a trivial algebraic reason across
all left-permutive rules.  Any route to Rule 30's prize sequence must retain
the selected seed ancestry or an equivalent orbit selector.  The next honest
finite question is whether one can find an ancestry-preserving observable
strictly weaker than two full adjacent columns yet strong enough to rule out
eventual one-column periodicity.  The present work does not supply it.
