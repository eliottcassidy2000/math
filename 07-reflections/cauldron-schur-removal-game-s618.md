# Cauldron Schur Removal Game (S618)

The user-defined cauldron game starts with natural numbers placed one at a
time into `k` cauldrons. A cauldron boils when it contains values satisfying an
additive equation, with the base rule `A+B=C`.

The first correction is definitional. Since each natural number is placed once,
the literal "three values `A,B,C`" rule uses distinct summands. That is the
weak Schur convention, not the classical Schur convention that allows `A=A`.
The script therefore reports both:

- literal distinct two-term, `k=3`: safe through `23`, forced at `24`;
- repeated-summand comparison, `k=3`: safe through `13`, forced at `14`.

The removal version is not the same problem. Once a cauldron boils and is
removed, the player can treat that color as a sacrificed resource. This turns
the question from "does a global sum-free coloring exist on `[1,n]`?" into
"how long can active sum-free resources be staged before every cauldron has
been spent?"

For the literal distinct two-term rule, exact active-state search gives last
boils at `3,10,27` for `k=1,2,3`. For the repeated-summand comparison, the
last boils are `2,7,20`. The optimal `k=3` literal strategy waits until boils
at `24`, `26`, and `27`; the intermediate placement at `25` is the
sacrifice-game move that is invisible in a first-boil Schur number.

The "more terms" variants tighten the first-boil base game in the computed
prefix: distinct two-or-three-term cauldrons force at `23`, and distinct
finite-sums cauldrons force at `22`. They also tighten the all-boiled removal
game: finite-sums removal has exact last boils `3,10,25` for `k=1,2,3`, so
the three-cauldron removal game loses two steps relative to the weak two-term
rule's `27`.

Literature anchor: weak Schur tables use the largest safe-prefix convention,
with `WS(3)=23`; generalized Schur/weak-Schur numbers for more summands are
nontrivial and should not be silently inferred from the two-term scout. See
`https://www.sciencedirect.com/science/article/pii/S0898122111009722` for the
weak Schur convention and `https://www.sciencedirect.com/science/article/pii/S0166218X18303846`
for the repeated-summand three-color generalized Schur formula.

Artifacts: `04-computation/cauldron_game_s618.py` and
`05-knowledge/results/cauldron_game_s618.out`.

## Assumption Challenge

Do not assume the vertices of the Tournament Analysis are the cauldrons. Other
candidate vertices are active-state signatures, residue classes, subset-sum
obligations, first-boil witnesses, removal order types, finite-sum arities, and
proof routes. The quotient should preserve the relevant survival predicate and
state honestly what witness information it forgets.

## Tournament Analysis

Vertices are proof/computation lenses, not cauldrons. The pairwise observable
is which lens better preserves the cauldron survival predicate while retaining
useful witnesses. The route tournament is transitive:

`removal active-state DP > finite-sums subset obligations > weak Schur
first-boil > fixed-arity generalized Schur > classical Schur comparison >
greedy residue coloring > raw cauldron labels`.

The score histogram is one vertex in every outscore `0..6`, with no directed
3-cycles and exactly one Hamiltonian path.
