# A typed reduction map for the open-problem portfolio

2026-10-26. **PROVED / standard elementary computability:** the reductions
and quantifier statements below. They are a classification framework, not
new equivalences between the named open conjectures. **OPEN:** any such
specific equivalence absent an explicit map and proof. No novelty is claimed
for the standard hierarchy or completeness results.

## First choose the object being classified

A fixed conjecture is one sentence. NP-completeness classifies a *language
of arbitrarily many finite instances*, equipped with an input encoding and
a polynomial cost model. A single still-unknown truth value does not acquire
NP-hardness from being written as a program. Extensionally the language
`{0 : sentence phi is true}` is either empty or `{0}`; both are decidable.
Our lack of a proof deciding which is an epistemic issue, not a language
lower bound. Formal implication between two fixed conjectures is instead
proved relative to a specified base theory.

Use a reduction record with these fields:

    source family and input encoding; target family and input encoding;
    effective map; exact yes/no direction; designated initial state;
    operations preserved; clock/size overhead; lost data and sidecar;
    proof system; finite hostile and source status.

Graph isomorphism, simulation, an arithmetic reparameterization, and
polynomial-time many-one reduction are four different entries in that
record. A common visual motif supplies none of their missing fields.

## Quantifier audit

Here R denotes a decidable predicate, not a heuristic numerical test.
The displayed levels are presentations/upper bounds for fixed sentences;
we do not assert they are minimal or complete instances.

| Statement / family | Quantifiers and class | What can be finitely witnessed |
|---|---|---|
| One designated machine halts | exists t: the t-step simulation halts; Sigma^0_1 | A literal finite run |
| One designated machine never halts | forall t: the t-step simulation has not halted; Pi^0_1 | No general finite witness for truth |
| Goldbach | forall even N>=4: bounded search finds p,q<=N; Pi^0_1 | A counterexample N, verified by exhaustive finite primality search |
| Integer-speed LRC at a fixed dimension | forall finite integer speed tuples: the exact finite piecewise-linear maximum meets the target; Pi^0_1 | A bad tuple and its rational event/max calculation |
| All-order finite black-core connectivity, with its exact metagraph definition fixed | forall n and admissible finite source: finite exhaustive component test succeeds; Pi^0_1 | A finite disconnected instance and exhaustive certificate |
| A fixed proposed total A000568 formula is correct at every n | forall n: its value equals finite tournament enumeration; Pi^0_1 | A mismatching n; totality of the proposed formula is required |
| Classical Collatz | forall n>=1 exists t: C^t(n)=1; Pi^0_2 presentation | Each successful orbit separately; a nontrivial finite cycle certifies only one kind of failure |
| Twin-prime infinitude | forall B exists p>B: p and p+2 are prime; Pi^0_2 presentation | Arbitrarily distant individual pairs, not the entire assertion |
| Totality of an arbitrary program | forall input x exists t: halt(x,t); Pi^0_2-complete family | A finite halting run for each fixed input, with no uniform finite coverage |
| P!=NP | forall polynomial-clocked deterministic SAT candidates exists finite input on which the candidate errs; Pi^0_2 presentation | An error against one candidate, not all candidates |

For LRC this row deliberately concerns the repository's integer-speed
formulation. It does not silently assert an additional real-to-integer
reduction. Finite exact testing comes from the finitely many breakpoints
and pairwise intersections of the distance-to-integer linear pieces on a
period. For the metagraph row every legality/port rule must itself be
decidable on finite data; an unspecified “black core” is not yet a language.
For P versus NP, enumerate all machines with explicit polynomial clocks,
including constants. Exact SAT truth is decidable by a finite exhaustive
search, so candidate correctness on one input is decidable even if slow.

RH has a Pi^0_1 formulation as well. [Lagarias's criterion](https://websites.umich.edu/~lagarias/doc/elementaryrh.pdf)
gives the equivalent universal inequalities
`sigma(n) <= H_n + exp(H_n) log(H_n)`. A strict violation is semidecidable:
rational interval bounds for the computable real right side eventually
certify it. Dovetail over n and accuracy. Thus failure has a finite encoded
rational certificate. One need not pretend that comparing arbitrary
computable reals for equality is decidable. This is a *derived quantifier
consequence* of the cited equivalence, not an RH proof or an equivalence
with Goldbach. The cited paper's strengthened equality clause can be
omitted for this standard universal-inequality criterion.

## Uniform completeness, and the precise obstruction

For arbitrary effective decidable R(e,n), compile the program that checks
n=0,1,2,... and halts at the first failure. Then

    forall n R(e,n) iff the compiled blank-tape machine never halts.

The compiler uses an effective class of total predicates (for example
primitive-recursive predicates), so it never hides an undecidable predicate
test. This proves the usual Pi^0_1 universality of nonhalting. It does not
show that a fixed Goldbach machine simulates every compiled machine.
Goldbach25 is *one* instance in this common target family. RH's machine is
another; two reductions into the same family do not compose into a
reduction between those two instances.

For a family `forall x exists t R(e,x,t)`, compile a machine on input x
that searches successive t. Its totality is equivalent to the source.
This proves Pi^0_2 hardness of program totality; membership follows from
the bounded native simulation predicate. Classical Collatz already has
the form of one totality assertion, with its specific map fixed.

**Uniform obstruction.** There is no computable compiler mapping every
program-totality instance e to one blank-tape machine M_e such that
`Total(e) iff M_e never halts`. Otherwise totality would be Pi^0_1 by
substitution, contradicting strictness of the arithmetical hierarchy.
Replacing “never halts” by “halts” fails for the corresponding Sigma^0_1
reason. This is an obstruction to a *uniform unrestricted compiler*, not
a theorem that a particular Collatz sentence lacks a special Pi^0_1
reformulation. Such a reformulation would need a new mathematical argument.
An outer loop cannot detect failure of an orbit merely by waiting longer.

The same distinction explains why a Busy Beaver oracle settles the
Goldbach machine's halting but a single numerical Busy Beaver value does
not automatically discharge all inputs of a totality problem. The
[Busy Beaver Frontier](https://www.scottaaronson.com/papers/bb.pdf)
states the exact native step convention and the halting-oracle equivalence.
No claim about a practical computable bound follows.

NP-completeness uses a different restriction: the map must be polynomial
in input length and instances have finite polynomially checkable witnesses.
A standard route uses a bounded nondeterministic computation tableau and
SAT, as in [Cook's theorem-proving paper](https://www.cs.cmu.edu/~15455/resources/Cook1971-complx-thm-proof.pdf).
An unbounded halting compilation is not such a reduction. Nor does expanding
a binary-encoded enormous time bound into that many tape rows give a
polynomial construction. Our machine/tournament carriers make no NP claim.

## What the machine clocks genuinely transfer

The inherited Session25 interface already proves the useful theorem. With
an absorbing raw halt, injective section encoding e, closed macro map G,
positive finite native durations tau, a verified initialization, and halt
reflection, require

    F^(tau(s))(e(s)) = e(G(s)).

Iteration preserves the orbit at the accumulated native clock. Positivity
makes an infinite macro orbit have unbounded native time. Absorption and
halt reflection prevent a hidden halt between two live checkpoints.
To quotient s further by q, require both qG and terminal status constant
on q-fibres; exact per-edge clocks additionally require tau constant there.
This is the original packet's criterion, cited rather than rediscovered.

Our bell calculations give complete, explicit examples with different
native clocks. Goldbach's native proof instead retains the tested even
number, partition search and division invariants. The “small control,
large arithmetic” theme is real, but the retained arithmetic state is what
makes its proof work. Removing the unbounded tape is the very operation
that must be justified, not a harmless graph simplification.

The ordered configuration tournament in `README.md` has exact one-step
conjugacy to the native dynamics, with all non-orientation information
explicitly colored. Root's program-table tournament is the complementary
finite program carrier. Neither supports a bare-tournament halting test.
In particular directed Hamiltonian-path *existence* cannot encode an
arbitrary yes/no problem on unrestricted tournaments: it is always true.
Any useful reduction must change the predicate or retain constraints.

## A concrete research interface, not a claimed solution

The useful common work unit is an **operation-preserving checkpoint proof**:
retain exact configuration parameters, prove all applicable native rules,
prove closure from the designated source, and attach a clock and target
predicate. Bouncer growth, the corrected bell, Goldbach trial division,
Collatz precision resets, and port-preserving metagraph repairs can share
this proof interface while retaining different unresolved obligations.

For an actual connection to settle two open problems, produce a map in at
least one direction that transports those remaining obligations. A match
of state counts, triangular numbers, parity frequencies, or dynamical
silhouettes is a candidate generator; it is not such a map. The earlier
Collatz extended-graph packet explicitly separates existence of a chosen
path from all-path/deterministic convergence
(`../../collatz_mod6_20260917_extended_collatz_scc.md`, inheritance and
section 2); that quantifier warning survives every encoding here.
