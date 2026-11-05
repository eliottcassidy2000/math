# Session 196: an exact two-apex low-cycle language

## Status and scope

This packet gives an independently derived, all-order finite-state classification
of the tournaments with exactly two Hamiltonian cycles **within the two-apex
class defined below**. It also gives exact Hamiltonian-path transfer matrices
for the class and a closed formula for its orientation-word inventory.
It does **not** classify all tournaments with two Hamiltonian cycles, and it does
**not** prove the full `{7,21}` Hamiltonian-path spectrum conjecture.

The motivating commit metadata about a general gluing classification was not
used as a theorem: no permitted precise source path for its proof was supplied.
The only repository source consulted was the explicit Session 195 guide at
`00-navigation/S195.md`, commit `ee07a5aac9cdfd8d3834fa834df4fb357fcaf769`.
All proofs and programs below are independent of the unverified metadata.
“New” here means derived in this session, not a claim of priority in the literature.

## 1. The four-letter model

Fix a transitive base `x_1 -> ... -> x_m` (every earlier vertex points to every
later vertex) and two distinguished vertices `a,b`, with `a -> b`.
A word `w=w_1...w_m` over `{0,1,2,3}` describes the remaining arcs:

- bit 0 of `w_i` is 1 exactly when `a -> x_i`;
- bit 1 of `w_i` is 1 exactly when `b -> x_i`.

Write `T(w)` for the resulting tournament, `H(w)` for its number of directed
Hamiltonian vertex words, and `C(w)` for its number of directed Hamiltonian
cycles counted modulo cyclic rotation. The named roots and the transitive base
order are part of this model. Counts of words are not unlabeled-tournament
counts and therefore are not new coefficients of A000568.

Every tournament that becomes transitive after deleting at most two vertices
has a presentation of this form, after choosing and ordering two roots, whenever
its order is at least two. Different presentations can describe isomorphic
unrooted tournaments.

## 2. Exact cycle transfer: nine coordinates

Start each Hamiltonian cycle at `a`. It has a unique form

`a, increasing(S), b, increasing(T), a`,

where `(S,T)` is an ordered partition of the base vertices. Indeed, every run
containing only base vertices must be increasing. Conversely, every such
partition satisfying the four boundary arcs gives a Hamiltonian cycle.

For either slot, retain one of three states:

- `E`: empty;
- `G`: nonempty and its current last vertex points to the following root;
- `B`: nonempty and its current last vertex does not point to the following root.

There are nine pairs of states. Begin with weight 1 at `(E,E)`, 0 elsewhere.
Read the base vertices in increasing order. For each vertex choose one slot.
If that slot is empty, require that its preceding root points to the vertex.
Then replace that slot's state by `G` or `B` according to whether the vertex
points to its following root; leave the other state unchanged. Add weights when
transitions coincide.

For `S`, the preceding/following roots are `(a,b)`. For `T`, they are `(b,a)`.
The accepting states are `(E,G)` and `(G,G)`: an empty `S` is permitted by
`a -> b`; an empty `T` would require the absent arc `b -> a`.
The sum of these two coordinates is exactly `C(w)`.

This is a bijective proof, not an extrapolation from finite enumeration.

## 3. Exact path transfer: two copies of 27 coordinates

A Hamiltonian path has one of the two root orders `ab` and `ba`. For order `ab`,
its unique slot form is

`increasing(S_0), a, increasing(S_1), b, increasing(S_2)`.

Use the same `E,G,B` states for each of the three slots. The preceding/following
boundaries are `(none,a)`, `(a,b)`, `(b,none)`. A missing preceding boundary
imposes no first-vertex condition. A missing following boundary makes every
nonempty last-vertex state good. All three slots may be empty.

For root order `ba`, use boundaries `(none,b)`, `(b,a)`, `(a,none)`. The middle
slot must now be nonempty because `b -> a` is absent; the endpoint slots may
still be empty. Sum the accepting coordinates for both orders to obtain `H(w)`.

Thus `H` is an explicit nonnegative integer weighted-automaton series of
dimension at most 54 on the same four-letter alphabet as the cycle filter.
This representation is uniform in base size. It counts actual tournaments, so
any word manufactured by an arithmetic construction is already a graph witness.

## 4. Saturation gives an exact finite cycle filter

Replace every coordinate by `min(3, coordinate)` after each cycle-transfer step.
This preserves the exact alternatives `C=0`, `C=1`, `C=2`, and `C>=3`, because
saturation commutes with sums of nonnegative integers. There are at most `4^9`
possible vectors. The complete reachable-state construction has **251 states**.
Its output-preserving minimal deterministic quotient has **32 states**.

`cycle_dfa.json` records the 32-state transition table and the output in
`{0,1,2,3}`; output 3 means `C>=3`. These are exact finite certificates:
`verify_language.py` reconstructs every reachable saturated vector, checks the
quotient's output and every transition, and checks that every pair of quotient
states has a distinguishing continuation.

Consequently, the words with exactly two Hamiltonian cycles form an explicit
regular language. The same construction works for any fixed upper cycle bound,
using saturation at one more than that bound.

## 5. A finite explicit regular expression for exactly two cycles

Juxtaposition denotes concatenation, `|` denotes disjoint alternatives in the
first-exit decomposition below, `epsilon` is the empty word, `*` permits zero or
more copies, and `+` requires one or more copies. Define

```
P = epsilon | (1|3)*1
Q = epsilon | (2|3)*2
R = epsilon | 3*1P | 3*2Q
```

Then `C(w)=2` if and only if `w` belongs to the following union:

```
1(0|1)* [ (20 | 310 | 32 3*1)P | 2 3+0 ]
|
2(0|2)* [ (1 3*2 | 30)Q | 33 3*0 ]
|
3 [ 0+ (2R | 3+2Q) | 10P | 20 | 2 3*1P ].
```

Spaces here are only for legibility, not alphabet symbols. This expression is
an exact all-order characterization of the stated subclass.

One can derive it by solving the right-language equations of the 32-state
quotient. For a second mechanically independent all-order check,
`verify_language.py` compiles the displayed expression using Thompson's
construction, then explores the complete product of its epsilon-NFA subsets
and the cycle DFA. The product has 66 reachable states, and every state agrees
on acceptance. This is a finite equivalence proof for all words, not a check
only through a chosen word length.

The expression has free binary prefix/suffix blocks and unary middle blocks.
That exposes the specific degrees of freedom available for a future inverse
arithmetic argument. Merely having the expression does not prove that its
weighted path-count image covers all admissible odd integers.

## 6. Exact word inventory

Let `A_m` count length-`m` words with `C(w)=2`. The expression is unambiguous:
the leading letter chooses its outer branch, the free binary prefix ends at a
letter outside its alphabet, unary runs have a unique first-exit position, and
`P`, `Q`, and `R` have unique indicated stopping boundaries.

Its ordinary generating function is

```
sum_{m>=0} A_m z^m
  = z^3 (7 - 14z + 8z^2 - 5z^3) / ((1-z)(1-2z)^2).
```

In particular,

```
A_0=A_1=A_2=0,  A_3=7,
A_m=(11m+6) 2^(m-5)-4  for every m>=4.
```

The first values from `m=3` are `7,21,57,140,328,748,1676,3708,8124`.
These happen to begin with 7 and 21, but they count decorated orientation words;
they have no direct bearing on the two forbidden Hamiltonian-path values.

For verification, `P(z)=Q(z)=(1-z)/(1-2z)` and `R(z)=1/(1-2z)`.
The three outer branches contribute respectively

```
z/(1-2z) [ (z^2+z^3+z^3/(1-z))P(z) + z^3/(1-z) ],
z/(1-2z) [ (z^2/(1-z)+z^2)P(z) + z^3/(1-z) ],
z [ z/(1-z)(zR(z)+z^2/(1-z)Q(z))
    + z^2P(z)+z^2+z^2/(1-z)P(z) ].
```

Their sum simplifies to the displayed generating function. For `m>=4`, extract
coefficients from

```
5z^3/4 + z^2/2 + 47z/16 + 61/16
 - 5/[32(1-2z)] + 11/[32(1-2z)^2] - 4/(1-z).
```

## 7. Verification and remaining scope

`apex_transfer.py` was checked against a separate adjacency/subset-DP counter
(`census.py`) on all 1,364 four-letter words of lengths 1 through 5, covering
all the modeled tournaments on 3 through 7 vertices. Both `H` and `C` agreed.
The empty word gives `H=1,C=0` directly.

`verify_language.py` additionally verifies the complete finite quotient,
all-order regular-expression equivalence, generating-function identity, and
the word inventory through length 40. Only the finite graph cross-check has a
size bound; the transfer bijection and finite-automaton equivalence have none.

This packet supplies a rigorous replacement for relying on the inaccessible
low-cycle gluing claim **only inside this explicitly stated subclass**. It
neither asserts that every two-cycle tournament has transitive two-vertex
deletion nor identifies its full scalar path-count image.

An appropriate next step is to combine the expression's free binary blocks
with the 54-coordinate path transfer and prove inverse-target rules, or prove
that some targets require leaving this two-apex subclass. Finite scalar coverage
alone would not answer that all-order question.

## 8. Exact bounded scalar image at arbitrary word length

There is an exact finite procedure to compute all values `H<=B` in this class
subject to `C` being 1 or 2, **without any word-length bound**. Saturate each of
the two 27-coordinate path vectors at `B+1`, and keep the 32-state cycle filter.
Nonnegative additivity again ensures that final path counts at most `B` are
unchanged. If an accepting-coordinate sum is greater than `B`, the true sum is
greater than `B`; if it is at most `B`, it is the exact true sum.

The following elementary monotonicity is needed to prune prefixes.

**Induced-extension lemma.** If `U` is an induced subtournament of `V`, then
`H(U)<=H(V)`.

**Proof.** It suffices to add one vertex `v`. For each Hamiltonian word
`x_1...x_n` of `U`, insert `v` as follows. If `v->x_1`, prepend it. Otherwise,
choose the first index `j>=2` with `v->x_j`, if one exists. All earlier
vertices, including `x_{j-1}`, point to `v`, so inserting `v` immediately before
`x_j` is valid. If no such `j` exists, append `v`, since `x_n->v`. This assigns
a Hamiltonian word of `V` to every Hamiltonian word of `U`; it is injective
because deleting `v` recovers the original word. Iterate for multiple added
vertices. QED.

Each prefix tournament `T(w_1...w_k)` is induced in any continuation, so once
its `H` exceeds `B`, no continuation can return to the target range. We may
also discard cycle-DFA states from which no continuation reaches output 1 or 2.
Explore all remaining reachable triples consisting of the cycle state and the
two saturated path vectors, merging equal triples. There are at most
`32*(B+2)^54` triples, so the procedure terminates. (Most are unreachable.)
Every kept transition and accepting scalar count is exact. Every admissible
word has all prefixes kept, by the induced-extension lemma and coaccessibility.
Thus the eventual closed reachable set proves both presence and absence for
all word lengths. The maximum discovery depth is only a diagnostic, not a
cutoff used in the algorithm.

The complete results produced by `bounded_image.py` are:

- `B=67`, saturation at 68: **3,388 states**, **3,999 retained transitions**,
  maximum discovery depth 15. The missing odd strong counts in `[3,67]` are
  exactly **7,21,67**.
- `B=501`, saturation at 502: **50,749 states**, **58,358 retained transitions**,
  maximum discovery depth 20. The missing odd strong counts in `[3,501]` are
  again exactly **7,21,67**.

The JSON output records one explicit orientation word for every positive result.
In particular, this proves an **all-order exclusion of 67 for the two-apex
class with one or two Hamiltonian cycles**. The same conclusion for all
low-cycle tournaments would require a further theorem not proved here.

This exclusion is genuinely about the low-cycle requirement: the two-apex word
`(1,2,3,0,1)` has **H=67 and C=4**, as checked by both the transfer and independent
adjacency/subset-DP. The global Hamiltonian-path value 1 is realized by a
transitive tournament; it is absent from the reported strong-count interval
only because `C=1 or 2` was required.
