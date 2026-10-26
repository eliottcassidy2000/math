# Rooted composition of actual black/blue complement lines

2026-10-26. **PROVED:** the joint two-port convolution and reflection-fold
rule for ordered sums, the run-forest reconstruction for arbitrary modular
substitution, and the marked Frobenius stabilizer statements in the companion
note. **FINITE-EXACT:** the explicit censuses and counterexamples below.
**OPEN:** a smaller universally closed response algebra for arbitrary
complement/re-marking dynamics, a complete prime atlas, and global Collatz
cycle classification. Local labels R1-R5 belong to this file. No external
novelty or independent human/agent audit is claimed.

The new obstruction is at the level of **joint roots**, not merely endpoint
degrees: two six-vertex presentations with equal separate entrance and exit
types yield eight-vertex complement targets with `H=289` and `H=253` after
the same singleton sandwich. For blue lines, even the correct root histogram
must be folded along reflection; independent multiplication overcounts.

## 1. Inheritance and conventions

The closest proved mechanism is
[marked composition M1-M5](../tournament_composition_classification_20261026/metagraph/RESULT.md):
ordered join, its complement-conjugate hinge product, and one-sided rooted
response histograms. Its canonical hostile is the order-four blue line with
half-edges producing `H=11/13` after singleton extension. The corrected near
miss is **THM-4804**,
[module shattering](../../../01-canon/theorems/THM-4804-complement-lines-shatter-every-module-and-the-carousel-fibre-as-frobenius-orbits-of-a-normal-basis.md):
old labelled module sets disappear, but an unlabelled tree shape can recur.
The least-used sidecar is the simultaneous orbit of all attachment ports.

Use the same carrier `X_n`: tournaments on positions `0,...,n-1` containing
the Hamilton path `0->1->...->n-1`. Let `tau` reverse all non-path arcs, and
let `sigma` reverse the positions and take full converse. They commute. For
`n>=3`, a literal line is `{x,tau x}`; it is blue iff `sigma x=x`, and black
otherwise. Ordinary class means isomorphism, without merging converse.
These conventions are those of **THM-643**,
[blue/black line parity](../../../01-canon/theorems/THM-643-blue-black-line-parity-structure.md).
Orders one and two are formal composition inputs; their `tau` does not form
two-element line pairs.

For a fixed labelled tournament `T`, ordered Hamilton paths modulo `Aut(T)`
are exactly its presentations in `X_n`; the action is free. Thus
`tc(T)=H(T)/|Aut(T)|`. A path count, a presentation count, a rooted response,
and a literal line count are different quantities.

Write `J(x,y)=x->y`. The complement of this join has interiors `tau x,tau y`,
all cross arcs pointing from the second block to the first, and one forward
hinge from the first block's exit to the second block's entrance.

## 2. R1: the joint two-port response and its convolution

A two-port state is an isomorphism class `h=[C;s,t]` of a tournament with
**ordered, jointly marked** entrance and exit vertices. The ports can coincide
for a singleton. Isomorphisms must preserve both roles simultaneously.
For any ordinary source class `A`, define

    F_A(h) = #{x in X_|A| : [x]=A, [tau x;0,|A|-1]=h},
    B_A(h) = #{such x with sigma x=x}.                              (1)

In particular `sum_h F_A(h)=tc(A)`. For `h=[C;s,t]`, `k=[D;u,v]`, define
`h star k` by putting every cross arc `D->C` except `t->u`, and retaining
the new ports `(s,v)`. This is well defined on joint rooted types.

**R1 (PROVED).** This rooted hinge operation is associative. If
`A=A_1->...->A_r`, then

    F_A(h) = sum_(h_1 star ... star h_r = h) product_i F_(A_i)(h_i).  (2)

Here the factors may themselves be nonstrong for the total-response formula.
Every Hamilton path of an ordered sum traverses each factor once, in order;
therefore its presentation fibre is the Cartesian product of the factor
fibres. Its complement is exactly the displayed rooted gluing. Grouping the
product by joint rooted target type proves (2). For associativity, either
parenthesization has unchanged interiors, backward arcs between every pair
of factors, exactly the forward hinges between consecutive factors, and the
same outer ports. This also proves the singleton boundary. ∎

For a source with at least two nonempty ordered factors and total order at
least three, its complement is strong: the displayed Hamilton path closes
using the last-to-first cross arc. The source is nonstrong. Thus each source
presentation counts exactly one literal line, with no division by two.
For an ordinary target class `Z`, sum (2) over states whose underlying class
is `Z`. The blue refinement below then gives black counts by subtraction.

## 3. R2: an exact blue-line rule from reflection folding

Use the **canonical ordered strong-component factorization**
`A=A_1->...->A_r`. Each factor is a strong component or a singleton.
Define

    rho([C;s,t])=[C^op;t,s].                                        (3)

This is the effect of `sigma` on the rooted complement target. Blue
presentations exist only if the component word is converse-palindromic:

    A_i is isomorphic to A_(r+1-i)^op for every i.                  (4)

**R2 (PROVED).** Subject to (4), form the blue target histogram as follows.
For each factor in the left half choose `h_i` with weight `F_(A_i)(h_i)`.
Its right mate is **forced** to be `rho(h_i)`, with no additional weight.
If `r` is odd, choose the central state with weight `B_(A_mid)(h_mid)`.
Glue this complete sequence by `star`, and sum the weights of equal output
states. The resulting histogram is exactly `B_A`.

*Proof.* Strong components are canonical and contiguous in every Hamilton
presentation. Reflection reverses their order and takes their converses.
Hence `sigma x=x` holds exactly when each right component presentation is
the literal reflection of the corresponding left presentation, and the
central presentation, if present, is reflection-fixed. Each chosen left
presentation has one such mate. Equation (3) follows from `sigma tau=tau sigma`.
Grouping these exact choices proves the rule. ∎

Canonical factors matter. A noncanonical cut can pass through the mirror
of another cut, so an arbitrary binary grouping is not a reflection test.
For the canonical factorization, (1)-(4) give both coloured line multiplicities
incident to every nonstrong source, given the atom response histograms.
This improves M5's full-presentation colour sidecar to a sufficient aggregate
interface for **ordered sums**. It does not assert that these histograms
alone close under every complement or legal re-marking of a strong target.

**FINITE-EXACT diagonal control.** The strong five-tournament with adjacency
row masks `(2,4,9,19,7)` has nine presentations. Its joint target-state fibre
sizes are `1,1,1,1,1,1,1,2`. Joining it with its converse gives exactly nine
blue presentations, one per reflected pair. Multiplying matching root-type
counts independently gives `7+2^2=11`, which is false. All 81 factor pairs
are checked independently in `audit.py`.

## 4. R3: separate entrance/exit types are insufficient

**FINITE-EXACT counterexample, independently checked by all permutations.**
On six vertices take free bits in lexicographic pair order `(i,j)`, `j>=i+2`,
least significant bit first. Presentations of masks `233` and `598` have rows

    x=(38,52,24,19,33,12),    y=(26,44,9,48,38,5).

They have the same strong source class. Their complement targets are

    C=(26,12,41,48,38,3),    D=(38,20,56,19,33,10).

Both lines are black; `C` and `D` are isomorphic. They are also isomorphic
with the entrance alone marked, and with the exit alone marked. Explicit
maps from `C` to `D` are respectively

    (0,5,3,1,2,4),    (3,1,2,4,0,5).

The first fixes vertex 0; the second fixes vertex 5. **No** isomorphism fixes
both 0 and 5. With `e` the singleton, direct counts give

    H(tau J(e,x,e))=289,    H(tau J(e,y,e))=253.                      (5)

Thus even the source class, target class, colour, and both individual rooted
types together do not determine iterated composition. The missing datum is
their **joint orbit**, not a finer endpoint degree. R1 retains it. This is
an explicit obstruction, not a claim of universal minimal state-space size.

## 5. R4: arbitrary prime substitution requires an ordered run interface

Let `T=Q[T_i]` use a specified modular partition, preferably the canonical
one when classifying. Given a Hamilton path `P`, cut it into maximal runs
inside children. Keep the quotient visitation word `omega`, whose adjacent
letters are distinct and follow arcs of `Q`. For child `i`, keep its ordered
path cover `F_i=(R_(i,1),...,R_(i,k_i))`, ordered by visits in `P`.

Define the **forest complement**

    C_i = T_i^op xor E(F_i).                                       (6)

It keeps the arcs used inside runs and reverses the other internal arcs.
Every run has named entrance and exit ports, with coincidences retained for
singleton runs. Mark all `2k_i` roles **jointly**, including their pairing
and their order of occurrence.

**R4 (PROVED).** The complete uncoloured partner `K_P(T)` is reconstructed by:

1. placing the rooted child targets `C_i` on disjoint vertex sets;
2. reversing every cross-child quotient arc;
3. restoring, for each adjacent pair of runs in `omega`, its exit-to-entrance
   path hinge.

Every pair of vertices is internal or crosses two children. In the first
case (6) is exactly `T^op xor E(P)`; in the second case the only exceptions
to the reversed quotient are the named hinges. This proves equality on
every arc. All three ingredients are therefore a sufficient one-step
composition interface, without assuming that a child remains a module
after the operation. ∎

To retain the **literal colour** and support legal re-marking, also retain
the ordered paths inside the runs, or equivalent data specifying their full
positional reflection. Then the data recover `P`, so blueness is the exact
predicate that its reversal is an anti-automorphism of `T`. Counts of path
covers alone do not supply these correlations. The full run-marked data
give a lossless representation; minimality among all encodings is not claimed.

Every global path corresponds uniquely to an admissible `omega` and these
ordered covers. One may therefore enumerate their responses and group by
target class and colour. Starting with a fixed labelled `T`, divide each
such path count by the actual `|Aut(T)|` to obtain presentation incidences:
the action on full ordered paths is free. If target and source are the same
ordinary class, divide its presentation incidences by two to count literal
self-lines (`n>=3`). This last division is unnecessary in the ordered-sum
case of R1. Apply any converse-node quotient only after fixing this convention.

Run arity is unbounded: in `C3[TT_k,TT_k,TT_k]`, visit the blocks in the word
`(0,1,2)^k`, taking successive vertices of each child. This is a Hamilton
path with `k` runs per child. Thus a universal one-entrance/one-exit block
description cannot serve all prime substitutions. The construction does not
rule out more efficient encodings of particular run patterns.

## 6. R5: how shattering and the symmetry threshold enter this interface

**Inherited PROVED result, THM-4804:** a complement line shares no proper
non-singleton module **vertex set** between its ends. At a path boundary
leaving a module `M`, the outside successor sees the departure vertex by a
retained path arc and every other vertex of `M` by a reversed arc. The next
path vertex cannot re-enter `M`, because uniform adjacency would contradict
the old path. The entering boundary is symmetric. Applying the involution
in reverse proves the two-sided statement.

Consequently, a fixed labelled substitution fibre cannot contain both ends
of a line. It does not follow from the same unlabelled tree shape alone.

**R5 (PROVED, unlabelled nonclosure of every carousel family).** For odd
`m>=3`, let `Q_m` have arcs `i->j` when `1<=(j-i mod m)<=(m-1)/2`, and
`A+=C3->1`, `A-=1->C3`, as in the
[carousel classification](../cyclictriple_composition_20261026/README.md).
Start with `T=Q_m[A+,...,A+]`. Traverse its blocks in order and each
block by its triangle followed by its sink. In any family member with `a`
plus blocks, degree multiplicities at `2m-2,2m-1,2m,2m+1` are
`a,3(m-a),3a,m-a`. For the chosen path the partner has multiplicities
`0,3m-1,2,m-1`, impossible for any `a`. To see this, complement sends each
interior degree `d` to `4m-1-d`; the first path vertex gains one further win
and the last loses one. The chosen endpoints are a triangle vertex and a
sink, producing the two degree-`2m` exceptions. This proves an actual line
leaves the family at every `m`, even after arbitrary relabelling. It does
not claim that every line leaves every unlabelled family member. At order
12 the incoming exhaustive census supplies that stronger finite statement.
The response construction in R4 works whether or not the target has a new
decomposition of the same shape; recompute its actual tree after gluing.

The incoming [Session 68 boundary proof, section 2](../session68_tournament_composition_classification_20261026/session68-root/METAGRAPH_BOUNDARY.md)
gives a stronger uniform-source conclusion, recovered from THM-785 and
THM-646: `c3(K_P(T))-c3(T)=d_T(first)-d_T(last)-1`. All degrees of a uniform
four-vertex decoration of any regular quotient have the same parity.
Hence **every** line from either uniform decoration changes the cyclic
triple count and exits that whole family. R5 is an alternative explicit
degree-profile witness, not a claim of priority for family nonclosure.

The [finite-field and sheet companion](ARITHMETIC.md) proves the relevant
marking law: a normal-basis port kills the residual cyclic stabilizer.
It connects this to the sharp `4p` threshold while distinguishing a mark
from an observer-preserving change of unmarked child type. The decorated
tree computes `|Aut(T)|`; the run response determines where its
`H(T)/|Aut(T)|` presentations actually go. Neither replaces the other.

## 7. Verification and research boundary

Run `python3 05-knowledge/results/rooted_lines_frobenius_sheets_20261026/verify.py`.
It verifies the manifest and replays three scripts normally and under `-O`
against frozen receipts. All use standard-library exact integer arithmetic.

`check.py` exhausts all path presentations for `n=1..6` (1,1,2,8,64,1024),
groups them by canonical strong-component word, and compares literal
coloured targets with both convolution formulas on all 32 nontrivial source
words. Eight such words have blue presentations. It reconstructs every
Hamilton-path partner in `C3[TT2,TT2,TT2]`, `C3[A+,1,1]`, and `C3[C3,C3,1]`
by the forest rule; each has paths re-entering a block. Its isomorphism
canonicalizer enumerates degree-compatible permutations, respecting all roots.
It checks the degree obstruction of R5 at odd `m=3,5,7,9,11,13`.

`audit.py` imports no producer code. It checks the six-vertex obstruction
using **all** 6! permutations, recounts both eight-vertex path totals using
**all** 8! orders, and counts the 81 reflection-pair candidates literally.
`arithmetic.py` has separate finite-field, modular CRT and bounded cycle-word
universes detailed in ARITHMETIC.md. Finite controls do not prove the
all-order statements; those follow from the displayed arguments.

The [research board](RESEARCH.md) records maps, losses, inherited results,
and stopping conditions. The completed result is a sufficient exact
interface with explicit obstructions to coarser choices. An effective
classification of all strong atoms and a universally small response algebra
closed under all black/blue dynamics remain **OPEN**.
