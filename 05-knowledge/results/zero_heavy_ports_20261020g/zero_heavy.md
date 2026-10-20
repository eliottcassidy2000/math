# Zero-heavy cores through exact feedback anchors

**PROVED:** the exact unlabelled count for tournaments made transitive by
one vertex deletion, its pair-free subfamily, the port-word clone criterion,
and a sharp family requiring linearly many separating anchors. The exact
finite-state canonical re-anchoring construction in section 5 also proves
rational class-count generating functions at every fixed feedback-vertex
bound, including the pair-free subclass.
**FINITE-EXACT:** independent class, word, clone, and re-anchoring controls.
No literature-priority claim is made.

## Inheritance and the changed parameter

The closest mechanism is the
[fixed feedback-deficiency kernel theorem](../kernels_games_detours_20261020f/general_fixed_deficiency.md).
Its hostile is an arbitrary pair-free core with all weights one: at large
feedback budget it has zero heavy coordinates, outside that theorem's
stable range. The corrected near miss is that an owner needs the full
ordered vertex data and cannot compare variable-length edge lists without
a guard. Here the parameter is **feedback vertex number**, the number of
vertices deleted to leave a transitive tournament. It is different from
feedback arc budget and from heavy-coordinate deficiency.

The niche incoming object is the
[full path-cover response](../strong_blue_games_20261020f/README.md), which
retains a profile rather than just one Hamilton count. The new incoming
[feedback-word packet](../feedback_word_hamilton_spectrum_20261020/proofs.md),
read at commit 422c69c14b before integration, already supplies the port-word
encoding and an exact finite matrix grammar for Hamilton-path counts. Its
[bounded-run barrier](../feedback_word_hamilton_spectrum_20261020/run_density.md)
concerns the numerical image of H, not the number of tournament classes.
We inherit the encoding and add the isomorphism/anchor quotient. We retain
arbitrary word lengths and arbitrary run counts.

The live concept board is: ordered transitive remainder; chosen feedback
anchors; complete adjacency letters; alternate-anchor owner; homogeneous
pairs; and the distinction between a word language and its H-value image.
Deleting anchor labels loses which presentation represents a class. The
restoration data are the anchor tournament, the full word, and the exact
re-anchoring map, rather than a scalar path count.

## 1. A closed scalar slice of the zero-heavy boundary

Let A_1(n) count isomorphism classes of tournaments on n vertices for which
some single vertex deletion leaves a transitive tournament. Include every
tournament on n<=1 by convention. Then

```text
A_1(0)=A_1(1)=1,
A_1(n)=2^(n-1) - (n-1) - binom(n-2,2),       n>=2.        (1)
```

The first nonempty values are

```text
n:       1  2  3  4  5   6   7    8    9    10
A_1(n):  1  1  2  4  9  21  48  106  227   475.
```

The ordinary generating function is therefore

```text
sum_(n>=0) A_1(n) z^n
 = 1 + z/(1-2z) - z^2/(1-z)^2 - z^4/(1-z)^3.             (2)
```

This is a genuine class-count recurrence without a growing core bank, on
the stated one-feedback-vertex slice. Explicitly, for n>=6,

```text
A_1(n)=5A_1(n-1)-9A_1(n-2)+7A_1(n-3)-2A_1(n-4).
```

It does not count all tournaments.

### Rooted words and the exact quotient correction

Distinguish a vertex x and order the remaining transitive vertices
`v_1 -> ... -> v_(n-1)`. Write b_i=1 when x->v_i. The word is a complete
rooted isomorphism invariant, because the transitive remainder has a
unique order. Hence there are exactly `2^(n-1)` rooted classes.

For an unrooted T, call a vertex eligible if deleting it leaves T
transitive. The number of rooted classes lying above T is the number of
automorphism orbits of eligible vertices, not their number. The following
classification computes that correction exactly.

**Lemma.** If a nontransitive tournament has at least two eligible vertices,
then it is exactly

```text
TT_a -> C3[1,1,m] -> TT_b,       a,b>=0, m>=1.              (3)
```

The brackets mean transitive inflation of a directed triangle. If m=1,
there are three eligible vertices and one eligible-vertex orbit. If m>=2,
there are exactly two eligible vertices and two eligible-vertex orbits.

**Proof.** Let u,v be two eligible vertices and suppose u->v. Their common
remainder X=T-{u,v} is transitive. In the unique order of X, both u and v
have a cut at which they can be inserted individually, since T-v and T-u
are transitive. If u's cut lies no later than v's, inserting both, with
u before v when the cuts agree, makes T transitive. Thus u's cut lies
strictly later. The nonempty interval W between the two cuts is transitive
and satisfies `v -> W -> u -> v`. Everything before both cuts dominates
this cyclic block; everything after both cuts is dominated by it. This
is (3).

When |W|=1 the only cycle is the displayed triangle, and its three rotations
are automorphisms. All three triangle vertices are eligible and equivalent.
When |W|>=2, deleting any one W vertex leaves a cyclic triangle, whereas
deleting u or v removes every cycle. No automorphism exchanges u and v:
the unique nontrivial canonical chain W within the cyclic strong component
is fixed, its two distinct incident directions distinguish u and v, and
the transitive outer chains are ordered. The outer chains themselves may
also be nontrivial canonical chains.
The two eligible roots are therefore inequivalent. QED.

There are consequently three cases contributing an excess over one rooted
class:

* `TT_n`: n inequivalent eligible roots, contributing excess n-1;
* (3) with m>=2: two root classes, contributing excess one;
* every other nontransitive class: one eligible-root orbit, no excess.

The distinct triples `a,b>=0,m>=2,a+b+m+2=n` give exactly
`binom(n-2,2)` classes. The strong cyclic block and its two exterior chain
lengths recover the triple intrinsically. Subtracting these excesses from
`2^(n-1)` proves (1). In particular dividing by n, or subtracting the number
of eligible vertices instead of their orbit count, would already fail on
the directed triangle.

### An explicit least-word language

With 0<1, the lexicographically least binary-word representatives are
exactly all binary words except the two disjoint regular languages

```text
0* 1+                 and                 0* 1{2,} 0 1*.
```

The first removes every transitive presentation except the all-zero word.
For (3) with m>=2, the two inequivalent root words are
`0^a 1 0^m 1^b` and `0^a 1^m 0 1^b`; the first is smaller, so the second
is removed. The m=1 presentations coincide and must not be removed.
Every other nontransitive class has one eligible-root orbit and thus one
word. This proves the language claim, rather than inferring regularity
from the rational count. It provides a concrete instance of section 5's
general canonical-owner construction.

## 2. Pair-free cores in the one-anchor slice

Let C_1(n) count the pair-free tournaments in the same slice. Then

```text
C_1(1)=1, C_1(2)=0, C_1(3)=1,
C_1(n)=2 for every n>=4,
sum_(n>=1) C_1(n)z^n = z+z^3+2z^4/(1-z).                 (4)
```

For base vertices, consecutive equal bits give a homogeneous pair.
Nonconsecutive base vertices are distinguished by an intervening base
vertex. Thus a pair-free one-anchor word must alternate. The anchor x
and base vertex v_i form a homogeneous pair exactly when

```text
b_1=...=b_(i-1)=0 and b_(i+1)=...=b_m=1,                 (5)
```

with b_i unrestricted. For m=0 the singleton survives; for m=1 neither
word survives; for m=2 only 10 survives, giving C3. At every m>=3 both
alternating words survive (5).

These two words give nonisomorphic tournaments. For even m>=4, one begins
with 1 and ends with 0 and is strong, while the other has a source and a
sink. For odd m>=3 one has a source but no sink, the other a sink but no
source. This proves (4), with no implicit multiplicity division.

The words retain exact clone information. Their long transitive stretches
must not be replaced by triangle modules: C3 is an internal cyclic atom,
whereas repeated letters give transitive homogeneous chains.

## 3. Exact ports and clone recognition for a anchors

Fix a ordered anchors X=(x_1,...,x_a) with induced tournament F. Their
complement is a uniquely ordered transitive base. Encode each base vertex
v by a letter b(v) in `{0,1}^a`, with b_i(v)=1 iff x_i->v. The pair
`(F,w)` reconstructs the entire tournament. This is the same exact encoding
used by the incoming Hamilton-count automaton, but here the target is an
unrooted isomorphism class.

Every homogeneous pair belongs to exactly one of three cases:

1. Two base vertices: they must be consecutive, and their letters agree.
2. Two anchors x_i,x_j: they agree against every other anchor, and
   `b_i=b_j` in every word letter.
3. Anchor x_i and base vertex at position p: all preceding letters have
   bit i equal to 0, all subsequent letters have bit i equal to 1, and
   the letter at p satisfies
   `b_j=1-F_ij` for every other anchor j. Bit i at p is unrestricted.

These are necessary and sufficient conditions obtained by comparing their
rows against every third vertex. The checker evaluates them separately
from the inherited canonical-chain routine. They are also finite-language
tests: equality of adjacent letters, a global bit equality, or existence
of one marked position satisfying fixed prefix/letter/suffix conditions.

No deletion-induced chain is silently contracted. All base vertices and
all their positions remain in w, including repeated letters. Passing to
pair-free cores means applying the complete three-case filter to the full
tournament, rather than assuming a base run is the final canonical chain.

## 4. Why a uniformly small separating anchor cannot suffice

Let T_m be the ordered join of m directed triangles. It is pair free.
Call S separating if the vectors `(1[s->v])_(s in S)` distinguish all
vertices outside S. Its minimum possible size is exactly

```text
sep(T_m)=m+floor(m/2)=floor(3m/2).                         (6)
```

Each triangle needs at least one selected vertex, since all its vertices
have identical comparisons to every vertex outside that triangle. If two
consecutive triangles each have only one selected vertex, the unselected
vertex losing to the first selected vertex and the unselected vertex
beating the second have the same complete signature. Thus the triangles
with at least two selected vertices form a vertex cover of the m-vertex
path. There must be at least floor(m/2) of them.

Conversely, select one vertex in every triangle and a second vertex in
every even-indexed triangle. Within each triangle the remaining signatures
are distinct. Two nonconsecutive triangles are distinguished by an anchor
between them. For consecutive triangles, one is doubly anchored, and its
remaining vertex has mixed comparisons to those two anchors, unlike every
vertex outside that triangle. This construction attains (6).

The same T_m has feedback vertex number m, since each of its disjoint
triangles must be hit, and one deletion per triangle suffices. Therefore
both the anchor size and the feedback-vertex parameter must grow on an
explicit zero-heavy family. Fixed-a results cannot be read as a uniform
solution for unrestricted tournaments.

## 5. All-fixed-a canonical regular language

**Theorem.** For each fixed a>=0, the numbers of ordinary
tournament classes with feedback vertex number at most a, and the numbers
of pair-free such classes, both have rational ordinary generating functions.
This is stronger than a finite Hamilton-count automaton: it must remove
different anchor choices exactly. It does not assert a bounded automaton
for unrestricted a.

For n>=a, any tournament in this class admits an anchor set of size exactly
a, padding a smaller feedback set if needed. Order its anchors and encode
it by `(F,w)`. Compare records lexicographically, first a fixed bit encoding
of F, then w. Both records for the same tournament have the same base-word
length n-a. Select the least numerical record over every ordered feedback
set of size a. This set of records is an isomorphism invariant, and any
record reconstructs the complete tournament. Thus there is exactly one
least record per unrooted class, regardless of how many anchor choices or
automorphisms yield that same record.

To prove the theorem, recognize the existence of a smaller record by a
finite nondeterministic automaton on w. The following specifies its finite
metadata and transition rules; no unbounded position counter is used.
There are finitely many anchor headers F. Treat each as a separate initial
branch and take their finite union; header comparisons below always use
the same fixed ordering on these finite bit strings.

### 5.1 Metadata for an alternative anchor set

Suppose the alternative set keeps a-s old anchors and selects s old base
vertices, where 0<=s<=a. Its new transitive base retains the old base in its
original order, except for those s deletions, and inserts the s displaced
old anchors. Guess at the start:

* the retained old anchors and their new anchor labels;
* the order in which the s selected base-anchor labels will occur;
* their s port letters, each from the finite alphabet `{0,1}^a`;
* the labels of the s displaced old anchors and their relative insertion
  order in the new transitive base.

This is finite metadata depending only on a. The guessed letters will be
verified when their selected positions are read. From this metadata one
already computes the new anchor tournament F': old-old comparisons come
from F; selected-base comparisons from their order; and old-base comparisons
from the guessed letters. There is no future adjacency information left
unrecorded.

### 5.2 Reading, deleting, and inserting

While reading w, nondeterministically mark exactly s positions as the new
base-origin anchors in the guessed order, verifying their guessed letters.
Those positions do not appear in the output base word. At word boundaries,
insert the displaced old anchors in their guessed order.
These insertions are epsilon transitions, of which at most a can occur
in an accepting run. Epsilon transitions do not change regularity and can
be removed by ordinary finite-state epsilon closure.

To certify that the new base is transitive, every ordinary input position
not selected as an anchor must have bit i equal to 0 if displaced anchor
x_i has not yet been inserted, and bit i equal to 1 if it has been inserted.
The guessed order of the displaced anchors must agree with their arcs in
F. These conditions check every pair in the new base: old base pairs were
already forward, the tests check old-base/inserted-anchor pairs, and the
constant F test checks inserted-anchor pairs. Selected input positions are
new anchors and are correctly exempt from these base tests.

The output letter of an unselected old base vertex is computed as follows.
For each retained old anchor, use its bit in the current letter. For each
selected base-origin anchor, use 1 if that selected position has already
been passed and 0 otherwise. The output letter of a displaced old anchor
uses F against retained anchors, and the complement of the corresponding
guessed input bit against every selected base-origin anchor. Thus the
entire alternate word is reconstructed without forgetting any port.

### 5.3 Bounded comparison delay and exact quotient

The alternative base word differs in its underlying vertex sequence by s
deletions and s insertions. At any stage, the difference between the number
of input letters consumed and output letters emitted is at most a in
absolute value. A FIFO buffer of at most a letters therefore suffices to
compare the two word streams lexicographically, retaining a three-valued
flag for equal/smaller/larger. The finite anchor-header comparison F' versus
F takes precedence. Require every guessed selection/insertion to occur and
the final stream lengths to agree.

Conversely, every alternate ordered feedback set produces exactly such
metadata, selected positions, and insertion events. Its remaining base is
transitive, and restriction to the surviving old base has the unique old
order. Hence this construction represents all competing records, not just
a convenient family of anchor moves.

The smaller-record language is therefore regular. Determinize and
complement it to obtain the owner
language; intersect with the finite clone filter in section 3 when counting
pair-free cores. Count accepted **words**, not nondeterministic accepting
paths: different metadata guesses can produce the same competitor and
must not introduce multiplicity. The deterministic transfer matrix has
integer entries counting letters, so its generating function is rational.
The finitely many n<a cases contribute a polynomial.

The construction's finite memory is of order a^2 bits before determinization
(anchor data, guessed letters, flags, and bounded queue), with constants
independent of word length. Determinization can make the state set very
large. This is an effective fixed-parameter route, not an efficient global
recurrence or a claim that all zero-heavy interfaces are bounded.

For clarity, if M is the adjacency matrix of the resulting deterministic
automaton, its entries count input letters producing transitions, and
the order-n count for a fixed header is `u M^(n-a) v`. Its generating
function is `z^a u(I-zM)^(-1)v`. Sum over the finite header branches and
add the finite n<a prefix. This supplies a constant-coefficient recurrence
for each fixed a. The effective construction can have a very large state
space; a full general deterministic automaton has not been materialized
in this packet.

## 6. Exact tests and stopping boundary

The checker generates all tournament classes through order seven by exact
one-vertex extension and canonicalization, independently checking the
known universe sizes `1,1,2,4,12,56,456`. It then directly tests deletion
to a transitive remainder. The one-feedback counts agree with (1).

Separately, all 1,023 binary words of lengths zero through nine are
constructed, quotient-counted, and tested for homogeneous pairs. The
ordinary counts agree with (1), and the core counts with (4). Every
one-anchor word is also checked against the explicit regular owner
language in section 1 and direct minimization over all eligible roots.
The order-four scalar recurrence passes on the available values. Every
eligible-root multiplicity in the full class census agrees with the
classification in section 1, including the essential three-vertex/one-orbit
triangle case.

The two-anchor least-word owner count, obtained by enumerating all alternate
valid anchor sets, agrees with the independent full-class census:

```text
n:                   2  3  4   5   6    7
feedback number<=2:  1  2  4  12  55  335.
```

For 3,431 actual re-anchoring operations the checker reconstructs the new
record using only old header/word and finite selected-position metadata;
it agrees with direct reconstruction from the tournament. The observed
maximum stream delay is two, within the claimed a bound. The three-case
clone predicate agrees with canonical-chain detection on every tested
representation. Exhaustive anchor subsets on joins of one through four
triangles give separating minima `1,3,4,6`, as in (6).

Reproduce from the repository root:

```text
python 05-knowledge/results/zero_heavy_ports_20261020g/zero_heavy_check.py
python -O 05-knowledge/results/zero_heavy_ports_20261020g/zero_heavy_check.py
```

The [receipt](zero_heavy_check.out) distinguishes the proved one-anchor
formula from bounded tests of the general owner. No polynomial
or recurrence is fitted to the two-anchor counts. All checks use explicit
exceptions rather than removable assertions.
Normal and `-O` runs passed in 1.58 and 1.57 seconds, respectively, with
byte-identical 899-byte output.

The parent and an independent agent audited the general finite-state proof,
including future selected-anchor incidences, displaced-anchor comparisons,
bounded stream delay, anchor relabeling, and nondeterministic duplicate
witnesses. Both accepted the construction. The only correction was to
scope W's canonical-chain uniqueness to its cyclic strong component in
section 1; no quotient or state condition was removed.

The achievement is a complete quotient correction on a nontrivial
unbounded zero-heavy slice, and a finite-state canonicalization theorem
for every fixed feedback-vertex bound. Unrestricted A000568 still needs
control as the anchor parameter grows. Neither the incoming H automaton
nor a rational fixed-a class count alone settles that remaining step.
