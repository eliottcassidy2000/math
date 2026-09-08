# Odd-braid support runs and the H4 single-cycle obstruction

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The main theorem is an abstract permutation statement. Its application to
mixed-cusp Keller support is conditional on a separately certified actual
marked H4 supplier. Numerical braid words alone are not that supplier.

## 1. Exact statement and inheritance

Let `a,b,c,d` be permutations of a finite set, each a single nontrivial
cycle of the same length `m>=2`, with arbitrary additional fixed points.
Suppose

    aba=bab,   bcb=cbc,   cdcdc=dcdcd,
    ac=ca,     ad=da,     bd=db.                         (1)

**Theorem.** The four permutations are equal, and their moved
supports are identical. In particular,
if the four permutations generate a transitive action of degree D, then
`D=m` and every generator is fixed-point-free. There is no such transitive
action in which a generator has a fixed point.

This is uniform in the cycle length and the ambient degree. It does not
assume the generators are involutions, use a Coxeter-group order, or classify
finite quotients of the Artin group. Equality follows from the remaining
relations after support coincidence is proved. The relations in (1) are the
standard Artin H4 chain with labels `3,3,5`; the proof only uses the displayed
relations, so no literature identification is needed as a dependency.

The closest proved mechanism is the arbitrary-cycle ordinary half-support
bound and the commuting-single-cycle support dichotomy in
[the marked D4 closure](planar_jc48_sep08_three_cusp_closure.md) and
[the earlier single-cycle theorem](planar_jc48_sep08_cycle_support.md).
The canonical hostile is `(123),(345)`: it satisfies the five-term braid
relation with moved supports of size three intersecting in one point.
Thus the ordinary half-support bound fails at a five-cusp. The corrected
near miss is a third-support bound, obtained by retaining consecutive
positions within a cycle rather than only the cardinality of its support.
The least-used sidecar is the equality case: an ordinary half-overlap forces
strict alternation along each involved moved cycle.

The live concepts are actual positive meridians, full moved supports,
consecutive cycle runs, commuting support components, and retained sheets.
The connection from (1) to support sets preserves pointwise fixedness and
cycle images. It forgets cycle order only after the required alternation
has been proved. Mere intersection cardinalities lose the final obstruction.
No tournament is intrinsic to these relations.

## 2. A general odd-braid run lemma

Let `sigma,tau` be arbitrary permutations satisfying the alternating relation
of length `2r+1`, where `r>=1`. Put `S=supp(sigma)`, `T=supp(tau)`.
Every sequence

    x, tau(x), ..., tau^r(x),   x in T,                (2)

meets S. Indeed, if sigma fixes every point in (2), apply the two alternating
words to x. The word beginning and ending in sigma gives `tau^r(x)`;
the other gives `tau^(r+1)(x)`. Equality would imply `tau(x)=x`, contrary
to `x in T`. Repetitions in a short tau-cycle cause no exception.

Therefore every tau-cycle in T meets S, and between successive points of
`S intersect T` there are at most r points outside S. Counting the cycles
separately gives

    |T| <= (r+1)|S intersect T|.                       (3)

By symmetry the same bound holds with S and T exchanged. Odd braiding also
conjugates the generators, but that fact is not needed for (2)--(3).
For length three this is the familiar half-support bound. For length five
it gives `|S intersect T|>=ceil(|T|/3)`. The five-letter hostile attains
this latter equality and violates the former bound.

A useful equality sidecar is independent of the five-term relation. If the
ordinary relation holds, `|S|=|T|=m`, and their overlap has size m/2, then
tau alternates between `T intersect S` and `T minus S`. The run lemma permits
no two successive points outside S. Each outside point has a distinct next
point inside; equal cardinalities make this a bijection, hence no two
successive inside points occur either. This works cycle by cycle after
summing their equalities, and in particular for a single m-cycle.

## 3. The complete support reduction

Write `A,B,C,D` for the supports of `a,b,c,d`. If two single m-cycles
commute, their supports are equal or disjoint. To see this, each support
is invariant under the other permutation. A nonempty intersection contains
an entire nontrivial cycle of the other permutation, hence its entire
m-element support. Equality follows from the common size. Thus all three
commuting edges of (1) carry an equality/disjointness dichotomy.

First suppose `A=C`. The five-term relation forces `C intersect D` nonempty
by (3). Since a commutes with d, this gives `D=A`. The ordinary relation
forces `A intersect B` nonempty, and b commutes with d, giving `B=D`.
All supports are equal, as claimed.

It remains to rule out `A intersect C` empty. The two ordinary relations
give

    |A intersect B|>=ceil(m/2),
    |B intersect C|>=ceil(m/2).

The intersections are disjoint subsets of the m-element set B. Hence m
is even, each intersection has size m/2, and

    B subset A union C.                              (4)

The equality sidecar of Section 2 says that c strictly alternates between
`C intersect B` and `C minus B`.

Since d commutes with a, D is equal to or disjoint from A. Equality is
impossible: the five-term relation requires `C intersect D` nonempty,
whereas C is disjoint from A. Thus D is disjoint from A. Since b commutes
with d, D is equal to or disjoint from B. Equality is again impossible,
because B meets A whereas D does not. Consequently

    D intersect A = D intersect B = empty.            (5)

Every point of `C intersect B` is therefore outside D. If any point of
`C minus B` were also outside D, its predecessor and successor under c
would both lie in `C intersect B` by strict alternation. Three consecutive
points of the c-cycle would be outside D, contradicting (2) for r=2.
We have proved

    I := C intersect D = C minus B,   |I|=m/2,        (6)

and c interchanges I and `C minus D`.

## 4. Five-term braiding forbids the remaining equality case

Take any `x in C minus D`, and put `v=c(x) in I`. Since d fixes x,
applying `dcdcd=cdcdc` to x shows that

    y=d c d(v)   satisfies   c(y)=y.                 (7)

If `d(v)` lay in I, then `c d(v)` would lie in `C minus D` by (6), so d
would fix it. Thus `y=c d(v)` would lie in C, where c has no fixed point,
contradicting (7). It follows that `d(v) in D minus C`. Since c fixes
that point, (7) now says `y=d²(v) in D minus C` as well.

The map c takes all of `C minus D` onto I. Thus the preceding conclusions
hold for every v in I:

    d(I) subset D minus C,
    d²(I) subset D minus C.                          (8)

Moreover d(I) and d²(I) are disjoint. An equality `d(v)=d²(w)` with v,w in I
would imply `v=d(w)`, contradicting the first inclusion in (8). Both images
are disjoint from I. The three pairwise disjoint m/2-element sets
`I,d(I),d²(I)` lie inside the m-element support D. This is impossible.
The disjoint case has been ruled out, proving support coincidence.

There is a further exact conclusion. On the common m-element support, the
centralizer of the single m-cycle d consists of its powers: a commuting
permutation is determined by the image of one point and then by following
the cycle. Both a and b commute with d, so they commute with each other.
The ordinary relation then gives a=b. Since a commutes with c and b braids
with c, it follows that c=a. Finally a commutes with d and c braids with d
in length five, forcing d=a. Every generator fixes the complement, so this
is equality on the entire ambient set. This proves the theorem.

If the action generated by the four permutations is transitive, their
common nonempty support must be the whole ambient set, since every point
outside it is fixed by the entire group. Hence D=m in the statement of
Section 1. The support letter D in Sections 3--4 denotes the fourth
support; this final D denotes the action degree, as declared there.

## 5. Scope, failure boundaries, and the mixed-cusp connection

The result requires one nontrivial cycle in each generator, with the same
length. Multiple moved cycles destroy the equal-or-disjoint commuting
support dichotomy; no exclusion of all H4 actions is claimed. The result
also needs the terminal five-term relation. The adjacent transpositions

    a=(12), b=(23), c=(34), d=(45)

on five letters satisfy both first ordinary relations and all three
commutations, and generate a transitive action with fixed letters. They
satisfy the terminal three-term relation but fail the terminal five-term
relation. This is a hostile to discarding or shortening that edge.
The equal-generator m-cycle action on m letters is a positive fixed-free
control and shows why the retained fixed-sheet hypothesis is essential
for a geometric application.

For the prospective mixed-cusp geometry, the source is an actual common
positive-meridian generating tuple for the complement of a whole
irreducible nonproperness curve. Its target is (1), with the actual cusp
pairs transported to the `3,3,5` edges and the node pairs to the three
commuting edges. The preserved predicate is the common cycle type and
pointwise fixedness of genuinely retained sheets. Conjugate meridians
have equal cycle type, and an actual retained subset has size k>=1.
The map loses the specific subsets and Euler counts; neither is needed
for this particular obstruction. The required sidecar is a certified
actual path/cluster supplier with generation, positive access, and the
correct local pairs. A numerical word list is insufficient.

**Conditional corollary.** Once that actual supplier is proved, the
corresponding whole-support Keller monodromy cannot have a single
nontrivial cycle in its meridian cycle type, in any degree. This note
alone does not supply the actual geometry, exclude arbitrary multi-cycle
monodromy, close that curve class, or prove JC(2).

## 6. Exact controls and reproducibility

The standalone source
[planar_jc48_sep08_h4_single_cycle.py](../../04-computation/planar_jc48_sep08_h4_single_cycle.py)
imports no inherited mathematical implementation. Its gates use explicit
runtime checks, so optimized Python cannot disable them.

The declared universes are:

1. Every ordered pair in each S_d for 2<=d<=5, tested for literal odd
   braid lengths three, five and seven, including identity and multiple
   cycles. Every accepted pair is checked for the full consecutive-run
   lemma, the support bound, and odd conjugate support cardinalities.
2. All relative positions with positive overlap of two single m-cycles
   for 2<=m<=7. Disjoint supports are already ruled out by the odd-braid
   run lemma. The first
   cycle is standard; every overlap subset and every cycle ordering on
   the second support is retained, with union degree 2m-j. This includes
   every positive overlap j and all fresh labels. Exact braid counts and
   the ordinary equality alternation are recorded.
3. A full H4 tuple census for 2<=m<=5 and m<D<=min(3m,10). After fixing
   the first cycle by simultaneous conjugacy, all remaining m-cycles on
   the declared D letters are enumerated; the six literal relations are
   the only filters. Every accepted tuple is checked for identical
   supports, and its actual generated orbit is computed separately.
4. Every cycle ordering on the remaining alternating equality support
   for m=2,4,6,8. The terminal braid-five relation always fails.
5. The sharp five-letter hostile, the transitive five-letter chain with
   the wrong final label, and same-support fixed-free positive controls.

The analytic proof supplies the unbounded quantifiers. No finite cutoff
is extrapolated, and no theorem about actual curve monodromy is inferred
from an abstract tuple census. All15,111 always-active gates pass with normal and optimized outputs
byte-identical to the frozen output. The [independent full analytic/source audit](planar_jc48_sep08_h4_single_cycle_audit.md)
accepts the equality theorem, all finite universes and both complete replays. Reproduce from the worktree root with:

```bash
python3 04-computation/planar_jc48_sep08_h4_single_cycle.py
python3 -O 04-computation/planar_jc48_sep08_h4_single_cycle.py
```

Frozen SHA256 pins:

- `planar_jc48_sep08_h4_single_cycle.py`, 7003 bytes: `c42119cd358e291f3adf1d1d6b3c049fb44f65470ad30ae60ff54616b92d54b2`.
- `planar_jc48_sep08_h4_single_cycle.out`, 3650 bytes: `968390414754d58917891eb2f864aa8d62ebbf9e7ff2d6d8182c6f1dfaaa34ae`.
