# No finite tournament has exactly 7 or 21 Hamiltonian paths

**PROVED, previously established. Independent finite recertification here.**
The result is inherited from [THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
which records the 2016 provenance. This write-up does not claim a new exclusion
or prove that every other odd count occurs. All graphs are Boolean tournaments;
H counts directed Hamiltonian vertex words, without division by automorphisms.

## 1. Insertion is an injection

Every directed Hamiltonian path of T-v extends to one of T by inserting v.
For a path x1,...,xm, prepend v if v->x1, or append v if xm->v.
Otherwise x1->v and v->xm, so some consecutive pair satisfies
xi->v->x(i+1). Choose the first feasible position deterministically.
Deleting v from the output recovers the input, proving

    H(T) >= H(T-v).

Iterating gives monotonicity for every induced subtournament. It also proves
existence of a Hamiltonian path, starting with the one-vertex tournament.
Consequently every isomorphism class has a labelling with spine
0->1->...->n-1. Enumerating its remaining (n-1)(n-2)/2 arcs is a complete
cover, although it contains duplicate isomorphism classes.

## 2. Elementary strong-deletion lemma

Every strongly connected tournament T of order n>=4 contains an induced
strong tournament on n-1 vertices. This is a special case of
[Moon's subtournament theorem](https://doi.org/10.4153/CMB-1966-038-7);
the following proof supplies the needed case directly.

Choose a proper strong vertex set S of maximum size, with |S|>=3. Such a
set exists: a strong tournament contains a directed triangle (a shortest
directed cycle has length three, since a chord shortens any longer one).
If |S|=n-1, we are done. Assume at least two vertices remain outside S.

An outside vertex with both an in-neighbour and an out-neighbour in S
could be added to S while retaining strongness, contradicting maximality.
Thus every outside vertex either dominates all of S (class A), or loses
to all of S (class B). Strongness of T forces both classes to be nonempty
and forces an edge b->a with b in B and a in A: otherwise A dominates
S union B. The set S union {a,b} is strong via a->S->b->a.
If at least three vertices were outside S, this would be a larger proper
strong set. Therefore exactly two lie outside S, namely a and b.

Remove any vertex s of S. Every remaining x in S has a->x->b->a;
these triangles connect all remaining vertices in both directions.
Hence T-s is strong. This completes the proof, including the n=4 case,
where the initial directed triangle already has size n-1.

In particular, every strong tournament of order at least seven contains
an induced strong seven-vertex tournament. No unbounded computational
search or external asymptotic estimate is used in this step.

## 3. The complete finite head

`census.cpp` tests every Hamilton-frame presentation through order seven.
Subset endpoint dynamic programming and literal path recursion agree on
every presentation. Independent Python permutation counts check every
canonical class. Their exact strong minima are:

| Order n | Frame presentations | Ordinary classes | Minimum strong H |
|---:|---:|---:|---:|
|3|2|2|3|
|4|8|4|5|
|5|64|12|9|
|6|1024|56|15|
|7|32768|456|25|

The complete strong spectra needed below are

    order3: {3}
    order4: {5}
    order5: {9,11,13,15}
    order6: {15,17,19,23,25,27,29,31,33,37,41,43,45}.

The singleton has H=1; no two-vertex tournament is strong. Neither7 nor21
occurs among these strong classes. A strong tournament of order >=7 has
H>=25 by Sections1-2 and the seven-vertex finite minimum. Thus neither
value occurs in a strong tournament of ANY order.

## 4. Pass from strong tournaments to all tournaments

The strongly connected components of a tournament have a unique linear
dominance order. Every Hamiltonian path visits them in that order, and
choosing a Hamiltonian path in each component gives exactly one path in T.
Therefore

    H(T) = product_C H(C).

If H(T)=7, some component must have H=7, because7 is prime. If H(T)=21,
every component count divides21, so a factor containing7 must be7 or21.
Section3 excludes both possibilities. Hence

    For every finite tournament T, H(T) is neither7 nor21.

More precisely, the complete attainable set with H<=21 is

    {1,3,5,9,11,13,15,17,19}.

The finite strong witnesses supply all these values; products of the strong
counts below21 introduce no others. This conclusion does not require an
oddness theorem imported from outside the packet.

## 5. The inherited 56-obstruction characterization, independently recovered

Define C={T:H(T)<=21}. Monotonicity makes C hereditary. A minimal forbidden
induced tournament is exactly one with H>21 and all one-vertex deletions
having H<=21. No such tournament has a singleton strong component.

A strong minimal obstruction has at most seven vertices: otherwise it
contains a proper strong seven-vertex subtournament with H>=25. The complete
finite census gives29 strong obstructions of order6 and13 of order7.

For a nonstrong minimal obstruction, each component belongs to C and so
has order3..6. There are at most three components, since four would contain
three components with product >=27 as a proper forbidden subtournament.
With three components, each must be a triangle: if one were larger, it
would contain a triangle, and three triangles would already be a proper
forbidden induced subtournament. This gives the one order9 obstruction.

With two components, if both have order>=4, both contain the unique strong
order4 tournament Q. Its ordered square Q->Q has H=25, so minimality forces
exactly that graph. Otherwise one component is a triangle. Its partner must
have H>7, hence order>=5. It cannot have order6: a strong five-vertex
subtournament gives a proper forbidden product >=3*9. Thus the partner
is one of the six strong order5 classes. Either component order is allowed,
giving twelve further graphs. Their every vertex deletion has H<=21,
since deleting from the five-vertex part leaves H<=5, and deleting from
the triangle leaves the partner's H<=15.

The resulting basis has counts29+13+13+1=56. Our regenerated adjacency
sets agree exactly, up to isomorphism, with the pre-existing
[Session133 basis](../session133_fractran_spectrum_refinement_20261101/spectrum/README.md).
This is an independent check of that theorem, not a discovery claim.
At arbitrarily large order, absence of these56 induced graphs is equivalent
to H<=21 (and therefore to H<=19). The bound9 is sharp: three ordered
triangles have H27 and every one-vertex deletion has H9.

## Boundary of the result

The forbidden values are excluded globally. Proving that they are the ONLY
missing positive odd values is a different, still-open existence problem.
It appears as Conjecture6.2 of
[Chan--Pak, Computational complexity of counting coincidences](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).
The inherited [two-prime-lane reduction, THM-4094](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
identifies the unclosed targets. A finite forbidden basis for one fixed
threshold does not give witnesses for all larger counts.
