# Exact pure-blue counts and branch geometry

Status: PROVED corollaries, with the structural dependencies stated below.
The fixed-path complement metagraph statements use n>=3. Orders one
and two are included in the counting series only, not silently treated
as free complement-line graphs.

## Dependencies and unique normal form

The uniform prefix-rerouting theorem in this packet classifies the
nontrivial strong pure-blue atoms as C3 and the regular C5.
Two additional proved inputs are used here:

1. The SCC center law, THM-2454,
   `THM-2454-pure-blue-alphabet-and-center-law.md`: a self-converse
   ordered stack is pure-blue exactly when its paired outer components
   are path-rigid and its optional central component is pure-blue.
2. The all-order rigid-strong classification in the companion
   `FINITE_ATOMS_AND_RIGID_CLOSURE.md`: the rigid strong alphabet is
   exactly {singleton,C3}. This is separate from pure-blue classification,
   and uses the cited tournament Hamiltonian-bypass theorem.

Consequently every pure-blue class has the unique normal form

    T = A -> Z -> A^op,

where A is an arbitrary ordered word in singleton and C3 components,
A^op is its converse, and Z is either empty, singleton, C3, or C5.
The empty center occurs exactly at even total order. This is a statement
about isomorphism classes; all these classes are self-converse, so
passing to converse-merged nodes introduces no further factor.

Let a_j count words in letters of weights 1 and 3. Set a_0=1,
a_j=0 for j<0, and

    a_j = a_(j-1) + a_(j-3),  j>=1.

## 1. All-order generating function and exact counts

Including the empty object, the pure-blue node counting series is

    P(z) = (1+z+z^3+z^5)/(1-z^2-z^6).

Writing p_n for the coefficient,

    p_(2m)   = a_m,
    p_(2m+1) = a_m + a_(m-1) + a_(m-2).

Thus p_n=p_(n-2)+p_(n-6) for n>=6, using the zero convention for
negative indices. The first counts at n=3,...,16 are

    2,1,3,2,4,3,6,4,9,6,13,9,19,13.

These are no longer bounded-census predictions.

## 2. Every pure-blue node is a leaf or a trivalent branch

The fibre mass of T is H(T)/|Aut(T)|. Both H and |Aut| multiply over
the canonical ordered strong components. The paired rigid components
have mass one; the center masses are

    empty: 1,  singleton: 1,  C3: 1,  C5: 3.

For n>=3, complement acts freely on tilings. Therefore the incident
half-edge mass of every pure-blue node is exactly one or three, and
all those half-edges are blue.

The mass-one node counts L_n and mass-three node counts J_n are

    L_(2m)=a_m,                     J_(2m)=0,
    L_(2m+1)=a_m+a_(m-1),          J_(2m+1)=a_(m-2).

The word "leaf" below refers to mass one and no loop. The word
"trivalent" refers initially to half-edge mass three; Section4 proves
that these are in fact three distinct neighbours, not repeated lines
to one or two neighbours.

## 3. All pure-blue neighbours are mixed for n>=4

The companion endpoint-transport theorem
`endpoint-transport-and-complement-leakage.md` proves:

- A tournament's complement-neighbour triangle count is constant
  exactly when its source and sink strong components are regular.
- For n>=4 the resulting triangle-current-rigid nodes form an
  independent set and have no complement loops.

Every pure-blue normal form has regular source and sink components:
they are singleton or C3 unless the whole tournament is C5, which is
regular as well. Thus every pure-blue node is triangle-current-rigid.
No complement line at a pure-blue node can end at another pure-blue
node or be a loop, for n>=4.

Because that line is blue, its other end is blue-touched. A
blue-touched node that is not pure-blue is mixed. Hence:

**At every order n>=4, each pure-blue node is a blue leaf or a
three-branch point, with every branch ending at a mixed node.**

The bound is sharp. At n=3 the transitive tournament and C3 are both
pure-blue mass-one nodes, and their unique complement line joins
them. This pure-blue-to-pure-blue exception must not be erased.

## 4. Three branches have three distinct mixed neighbours

Consider a trivalent pure-blue node

    T=A -> C5 -> A^op,

where |A|=r, so n=2r+5. If r=0, the three normalized Hamiltonian-path
orbits of C5 are represented, for S={1,2}, by

    (0,1,2,3,4), (0,2,3,4,1), (0,2,4,1,3).

Flip all arcs except the four path arcs. The resulting tournaments
have respectively

    (H,|Aut|)=(11,1), (15,3), (13,1).

In particular they are pairwise nonisomorphic, even after converse
merging since H is converse invariant.

Now assume r>=1. In the original T:

- each left vertex has score at least r+5;
- each right vertex has score at most r-1;
- each of the five central vertices has score exactly r+2.

Along any Hamiltonian path, all left vertices precede all center
vertices, which precede all right vertices. In the complement partner,
an interior vertex of original score s has new score n-1-s.
The first path vertex instead has score n-s; the last has n-2-s.
Thus the five central vertices retain score r+2, while every left
vertex has score at most r and every right vertex has score at least
r+4.

It follows that the central five vertices are intrinsically identified
in the partner as **exactly the vertices of score r+2**. This
identification survives converse: r+2 is the middle score because
n-1=2r+4.

The induced tournament on those five vertices is exactly the
size-five complement partner associated with the central Hamiltonian
path. The three central path orbits have the three distinct H-values
11,15,13 computed above. Consequently their full complement partners
cannot be isomorphic or converse-isomorphic.

The outer rigid components provide no further path-orbit choices,
by the center factorization. Therefore every trivalent pure-blue
node has **exactly three distinct mixed neighbours**, connected by
one blue line each.

This is stronger than a degree count. The central five-vertex
score shell retains information that the full score multiset alone
would lose: the three induced center types have the same score
multiset {1,2,2,2,3} but different Hamiltonian-path counts.

## 5. Exact line budget, without claiming the whole arrangement

For n>=3, the number of blue tilings is

    2^floor((n-1)^2/4),

because each orbit of the staircase reflection gives one free blue
bit. Hence the total number of blue complement lines is

    beta_n=2^(floor((n-1)^2/4)-1).

At n>=4, precisely

    L_n + 3 J_n

blue lines have a pure-blue endpoint, and all of them are
pure-blue-to-mixed. The remaining

    beta_n - L_n - 3 J_n

blue lines are mixed-to-mixed, allowing loops and multiplicities
inside the mixed locus.

At n=3, instead, the single blue line has two pure-blue endpoints.
The expression L_3+3J_3=2 counts its two half-edges, not two lines.

These statements fix the pure-blue periphery exactly but do not
claim that counts, degree sequences, endpoint inventories, or
transport moments uniquely determine the entire merged metagraph.
The mixed-to-mixed wiring and further realizability constraints remain
separate.

## 6. A size-three count bijection

The identity

    p_(2m+1)=a_(m+2)=p_(2m+4)

has a direct normal-form bijection. Encode an odd pure-blue class
by its left word A and its center. Replace them by an even pure-blue
class whose left word is

    center singleton: A11,
    center C3:        A3,
    center C5:        A31.

Here 1 and 3 denote the singleton and C3 letters, not vertex labels.
Each new left word has weight m+2. Conversely, every word in 1 and 3
of weight at least two ends uniquely in one of the disjoint suffix
classes 11, 3, or 31. Removing the suffix recovers the center and A.

This is a size-changing bijection of pure-blue classes. It is not
asserted to be a graph isomorphism or to preserve complement
neighbours across different orders.
