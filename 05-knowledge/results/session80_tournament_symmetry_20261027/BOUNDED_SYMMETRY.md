# Stronger support theorem: sharp at every even excess

**Status: PROVED, independently audited.** No literature novelty or optimality at odd excess is claimed. The
structural input remains canonical tournament modular decomposition and the
[decorated symmetry law](GROUP_LAW.md).

For a nonempty tournament T with prime-node height h>=1 and n vertices, put
r=n-2h-1>=0. Then the union of the supports of Aut(T) has size

    s(T) <= floor(3(n-2h+1)/2) = floor(3r/2)+3.             (1)

Consequently Aut(T) embeds in S_(floor(3r/2)+3), and at fixed r there are finitely
many possible abstract groups and moved permutation actions up to isomorphism.
The full vertex action may have arbitrarily many fixed points. If h=0 the
tournament is transitive and its automorphism group is trivial. The earlier
bound s(T)<=3(n-2h) follows as a weaker corollary (also at h=0).

## Prime quotient lemma

If an automorphism of a prime tournament fixes a vertex and moves at most five
vertices, it is the identity. Every cycle of a tournament automorphism has odd
length, since an even cycle has an antipodal pair reversed by a power. Thus a
nonidentity permutation moving at most five vertices consists of exactly one
3-cycle or one 5-cycle, together with fixed points. Every fixed point sees that
single cycle uniformly by arc preservation. Its support is therefore a proper
nontrivial module, contradicting primeness. The fixed-vertex hypothesis is
essential: a prime C3 or regular cyclic C5 can admit a full cycle.

## Induction with the small-outside cases explicit

It is enough to establish the real-valued bound 3(n-2h+1)/2: s(T) is an integer.
Choose a maximal-height canonical child A_j, of size n_j and height h_j;
let q=n-n_j. Use induction on n whenever the child's height is positive.

At a linear root h_j=h>=1. Its quotient has no automorphisms, so

    s(T)<=s(A_j)+q<=3(n_j-2h+1)/2+q<=3(n-2h+1)/2.

At a prime root h=h_j+1 and q>=2. If h_j=0, then h=1 and n>=3; hence the
trivial bound s(T)<=n<=3(n-1)/2 is sufficient. Now assume h_j>=1 and write H
for the decorated quotient stabilizer.

If H moves j, its orbit has odd size at least three, and all children on that
orbit are isomorphic. Using n_j>=2h_j+1 gives

    n>=3n_j>=6h_j+3,  hence  n<=3(n-2h_j-1)/2.

This again proves the desired bound from s(T)<=n.

Suppose instead H fixes j. Restriction to A_j lies in Aut(A_j), and the
inductive bound leaves the outside support allowance

    3(n-2h+1)/2 - 3(n_j-2h_j+1)/2 = 3(q-2)/2.

For q>=6, all q outside vertices fit this allowance. For 2<=q<=5, H is
trivial by the prime quotient lemma: it fixes j and can move at most q other
slots, since every outside child is nonempty. Thus outside support arises
only internally in outside children. There are at least two such children,
so each has size at most four. On one or two vertices all automorphisms are
trivial; on three or four vertices at most three vertices move. The latter
follows from odd group order and its possible orbit sizes: a nontrivial orbit
has odd size at least three, and there is no room for two such orbits. Thus:

| q | Maximal outside support | Allowance 3(q-2)/2 |
|---:|---:|---:|
| 2 | 0 | 0 |
| 3 | 0 | 3/2 |
| 4 | 3 | 3 |
| 5 | 3 | 9/2 |

For q=5 two nontrivial child supports are impossible because each would need
at least three vertices. This exhausts the cases and proves (1).

## Sharp every even excess, with arbitrary added depth

Start U=C3, and apply U <- C3[U,C3,1] exactly k times. Each step adds four
vertices, increases prime height by one, and adds an independent C3 action
on three moved vertices. The root quotient stabilizer is trivial: the singleton
slot is distinguished. This remains true at the first step when the other two
children are both C3; a C3 quotient has no nonidentity rotation fixing the
singleton slot. Therefore

    h=k+1,   n=4k+3=2h+1+2k,
    Aut(U) is C3^(k+1),   s(U)=3k+3.

Now apply any number of ordinary lifts U <- C3[U,1,1]. The non-singleton slot
is distinguished, so each lift adds two vertices and one unit of height,
leaving group and support unchanged. For every k>=0 and every h>=k+1 this
yields r=2k and equality s=3k+3=floor(3r/2)+3. This establishes sharpness for
every even r in the stated depth range, not a classification of all equality
cases or an optimal odd-r bound.

## Independent finite verification and freeze boundary

The augmented `check_support.py` checks both support bounds on all 33,867
labelled tournaments of orders 1 through 6 using actual full automorphism
permutations and directly extracted canonical trees. It retains the zero-excess
spines at heights 1..8 and adds equality constructions for k=0..3 and zero,
one, or two subsequent ordinary lifts (12 specified examples, up to order 19).
For these it verifies height, excess, support, and exact automorphism order
3^(k+1). The abstract direct-product identification is proved structurally above,
not inferred from group order alone. Normal/-O receipts match byte for byte.
No larger-order exhaustive census or optimal odd-excess bound is claimed.

Run `python3 check_support.py` and `python3 -O check_support.py` from this
directory. The checker is standalone and reads no saved data. Its larger
constructed examples use direct subset-module extraction, so their runtime
is intentionally higher than the small-order checks.

This result supplies a finite possible bank of moved symmetry actions at fixed
excess. It does not say that the aggregate prime cycle index determines which
action occurs, nor that the whole fixed-excess tournament family has finitely
many classes. Depth and the number of globally fixed vertices remain unbounded.
The [prime collision](SYMMETRY_INFORMATION.md) remains a separate information
barrier; the [polynomial family theorem](COMPOSITE_FAMILIES.md) gives a separate
counting simplification under the same fixed-excess restriction.
