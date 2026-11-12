# Exact crossing families and nonunique seed ancestry

**PROVED scoped.** The physical obstruction in README J2 belongs to an
infinite strong family with nonunique ordinary birth ancestry at every
order at least eleven. Its two birth rays intersect at a real seed.
The intersection must be counted and its endpoint data retained.

## C1 Classification of two run crossing tensors

Suppose the complete two-port defect has at most two runs per port. Set
x=2^a and y=2^b. Its unique exponential form is

```
delta(a,b)=k*x*y+l*x+m*y+d.
```

Suppose both directions are individually stable at positive integer block
lengths (A,B). The identity delta(a,B)=1 for every positive a forces
k*2^B+l=0 and m*2^B+d=1. The other identity forces k*2^A+m=0.
Consequently either k=0 and delta is identically1, or

```
delta(a,b)=1+k*(2^a-2^A)*(2^b-2^B),     k != 0.         (C1)
```

Conversely this formula gives the two stable directions at (A,B).
Its entire positive integer unit locus is a=A or b=B. This is a complete
classification within the declared two-run class, not a classification
of general multi-mode tensors.

For the physical J2 quotient, port0 has a unique out-neighbour, selected
vertex1, and port8 has a unique in-neighbour, selected vertex5. Both remain
single vertices while the ports vary. The selected-endpoint bound is
therefore uniformly R<=2 in each coordinate, even when the other module
grows. Thus four literal corner defects determine the entire function.
The corner values at (1,1),(1,2),(2,1),(2,2) are5,1,1,1, giving
k=1 and A=B=2 in (C1). The independent full skeleton supplies the same
identity without using this four-corner shortcut.

This is also a sharp boundary for small screens: when the uniform two-run
bound is absent, three unit corners do not suffice. The complete order-nine
quotient probe retains two converse hostile presentations whose later
single-axis defects are1,1,-35,-215.

## C2 An infinite strong nonconfluent family

Let Q and O be the J2 quotient and selection. Write T(a,b) for substitution
of TT_a and TT_b at its ports0 and8. Q is strong, has no unselected twins,
has trivial marked automorphism group, and is marked self-converse. These
claims are checked directly from its nine rows. Substitution of nonempty
modules preserves strength. Its canonical marked unselected-twin quotient
is Q with weights a,b and all other weights1. Hence distinct ordered pairs
(a,b) yield nonisomorphic marked tournaments.

By (C1), T(a,b) is a unit seed exactly when a=2 or b=2. The only nonsingleton
unselected twin blocks are the two named modules. At (a,2), the first block
is stable; the second is stable exactly when a=2. The symmetric statement
holds on the other line. No singleton exterior vertex can be contracted,
regardless of whether it happens to acquire stability.

Therefore the complete eligible-contraction graph inside this family is
known at every size. Its terminal births are exactly T(1,2) and T(2,1).
T(2,2) contracts to either. Every other unit state of order at least eleven
contracts along its variable block to T(2,2), and thus reaches both births.
This proves nonunique ordinary birth ancestry in every order n>=11.

Since |T(a,b)|=7+a+b, the exact marked-class series is

```
2*z^10/(1-z)-z^11.
```

The two infinite birth cylinders intersect once, at T(2,2). Their
nonconfluent subfamily has series z^11+2*z^12/(1-z), and their contraction
edges have series2*z^11/(1-z). The ordinary state counts are2 at order10,
1 at order11, and2 at every order>=12. The quotient involution exchanges
(a,b) and(b,a), giving one converse orbit at every order>=10.

The previous conditional local-diamond theorem is still useful on charts
where the chosen pair is jointly stable. It cannot make this crossing
confluent. A correct general implementation retains all terminal births
and all cylinder intersections; it must not choose one root as an invariant.

## C3 The scalar count hides the continuing interface

The exact Hamiltonian count and selected slope are

```
H(a,b)=175*2^(a+b)-9*2^a-9*2^b+1,
U(a,b)=-8-49*2^b-95*2^a+334*2^(a+b)
                    -64*3^a+506*3^a*2^b.               (C2)
```

The H formula follows from the uniform two-run bound and four endpoint
counts; U comes from the complete independent endpoint skeleton. Both are
checked on literal inflated graphs. On either unit line, H=691*2^m-35.
Thus the physical growth step has the exact scalar recurrence H_next=2H+35.
But the two unit lines have different slope profiles:

```
U(a,2)=-204+1241*2^a+1960*3^a;
U(2,b)=5841*2^b-964.
```

The mode asymmetry has a concrete endpoint cause. At port0 the single
out-neighbor bounds H and delta by two runs, but U allows both paths to
end in the port, permitting three runs. At port8 the single in-neighbor
still bounds U by two runs because its second path starts selected.
Thus the respective U mode budgets are(3,2), explaining the3^a*2^b term
and the different continuing profiles under converse. These integer
bases count run allocations; primality alone is not the mechanism.

At the terminal births, both H values are1347, while U is8158 or10718.
Their odd parts are4079 and5359. At the shared state T(2,2), U=22400,
whose odd part is175 and whose dyadic valuation is7. Neither H nor a
converse-merged state determines the selected tail coefficient.

For b>=3 the second line gives the seed-preserving coefficient family

```
c_k=5841*2^k-241,  k=b-2>=1,
U=4*c_k,          c_(k+1)=2*c_k+241,
source order=k+11.
```

The inherited selected-tail construction realizes1+c_k*2^s for every
s>=3. The two lower heads are now closed uniformly by explicit nonunit
endpoint interfaces and an H=9 factor; see [FORBIDDEN.md F5](FORBIDDEN.md).
Consequently every s>=1 is realized for every k>=1. This is a new explicit
recursive family from the obstruction, not coverage of all odd c and not
a minimum-order claim. It has an exact
inverse membership test: c+241 must equal5841 times a positive power of2.
The special crossing value c=175 has a later tail start s>=8 because of
its larger dyadic valuation. A small odd part does not imply a small
below-ray gap.

## C4 What survives under converse folding

The two births above are converse, so their scalar ancestry collapses to
one orbit after forgetting orientation. That fold loses the U values in
C3 and cannot replace the ordinary marked source for continuing tail
construction. It is a valid observer quotient only for predicates that
are converse invariant.

Even confluence modulo converse is false without the strong restriction.
Prepend an unselected singleton source to every state in C2. H and W are
unchanged: the source must begin the first path of every taxed cover, since
the second path starts at a selected vertex. Thus the defect and old
contractions survive. The singleton prefix is a stable but noncontractible
SCC and has no twin in the strong suffix. The two terminal states now have
SCC order (singleton,strong); conversing either reverses that order, so
they are not converses. This gives nonconfluence even after converse
folding at every order>=12, in the unrestricted tournament class.

Strong confluence modulo converse is a separate OPEN statement; it is not
needed for enumeration by exact cylinder unions. The complete fixed-mark
one-vertex extension probe of Q finds no strong extension preserving all
three necessary unit corners. That failed local repair does not classify
other cores or more general operations.

## C5 The exact local transport test is weaker than full joint stability

Suppose two eligible blocks have current lengths(A,B), both at least2,
and let f(a,b) be the defect with the rest of the marked graph fixed.
Individual stability gives f(a,B)=f(A,b)=1 at every positive length.
After contracting the second block, the first remains stable **iff**
f(a,B−1)=1 for every a>=1. The corresponding condition after the other
contraction is f(A−1,b)=1 for every b>=1. These two univariate identities
are exactly the test that the two elementary contractions form a legal
commuting diamond. Full joint stability implies them, but is stronger.

If this local test holds at every fork in a reduction-closed class,
termination and the local-diamond argument give ordinary confluence.
Failure of the test at a particular fork does not by itself exclude a
later common descendant through other blocks; the explicit J2 roots
prove nonconfluence because they are terminal and distinct.

A unit double-contracted corner alone is insufficient algebraically.
For example the **formal**, not asserted physical, function

```
f(a,b)=1+(2^a-4)*(2^b-4)*(2^a-2)
```

has stable axes at(2,2) and f(1,1)=1, but f(a,1) is not identically1.
Thus the missing condition is stable-line transport, not just the four
corner values, once the uniform two-run hypothesis is removed. A different
formal example1+F_3(a)F_3(b), with F_3(m)=3^m−3*2^m+3, has the legal
diamond at(2,2) but fails global joint stability. These separate the
quantifiers a repaired proof must keep.

For physical two-run crossings(C1), nonzero k and A,B>=2 always fail
both adjacent transport identities. For higher-mode charts use the full
univariate run bound at each contracted child. This is a cheaper sufficient
confluence certificate than computing every multivariate coefficient.

## Map and verification contract

Source: a marked strong quotient with two transitive modules. Target:
literal tournaments, their defect, H, U, and contraction graph. Map:
substitution of positive chains. Preserved predicates: strength and the
stated exact unit-locus equation. Lost data under H or converse folding:
the selected slope and ordered block weights. Required sidecar: the full
marked quotient, weights and endpoint tensor. Cheapest hostile: the four
corners5,1,1,1; the nonunit corner is essential.

Reproduction: contractible_obstruction.py and crossing_family.py. The
two-run theorem is algebraic; the physical constants are exact finite
certificates, checked by independent skeleton and literal algorithms.
