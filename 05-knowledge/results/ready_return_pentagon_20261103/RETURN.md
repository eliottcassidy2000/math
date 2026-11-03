# Both singleton-ladder branches have literal graded returns

**PROVED for every k>=1**, with exact independent controls. These are
returns through births, global converse and a retained root, not a
closed construction of every ternary descendant.

## R1. Exact maps on the actual tournaments

Let E_k be the word11^k, with ports a=0,b=1, chain2,...,k+1 and root
r=k+2. It has h=3*2^k-1. The inherited singleton-core calculation gives
two ready sources: E_k^op at b is1-ready, and E_k at r is2-ready.

In either source add x beaten only by the selected root, then y beaten
only by x. Take global converse. The following retained old root and
literal vertex order identify the returned tournament:

|Input / selected root|Retained return root|Order after converse|Returned word|H|
|---|---|---|---|---|
|E_k^op / b|r|a,b,u1,...,uk,y,x,r|11^k,00,01|9*2^k-1|
|E_k / r|a|r,b,uk,...,u1,y,x,a|10^k,00,10|9*2^k+1|

To prove the identities, inspect the birth arcs and reverse them. In
the first row the old prefix remains11, y loses to both in-ports, and
x has type01. In the second row the ordered in-ports are old r->old b;
each reversed chain vertex has type10, y has type00 and x type10.
The retained root beats the chain and loses precisely to the displayed
in-ports. Every old pairwise relation is accounted for, including the
chain reversal. These are isomorphisms with named vertex transport,
not equalities inferred from Hamilton counts.

The two input phases agree exactly when k=2 mod6, as in the inherited
ladder. Other k retain their different phases. Neither return is claimed
to preserve the input phase; the receipt records input, sealed and
returned phases separately.

The complete returned profiles are

    minus: (h,s,q,F)=(9*2^k-1,3,3,9*2^k+3),
    plus:  (h,s,q,F)=(9*2^k+1,3*2^k-1,3*2^k+3,12*2^k+5).

Thus the minus return lies in the complete neutral C4 atlas. The plus
return has q-s=4, F-h=3*2^k+4 and joint defect delta=6*2^k. Its growing
defect is explicit, rather than hidden in a scalar source reset.

## R2. The joint transport is a two-live-port program

Start with the retained spare root in slot0 and selected root in slot1:
(r,b) for the minus input, (a,r) for the plus input. Use the inherited
eighteen endpoint-skeleton polynomials, whose degree counts anonymous
path components. A command is `(polarity,mask,discard)` as defined in
[graded_carrier.py](../source_selection_closure_20261103/graded_carrier.py).

The two births are exactly `(0,2,1)` and `(0,2,2)`: replace slot1 by x,
then forget the newborn y while keeping (spare,x). Converse swaps each
component's start/end slots, retaining its multiplicity and anonymous
degree. Select the **already retained** slot0. No forgotten old vertex
is freely renamed back into the interface.

Initial E_k construction is also literal: start a->b, append k chain
vertices using `(1,0,2)`, then add r using `(0,3,discard)`. The last
discard selects the required spare/active pair. Converse the source
for the minus row. This proves an exact program on full joint covers.

The two lanes start with different two-port markings. A **common** input
that can select either lane is supplied by retaining all three physical
ports(a,b,r) in E_k. Beginning with E_1, each prepending11 birth uses
the three-port command `(0,4,3)`. Project its joint covers to(r,b), then
converse for the minus lane, or project to(a,r) for the plus lane.
Projection forgets endpoint names and merges anonymous components with
their multiplicities; it never combines independent marginals. The
script verifies both full joint projections. Three common ports are
sufficient here; no minimality or two-port common-selector theorem is
claimed.

The integer polynomial degree is unbounded. The two slots and eighteen
skeleton types do not mean eighteen scalar states. We retain the actual
tournament as well; the carrier does not certify arbitrary deletion or
an unnamed native re-marking. Global converse is explicitly an operation
here. The inherited stabilized arc compiler can realize it through native
complements with temporary vertices, but no claim is made that converse
is one native line or that its intermediate H values stay fixed.

## R3. Independent run skeletons prove the uniform observer formulas

A transitive k-clone block with fixed external neighbourhood can occur
in j ordered nonempty runs among m ordered paths. Once the external
vertices and run slots are specified, assigning the k ordered clone
vertices surjectively to the runs gives `j! S(k,j)` realizations; each
run is internally increasing. With e external vertices, j<=e+m.
Thus each fixed-bin observation has an integer expansion in these
surjection counts. This is an exact skeleton decomposition, not an
interpolation assumption.

For our returned sources e=5. Enumerating their finite external/run
skeletons gives the following nonzero coefficients at j=1,2,3; j0
terms do not contribute for k>=1:

|Return|H coefficients|E(1) coefficients|S(2) coefficients|
|---|---|---|---|
|minus|17,9|1|33,24|
|plus|19,9|2|52,50,14|

Using `1!S(k,1)=1`, `2!S(k,2)=2^k-2` and
`3!S(k,3)=3^k-3*2^k+3`, the golden observer v=(h,X,Y) is

    v_minus=(9*2^k-1,12*2^k-8,12*2^k-7),
    v_plus =(9*2^k+1,7*3^k+4*2^k-4,7*3^k+4*2^k-2).    (R1)

The coefficient7 is explained by the fourteen admissible three-run
S(2) skeletons and the observer's division by two. The minus return
has no three-run term. Thus the direction erased at characteristic7
comes from an actual incidence multiplicity, not a fitted residue law.

The script counts skeletons directly and independently compares the
joint polynomial program and ordinary subset path DP. All five Q4
presentations are applied both by their inherited integer matrices and
by literal tournament attachments on the smallest returned sources.

## R4. Exact field phase, Frobenius, and the extra precision bit

Keep the established convention F49=F7[theta]/(theta^2-4theta-2) and

    alpha=a+b theta,  a=4(h+2X-Y), b=Y, g=h-Y mod7.

The pair(alpha,g) is equivalent to v modulo7; it is not reduction
modulo49. Writing t=2^k mod7, R1 gives

    minus: (a,b,g)=(2,5t,4t+6),
    plus:  (a,b,g)=(3t+1,4t-2,5t+3).                    (R2)

Both have exact period3. The term7*3^k is invisible here. All five Q4
phase actions are retained; in particular q2 multiplies alpha by theta
and fixes g. Frobenius is `(a,b,g)->(a+4b,-b,g)` and conjugates theta
transport to multiplication by4-theta. The script tests this relation.
Advancing k is a separate clock, not Frobenius.

There is an all-precision refinement. For the exact integer observer
v modulo7^a, the minimal periods are

    minus: 3*7^(a-1),                         a>=1;
    plus:  3 if a=1; 6*7^(a-1) if a>=2.        (R3)

Proof: ord_(7^a)(2)=3*7^(a-1) and ord_(7^a)(3)=6*7^(a-1), by the
base orders and valuations v7(2^3-1)=v7(3^6-1)=1. In the plus row the
H coordinate forces the full2-power period; the remaining term forces
the3-power period modulo7^(a-1). Their least common multiple is R3.
The minus row is affine in2^k with a unit H coefficient.

For a>=2, fixing the dyadic phase moduloL=3*7^(a-1) leaves exactly two
plus-observer values. Advancing by L changes3^k to its negative modulo
7^(a-1), so one binary sign distinguishes them. Thus one additional bit
is necessary and sufficient beyond that dyadic clock for this integer
observer, at every7-adic precision. It is not an identification of this
sign with Frobenius, nor a reconstruction theorem for full joint covers.

## R5. A synthetic phase closure and its trace-zero correction

There is a useful separately scoped operation model. Extend the growth
maps seen on R2 affinely to the entire phase space, and adjoin field
Frobenius sigma:

    G_minus(a,b,g)=(2a-2,2b,2g-6),
    G_plus (a,b,g)=(2a-1,2b+2,2g-3),
    sigma(a,b,g)=(a+4b,-b,g).                           (R4)

G agrees with prefix growth on the actual returned-family phases.
We do **not** have a tournament realization of sigma or G on every
new synthetic state, nor closure under all five Q4 operations here.
With this precise scope, the generated phase orbits have6 and21
states, and their permutation groups have orders6 and42, respectively.

For the minus map G and sigma commute. Its three growth phases have
nonzero theta coefficient, so adjoining conjugation doubles the orbit.
For the plus map the commutator is

    sigma G sigma G^(-1): (a,b,g) -> (a+1,b+3,g).

The translation delta=1+3theta has trace zero, sigma(delta)=-delta,
and its seven multiples form the entire trace-zero line. G scales
that line by2; sigma negates it. These generate the faithful group
C7 semidirect(C3 times C2), of order42. The g coordinate cycles through
three different values, and each fibre contains seven translated
phases. Hence its orbit has exactly21 states. This is a literal
mechanism for7*3 in this phase language, not an explanation for the
forbidden Hamilton count21.

This controller has a literal three-to-one Paley quotient. Define

    phi(a,b,g)=5b+3 in F7.

Then phi(Gp)=2phi(p), phi(sigma p)=-phi(p)-1, and the commutator translates
phi by1. Every quotient point has three phase states. On F7 put the
Paley arcs x->y when y-x belongs to{1,2,4}. Growth and the commutator
preserve these arcs; Frobenius sends them to their converses. The odd
subgroup is exactly the classical Paley7 affine automorphism group of
order21. The full controller acts as all42 affine maps. This supplies
an equivariant map, rather than matching equal cardinalities.

Pairs in one phi fibre are ties; no arbitrary tournament orientation
is assigned to them. The map forgets a three-valued phase coordinate.
It is a quotient of the **synthetic controller**, not a map proving
anything about the Hamilton count of the original returned source.

Finally {theta,theta^7} is a normal basis. Coordinates
`(b+2a,2a)` of alpha in that basis are swapped by sigma. The model
therefore connects rotation, Frobenius and a nontrivial loop correction
without identifying any of them with an unproved source operation.
An exact permutation closure in `return_program.py` checks both groups.
Lifting these trace-zero corrections to actual returned sources is a
concrete additional construction target.

## R6. What remains open

These maps return **both** branches of an infinite ready ladder to a
uniformly described actual source class. They do not return to the old
ready interface; O3 proves that impossible at H35. A repeatable source
machine must act on the enlarged grade and still produce exact offsets2
and4, or use another inverse arithmetic language. No such all-size
closure, complete prime/7p supply, or full {7,21} classification is
claimed here.
