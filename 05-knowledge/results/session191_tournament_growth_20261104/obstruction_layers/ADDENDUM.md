# The complete layer one vertex below the sharp obstruction ceiling

**PROVED structural reduction + FINITE-EXACT small core classification;
independently reviewed.** This addendum uses Theorem A of `PROOF.md`.
It does not change that file or claim a full Hamiltonian-path spectrum proof.

Let 3^(r-1)<=B<3^r and theta=B/3^r. Let A_1(B) be the number of ordinary
isomorphism classes of induced-minimal H>B obstructions of order 3r-1.
For binomial coefficients below, set binom(m,k)=0 whenever m<k or m<0.

**Theorem.** For every integer B>=1,

    A_1(B) = 6*max(r-1,0) * 1_{theta>=5/9}
           + binom(r-1,2) * 1_{5/9<=theta<25/27}
           + binom(r-2,5) * 1_{theta>=625/729}
           + [2*binom(r-2,2) + 6*binom(r-2,3)
              + 8*binom(r-2,4)] * 1_{theta>=25/27}.    (1)

All threshold comparisons are rational and exact. The formula includes
r=1, when the layer has order two and is empty. It includes B=21, where
r=3 and A_1=12+1=13, agreeing with the independently certified old basis.

## 1. Complete template table

As before Q denotes the unique strong four-vertex class. Let:

* S5 range over all six strong five-vertex classes: H is9,9,11,13,15,15,
  and each has maximum vertex-deletion count5.
* M5 range over the two of these with H9.
* M6 range over the two strong six-vertex classes with H15; their maximum
  vertex-deletion count is9.
* M7 denote the unique strong seven-vertex class with H25; its maximum
  vertex-deletion count is15.

These are finite complete classifications checked by the supplied code,
not assumptions about generic strong tournaments. Every obstruction in
the layer has one of these exceptional SCC multisets; all other SCCs are
triangles, in any positions.

| exceptional SCC multiset | triangles a | threshold theta | number of classes |
|---|---:|---|---:|
| one S5 | r-2 | [5/9,1) | 6(r-1) |
| Q,Q | r-3 | [5/9,25/27) | binom(r-1,2) |
| M7,Q | r-4 | [25/27,1) | 2 binom(r-2,2) |
| M6,Q,Q | r-5 | [25/27,1) | 6 binom(r-2,3) |
| M5,Q,Q,Q | r-6 | [25/27,1) | 8 binom(r-2,4) |
| Q,Q,Q,Q,Q | r-7 | [625/729,1) | binom(r-2,5) |

A row is absent when a<0. The indicated binomial convention enforces this
automatically. Sum the rows to obtain (1).

## 2. Why the table is exhaustive

For deficit d=1, Theorem A gives kappa^(E-1)<9. Since kappa^6>9,
E<=6; the order congruence 3r-1=3*(number of SCCs)+E gives E=2 or5.
There is no pure-triangle obstruction at positive deficit.

The possibilities for the orders of exceptional SCCs are exactly

    E=2: (5), (4,4);
    E=5: (8), (7,4), (6,5), (6,4,4),
         (5,5,4), (5,4,4,4), (4,4,4,4,4).

For an exceptional list of total order s, put t=(s+1)/3 and a=r-t.
Write D(C)=max_v H(C-v). Minimality requires the maximum deletion count
divided by3^a to be less than3^t, since B<3^r.

* **(5).** D(S5)=5, and H(S5)<=15 makes triangle deletion no larger than
  5*3^a. Also H(S5)>=9. Thus the row begins at5/9 and persists to1.
* **(4,4).** H=25*3^a, maximum deletion15*3^a. Thus the exact threshold
  interval is[15/27,25/27).
* **(8).** This would require D(C)<27. No strong eight-vertex tournament
  has that property, by the finite certificate in Section3 below.
* **(7,4).** Deleting in Q gives3H(C7)<81, so H(C7)<=26. The complete
  order-seven census leaves only M7. Its product count is125 and its
  maximum deletion is75, giving75/81=25/27.
* **(6,5).** Strong deletion inside C6 leaves a strong five-vertex class,
  so D(C6)>=9. Multiplication by H(C5)>=9 yields a deletion count at
  least81, contradicting the strict bound3^t=81.
* **(6,4,4).** Deleting in Q forces15H(C6)<243, hence H(C6)<=16. The
  complete six-vertex census leaves exactly the two M6 classes. Their
  product count is375 and maximum deletion225, giving225/243=25/27.
* **(5,5,4).** Deleting in Q leaves a count at least3*9*9=243, impossible.
* **(5,4,4,4).** Deleting in Q forces75H(C5)<729, hence H(C5)=9. The
  two M5 classes have product count1125 and maximum deletion675, giving
  675/729=25/27.
* **(4,4,4,4,4).** The product count is3125 and maximum deletion1875,
  giving1875/2187=625/729.

For every surviving row except(4,4), the product count divided by3^t is
at least1, so it remains outside H<=B throughout its indicated interval.
For every row, triangle deletion is dominated by the stated maximum
deletion of an exceptional component, including at a=0.

The exceptional skeletons are intrinsically ordered. For M7,Q there are
two orders; for M6,Q,Q there are two choices of M6 and three positions;
for M5,Q,Q,Q there are two choices of M5 and four positions. Triangle
gap counts then give the binomial factors in the table. No converse
classes are merged.

## 3. Small complete certificates

`verify_next_layer.py` independently generates the full fixed-Hamilton-path
cover through order seven, obtaining ordinary class counts
1,1,2,4,12,56,456. Degree-cell canonicalization checks every permutation
that could give an isomorphism. It counts H and every one-vertex deletion
by endpoint subset DP. It checks the required M5/M6/M7 data and literally
counts all permutations for the unique M7 witness and all its deletions.

For the missing order-eight template, let C be strong with D(C)<=26.
Strong deletion gives a strong seven-vertex C-v with H<=26. The complete
census shows this must be M7. Therefore C occurs among the128 one-vertex
extensions of the displayed M7 representative. The code tests all128
literal adjacency extensions and finds none strong with D<=26. This
proves the exclusion at order eight, not merely a sample assertion.

Finally the verifier exhausts every exceptional-class choice in the nine
size patterns, calculates its exact admissible threshold interval, and
compares it with the table. It checks the count formula at every rational
endpoint and in every intervening chamber for r through40, with hostiles
on the disappearing Q,Q family at theta=25/27.

## 4. The change with size

The largest-order obstruction layer has only the two families of
`PROOF.md`; one vertex lower already has six families and three rational
transition points. The threshold25/27 is an actual replacement boundary:
Q,Q ceases to be an obstruction there because its H is now allowed, while
three different mixed families become minimal obstructions. Consequently
obstruction counts need not increase at each threshold, even though the
hereditary classes H<=B increase. Equal cardinalities of two obstruction
bases still say nothing about equality of those bases.

These results describe a complete size stratum at every threshold, while
leaving the global existence problem for prescribed H untouched.
