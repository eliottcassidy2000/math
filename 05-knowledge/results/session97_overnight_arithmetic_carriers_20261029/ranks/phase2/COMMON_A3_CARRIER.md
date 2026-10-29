# One exact A3(2) carrier inside three different rational K3 surfaces

**PROVED conditional only on the standard no-reducible-fibre height formula;
all model hypotheses checked exactly in the companion geometry audit.**

Let a rootless ellipticK3 fibration be written

    y²+(Lx+B)y=x(x+D)(x+E), B=pE+qD+pq(L−p−q).

Take the integral sections P=(0,0), Q=(−D,0), R=(−pq,p(pq−E)).
Suppose the homogeneous intersection counts are

    P·Q=4,  P·R=2,  Q·R=0.

Then their height Gram matrix is

    [ 4 −2  0 ]
    [−2  4  2 ]
    [ 0  2  4 ].

Indeed every integral section has height4, and the off-diagonal pairing is
2 minus its section-intersection number. NegatingR identifies this matrix
with twice the A3 Cartan matrix. Thus P,Q,R span an A3(2) lattice of
determinant32 and are independent. No specialized high-rank fibre is needed
to prove that independence.

The hypotheses are algebraic guards, not generic handwaving. In affine
coordinates they are checked by polynomial gcds, plus the projective
intersection term

    min(4−deg(xP−xQ), 6−deg(yP−yQ)).

For example gcd(D,0)=D gives the first intersection4; gcd(pq,p(pq−E)) has
degree2 in the named models, while gcd(D−pq,p(pq−E)) has degree0. Every
infinity contribution is0. `check_common_a3_carrier.py` verifies all three
sets of guards exactly for ICARM273,302,4239 and obtains the same matrix.
If the guards are dropped, even p=0 collapsesR=P and destroys rank3.

## What the carrier does not retain

The three ambient K3 surfaces have different F233 point counts
58509,58404,58374 at common good reduction, so they are pairwise
non-isomorphic overQ. Therefore this same marked abstract height lattice
cannot determine the ambient rational K3 surface. Its local intersection
information is exact, but it loses the Frobenius-sensitive global geometry.

This supplies an explicit source/target/loss ledger:

* source: a named ellipticK3 plus the three polynomial sections;
* target: their integral height lattice;
* preserved: all pairwise canonical heights, rank3, determinant32;
* lost: the ambient Frobenius representation and hence rational surface
  isomorphism data;
* hostile control: the three audited record-fibration formulas;
* sufficient restoration sidecar for this particular distinction: one named
  good-prime point count; not a universal reconstruction theorem.

The A3 occurrence is structural and verified, unlike a comparison based only
on a matching integer. It does not identify these lattices with A1/A2 cost
rules, tournament forbidden values, or Selmer quadratic sample spaces.
