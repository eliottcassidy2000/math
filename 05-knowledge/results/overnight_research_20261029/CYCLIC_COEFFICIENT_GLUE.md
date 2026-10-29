# From actual Paley7 coefficients to a clocked Eisenstein lattice

2026-10-29. **PROVED elementary integral-coordinate identities;
FINITE-EXACT native inputs replayed from incoming commit4ff64c524c.**
The decomposition is classical C3 representation arithmetic; no historical
novelty is claimed. Its value here is a precise source-preserving bridge.

## C1. The integral guard and clock

Start with an ordered integer coefficient triple(c0,c1,c2), retaining
which group element of C3 each coefficient belongs to. Define

    h=c0+c1+c2,  a=c0-c2,  b=c2-c1.                      (C1)

This is an isomorphism of Z^3 onto the index-three sublattice

    L={(h,a,b) in Z^3: h-a+b=0 modulo3}.                 (C2)

Indeed c2=(h-a+b)/3, c0=c2+a and c1=c2-b. The guard is both necessary
and sufficient. Over Q the division is allowed everywhere; over Z the
triple(1,0,0) is the smallest positive-sum unrealizable split state.

Cyclic re-marking(c0,c1,c2)->(c1,c2,c0) acts by

    h->h,  (a,b)->(-a-b,a).                              (C3)

The phase matrix is [[-1,-1],[1,0]]=A^(-T), with A the clock used in
the matched cross-join family. It preserves both Q(a,b)=a^2+ab+b^2
and(C2). This is an explicit intertwining map of labelled C3 modules,
not a scalar resemblance or a claim that two attachment languages agree.

If omega^2+omega+1=0 and eta=-omega, then

    c0+c1*omega+c2*omega^2 = a+b*eta,
    Norm(a+b*eta)=Q(a,b)
      =c0^2+c1^2+c2^2-c0*c1-c0*c2-c1*c2,
    h^2+2Q(a,b)=3(c0^2+c1^2+c2^2).                      (C4)

Keeping only Q destroys phase and h. Keeping h and phase but omitting
the guard invents integer source triples. Keeping(h,a,b) with(C2)
recovers every original coefficient and the clock for free.

## C2. The actual15/12 and11/13 sources

The incoming [two-run Paley7 module, R4](../overnight_rank_return_20261029/TWO_RUN.md)
has a regular order21 affine action. Summing its convolution kernel over
the seven translations produces the following C3 quotient coefficients:

| Run-exchange sector | (c0,c1,c2) | h | (a,b) | Q |
| --- | --- | ---: | --- | ---: |
| plus | (15,12,12) | 39 | (3,0) | 9 |
| minus | (-11,2,2) | -7 | (-13,0) | 169 |

Both sources obey(C2); the inverse gives c2=12 and2. The trivial
eigenvalues39 and-7 must remain alongside the norms9 and169 in the
native determinant formula. The roles of7 and21 come from the actual
C7 translation group and C7 semidirect C3 transporter, not from an
elliptic rank heuristic.

This map is the direct connection to our
[matched-clock norm](EISENSTEIN_CLOCK_JOINS.md). It does NOT mean that
every integer phase is attainable by native strong-atom attachments.
For example, fixing h=39 restricts a-b to0 modulo3; the unconstrained
cross-join classification did not fix such a trivial-character sum.
For the native pair of run sectors there is an additional parity glue,
documented in TWO_RUN R3. A full integral construction must satisfy
all of its realization guards simultaneously.

## C3. Complete modular recurrence for this linear bridge

Modulo any m, the kernel of(C1) consists of equal triples(w,w,w)
with3w=0. It has size g=gcd(m,3). Thus

    |image|=m^3/g,
    image={(h,a,b) modulo m: h-a+b=0 modulo g}.           (C5)

In particular the image count at3^k, k>=1, is3^(3k-1), multiplying
by27 at each further digit. This is a different observer from the
nonlinear determinant quotient, whose count is7*27^(k-1). Equal
eventual branching factors do not identify their kernels or source data.

One extra ternary digit of(h,a,b) suffices to recover the original
triple modulo3^k, by the displayed division by3. It is sharp: adding
3^(k-1) simultaneously to the three coefficients changes the source
modulo3^k but not its image there. The compatible infinite tower is
faithful. Again, this is exact precision control rather than finiteness
of an arithmetic inverse limit.

## Inheritance, hostiles and next decisive test

The underused ancestor is the
[Session35 augmentation lattice](../session35_dimensions_dualities_20261023/recursions/AUGMENTATION-LATTICE.md).
The closest incoming mechanism is the actual sum/difference parity
glue in TWO_RUN. Their common warning is that a rational decomposition
can erase an integral extension. Here(C2) gives that extension exactly.
The source is the coefficient triple; the target is L; the map is(C1);
the preserved data are all three coefficients and cyclic clock; deleting
the guard or source sum loses realization or uniqueness respectively.

The81/27 joint-run control passes its replay: both have9,011,898 ordered
two-run covers, but their response ranks modulo3 are24 and22. Thus a
total coefficient sum does not replace joint incidence. The prime-Paley
pencil checker also replayed all556 declared cases successfully.

Next inspect a native three-run response, retaining both its S3-sector
integer gluing at2 and3 and the original complemented target relation.
Test whether a proposed phase quotient can distinguish the81/27 control
before claiming native minimality. These statements do not settle
A000568 generation, LRC(14), Collatz or Sha finiteness.

## Reproduction

`python -B cyclic_coefficient_glue.py` and the same with `-O` check
every triple in{-12,...,15}^3, the complete modular source/image sets
for1,2,3,4,7,9,11,12,13,27, the literal clock map, and both actual
Paley7 coefficient rows read from the incoming two_run.json. No scalar
count is substituted for the full source-to-image comparison.
