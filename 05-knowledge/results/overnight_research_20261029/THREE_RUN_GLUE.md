# Native three-run responses: the exact S3 lattice and a Sha obstruction

2026-10-29, third hourly continuation. **PROVED integral identities;
FINITE-EXACT native matrices and primary certificates; CITED arithmetic
boundary. Same-agent alternative-algorithm audit.** No general Sha,
Collatz, LRC(14), A000568 or metagraph-minimality conclusion is claimed.

**Concurrent inheritance:** incoming commit e9e53cf3ed independently
obtained the same S3 lattice and native2/3-primary layers in
[THREE_RUN](../overnight_rank_return_20261029/THREE_RUN.md), and additionally
proved a sharp T5 native-complement reconstruction. Those shared results
are corroborating derivations, not separate discoveries. Our additional
outputs here are the complete Paley7 Smith form, the order4 minimality
census, singleton-run hostile decks and the pairing/clock obstruction.

## R1. The source and its six run-order coordinates

For a labelled tournament T on n vertices, let M3 have ordered triples
of distinct vertices as row and column indices. Its(s,t) entry counts
ordered spanning covers by three vertex-disjoint directed paths, whose
starts and ends are s and t in the given run order. Singleton paths are
allowed. As in the inherited square-zero exclusion construction,

    M3[s,t]=[x_V] P_(s1,t1) P_(s2,t2) P_(s3,t3),
    P=X(I-AX)^(-1),  x_v^2=0.                            (R1)

The three factors retain their names; multiplication enforces disjoint
support and the coefficient extraction enforces spanning. Simultaneous
permutation of the runs leaves M3 invariant. Each unordered three-set
of ports therefore supplies one regular S3 orbit of six coordinates.
The S3 action permutes run names. It is NOT an even-order automorphism
group of the tournament, whose actual vertex symmetries are kept separately.

The closest proved mechanism is [TWO_RUN](../overnight_rank_return_20261029/TWO_RUN.md),
especially its parity fiber product. The anchor is native integral
joint response; the niche is a source-preserving group-ring lattice;
the wildcard is whether its actual cokernel can even have a Sha pairing.
The live board is exclusion / named ports / symmetry sectors / integral
precision / alternating pairing / arithmetic realization.

## R2. A complete integral Fourier chart, with its clock retained

Order S3 lexicographically as

    (012),(021),(102),(120),(201),(210).

Let rho(g) be the action of g on the augmentation basis(e0-e2,e1-e2):
rho(g)_ij=[g(j)=i]-[g(2)=i], for i,j=0,1. For six integer coefficients
c_g define

    u=sum c_g, v=sum sign(g)c_g, W=sum c_g rho(g).

With W=[[a,b],[c,d]], the coordinate transform is the following matrix F:

    [1  1  1  1  1  1]
    [1 -1 -1  1  1 -1]
    [1  1  0 -1  0 -1]
    [0  0  1 -1  1 -1]
    [0 -1  1  1 -1  0]
    [1 -1  0  0 -1  1].                                 (R2)

Its determinant is54 and its Smith factors are1,1,1,3,3,6. Its image
is EXACTLY the following sublattice of Z^2 plus Mat2(Z):

    u=v modulo2;
    a+b=u, c+d=u, a-c=v modulo3.                          (R3)

For a constructive inverse, H=6F^(-1) is

    [1  1  2  0  0  2]
    [1 -1  2 -2  0 -2]
    [1 -1  0  2  2  0]
    [1  1  0 -2  2 -2]
    [1  1 -2  2 -2  0]
    [1 -1 -2  0 -2  2].                                 (R4)

Multiplication verifies HF=6I. Reducing its six divisibility conditions
modulo2 gives the single parity condition; reducing modulo3 gives the
three conditions in(R3). Hence(R3) is both necessary and sufficient,
and c=H(u,v,a,b,c,d)^T/6 recovers every original coefficient. The
reused letter c in the matrix entry is a scalar; the inverse yields
the six-component coefficient vector.

The conditions have a representation-theoretic explanation: modulo2,
the trivial and sign characters coincide. Modulo3 the standard module
has the invariant line spanned by(1,1), with trivial character on that
line and sign character on the quotient. The three ternary conditions
retain precisely this extension. The cycle(120) acts by
rho(120)=[[-1,-1],[1,0]], the same explicit Eisenstein phase clock in
[CYCLIC_COEFFICIENT_GLUE](CYCLIC_COEFFICIENT_GLUE.md). Reflections have
determinant-1 and retain the run-permutation sign.

### Closure and what it does not assert

The map c -> (u,v,W) is the integral group-ring representation. Thus
its image is a ring with componentwise scalar multiplication and matrix
multiplication on W. For any numbers of regular port orbits, matrices
whose entries satisfy(R3) classify ALL S3-equivariant integer linear
maps between their coefficient spaces. The inverse(R4) recovers each
six-entry convolution kernel. This is an exact closed algebra for those
linear operations, with the source lattice and run clocks retained.

For a native M3 block indexed by source/target three-sets i,j, put
k_ij(g)=M3[(i,identity),(j,g)]. Its standard block is

    L_ij=sum_g k_ij(g) rho(g^(-1)).                       (R5)

The two scalar blocks are the unsigned and signed sums. The inverse
inside rho is essential: the original block has entries k_ij(hg^(-1)).
In Fourier coordinates it acts by left multiplication L_ij on W_j,
so the standard sector occurs twice, once for each matrix column.
Every block satisfies(R3), and the checker verifies the full6-by6
conjugation identity, not only the sector determinants.

This closure is for equivariant linear maps. Multiplying endpoint
matrices without the square-zero support algebra can revisit vertices.
It also omits the full complemented target relation needed by native
metagraph composition. Those data are not supplied by a Fourier chart.

### All-modulus observer recurrence

From the Smith factors of F,

    |image(F modulo m)|=m^6/[gcd(m,3)^2 gcd(m,6)].         (R6)

For b named regular orbits take the b-th power. At2^k,k>=1 the single-
orbit exponent is6k-1; at3^k it is6k-3. Subsequent digits multiply
the count by64 and729 respectively. One extra digit at2 or3 is exactly
the uniform precision needed to reconstruct original coefficients;
all other primes have no loss. These counts concern this linear chart,
not the earlier nonlinear K(m) determinant quotient.

## R3. The first native Smith discrepancy is already at five vertices

For T5, with edges i->i+1,i+2 modulo5, M3 has order60, total mass210
and exact Smith form

    1^[52], 3^[2], 6^[2], 0^[4].                        (R7)

The separate sector Smith forms are

    plus:     1^[6], 0^[4];
    minus:    1^[6], 2^[4];
    standard: 1^[16], 3^[3], 6,   occurring twice.

Ignoring(R3) therefore yields1^[48],3^[2],6^[6],0^[4]. Its cokernel
has four extra Z/6 factors. The rational rank56 and free rank4 agree,
and both determinants vanish, so neither statistic detects the error.

For EVERY modulus m the correct response-image count and discrepancy are

    |image(M3 modulo m)|=m^56/[gcd(m,3)^2 gcd(m,6)^2],
    correct image size / unglued image size = gcd(m,6)^4. (R8)

Thus the missing factors are16 at powers of2,81 at powers of3 and1296
when6 divides m. At all positive precisions the dyadic image exponent
is56k-2 and the ternary exponent56k-4.

This five-vertex boundary is FINITE-EXACT minimal for the specified
root-lattice sector split: at order3 all three paths are singletons,
so M3=I6 regardless of orientation. A complete64-labelled-tournament
census at order4 reaches all four isomorphism classes and finds no
Smith discrepancy. Their ranks are18,18,24,18, with all nonzero Smith
factors1. Orders below3 have no ordered three-run port space.

## R4. Paley7: complete primary data and a small exact verifier

For edges i->i+s, s in{1,2,4} modulo7, M3 has order210, full rational
rank210 and total mass6678. Its complete exceptional-prime multiplicities
(the entry at position e counts factors with valuation e) are

| prime | actual | unglued |
| --- | --- | --- |
| 2 | (165,31,2,3,9) | (160,35,6,3,6) |
| 3 | (201,2,6,0,1) | (197,9,3,1,0) |
| 5 | (198,12) | (198,12) |

The full actual Smith form, using x^[r] for multiplicity r, is

    1^[165], 2^[31], 4^[2], 40^[3], 240^[2], 720,
    46419120^[2], 789125040,
    24488058009151440^[2], 220392522082362960.             (R9)

Its dyadic response exponent is
165k+31(k-1)_++2(k-2)_++3(k-3)_++9(k-4)_+;
its ternary exponent is
201k+2(k-1)_++6(k-2)_++(k-4)_+.
Both become linear after depth4, with slopes210 and offsets80 and18.
These give exact all-depth recurrences for this specified native observer.

The reference sector Smith forms were first obtained by direct integer
Smith reduction. A much cheaper replay computes exact Bareiss sector
determinants, divides them completely over the independently prime-checked
bank{2,3,5,11,17,79,5861,392809}, and certifies every local valuation by
unit-pivot elimination. The remaining determinant cofactor must equal1;
the bank is a certificate, not an assumed search cutoff. The resulting
sector forms agree with the direct reference. Outside2 and3, F is an
integral local isomorphism. At2 and3 the complete actual210-dimensional
layers replace the unglued ones. This reconstructs every factor in(R9).

The initial general Smith runs took minutes; the final declared verifier
takes a few seconds on this workspace. These are observed runtimes, not
an asymptotic improvement or a benchmark claim. `--direct-smith` retains
the slower independent path. No enormous-prime-factor assumption is hidden.

## R5. The81/27 control with a designated singleton run

Fix the third run to be the singleton vertex v. Then its corresponding
submatrix is exactly the two-run matrix of T-v:

    M3[(a,b,v),(c,d,v)]=M2_(T-v)[(a,b),(c,d)].             (R10)

No approximation is involved: deleting that singleton gives a bijection
of the counted path covers. The complete third-run source and target
are retained. For the original twelve-vertex pair, the unordered decks
of these slice ranks modulo3 are

    order81 control: 22^[3],28^[9];
    order27 control: 22^[3],28^[3],32^[3],34^[3].          (R11)

The parity decks agree:56^[3],58^[9]. The scalar mass decks also agree:
1268622^[3],1243890^[9]. Modulo5 the decks are106^[12] versus
105,106^[11]. The original unsliced two-run matrices already separate
modulo5 by ranks132 and130, so this is not a claim of first separation
there. Their ternary ranks24/22 remain the preceding stronger observer.

The script enumerates the actual vertex automorphisms and uses one
representative per orbit, with exact multiplicities. It checks the same
67 labelled cyclic triples and all inherited source rows. Selected
entries receive an independent walk inclusion-exclusion verification.
It does NOT claim the full1320-by1320 three-run Smith form was computed.

## R6. A genuine obstruction to a direct Sha identification

The order of the finite Paley7 group in(R9) factors as

    2^80 3^18 5^12 11^6 17^4 79^3 5861^6 392809^3.       (R12)

Its squarefree part is79*392809=31031911. Thus this native cokernel
is neither of square order nor twice-square order.

**CITED arithmetic boundary.** Poonen--Stoll,
[The Cassels--Tate pairing on polarized abelian varieties, Corollary9](https://math.mit.edu/~poonen/papers/sha.pdf),
states that a finite Sha for a principally polarized abelian variety
is abstractly T times T, or Z/2 times T times T. Consequently(R12)
excludes identifying this entire finite native response cokernel with
such a Sha group. The source was checked2026-10-29. This is an obstruction
to the proposed direct bridge, not evidence against Sha finiteness.

There is a sharper elementary group test than square order. For each
prime p and exponent e, count the cyclic summands Z/p^e. A perfect
alternating pairing exists iff every such multiplicity is even. To
prove necessity, choose an element of largest order and its primitive
pairing partner, split off their nondegenerate hyperbolic plane, and
iterate. Pairing equal cyclic summands proves sufficiency.

Therefore the smallest DIRECT-SUM addition making an abstract finite
abelian group alternating adds one Z/p^e for each odd multiplicity.
Applied to(R9), the required additional cyclic orders are exactly

    2, 8, 16, 81, 79, 392809.                            (R13)

Allowing the exceptional single Z/2 in the cited Sha group shape removes
the first2 from this minimum abstract repair. This does NOT assert that
the remaining completion is realizable arithmetically or extends a given
clock. For comparison, even the square-order cyclic group Z/25 fails
the alternating test and needs another Z/25. The ternary81 in(R13)
is required despite the original group's even total3-adic exponent18:
the multiplicities at each exponent, not just their weighted sum, matter.

### The canonical clock-preserving doubling and its limitation

For ANY nonsingular integer matrix M there is a canonical alternating
presentation

    H(M)=[[0,M],[-M^T,0]].                              (R14)

Its Smith factors duplicate those of M; its discriminant group is the
hyperbolic pairing on coker(M) and its dual, with squared order.
If R_out M=M R_in are specified unimodular source/target clocks, then
diag(R_out^(-T),R_in) preserves H(M). Literal multiplication proves this
identity. The induced cokernel clock is its contragredient
diag(R_out,R_in^(-T)), preserving the source and dual names.
For native permutation transporters this keeps the named
ports and their existing clock action. All496 nonsingular2-by2 matrices
in the declared entry box{-2,...,2} verify the duplicated Smith law.

This makes the doubling/squaring connection constructive. It is also
available for every finite presentation, so it supplies no arithmetic
Selmer realization or finiteness theorem by itself. The smaller abstract
repair(R13) lacks a supplied clock extension; finding its least compatible
replacement was the next question. The subsequent
[clock-minimum calculation](CLOCK_MINIMAL_PAIRING.md) now resolves all
six elementary odd-primary translation modules. At11,79,5861,392809,
every invariant bilinear form vanishes and hyperbolic doubling is optimal
even with the full source automorphism action. The higher-exponent2/3
parts remain open.

## Inheritance and next decisive tests

The corrected near miss is replacing a rational sector decomposition
by an integer direct sum. The least-used ancestor recovered here is
[THM-2785, minuscule carry Fourier/common-atom boundary](../../../01-canon/theorems/THM-2785-a12-minuscule-carry-fourier-character-and-common-atom-boundary.md):
matching Fourier information does not supply the fully marked physical
object. Our exact source lattice repairs an analogous algebraic loss;
it does not transfer that LRC mechanism or its source clocks.
The existing META-PATTERNS rule to integralize after descent applies.

The incoming complement result already resolves one concrete native test:
T5's four hidden directions are tetrahedron boundaries, and every one of
its15 Hamilton-path complements exposes all four in the named weight
language. The fully vertex-invariant source instead has dimension12
and is already injectively observed. Preserve these distinct scopes.
The [primary-clock continuation](PRIMARY_CLOCK_COMPLETION.md) now
classifies the higher-exponent2/3 repair for C7 and the exact difference
between adjoining roots and adding independent summands. Next extend
the remaining2/3 affine/run clocks and the native complemented-relation
test beyond this source.
A sector-image count alone cannot replace that relation or recreate
A000568.

## Reproduction

Run `python -X utf8 -B three_run_glue.py` and the same with `-O`.
The receipt records all labelled universes, matrices by canonical hash,
full small Smith forms, sector references, complete primary layers and
hostile slice orbits. `--direct-smith` additionally repeats the slower
general integer Smith calculation. Every native matrix through7 is
constructed both by disjoint subset-path contraction and by complete
permutations with two cuts. Exceptions remain active under optimization.
