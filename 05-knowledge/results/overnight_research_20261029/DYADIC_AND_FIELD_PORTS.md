# Complete modular carriers: dyadic branching and finite-field recovery

2026-10-29, second hourly continuation. **PROVED in the stated query language;
FINITE-EXACT supporting checks; same-agent audit, no external priority claim.**

This closes the dyadic obligation left in
[Eisenstein clock joins, E5](EISENSTEIN_CLOCK_JOINS.md). It does not identify
integer residue rings with finite fields of the same size.

## D1. Source, queries, and exactly what is observed

Keep the matched order-three clock and the labelled cross chart

    B(a,b)=[[a,b],[-a-b,a]],  Q(a,b)=a^2+ab+b^2,  z=de.

The source state modulo m is any (z,a,b) in (Z/m)^3. A named integer
attachment (u,v) produces the determinant observation

    R_(z,a,b)(u,v)=[z-Q(a+u,b+v)]^2 modulo m.              (D1)

Translations, global sign reversal, and the clock-centralizing unit
re-markings from E4 are allowed. They reindex the complete response
function, so equality of those functions is a closed equivalence.
K(m) counts its classes. These are full functions, not merely equal
initial determinants. Observations do NOT include exact singularity,
integral cokernels, or named classes in those cokernels. Formal modular
updates through singular residues are allowed. The clock generator,
port names and query coordinates are fixed and known throughout.

The 5-by-5 grid of integer queries determines the complete function:
fifth forward differences in either argument vanish. The resulting
integer recurrence propagates forward and backward without division.
It is not valid to infer coefficient equality over every residue ring.

## D2. Exact dyadic equivalence and the complete recurrence

For k>=4 put m=2^k. Two states have equal complete responses iff

    z'=z modulo m/2;
    if z is even: a'=a and b'=b modulo m/2;
    if z is odd:  a'=a and b'=b modulo m/4.                (D2)

The fibers have size8 in the even-source stratum and32 in the odd-source
stratum. At k=3 the phase conditions reverse: even z retains a,b modulo2,
while odd z retains them modulo4; z is still retained modulo4.
At k=1 the response is faithful. At k=2 it retains exactly the three
parities. Therefore

    K(1)=1, K(2)=8, K(4)=8;
    K(2^k)=5*2^(3k-6), k>=3.                             (D3)

In particular the sequence at moduli1,2,4,8,16,32,64,128,256 is
1,8,8,40,320,2560,20480,163840,1310720.

This has a full refinement-tree law, not just a count recurrence:

| Precision step | Children per old response class |
| --- | --- |
| 1 to 2 | 8 |
| 2 to 4 | 1 |
| 4 to 8 | 2 for even z, 8 for odd z |
| 8 to 16 | 32 for even z, 2 for odd z |
| 2^k to 2^(k+1), k>=4 | 8 for every class |

The numbers in the first column are moduli. In particular the factor8
in K(16)/K(8) is an average over very unequal branch sizes. Only the
subsequent steps are uniformly eight-branching. The keys in(D2) give
the actual parent maps and closed attachment transitions.

For fixed source z=9 the phase counts at k=0,1,2,3 are1,4,4,16;
for k>=4 they are2^(2k-4). Thus even this source has a plateau at the
8-to16 step, followed by uniform factor4 refinement.

### Proof

Put L1=2a+b, L2=a+2b and t=z-Q(a,b). The difference of two response
polynomials has degree at most3. Its third differences along u,u,u
and u,u,v are respectively12 delta(L1) and4 delta(L1+L2).
For k>=3 these force L1,L2, hence a,b, to agree modulo q=2^(k-2);
the linear coordinate change has determinant3, a dyadic unit.

The second u-difference is

    2 delta(L1^2-2t)+12 delta(L1).

Writing delta(L1)=q eta, its terms other than -4 delta(t) are divisible
by4q: they are4q L1 eta+2q^2 eta^2+12q eta, with q even.
It follows that t and then z also agree modulo q.

Write the remaining increments as q(epsilon_z,epsilon_a,epsilon_b),
with each epsilon modulo4. For r=a+u, s=b+v put

    W=z-Q(r,s),
    H=epsilon_z-(2r+s)epsilon_a-(r+2s)epsilon_b,
    C=Q(epsilon_a,epsilon_b).

The new inner expression is W+qH-q^2 C. If k>=4, its squared difference
modulo4q is2qWH. Equality is therefore exactly WH=0 over F2.
Here Q(r,s)=r+s+rs vanishes only at(0,0), and
H=epsilon_z-s epsilon_a-r epsilon_b.

For even z, W=1 at all three nonzero pairs. These three equations force
all epsilons even. For odd z, W=1 only at(0,0), forcing just epsilon_z
even. This proves both directions of(D2).

At k=3, q=2 and the q^2 H^2 term survives: the difference is
4(WH+H^2). Equality now says H(W+H)=0 over F2, equivalently H=0
where W=0. This exchanges the two strata. It explains the exceptional
branching; it is not an extrapolation from a numerical fit.

At k=1 the minority position of the binary function identifies a,b
and its value identifies z. Squaring modulo4 only sees parity, giving
the k=2 statement. Counting the displayed fibers proves(D3). Reducing
the coordinate keys proves every transition count in the table.

## D3. Every integer modulus, and precision rather than finiteness

Write m=2^alpha 3^beta c with gcd(c,6)=1. Combining(D3) with the
proved ternary classification and CRT in E5 gives the complete formula

    K(m)=kappa2(alpha) kappa3(beta) c^3,
    kappa2(0)=1, kappa2(1)=kappa2(2)=8,
    kappa2(alpha)=5*2^(3alpha-6) for alpha>=3,
    kappa3(0)=1, kappa3(beta)=7*27^(beta-1) for beta>=1.    (D4)

Both states and named query arguments split by CRT; the product rule
does not presume independent ports where none were supplied.
For alpha>=3,beta>=1 the ratio K(m)/m^3 is35/1728.
Examples are K(11)=1331, K(12)=56, K(13)=2197, K(76)=54872,
K(228)=384104 and K(768)=9175040. The earlier5,10,25,50,100 sequence
has response counts125,1000,15625,125000,125000. The last equality is
the general plateau K(2c)=K(4c) for odd c, not an aliquot identity.

All three source coordinates modulo2^k are determined by responses
modulo2^(k+2). Two extra bits are a sharp uniform bound for k>=2:
at k=2 the even-source states(0,0,0),(0,2,0) coincide modulo8; for
k>=3 the odd-source states(1,0,0),(1,2^(k-1),0) coincide modulo2^(k+1)
but differ in their k-bit phase. Consequently the complete compatible
2-adic tower is faithful, despite every sufficiently deep finite level
having nontrivial fibers. The ternary classification likewise needs at
most one extra digit to recover all source digits. None of this bounds
an arithmetic Selmer tower or implies Sha finiteness.

## D4. Over F64 and F256, three queries suffice and are necessary

Replace the base ring by any finite field F_q of characteristic2, and
keep exactly the polynomial(D1), all q^3 states and all field queries.
Squaring is bijective: its inverse is x -> x^(q/2). Let w00,w10,w01
be the unique square roots of the observations at(0,0),(1,0),(0,1).
Then

    b=1+w00+w10,
    a=1+w00+w01,
    z=w00+Q(a,b).                                        (D5)

Indeed Q(a+1,b)+Q(a,b)=1+b and Q(a,b+1)+Q(a,b)=1+a.
Thus all q^3 states are separated. Three q-valued queries are also
necessary in the worst case by cardinality, even if chosen adaptively.
This includes zero and singular responses; restricting the source
universe would be a different problem.

| Base ring/field | Complete response classes |
| --- | ---: |
| Z/64 | 20,480 |
| F64 | 262,144 |
| Z/256 | 1,310,720 |
| F256 | 16,777,216 |

The equality of cardinalities of the base sets does not preserve their
multiplication, squaring, nilpotents, or polarization. Over the fields,
Frobenius removes the outer square and the nondegenerate polar form
recovers phase. Over the dyadic rings, the surviving quadratic error
at modulus8 creates the stratum reversal in D2. The field result is
an elementary polarization calculation, not a new finite-field theorem.

The irreducibility descriptor of the rational C3 clock must also change
on base extension. In F64 and F256, X^2+X+1 has two roots; squaring
exchanges them. Thus the clock splits there, and Frobenius retains a
nontrivial two-root permutation. The checker explicitly verifies these
roots and their exchange. The decoder does not require Q to be anisotropic.

## Inheritance, connections and limits

The closest mechanism is E5's ternary finite-difference proof. Its
corrected near miss is division by a Vandermonde determinant at a bad
prime: values determine the function but not all coefficients.
Canonical hostiles here are the opposite parity strata at8 and16,
and F4 versus Z/4 (64 versus8 classes).

The least-used recovered sidecars are
[THM-4435, universal Hermite precision](../../../01-canon/theorems/THM-4435-four-node-metric-blindness-and-universal-hermite-precision.md)
and [THM-4443, higher-jet unit boundary](../../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md).
The former computes sharp observation precision via an inverse matrix;
the latter shows that an LCA/distance tree can miss unit residues. Here
the finite-difference inverse has its own sharp extra-bit bound, and the
source parity controls the fiber. These are different operators; no
Smith-matrix or tree-isomorphism identification is asserted.

The recovered [11/12/13 sampled-unit clock](../weighted_crt_clocks_20261028/CLOCKS.md)
has losses at11 and13 because its observation is a sampled trace.
Our full cross-attachment observations are faithful at both primes.
The source and allowed tests, not the repeated prime alone, determine
the quotient. The Zagier/Fricke normalization in that earlier note is
an inherited source claim, not re-derived here.

The live board is: joint labelled ports / clock generator / finite
response quotient / inverse precision tower / arithmetic realization.
The new result closes the fourth item for this formal quadratic family.
It leaves the marked arithmetic realization open. It also distinguishes
an eventually regular refinement tree from a critical ternary return
process; a common recurrence metaphor does not transport a Collatz proof.

The81/27 tournament control remains separate: identical labelled cyclic
triples do not retain its stabilizers or ordered attachments. Incoming
[endpoint precision](../overnight_rank_return_20261029/ATOM_PRECISION.md)
separates that pair at modulus9 in the single-run channel. The newer
[two-run module](../overnight_rank_return_20261029/TWO_RUN.md), incoming
commit4ff64c524c, strengthens the observer and separates it already modulo3
by ranks24 and22. These statements have different observation languages.
Nothing in(D4) recreates A000568 or replaces the native joint-run carrier.

There is also a literal norm map in that incoming Paley7 calculation.
The two C3 quotient coefficient triples are(15,12,12) and(-11,2,2).
For ANY triple(c0,c1,c2), put a=c0-c2,b=c2-c1. Then

    c0^2+c1^2+c2^2-c0*c1-c0*c2-c1*c2 = Q(a,b).

This is the norm of c0+c1*omega+c2*omega^2 in the nontrivial C3
constituent, using eta=-omega as the Eisenstein basis. The two triples
give(a,b)=(3,0),(-13,0), hence norms9 and169. This is a specified
source/quotient/basis map connecting the user's15,11,12,13 motifs.
It loses the trivial-character coordinate c0+c1+c2; retain that sum
and the integral realization guard. The incoming sum/difference sectors
are glued modulo2, so their rational spectral factors alone do not
determine the integer response. No isomorphism of their operation sets
with our arbitrary cross-attachment language is inferred.
The [cyclic coefficient bridge](CYCLIC_COEFFICIENT_GLUE.md) now spells
out the exact index-three guard h-a+b=0 modulo3, its source inverse,
its dual C3 clock and its all-modulus recurrence.

## Reproduction

`python -B dyadic_field_ports.py` and the same with `-O` exhaust every
state modulo1,2,4,8,16,32,64. Every residue query is compared through8;
the complete 5-by-5 function certificate is used thereafter. The checker
also tests selected entire lift fibers through2^18 and the exact
parent-map branching, rather than just aggregate counts.

All64 states over F4, all262,144 over F64 and all16,777,216 over F256
are decoded from three directly evaluated outputs. The polynomials
x^2+x+1, x^6+x+1 and x^8+x^4+x^3+x+1 pass the finite-field
irreducibility criterion: x^(2^r)=x and gcd(x^(2^(r/p))-x,f)=1 for
every prime p|r. This works because every irreducible factor degree
must divide r, and a proper divisor divides some r/p. A separate
carryless-product/long-division path verifies every multiplication pair.
Formula-derived CRT values are labelled as such, not claimed censuses.
