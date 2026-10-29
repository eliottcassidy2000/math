# Four vertices, six run markings, and the exact T5 gluing recurrence

2026-10-29. **PROVED scoped identities; FINITE-EXACT matrix certificates;
CITED classical lattice input. Same-agent alternative-algorithm verification.**
Fermat-prime completeness, Sha finiteness, LRC(14), and a carrier closed
under every native tournament operation remain OPEN.

## Inheritance and the actual advance

The closest mechanism is the exact index54 S3 chart in
[THREE_RUN_GLUE](../overnight_research_20261029/THREE_RUN_GLUE.md).
The tetrahedron-boundary description of T5's kernel was already proved in
[THREE_RUN, section T4](../overnight_rank_return_20261029/THREE_RUN.md).
It is recovered here, not claimed as a new discovery. The new deductions
are the explanation of the four spurious Z/6 factors by kernel volume,
an exact all-power cokernel sequence with its order-five clock, complete
Smith forms at every positive power, and a minimal matrix recurrence.

The least-used sidecars recovered are
[THM-3045: K4 edge isotypic binary-ternary integral clutch](../../../01-canon/theorems/THM-3045-k4-edge-isotypic-binary-ternary-integral-clutch.md),
the [August tetrahedral representation atlas](../../../07-reflections/tetrahedral-k4-representation-atlas-tournaments-xor-berggren-h1-codex-20260816.md),
and [Session39's signed linking forms](../session39_signed_phases_and_prime_cylinders_20261023/linking/README.md).
The corrected near miss is identifying rational direct sums with integral
ones, or treating unsigned edge coordinates as oriented tournament data.
The canonical hostile is the 81/27 pair, retained with the same labels.

Portfolio: anchor = exact native three-run glue; niche = the old K4 lattice
atlas; wildcard = Fermat squaring dynamics, developed in [FERMAT](FERMAT.md).
The concept board is oriented edges / named path ports / primitive kernels /
integral gluing / source clocks / CRT relative phase. Kernel volume connects
the first four; the irreducible C5 defect and the arithmetic phase example
show why the last two cannot be discarded.

Concurrent integration: incoming Session98's
[growth synthesis](../session98_tournament_growth_20261029/GROWTH_THEORY.md)
and [projective A000568 recurrence](../session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md)
were read at d9d7e33349. They remove the conserved vertex-mass direction
from a weighted counting response, while retaining exact scalar weights.
The present quotient instead compares two integral representations of a
named linear response. Both require specifying future tests, but the
former's counting equivalence does not imply our integral equivalence.
Their N=100 checks and global generating-function conclusions are inherited
from that packet, not independently replayed or claimed here. This preserves
the separation between an A000568 coefficient recurrence and an individual
metagraph attachment carrier.

## 1. Three different six-coordinate objects

A tournament on four vertices has six directed edges. After forgetting
direction its support is K4. That does not identify its edge coordinates
with the six orders of three named runs.

| Carrier and decomposition | Integral quotient | Index |
| --- | --- | --- |
| Unsigned K4 edges, rational S4 types 1+[22]+[31] | (Z/2)^3 + Z/3 | 24 |
| Oriented K4 edge cochains, cut plus cycle lattices | (Z/4)^2 | 16 |
| Regular S3 run coefficients, trivial/sign/standard Fourier chart | Z/2 + (Z/3)^3 | 54 |

The first two rows are inherited mechanisms, independently replayed here.
Their Smith forms are respectively 1,1,1,2,2,6 and 1,1,1,1,4,4; the third
has Smith form 1,1,1,3,3,6. For the second row reorienting an edge changes
its coordinate sign, preserving the cut/cycle lattice calculation. For the
first row the S4 action simply permutes unordered pairs. Even restricting
the first action to a vertex stabilizer S3 does not give its regular
representation: a transposition fixes two edges, whereas a nonidentity
element fixes no point in the regular action.

There is an exact map from a four-vertex tournament to its three-run
response M3: a spanning three-run cover has one directed edge and two
singleton runs. Thus each directed edge gives 3 choices of its run name
and 2 orders for the other runs, for 6 marked covers per edge and mass36.

More generally, for n vertices and n-1 named runs there is one two-vertex
path. A port state is an ordered (n-1)-tuple, completed by its missing
vertex. Traversing a response edge replaces one port by that missing
vertex. This identifies the n!-state response with a directed vacancy
walk on permutations. Forgetting direction gives the star-transposition
Cayley graph; projection to the missing vertex is a labelled covering of
the converse tournament with fibre size (n-1)!. Therefore

    total marked covers = n! (n-1)/2.                         (1)

The supplied tournament orientation chooses the directed lifts. The
vacancy model is checked against independent path-cover enumeration on
all four order-four classes. Their ranks are18,18,24,18 and every nonzero
Smith factor is1. In particular the unique strong class has an integral
inverse, so every positive power of that response remains unimodular.
The already completed all64-labelled census is inherited, not recast as
a fresh enumeration. No order-four class has the T5 discrepancy.

## 2. What the four directions really are

Use T5 on vertices0,...,4 with i beating i+1 and i+2 modulo5. Let M be its
60-by-60 three-run matrix, retaining every ordered endpoint triple. Let
F be the block-diagonal Fourier chart, one block per unordered triple,
and put A=F M F^(-1). The ambient sector lattice is E=Z^60 and the actual
source lattice is Lambda=F Z^60. Both matrices are integral in their own
bases, but Lambda is a proper sublattice of E.

Let P be the 10-by-10 unsigned block of A. For the oriented tetrahedron
obtained by deleting vertex i, define

    b_i = (-1)^i partial[0,...,omit i,...,4].

The five b_i sum to zero. Any four form a primitive integral basis of
K=ker P=ker P^T, of rank4. These are the inherited tetrahedral directions.
The checker verifies the actual matrices, not just the dimensions. In
the basis b_0,...,b_3 their Gram matrix is

    G = 5 I_4 - J_4,       Smith(G)=(1,5,5,5), det G=125.   (2)

The unsigned sector contains K, while the sign and standard sectors are
nonsingular. On a vector whose other sectors vanish, the chart's integral
guards reduce to u divisible by2 and by3. Consequently

    Lambda intersect K = 6K.                               (3)

Equivalently, the native kernel basis repeats every triangle coefficient
on its six run orders. Its Gram matrix is 6G, of determinant6^4*5^3.
The corresponding unglued kernel basis has Gram matrix G. This is where
the four and the six enter: four independent tetrahedron relations,
each repeated on six named orders. It is not the assertion that four
vertices intrinsically contribute four Z/6 factors.

## 3. Kernel volume gives an all-power formula

**General integral lemma.** Let B be an integral square matrix with the
same right and left rational kernel. Let K be its primitive integral
kernel, with Gram determinant g, and let D be the absolute product of
its nonzero eigenvalues. Then for every r>=1,

    |tors coker(B^r)| = D^r/g.                             (4)

Proof: put S=K^perp intersect Z^N. The common-kernel condition makes B
invertible on S_Q. The standard unimodular-lattice gluing identity gives
[Z^N:K+S]=g. B maps Z^N/K injectively into S, and
[B Z^N:B S]=[Z^N:K+S]=g. Since [S:B S]=D, the torsion order is D/g.
Replacing B by B^r preserves K and replaces D by D^r. The gluing identity
is classical; see Bacher--de la Harpe--Nagnibeda below. This proof also
shows why rank stabilization and the common left/right kernel matter.

For the literal T5 matrix the nonzero eigenvalue product is

    D = 52,488,000 = 2^6 * 3^8 * 5^3.

Thus its native torsion order tau_r obeys the exact recurrence

    tau_1 = 324,
    tau_(r+1) = 52,488,000 tau_r,
    tau_r = 2^(6r-4) 3^(8r-4) 5^(3r-3),       r>=1.        (5)

The unglued torsion order is 6^4 tau_r at every positive power. Equation
(5) is an all-r proof from the kernel-volume identity, not a fit to the
sampled powers. The appearance of prime5 after the first power is
explained by the independent Gram determinant5^3 in (2).

## 4. Stronger than an order ratio: the exact gluing sequence

The quotient E/Lambda is (Z/2)^10 + (Z/3)^30. Its induced A-operators
are supplied explicitly in certificate.json. At p=2 and3 each has a
four-dimensional kernel, rank unchanged after squaring. The map from K
onto that kernel is surjective and has kernel pK. These finite exact
identities prove the same facts for every positive power.

Apply the snake lemma to Lambda -> E -> E/Lambda and the vertical maps
A^r. The surjection from ker(A^r on E) to the quotient kernel kills the
connecting map, giving

    0 -> coker(M^r) -> coker(A^r) -> (Z/6)^4 -> 0, r>=1.   (6)

The free-part map is an isomorphism: its observer is K^T u, which is
surjective since K is primitive and every integral u occurs as the
unsigned coordinate of a native coefficient vector. Hence (6) also
restricts to an exact sequence of torsion groups.

At r=1 the native Smith form is 1^52,3^2,6^2,0^4 and the unglued form
is 1^48,3^2,6^6,0^4. At r=2 the dyadic groups are C2^4+C4^2 and C4^6.
So (6) should not be replaced by a direct sum at every r: four of its
dyadic extensions are already nonsplit at r=2.

The natural vertex rotation i->i+1 cycles b_i. Its matrix C on four
generators satisfies Phi_5(C)=0. The quotient in (6) has this same
character: Phi_5 is irreducible over F2 and F3. Thus the defect is one
four-dimensional irreducible module at each prime, realizable additively
as F16 and F81 with an order-five multiplicative clock. It is not four
clock-independent scalar counters. Field Frobenius has order4; it is a
different operation from the native vertex rotation of order5. No native
Frobenius operation is inferred from this field presentation. The number
81 here is a state count, not the automorphism order of the hostile source.

## 5. The minimal recurrence for all named linear responses

Write

    f6=(x^2-5x+5)(x^4+5x^3+15x^2+25x+25),
    f10=(x^2-3x+1)(x^4-x^3+x^2-x+1)
                     (x^4+4x^3+16x^2+24x+16),
    f20=(x^4+4x^3+9x^2+7x+2) h16,
    h16=x^16-4x^15+7x^14-19x^13+49x^12-62x^11+66x^10
        -101x^9+126x^8-93x^7-83x^6+144x^5+115x^4
        -84x^3+99x^2+189x+81.

The sector characteristic polynomials are x^4 f6, f10, f20. Their
minimal polynomials are x f6, f10, f20: the nonzero factors are squarefree
and the zero kernel is semisimple. They are pairwise coprime. Therefore

    minimal_polynomial(M) = x f6 f10 f20, degree37.          (7)

Every fixed named linear response a^T M^r b obeys its coefficient
recurrence, over the integers and after reduction modulo every m. Some
individual observers can have smaller order;37 is the minimum common
matrix recurrence. After the initial zero mode, f6 f10 f20 gives order36
for the matrix sequence beginning at M^1. All coefficients are saved.

**Operation boundary:** powers of M compose this fixed weighted linear
response. They do not count disjoint spanning covers of successively
larger tournaments. The inherited legal-complement test exposes all four
T5 kernel directions in one move, so neither K nor the defect quotient
licenses a smaller carrier closed under native complement.

## 6. Complete Smith forms at every power and every modulus

The free rank is always4. To obtain the56 nonzero invariant factors,
sort each of the following lists of valuations and multiply the
corresponding prime powers coordinate by coordinate:

    p=2:  0^50, (r-1)^4, r^2;
    p=3:  0^48, (r-1)^4, r^4.                             (8)

Here superscripts denote multiplicities, not exponentiation of r. At
p=5 there are50 zeros and six more entries. Write r=4q+s with1<=s<=4.
Add2q to each entry of the following row:

| s | six five-primary valuations before adding2q |
| --- | --- |
| 1 | 0,0,0,0,0,0 |
| 2 | 0,0,0,1,1,1 |
| 3 | 0,1,1,1,1,2 |
| 4 | 1,1,1,2,2,2 |

No other prime occurs. These are all-r results, not interpolations.

**Proof at5.** F is invertible over Z5; the sign and standard sectors
are also invertible there. Put S=ker(P^T)^perp intersect Z^10. Smith(P)
is1^6,0^4, so P Z^10=S. Direct exact integer division shows P^5/25 is
integral with the same saturated image S. Hence

    P^4 S=25 S,
    P^(r+4) Z^10=25 P^r Z^10,                  r>=1.       (9)

The four initial Smith profiles prove the displayed four-step rule.
The tempting two-step extrapolation fails at r=5 precisely because the
correct return is four steps: all six valuations are then2.

**Proof at2 and3.** Use the integral p-adic Fitting decomposition into
unit and nonunit eigenvalue parts, in both E and Lambda. In the unsigned
sector the nonunit part is precisely K, with operator zero. In each
other nonunit sector, (i) the multiplicity of x in its characteristic
polynomial modulo p, (ii) its nullity modulo p, and (iii) the valuation
of its determinant are equal. The relevant triples are (4,4,4) for the
sign sector at2, (1,1,1) for the standard sector at2, and (4,4,4) for
the standard sector at3. Consequently that nonunit block is p times an
integral automorphism. There are two standard copies.

Thus E_bad=K+W_p with A zero on K and p times an automorphism on W_p;
rank W_2=6 and rank W_3=8. The quotient bad part has exponent p and
dimension4, and K maps onto all of it, by section4. It follows that
Lambda_bad projects onto W_p. Therefore

    A^r Lambda_bad = p^r W_p.

Each positive step multiplies the image lattice by p. The original
Smith profiles at r=1 give (8); the unit parts contribute50 and48 zeros.
This proves every primary factor. Absence of further primes follows
independently from (5).

For precision h>=1, the exact native response-image logarithms are

    L2(r,h)=50h+4(h-r+1)_+ +2(h-r)_+,
    L3(r,h)=48h+4(h-r+1)_+ +4(h-r)_+,
    L5(r,h)=50h+sum_(j=1)^6 (h-e_j(r))_+,                (10)

where e_j(r) is the five-primary list. Other primes give56h. Multiply
the prime-power image counts by CRT for an arbitrary modulus. In the
unglued model replace the last two terms of L2 by6(h-r)_+, and of L3
by8(h-r)_+. Consequently

    native image / unglued image
      = 2^(4[valuation_2(m)>=r]) 3^(4[valuation_3(m)>=r]). (11)

At r=1 this recovers gcd(m,6)^4. At higher powers a coarse observer can
miss the discrepancy completely: it first appears at digit r, although
the integral quotient (6) always has order1296. This connects composition
depth with the exact precision required to see the missing coordinate.

## 7. Verification, source map, and next decisive test

Run from the repository root:

    python -B 05-knowledge/results/four_vertex_glue_fermat_20261029/verify.py
    python -B -O 05-knowledge/results/four_vertex_glue_fermat_20261029/verify.py --output <temporary-json>

Both executions reproduce identical JSON. Exact paths are: independent
subset convolution and permutation/cut enumeration for M; integral Smith
form for the original matrix and small kernels; separate local elimination
for powers1..9,13,17,25 at2,3,5; and literal unit squaring walks for every odd
modulus3..151, including prime powers and nonsquarefree cases.

The 81/27 automorphism counts and67 common labelled cyclic triples are
recomputed. Their saved, label-checked two-run matrices still have ranks
24/22 modulo3. Those matrices are reused, not claimed freshly generated.

The connection maps and losses are explicit: unsigned K4 support forgets
tournament directions; the cut/cycle map keeps them with index16 gluing;
F keeps all named run coefficients only with its four congruence guards;
the quotient in (6) records integral error with a C5 clock; the arithmetic
CRT map in FERMAT keeps relative cycle phase. These are different objects
sharing a lattice/phase preservation mechanism, not a single unnamed six.

Next decisive test: transport the exact sequence through an actual
path-complement pair with the same source coordinates. The inherited
rank60 stacked-observer result predicts that a common four-dimensional
zero-mode explanation must break. Determine the replacement gluing
sequence for a mixed word M_U M_T, rather than extrapolating the fixed-T
power theorem to a changing native source.

Primary lattice source: [Bacher, de la Harpe and Nagnibeda (1997), The
lattice of integral flows and the lattice of integral cuts on a finite
graph](https://www.numdam.org/item/BSMF_1997__125_2_167_0/). The general
unimodular gluing mechanism is classical; the identified T5 matrices,
their clocks, and the deductions (5)--(11) are the present scoped work.
