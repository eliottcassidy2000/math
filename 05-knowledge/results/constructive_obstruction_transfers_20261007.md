# Constructive obstruction transfers: clocks, algebraic witnesses, and inverse support

**Status:** PROVED elementary transfer rules and FINITE-EXACT implementations;
the named OpenAI paper statements are CITED and accepted as premises for
this session, without auditing their proofs. Existing classical mechanisms
are identified below rather than claimed as new. Universal positive-integer
Collatz coverage, a displayed six-chromatic unit-distance graph, and JC(2)
remain OPEN in this work.

## 1. Inheritance and the useful change of representation

The anchor is source-specific Collatz coverage. The niche is extracting a
finite algebraic colouring obstruction. The wildcard is exact compression
and polynomial termination in the planar Jacobian program. The live board
is **source / guard / time pair / finite witness / realization / degree or
height**. These are compared by actual maps and preserved predicates, not
by matching numerals or theorem names.

The closest Collatz mechanism is the [finite-bank compiler](collatz_finite_bank_certificate_20261007.md),
integrated with [THM-4581, Haar coalescence](../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md).
Its hostile pair is (1,2): both positive sources reach ROOT but never merge
at equal Terras time. The repaired near miss is therefore an all-integer
equal-clock conclusion. Its neglected sidecar is the pair of actual clocks.

For geometry the inherited mechanism is compactness plus real-closed-field
transfer, already in [THM-4558, field-plane colourings](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md)
and Madore section 5.4. Its hostile is a complete search inside a field known
to be at most five-colourable. For JC the [audited planar board](planar_jc_long_20260906_board.md)
already separates compatible finite response rows from polynomial
termination. An arbitrarily deep identity jet can hide nontrivial global
composition; the neglected sidecar is a degree-support bound.

The useful common operation is to compile a proposed transfer into a finite
object with an independent checker. This does not settle the existence of
such an object for every input. It does make the next missing assertion
precise, and allows successful certificates to be reused without losing
their source, domain, or terminal condition.

## 2. Clock pairs make asynchronous Collatz certificates composable

For any total map F, a receipt (x,y;a,b) means the checked equality
F^a(x)=F^b(y). Across a matching middle source, composition is

\[
(a,b)*(c,d)=(a+\max(c-b,0),\ d+\max(b-c,0)).
\]

This is the classical bicyclic monoid: (a,b) acts as translation by a-b
only on the tail t>=b. Storing only a-b loses the entry barrier. For example,
the valid Terras receipt (5,4;3,3) has zero charge but does not say 5=4.
The [clock proof and checker](collatz_clock_holonomy_20261007.md) preserve
this barrier and verify source seams before composition.

A closed receipt (x,x;a,b), a!=b, proves that F^min(a,b)(x) is periodic
with period dividing |a-b|. For finitely many such loops at one source,
let M be the largest entry time and g the gcd of their positive charges.
Then F^g(F^M(x))=F^M(x), by Euclidean subtraction of periods.

For the positive Terras map T, the only solutions of T^2(z)=z are 1 and 2.
Thus an actual loop bank with charge gcd 2 supplies a ROOT certificate.
This turns some closed proof dependencies into a sound finite inference.
Zero-charge loops do not qualify; universal existence of the required
nonzero loops is unproved. The demonstration loops are deliberately built
from already known finite ROOT paths, not presented as new coverage.

## 3. The numbers 223, 233, and 322 give an exact family, with its limits

The recovered arithmetic roles are distinct: 223 and 233 are ordered
odd-word carries, with binary orders 37 and 29 at those primes; 233 is F13,
while 322=L12=89+233. The latter is even and cannot be the unreduced carry
of a nonempty odd Collatz word. The prime-three Kervaire paper also has
stem 322, but no predicate-preserving map between that homotopy stem and
these integers has been supplied.

Their actual common future is

\[
T^7(223)=T^{15}(233)=T^{25}(322)=425.
\]

The triangle formed by these three receipts has zero total charge. It
proves a common future and is not, alone, the periodicity certificate above.
Retaining the full words yields the exact simultaneous family

\[
m=223+62208t,\qquad n=233+65536t,\qquad H=425+118098t,
\quad t\ge0.
\]

Both sources reach H with the same respective valuation words as their
seeds. Moreover m=(243n+469)/256<n. A supplied authenticated ROOT word
for m can therefore be spliced into one for n. Every n in this particular
family already has ordinary first-step descent, so this is an exact
transport interface, not added generic stopping coverage.
The [prime/partition package](collatz_prime_partition_223_20261007.md)
proves the general simultaneous-fibre formula and the full 322 lift.

It also obtains a conditional transfer from the accepted prime-predecessor
Poisson–Dirichlet theorem. If p-1=2^a t with t odd and n0=p-2, then during
the initial a-1 one-valuations,

\[
n_j+1=3^j2^{a-j}t.
\]

Deleting all factors 2 and 3 leaves the same integer at every point of
this tube, even at an adaptively selected index inside it. Its ranked
logarithmic prime factors, with the original denominator log(p-1), inherit
the accepted PD(1) limit over sampled primes. The invariant ends at the
reset, and neither its value nor its distribution implies downward drift.

The Partition Principle paper sharpens the selector requirement. An
injection back into a set of certificates need not choose a certificate
for the correct source. The required section satisfies pi(s(n))=n.
Our explicit fibre has a computable parameter inverse and checks exactly
this condition. Choosing the least certificate code is effective when a
certificate exists; uniqueness does not prove universal termination.

## 4. The plane theorem does give a constructive finite search

Accepting unrestricted non-five-colourability, compactness gives a finite
non-five-colourable unit-distance graph. A vertex-critical subgraph has
chromatic number exactly six. Real-closed-field transfer realizes its exact
edge and nonedge pattern with algebraic coordinates. Enumerating finite
graphs and deciding their real-algebraic realization consequently gives a
search with a halting proof under that premise. This is the classical
Madore consequence, not a newly discovered implication of ZFC.

What is missing is a practical bound on vertices, algebraic degree, or
description height. Searching one preselected field is insufficient.
The [finite-algebraic package](finite_algebraic_witness_20261007.md)
implements exact chosen-root geometry, a colouring-refutation checker,
and a Moser-spindle positive control. It does not display a six-chromatic
graph. Root isolation matters: a nonzero polynomial remainder can vanish
at the selected root of a reducible defining polynomial.

A separate useful reduction emerged from the triangle-removal paper.
For every geometric unit-distance graph in the plane, form a graph whose
vertices are its equilateral triangles, adjacent when they share an edge.
This conflict graph is bipartite. For triangle centroid c, the nonzero
complex number product(z_i-c) changes sign across a shared edge, providing
the bipartition on each component. Hence maximum edge-disjoint triangle
packing reduces to bipartite matching, with an exact matching/cover
certificate. This holds even when edges of the original graph cross.

The accepted random-removal law for complete graphs does not apply to this
geometric optimization problem. Nor does deleting triangles preserve
chromatic number. The Mahler paper supplies an affine volume invariant,
but equal convex hulls can contain unit-distance graphs of different
chromatic numbers; the incidence information must remain in the packet.

## 5. Superstring compression and a finite Jacobian obstruction block

Actual invertible operation words can share one superstring while retaining
their occurrence intervals and inverse dictionaries. If P_j is the map of
a prefix of that tape, the operation on interval [l,r) is P_r P_l^{-1}.
The accepted superstring theorem optimizes explicit token length; interval
addresses, token payload, polynomial expansion, and evaluation complexity
are separate costs. Partial Collatz receipts still require their clocks
and domains. The finite demonstration stores twelve tokens on a six-token
tape and independently verifies every chart cocycle.

The group Kervaire theorem preserves coefficients in an abstract overgroup;
it does not realize its new element as a polynomial automorphism of the
same plane. The prime-three theorem retains an actual order-three lift,
an extra obligation beyond associated-graded survival. In the JC compiler
the corresponding elementary repair is relative: corrections must lie
in the kernel of the data already fixed, as well as kill the new defect.

There is also a finite termination target. For F in k[x,y]^2 of degree d,
F(0)=0 and invertible linear part, let G be its formal inverse. The classical
planar inverse-degree bound gives

\[
F\text{ is a polynomial automorphism}
\iff [G]_{d+1}=[G]_{d+2}=\cdots=[G]_{d^2}=0.
\]

Indeed, if this block vanishes, H=jet_d G makes F(H)-id a polynomial of
degree at most d^2 with every possible coefficient zero. This is a standard
degree-bound consequence, not a proof of JC(2). The
[Kervaire/superstring/JC package](kervaire_superstring_jacobian_20261007.md)
keeps the whole input polynomial and checks both inverse compositions.
Its quadratic control treats all normalized quadratic Keller coefficients,
not just numeric samples. The known quadratic case is therefore a positive
control for compiling universal coefficient obligations into exact identities.

More precisely, fix d, let I_d be the ideal of every coefficient of
det JF-1 in the full normalized input coefficients, and let p_j range over
the inverse block. The degree-at-most-d complex Jacobian statement is
equivalent to p_j in radical(I_d) for every j. Each obligation has the
finite polynomial certificate

\[
1\in(I_d,1-zp_j).
\]

The quadratic compiler proves the stronger ordinary ideal membership for
all 18 cubic/quartic inverse coefficients, with exact coefficient identities.
For higher degrees, a nonzero ordinary-ideal remainder would not refute
radical membership. This makes the next computation a well-typed algebraic
question and an independently checkable identity; it does not establish
those identities for all degrees.

## 6. Next decisive targets and evidence

| Lane | Finite object now checkable | Remaining precise target |
|---|---|---|
| Collatz induction | Actual asynchronous join to a smaller certified source, with source guards | Construct such joins for uncovered sources, with a well-founded dependency |
| Collatz closed dependencies | Verified loops and charge gcd, retaining entry times | Obtain useful charged loops without assuming the source is already rooted |
| Algebraic plane | Chosen-root coordinates, exact incidence and colouring proof | Find a six-chromatic packet while escaping known low-chromatic fields |
| Triangle packing | Complete triangle list, matching and complementary cover | Use packing as an auxiliary graph-search coordinate without discarding colour constraints |
| Planar JC | Full coefficient constraints and finite inverse block | Prove block vanishing for a larger complete polynomial family |

All four packages include proofs of their infinite-quantifier elementary
claims, finite universes, positive controls, hostile inputs, scripts and
deterministic outputs. Normal and optimized Python outputs agree, and each
package has an independent peer audit. The packages perform **308,328 exact
checks**: clocks 238,012; prime/partition 68,370; geometry 1,667; and
Kervaire/superstring/JC 279, including the symbolic quadratic extension.
These tests verify the new interfaces, not the accepted external theorems.

The concurrent audit commits preserve the THM-4581 coalescence core while
retaining its time/count conventions and separating the rate sketch.
They resolve HYP-9213/9214/9220 in the Haar model. The weighted fractional
moment controls the affine translation; the proposed T^(-1/2) rate has a
sketch-level upper bound up to logarithms, and its leading constant
(|k0|+E[J])sqrt(4/pi) remains heuristic. The
[concurrent coalescence note](oai3_two_orbits_twos_and_threes_20261007.md)
retains those distinct statuses and the corrected timing conventions.
They do not turn Haar-almost-sure hitting into prescribed-integer hitting.
The prior [orbit-information synthesis](collatz_orbit_information_openai_synthesis_20261007.md)
remains the route for the long-gap exponent laws and decimal digit bounds.
The narrow next Collatz test is a new guarded family that actually reaches
an uncovered source and yields a paid smaller dependency, rather than
another invariant of an already descending family.
