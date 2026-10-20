# Independent genus-six frontier audit

20 October 2026. **PASS in the stated primitive monic generic-source scope.**

## Verdict

The main frontier classification, all six explicit polynomial realizations,
the supplemental parameter families, and the sharp quartic/sextic
codimension-five-to-six boundary pass independent audit. No mathematical
correction to a theorem, identity, finite universe, or claimed consequence
was required.

The complete symmetric genus-six list has six members. Exactly four are
planar degree semigroups at infinity. The nonprincipal characteristic
presentations are exactly

    (2,13), (3,7), (4,5), (4,10,5), (4,6,9), (6,9,4),

up to reversal of the first two entries. Both alternate presentations are
genuine primitive polynomial realizations, not consequences of an assumed
realizability converse.

The inherited all-quadratic-semigroup theorem removes the first
presentation in the actual Keller setting. The five remaining unoriented
terminal degree pairs are exactly

    (3,7), (4,5), (4,6), (4,10), (6,9).

This is a necessary frontier map. It does not prove any of these terminal
pairs possible for a noninvertible Keller map, exclude generic
codimension six, improve the conductor bound to fourteen, classify all
embeddings, or prove the planar Jacobian conjecture.

## Independent computation and finite universes

`independent_audit.py` imports neither supplied checker. Its **119**
explicit checks pass normally and with optimization; the JSON outputs are
byte-identical.

The semigroup census is genuinely different from the producer's subset
censuses. If the multiplicity is m, then 1,...,m-1 are gaps, so m<=7.
For each m=2,...,7 the audit enumerates positive Kunz-coordinate tuples
of length m-1 with sum six. There are 32 tuples in the entire universe.
The Kunz inequalities select numerical semigroups; conductor twelve
selects exactly the six symmetric genus-six cases. This is exhaustive,
not an arbitrary numerical cutoff.

For each generating set, a shortest-path computation on residue classes
calculates the exact Apéry set. Membership, primitive generators,
conductor, and gaps then follow for all nonnegative degrees. This avoids
using a finite monoid horizon as the independence check.

The characteristic reconstruction uses the cited initial bound thirteen
and later bound eleven, recursing until the gcd reaches one. It imposes
no sequence-length cutoff. The maximum possible initial gcd is four;
strict later gcd drops allow at most four entries. The length-four
positive control (8,12,10,5) has conductor twenty. No length-four
sequence has conductor twelve. Selecting by the exact Apéry conductor,
rather than the characteristic conductor formula, recovers the twelve
oriented sequences. The formula is checked afterwards.

The producer's full Cartesian counts 108, 1,188, and 13,068 also
reproduce. Its 129-check main receipt and 48-check supplementary receipt
were replayed normally and with optimization. Each replay matches its
provided receipt byte for byte.

## Algebra and realization audit

The field-index argument is valid over every characteristic-zero ground
field: for every nonconstant f in the algebra, the intermediate index
divides [k(t):k(f)]=deg(f). The displayed degrees have gcd one, proving
primitivity.

Monic p makes k[t] finite integral over k[p]. Primitivity gives minimal
polynomial degree deg(p) for q over k(p); integral closedness puts its
coefficients in k[p]. Thus the power basis is a basis of the entire
algebra, not just a subalgebra.

The triangular basis change also works for
r=q²-U(p)q-V(p), not only r=q²-H(p). This was already the correct
argument in the beta-defect section; the author explicitly generalized
the preliminary lemma during review. The independent checker verifies
the universal determinant-one transitions for four and six entries.

Distinct residue degrees prevent leading cancellation and identify the
entire degree set and codimension. The two alternate degree lists
0,10,5,15 and 0,9,4,13,8,17 have the claimed exact Apéry minima.
Every concrete derivative gcd is one.

Here immersion means differential immersion: no simultaneous zero of
p′ and q′ over an algebraic closure. It does not mean injectivity or a
locally closed scheme embedding. Likewise, the producer's resultant
degree check is only a consistency control; resultant degree alone
would not prove primitivity. The field-index proof above is load-bearing.

## Supplemental family audit

All parameter identities were expanded independently. With D=a²-4b,
the sextic remainder has degree four precisely on aD!=0. Its full
module basis therefore gives semigroup <4,6,9> and codimension six.

The derivative resultant is exactly 6561 D⁴/16384. Together with the
value of q′ at zero, this proves immersion precisely when
D(a²+4b)!=0. On the genus-six stratum, the divisor a²+4b=0 has the
single common derivative zero t=0. The independent symbolic gcd is t.

On a=0,b!=0, the remainder is b⁴t²/64. Thus the algebra contains t²,
and the even/odd decomposition gives exactly <2,9>, codimension four,
and conductor eight. On D=0, p=h² and q=h³ with deg(h)=3, so the
field index is three, not one. The quartic b=0 boundary likewise has
field k(t²), not k(t).

For the narrow even/odd ansatz, the audit starts from all four unknown
odd coefficients of q and independently solves the degree-16, 14, 12,
and 10 equations. The degree-eight obstruction is exactly
D(5a³-36ab+32u)/128. For D!=0 this forces the advertised family; the
degree-six equation determines B. For D=0 the whole exceptional branch
has p=h², q=h³+c h, and zero remainder. Since
h=q/(p+c), its field is exactly k(h). No primitive exceptional family
was lost by division. The completeness statement is correctly confined
to this normalized parity ansatz.

## Sharp beta boundary and inherited scope

The first odd cancellation degree eleven has coefficient 2A and
codimension seven. After A=0 and the degree-ten subtraction, the next
coefficient is beta=2B-3b. For beta!=0 the whole-module basis has degrees
0,6,9,15 and codimension exactly six. Consequently delta<=5 forces
beta=0, whereas delta=6 does not.

The full scalar-retaining Jacobian was differentiated independently
using symbolic coefficient functions. Its degree-seven coefficient is
6a′-4c′ and its degree-six coefficient is -2 beta′. Hence the Keller
condition forces beta constant, not necessarily zero. The displayed
one-variable hostile has beta=2 and is not a Keller counterexample.

The exact predecessor proof, predecessor independent review, its source
note, the all-quadratic-semigroup proof, and that proof's audit were read
and hash-checked. The bound-twelve result is inherited with its stated
primitive monic generic-source hypotheses, not re-proved or strengthened.
Constant-leading terminal reduction, scalar extension, and the single
rational point at infinity have the required hypotheses. No low-degree
generic generator is silently globalized.

## Primary-source audit

The author-hosted [Sathaye–Stenerson text](https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf)
was freshly opened. Section 3.1 supplies the characteristic axioms,
conductor expression, structural bounds, and primitive-entry property.
Its characteristic-zero and one-place-at-infinity hypotheses match the
application after scalar extension. Its distinction between plane curves
and polynomially parametrized curves is respected.

The [original Abhyankar–Sathaye article, page 1062](https://scispace.com/pdf/uniqueness-of-plane-embeddings-of-special-curves-omujaliiec.pdf)
was freshly retrieved as indexed primary-PDF text. Section 2 identifies
the first entries with actual coordinate pole orders and describes
principal-pair reduction. Direct DOI/PDF opens failed; no downloaded
copy is claimed. No uniqueness theorem from that article is used.

## Reproduction and reviewed bytes

From this audit directory:

    python3 independent_audit.py > independent-results.json
    python3 -O independent_audit.py > independent-results-optimized.json
    cmp independent-results.json independent-results-optimized.json

The initial main proof received for audit had SHA-256

    01cc72919f4184f032dc9686540eda9efe5fd446b7c3032beef6faebd72ef14e

The reviewed main proof after the harmless monic-quadratic clarification
had SHA-256

    09015a1d536a1fcd68e4be94d9c26f0b27d14a0f737964f902f5e9cc90fbeb37

The initially inspected supplementary proof had SHA-256

    f652fdc72c50588d96bb8c86d450e5df8d1e7d0b0ab7e194895566d9dd65c6fe

The final inspected main proof, after explicit immersion and resultant
scope clarifications, has SHA-256

    709dc8bd7e95bef496c76fca824e30918e2a4a2a6747f0594fbea0438af34323

The final inspected supplementary proof, after its matching immersion
clarification, has SHA-256

    645f78de81ad90e8840dc11227885f13e1cfe81a564e6b4414a0cf491e51d107

These review changes clarify valid arguments and terminology. They do
not change the formulas, universes, hypotheses, or conclusions.

The exact final inspected packet and audit artifacts are also recorded
in `AUDITED_HASHES.json`, which preserves the inspected source-version
hashes. Publication-copy changes are recorded in
`../PUBLICATION_COPY_LINEAGE.json`.

For a publication-layout copy, dependency-map keys may be changed from
staging paths to the corresponding sibling publication-relative paths,
provided every pinned dependency hash is unchanged. The script derives
its dependency root from the containing packet's parent, so no
mathematical code change is needed. Check labels containing those paths
will change, requiring a fresh normal/optimized replay and a recorded
original-to-publication hash lineage. A replay through local sibling
copies or links verifies that copied layout only; it is not remote
publication verification.
