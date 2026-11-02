# Independent audit: width-three local universality

**Verdict: ACCEPT the all-modulus theorem in its stated scalar scope.**
Reviewed `width_phase2/PROOF.md`, final SHA256
`b7da98667ef7547dd8814d522363cde9c7209f39283bfe4bce7a41cb434d4bf2`.
The final convention indexes local width by positive integers, so K1 has
local width one; this avoids a width-zero convention ambiguity.

Let the grammar start at K1 and use only ordinal join J(T,U) and
B(T,U)=C3[T,U,K1]. For every M>=1 and every parity-compatible residue z,
this grammar supplies arbitrarily large distinct H values of strong
tournaments with H=z modulo M. Its strong scalar spectrum therefore has
profinite closure exactly the odd profinite integers. This is **not** an
exact integer-spectrum theorem, an all-context profile-universality theorem,
or a global bound on exact spectral substitution width.

## 1. Polynomial carrier and actual contexts

The falling-factorial path-cover polynomial F is an integral carrier. The
unordered path-cover substitution kernel is bilinear in the two module
profiles, and those coordinates are linearly equivalent to F. Transitive
factors give the monomial basis x^r, r>=1. In a two-transitive-factor cyclic
substitution, every cycle uses the singleton; selecting nonempty increasing
subsets of the two factors and retaining even total size gives exactly the
two finite-difference products in the displayed binary law. Thus its
monomial verification does prove the full bilinear identity.

The source-padding law F->xF and ordinal multiplication are correctly used.
The unary C,D formulas, including H(D(T))=3F(2)-H(T), are evaluation
identities of actual positive contexts, not permission to subtract graphs.

## 2. Reset, sidecars, and positive inverse scope

For M=2^a 3^b n with gcd(n,6)=1, the required N exists by finite unit-group
periods. Its size ensures both reset conditions g=0 modulo 2^a and f=0
modulo 3^b. At the auxiliary x=3 evaluation on the n-part, the equality of
the auxiliary values at 2,3,4 removes every uncontrolled higher evaluation.
This validates the U2,U4 triples and all four gate-table rows.

Independently multiplying the two-by-two matrices gives U's shear 15/4,
its primary actions (-27;16,-2), and the R,E matrices. The stated powering
indeed fixes both primary anchors. The shear coefficient has exactly one
factor of five on the coprime-to-six part, so gcd(tau,n)=gcd(5,n).

D1,A,V are permutations of the displayed finite ambient controlled state.
Their actual contexts preserve the reset class. Positive powers representing
their inverses therefore suffice for every future controlled response,
without returning the full polynomial or graph to an earlier value. E has
determinant -13/8 and is never inverted. This is an essential valid scope
distinction, not a hidden full-profile periodicity assumption.

## 3. Primary anchors and exceptional primes

The four literal seed pairs are correct. The unit-generation argument from
the full 24-unit group modulo 72 lifts: any proper subgroup at a higher
2/3-primary modulus has a quotient character to C2 or C3, and every such
character factors through the appropriate divisor of 72. This uses the
elementary finite prime-power unit structures, not a primitive-root conjecture.

Dividing off the prescribed three-adic valuation and solving the remaining
unit anchors is valid even when z is zero modulo 3^b. Adding a period of
five forces a positive H=5 seed exponent without changing the anchors.
The H and g seed factors then have no common prime on n, and at five the
state is deliberately (0,unit). The arbitrary-nonzero-vector replacement
would be false because of the fixed diagonal mod-five hostile.

The first shear makes h a unit on n0=n/5^c and R transfers it to g. The
ten-state mod-five R/E orbit has a state over every first coordinate with
nonzero second coordinate. Before E, setting h=0 on n0 makes the new g=-2g;
before R, setting h=g makes the new pair (g,g). These shears leave mod five
unchanged. Consequently the prime-thirteen singularity is avoided at every
actual E application, while the finite mod-five word is still implemented.
Once g is a unit everywhere, the final shear congruence has exactly the
already-enforced mod-five compatibility condition and no other obstruction.

All these controls act on one actual evolving tournament, so the CRT
conclusion does not splice together independently realized residues.

## 4. Strong extraction and unbounded fibres

Adding n to the final S exponent preserves the controlled state and makes
the last S block positive. Its final chronological primitive gate is A,
whose output is source padding followed by B(T,U2). Removing only that
final padding preserves the exact H and leaves a strong cyclic substitution.
No removal is performed before later gates that might detect the lost
profile coordinates.

For actual nonempty T,U, the three cyclic rotations of whole-module
Hamilton paths in B(T,U) give at least 3H(T)H(U) distinct Hamilton paths.
Thus every U=A D1^2 A^(-1) V includes three strictly growing cyclic contexts
and increases H by at least a factor 27; positive inverse powers, joins,
and padding cannot undo that growth. Adding further positive shear periods
therefore yields arbitrarily large distinct exact counts, not just larger
graphs with the same count.

Width at most two gives only transitive tournaments. Together with parity,
the theorem yields the claimed positive-index local-width classification.
The inherited exact-width witnesses at 13,19,31 remain inherited results;
no new exact-width proof or external priority claim is asserted here.

## 5. Independent replay and hostile controls

`audit_width_universality.py` imports no producer code. Its direct foundation
uses the independently written unordered-set-partition/Hamilton-subset
kernel in `audit_width_states.py`, not the producer's literal cut engine.
It checks all 121 ordered pairs of labelled factors of orders one through
three, at each integer x=-2..6, for both ordinal multiplication and the
binary cyclic law. All four seed pairs are separately checked on literal
graphs.

The independent controller computes primary periods and positive inverse
powers from its own finite matrices, finds anchor exponents by unit-group
BFS, and finds safe shears by direct finite congruence search. It discovers
the mod-five steering words by BFS rather than using the producer's table.
Every compatible residue modulo 1..96 and selected mixed/prime-power cases
is landed simultaneously. Small complete recipes are also executed as
full polynomial functions modulo M, with every inverse expanded into
positive primitive gate calls. The larger cases are explicitly controlled-
state checks, not claimed full graph expansions.

The fixed mod-five diagonal and E(2,3)=(0,0) modulo 13 are retained as
hostile controls. Every full replay checks both reset sidecars after every
primitive gate, and every recipe ends with a positive A block for strong
extraction. The finite receipts support the independently reviewed proof;
the all-modulus conclusion does not extrapolate from their tested range.

The terminal receipts contain 3,562 simultaneous landings and 132 complete
full-polynomial positive recipes, including six mixed cases at moduli
35,65,120. Those full recipes expand 33,705 primitive gate applications and
2,566 inverse requests into positive powers. Normal/-O receipts are
byte-identical with SHA256
`523f754eeb2ed25a23918d7c064fec66c12eeedf78f531bc3a649768c0ea7987`.
Minimal local script dependencies are `audit_width_states.py` and
`tournament_kernel.py`; no producer files or witness snapshots are needed
to run this checker.

    python3 -B audit_width_universality.py > width_universality.normal.json
    python3 -B -O audit_width_universality.py > width_universality.optimized.json
    cmp width_universality.normal.json width_universality.optimized.json
