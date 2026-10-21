# Independent audit of the dyadic conductor and blue-fibre bridge

Date: 2026-10-21. Final audit target: the approved current source files in `current_inputs/`. Exact input hashes and source paths are in `FINAL_INPUT_PINS.json`. This package contains only current results and current verification receipts. No producer or repository files were edited by the auditor.

## Verdict

**FINAL PASS for the current dyadic packet.** The all-level dyadic class group, literal reduced-form labels, extension maps, square/genus coordinates, positive finite-blue-parent bridge, two exceptional stabilizers, conductor-16 nonsplitting, and nonliftability of shore exchange are sound. The restricted idoneal conclusion is unconditional and does not use global completeness of the 65-number list.

The theorem is explicitly restricted to dyadic conductors f=2^a, a>=3. The general residue-to-class map requires gcd(N(alpha),f)=1; odd norm is equivalent only in the dyadic case. The current proof and checker retain the exact nondyadic boundary below. All current input pins and both verification modes pass.

An exact hostile is f=40. Its residue quotient has 24 classes. The norm-mod-8 kernel has 6 classes, while the square subgroup has 3. The class represented by 1+8 omega has norm 57 congruent to 1 mod 8 but is not a square. Analytically, the quadratic character of its norm modulo 5 is -1: this character is well-defined on G40 because scalar norms are squares and Eisenstein units have norm 1, whereas every square class has character +1. Its contraction form reduces from (57,30,25), via (25,-30,57), to (25,20,52). The independent checker retains the exact quotient hostile. The valid dyadic theorem is unaffected.

## 1. Proper ideal classes and the unit/scalar quotient

The relevant group is the Picard group of proper invertible fractional ideals. For f=2^a, a residue has odd norm precisely when it is a unit modulo fR. An integral ideal alpha R of odd norm is prime to the conductor, so its contraction is invertible. This excludes the noninvertible conductor ideal.

Conrad's exact sequence (5.7) gives the quotient of (R/fR)^* by scalar residues and the image of R^*. The maximal order R is norm-Euclidean and has trivial class group. For f>=2 there are 3 f^2/4 residue units. The scalar subgroup has f/2 elements, and omega contributes a disjoint factor of order 3, since its three powers are distinct in F4^*. Thus the denominator has 3f/2 elements, not merely f/2, and the quotient has f/2 classes.

The map is genuinely alpha mod fR -> [alpha R intersect O_f], with that direction verified directly in Conrad's equation (5.4). Multiplying alpha by an Eisenstein unit does not change alpha R. Multiplying by an odd scalar changes the contraction class by a principal ideal. For alpha congruent to 1 mod 2R, a denominator representation must use a unit reducing to 1 in F4^*, hence only plus/minus 1. Therefore alpha=A+B omega is trivial exactly when f divides B. This addresses both scalar and unit ambiguities, including the factor of 3.

For g|f, extension of prime-to-f contracted ideals gives (alpha R intersect O_f) O_g = alpha R intersect O_g. Both contractions are the unique ideals prime to the conductor with the same extension to R. Thus residue reduction models the actual ideal-extension map, not a homomorphism chosen from matching cardinalities.

## 2. All-order valuation and group structure

Put w=1+4 omega and j=1+2 omega. Since j^2=-3 and j is not scalar modulo 4, its class has exact order 2 at every f>=4. If w^(2^r)=A_r+B_r omega, squaring gives A_(r+1)=A_r^2-B_r^2 and B_(r+1)=B_r(2A_r-B_r). Starting from odd A_0=1 and v2(B_0)=2, induction keeps A_r odd. Since B_r is divisible by 4, v2(2A_r-B_r)=1 exactly. Hence v2(B_r)=r+2 for every r>=0.

The quotient order is a power of 2. The scalar criterion therefore proves the exact order 2^(a-2) of [w] at conductor 2^a, including its triviality at a=2. Every power of w is scalar modulo 4, whereas j is not. Their subgroups intersect trivially. Their product has the full quotient size, proving C_(2^(a-2)) x C2 for every a>=2, not just the checked range.

Reduction to conductor 2^(a-1) reduces the w exponent modulo 2^(a-3) and preserves the j bit for a>=3. The final G4->G2 map is trivial. No finite computation is being used to justify the universal quantifier.

## 3. Actual forms and extension labels

The sign convention is important. The paper's formula matches the conventional form-to-ideal map

    (A,B,C) -> <A, (-B+sqrt(D))/2>,  D=-3f^2.

For alpha=a+b omega of norm N and gcd(b,N)=1, reducing R modulo alpha identifies omega with -a/b modulo N. Contracting to O_f gives the lattice <N, f omega-r>, with r=-fa/b modulo N. This is exactly the displayed ideal convention when B=f+2r modulo 2N. Thus the signs in the conductor-16 labels are not an accidental inversion of the chosen ideal convention.

The independent checker does NOT reuse contraction_form or reduce_form. It enumerates primitive reduced forms directly. After the proper swap (A,B,C)->(C,-B,A) when needed, it takes the Eisenstein Euclidean gcd of the two actual ideal generators A and (f-B)/2+f omega. The gcd is a generator alpha of the extended R-ideal; its norm is checked to equal A. Its normalized residue gives the class coordinate. This recovers every printed form at f=2,4,8,16,32,64,128, including:

- f=8: 1 -> (1,0,48); w -> (4,4,13); j -> (3,0,16); wj -> (7,2,7).
- f=16: 1 -> (1,0,192); j -> (3,0,64); w -> (13,8,16); wj -> (7,4,28); w^2 -> (4,4,49); w^2j -> (12,12,19); w^3 -> (13,-8,16); w^3j -> (7,-4,28).

All printed extension targets agree with the independent ideal-generator computation. For the written form-extension shortcut, primitive forms of even discriminant have even B, and A,C cannot both be even. Hence a proper swap always provides odd A. Extension from f=2g replaces the square-root congruence B by B/2 modulo A; the unique even lift B' modulo 2A yields the stated form. This proves the shortcut's integrality and correct orientation.

## 4. Dyadic norm, genus, and idoneality

For dyadic f divisible by 8, changing a residue representative modulo f preserves its norm modulo 8. Odd scalars have square 1 modulo 8, and Eisenstein units have norm 1. Therefore the norm is a well-defined character of the quotient. N(w)=13 congruent to 5 and N(j)=3 give 5^t 3^e. These two residues independently generate (Z/8Z)^*, so the kernel is t even and e=0. This is simultaneously the square subgroup and the kernel of extension to G8.

Calling the square quotient the genus quotient is the classical principal-genus theorem. The character table and all signs in the packet are correct. In particular chi_-4 kills A and detects the j exponent; chi_-8 kills B; chi_8 kills C. The ramified prime above 3 is jR, and its norm is odd, so its contraction is the proper invertible class B. This decoration is genuinely arithmetic.

The inherited classification d=3m^2 translates to even f=2m and discriminant -4d=-3f^2. Its odd-factor argument is correct: the class-number formula for even f has a factor 3 if f contains 3 or any odd prime >=5, because every such prime ell has 3 dividing ell-chi_-3(ell). Exponent dividing 2 therefore forces f dyadic. The audited tower then gives exactly f=2,4,8, hence d=3,12,48. Grube's equivalence, via one-class genera, is correctly invoked. Neither this argument nor the dyadic proof uses the global 65-idoneal completeness question.

## 5. The positive blue-fibre bridge and its limits

The pinned balanced companion uses exactly the added membership bits in its left and reflected-right subset coordinates. Above a finite blue parent the old coordinates agree, so the child is blue precisely when epsilon=eta. The zero-bit section specifies the origin. The map

    theta(epsilon,eta)=B^epsilon C^eta=w^eta j^(epsilon+eta)

is a group isomorphism. It sends the two toggles to B,C and their diagonal product to A. Actual conductor extension kills A and retains epsilon XOR eta. Hence the two blue children map into the true conductor kernel {E,A}. The three character identities are exact, including their left/right labeling.

Above a finite black parent the old subsets remain unequal after inserting the new bit. Every child is black, including parity-zero children. Thus the correct predicate is parent-blue AND new-bit equality, and parent coincidence is necessary. A four-element quotient label cannot silently replace this retained predicate.

At a left-infinite parent the left toggle fixes the fibre, giving stabilizer <B>; at a right-infinite parent the right toggle fixes it, giving <C>. These distinct order-2 stabilizers both occur. The conductor quotient instead has stabilizer <A>. No single group isomorphism can identify both exceptional quotient types with the one conductor quotient. The complete local G8-set description using G8, G8/<B>, G8/<C>, and a singleton is valid, provided the inactive-port type is retained. A free fibre without a chosen origin is a torsor, not a canonically pointed group.

Swapping the two shores sends (t,e) to (t+e,e), fixing A and swapping B,C. This is the nonidentity automorphism preserving the conductor flag. It is not complex conjugation: ideal conjugation is inversion, which is trivial on G8. It also fails to preserve the ramified-prime decoration B.

## 6. Nonsplitting, carry, and nonliftability

At conductor 16 the class w has order 4, and the extension kernel is <w^2>. Both lifts of A are w and w^3; both square to the same nontrivial kernel class. A homomorphic section would have to lift the target involution A to an involution, so no section exists.

The stated carry rule follows by adding 1 modulo 4 to t=t0+2t1. It is not the rule for adding an independent commuting involution. At every later adjacent dyadic extension, a lift of the target's cyclic generator has twice its order, proving nonsplitting for every a>=3.

The 2-torsion subgroup of G16 is characteristic and projects to {E,B}. A commuting lift of shore exchange would preserve this image; it does not, because B maps to C. Thus shore exchange cannot lift even as an abstract group automorphism commuting with extension. Both the producer and the independent checker also exhaust the eight automorphisms of C4 x C2 and find none with the required action.

## 7. Scope firewalls retained

Negation fixes any ideal as a set, and conjugation acts by inversion. Neither is translation by a nonidentity class. A nonzero z in O_f generates a principal ideal, so norm-ball points cannot become the residue classes used here merely by renaming them. The residue generators w and j need not belong to the smaller order.

The local class-group bridge does not provide a map of global tournament metagraph incidence. The balanced companion retains the bypass-edge correction; nothing in the class-group quotient records those labelled cut ports. No electrical weights, Laplacian, current, response map, global metagraph classification, native-source identity, or global idoneal completeness follows from this audit.

## 8. Reproduction and independence

The final producer checker succeeds in normal and optimized Python and reproduces both final frozen receipts byte-for-byte: 28,754 live checks. All seven final listed producer SHA256 checks pass. The inherited divisor-norm-ball file hash agrees with SOURCE_PROVENANCE.json. The final companion text was separately copied and hashed.

The independent checker uses a second finite group model. Each dyadic residue class has a unique representative 1+s omega with s even: first rotate by omega until its F4 residue is 1, then scale by the inverse of its odd real coefficient. Its multiplication is

    s * t = (s+t-st)/(1-st) mod f.

It exhausts all residue units and verifies the normalization fibre sizes, norm ambiguity, multiplication, coordinate bijections, all reduced-form ideal generators and extensions, square kernels, port predicates/stabilizers, nonsplitting, and nonliftability. It imports no producer code. Both modes pass **46,719** explicit checks against the final pinned packet and have byte-identical JSON output. The all-order mathematical conclusions come from the written proof audit, not the finite computations.

Run from this audit directory:

    python3 -B independent_checker.py > /tmp/independent.normal.json
    python3 -B -O independent_checker.py > /tmp/independent.optimized.json
    cmp /tmp/independent.normal.json /tmp/independent.optimized.json
    cmp /tmp/independent.normal.json independent.normal.json

## 9. Primary-source verification

Sources were independently opened on 2026-10-21. SOURCE_CHECKS.json records the verified locations and the exact local dependency hashes are in FINAL_INPUT_PINS.json.

- Keith Conrad, The conductor ideal of an order, Theorem 3.8 (printed pp. 6-7), equation (5.4) and its discussion (printed p. 15), equation (5.7) and Theorem 5.3 (printed p. 16): https://kconrad.math.uconn.edu/blurbs/gradnumthy/conductor.pdf
- Ernst Kani, Idoneal numbers and some generalizations, section 2.2 equation (2.2) (printed p. 200), Theorem 6 and Remark 7(b) (printed p. 202): https://www.labmath.uqam.ca/~annales/volumes/35-2/PDF/197-227.pdf
- MIT 18.785, Fall 2018 Problem Set 9, Problem 5 introduction (printed p. 6) and parts (b)-(c) (printed p. 7): ideal-sign convention, nonmaximal-order scope, inverse form, and ambiguous boundaries. https://math.mit.edu/classes/18.785/2018fa/ProblemSet9.pdf

No extra research claim or literature-priority claim was added.

## 10. Approved current hashes and self-contained replay

- CONDUCTOR_TOWER.md: b6911f6479f9cb75a26b168573c97d6afa15843a61e535d63ca72a71088fd9f3
- verify_conductor_tower.py: 04cbac2aaed5262ead2226da02713b94dc8ce0d8e1412295eaac2e3ea7311a2a
- verification.json and verification.optimized.json: 55029b081f5f387aa6d86cd356117b3f85046f5c4237b80deb217222aa60b1e4
- BALANCED_FIBRE_GROUP_AND_PORT_DEFECT.md: 4ce1d99cfff02b97997551fa1edf35a3f6baf64755d45d654e6d9da0afd88603

The independent checker is byte-identical to the approved audit program. FINAL_INPUT_PINS.json selects current_inputs/, so every active source dependency is contained in this package. Normal and optimized execution was also repeated from an isolated temporary copy containing only these files; both outputs match the retained 46,719-check receipts. The producer program and its current 28,754-check receipts are likewise unchanged and reproduced in both modes. No superseded source copy or earlier execution receipt is included.

PUBLIC_APPROVAL.json is the machine-readable approval manifest. SHA256SUMS inventories this complete current-results package.
