# A clock-optimal eighteen-source circuit in the full twenty-eight-option supply

Status: FINITE-EXACT. In the twenty-eight specific authenticated native options on targets U=1965 and V=2234, with one slot per source, eighteen-source primitive circuits exist. Their minimum endpoint owner-clock cost is thirty-seven. Up to reversing every sign, exactly one primitive eighteen-source relation has this minimum cost.

This is a theorem about the supplied fixed options and their complete pinned native domains, not a minimum over every order-eight target pair, owner frame, path menu or full-atlas rewiring.

## 1. The unique minimum-cost relation

The positive and negative source sets are

    A = {1394,1425,1457,1480,2176,2185,2581,2835,2945},
    B = {1477,1549,1691,1921,1974,2837,2905,3000,3369}.

Move each source in A from U to V and each source in B from V to U. There are nine sources on each side. All eighteen sources are distinct and lie outside the target pair.

Seventeen of these options are the old two-clock transpositions. The sole added source, 2835, changes three clocks. Thus the new circuit costs 17*2+3=37 owner-clock changes. It escapes the old twenty-option sixteen-source ceiling by adding one genuinely different source option, not by reusing a normalized port or disguising edge multiplicity as another literal path.

For the original primitive potential X, both unsigned source sums are exactly 1161. Hence this circuit's chosen integer X current debt is zero, unlike the separate forty-three-clock circuit's modular debt of 512. The remaining four binary charge coordinates also cancel.

## 2. Literal native realization and preserved information

For each selected source the checker reconstructs its actual tournament matrix from the original canonical code. It fixes every common indexed owner clock of the two paths, then enumerates all permutations of the remaining owners. Exactly TWO actual directed Hamiltonian paths remain in every selected slot. No extra-path restriction is smuggled into the binary model.

The shared roles are the actual first owner, last owner and two ordered owners of an internal indexed arc. Their local occurrence clocks, all additional common clocks, source matrices, full owner-to-clock maps, canonical source codes and trivial source automorphism groups are retained. Both literal complementary matrices and both reflected-converse matrices and paths are checked. These are actual owner/kernel lifts, not extra resource capacity.

The F/R resource orbit of each selected path has four distinct normalized words. Both endpoints contain eighteen different normalized black resources, or seventy-two closed ports. The old and new alternatives are disjoint, giving thirty-six normalized resources and 144 closed ports across their union. The full supplied inventory contains fifty-seven distinct normalized resources, or 228 closed ports, including the one additional raw pinned path at unused source 1923.

The declared shared owner coordinates do not assert a common-host tournament amalgam. Different source classes may disagree on induced arcs not included in the pin record.

## 3. Full nodewise dyadic-current certificate

Let

    u = sum_(s in A) e_s - sum_(s in B) e_s,
    v = e_U - e_V.

The actual paired black Laplacian perturbation is

    Delta L = 2(uv^t + vu^t).

It has zero row sums and preserves the source and target degree margins. Two complete integer vectors, with all integer corrections retained, verify by direct multiplication over the original 438216 edge rows that

    L z_u = 2u + 256 w_u,
    L z_v = 2v + 256 w_v.

Every harmonic character modulo 256 therefore annihilates both charges, and Delta L h=0 at every physical owner. The checker does not trust the optional solver's elimination algorithm; it multiplies the saved equations anew.

The independently audited original critical group's two-primary part is

    (Z/2)^3358 + (Z/4)^4 + Z/256.

Its exponent is 256. Therefore the mod-256 certificate preserves every genuine dyadic harmonic bank, with the correct scaling of higher characters and the constant gauge killed by zero row sums. Odd-primary guards are not claimed.

The effective paired-charge coordinates are the four Y parities and X modulo 128, giving Q=(Z/2)^4+Z/128. The original targets have X_U=X_V=146, so each selected ordinary per-slot X phase modulo 8 agrees as well. This is not preservation of every ordinary per-port character: the genuine ordinary mod-8 harmonic 4e_U distinguishes the targets, as checked directly.

The selected packet preserves its baseline owned current signature. It need not be harmonic as a packet in isolation. It preserves a compatible remainder's guard if that remainder supplies the original cancelling signature. The fully occupied canonical atlas has not thereby acquired spare replacement ports.

## 4. Exact primitive-menu check

All 262144 assignments to the complete pinned binary domains are checked. Exactly 48620 preserve the target margins. Exactly 754 preserve the full original mod-8 bank. Precisely two preserve the full mod-256 bank and hence all dyadic guards: all-before and all-after.

Every other balanced assignment has a directly checked authentic harmonic separator among 64Y_0, 64Y_1, 64Y_2, 64Y_3 and X modulo 256. The respective separator counts are 24380, 12190, 5990, 2990 and 3068, totalling 48618. Their target-owned current is nonzero.

Thus the eighteen-slot circuit is primitive in its complete pinned menu. Its exact simultaneous changed-slot budget is eighteen, and its exact changed-owner-clock budget is thirty-seven. These lower bounds can change after additional charged slots or literal path alternatives are authenticated; they do not forbid a decomposition in a larger resource inventory with auxiliary slots eventually restored.

## 5. Why thirty-seven is sharp, and why the minimizer is unique

Every nonidentity permutation changes at least two indexed owner clocks. Thus eighteen changed slots cost at least thirty-six.

The supplied options have cost distribution

    twenty cost-two options, five cost-three options, three cost-four options.

A cost-thirty-six eighteen-source relation must select eighteen of the twenty cost-two options. A cost-thirty-seven relation must select seventeen cost-two options and exactly one cost-three option. These observations enumerate ALL supports at those two boundary costs.

The five-bit image of Q/2Q is a useful sign-independent support constraint: the four binary source coordinates must xor to zero, and the total cyclic coordinate must be even. For each surviving support, put its smallest source on the positive side to quotient by global sign reversal, then test every balanced partition.

The complete counts are:

- Cost 36: 190 support sets; 3 parity-admissible supports; 72930 balanced partitions; 1129 equal-Q relations; zero primitive relations.
- Cost 37: 5700 support sets; 183 parity-admissible supports; 4448730 balanced partitions; 69358 equal-Q relations; exactly one primitive relation, the displayed A/B pair.

For each equal-Q relation, a proper balanced relation exists if and only if one exists with at most four sources per side: a larger proper subrelation can be replaced by its smaller complement. A one-per-side relation is excluded by the distinct source labels. The checker tests all cardinalities two, three and four.

At cost 36 the first proper witnesses have half-sizes two, three and four with counts 288, 832 and 9. At cost 37 the corresponding counts are 22680, 45657 and 1020; the remaining one relation is primitive. A separate pure-Python boundary checker agrees with the C++ discovery census, and its normal and python -O receipts are identical.

Consequently thirty-seven is the exact minimum in this supplied twenty-eight-option family. The single three-cycle at source 2835 is forced by the unique minimizing relation in this scope. No claim is made that every possible native eighteen-source circuit requires that source or has cost at least thirty-seven.

## 6. The precise S2/S3 connection

For each selected slot, compare the old and new owner occupying each clock. This defines a relative permutation of THAT source's eight owner labels. Seventeen relative permutations are transpositions. At source 2835 the paths are

    old: 2,1,7,5,6,4,0,3
    new: 2,1,7,5,0,6,4,3.

The free owners (6,4,0) at clocks 4,5,6 become (0,6,4). The relative owner map is the three-cycle (0 4 6). This is a concrete distinction between an order-two relative permutation and an order-three one, not numerology based on slot count.

It is NOT a legal S3 action on the native domain. Applying this same three-cycle again gives the ambient permutation

    2,1,7,5,4,0,6,3,

which is not a Hamiltonian path in the actual source matrix. The source's pinned native domain still has only two states. Nor are these nontrivial relative permutations source automorphisms: the checked source automorphism groups are trivial. They are not assembled into a permutation of one globally shared host.

Let z mark changed slots and t mark changed owner clocks. Before the global current and target-margin predicates, the complete binary menu has weighted enumerator

    (1 + z t^2)^17 (1 + z t^3).

After the full target-margin and dyadic-current guards, it becomes exactly

    1 + z^18 t^37.

check_cycles.py verifies all eighteen literal relative permutations, the rejected third ambient iterate, every coefficient of the unguarded enumerator, and the full guarded enumerator by direct assignment enumeration. Its normal and optimized receipts agree.

### Complement reciprocity and the odd-clock distinction

At each dyadic precision, complementing every flip bit takes a legal state to another legal state. The full circuit difference vanishes at every such precision, and target margins are balanced. If P_h(z,t) records the legal states relative to the old endpoint, then

    P_h(z,t) = z^18 t^37 P_h(z^-1,t^-1).

The exact accepted-state counts at moduli 2,4,8,16,32,64,128,256 are respectively 48620,1526,754,350,152,72,34,2. The complete coefficient lists and reciprocity checks are saved in the cycle receipts.

Since the total clock cost 37 is odd, complement exchanges even and odd clock cost. Exactly half the legal states have each parity at every precision. The only odd-cost slot is source 2835, so clock-cost parity is exactly its flip bit. For example the full mod-8 bank has 377 states of each clock parity and exactly 377 with source 2835 flipped.

Equivalently, the univariate clock enumerator P_h(1,t) is divisible by 1+t over the integers. This is a consequence of a literal complement involution and odd total cost. In an all-transposition menu, all clock costs are even instead, so this parity-exchanging conclusion does not apply. No modular-form or broader representation-theoretic identification is being asserted.

## 7. Reproducibility and boundaries

Run check.py against the original repository exports, check_boundary.py against the frozen supplied inventory, and check_cycles.py against the original repository exports. All use explicit failure checks, and normal versus python -O receipts are byte-identical. SOURCE_PINS.json authenticates the original graph and inherited source inventory. MANIFEST.json records the deliverable hashes. Full integer-image witnesses are retained independently of the solver.

This result separates three genuinely different obstructions: finite-group cancellation can make eighteen changes atomic; native permutation costs control how expensive those changes are in actual owner clocks; and normalized resource capacity decides whether the path alternatives can be used without repeating a port. Algebra alone does not supply additional native options or a compatible global remainder.
