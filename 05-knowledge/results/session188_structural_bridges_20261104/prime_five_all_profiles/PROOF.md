# Every prime five-name quotient converse has a child-uniform separator

## 1. Local theorem and fixed names

Let Q be any prime tournament on five quotient names, and let A_0,...,A_4 be arbitrary nonempty tournaments on fixed disjoint physical name sets. Then the full named two-port responses of Q[A_i] and Q^op[A_i] differ. A child-uniform physical query already separates them. Arbitrary target interior replacements preserving each complete source child profile give the same conclusion.

All child orders, orientations and unequal profiles are allowed. The original physical names remain in the graph. A quotient-mask bit means the chosen port has that incidence to every name of the corresponding actual child. An isomorphism used to select a quotient representative simultaneously relabels its child sets and query coordinates; it never independently forgets source and target names.

This packet proves local baseline quotient-converse separation. It does not by itself claim arbitrary-ambient transfer or noncancellation for arbitrary families of original modules. Those require separate arguments and audits.

## 2. Five normalized profile types and a complete small-run functional

The independently reviewed five-type theorem in `../five_profile_types/PROOF.md` gives the normalized profile q/H modulo8, with H=q1 odd:

    S=(0,1,0,0,...), P=(0,1,2,0,...),
    E=(0,1,6,4,...), O=(0,1,6,6,...), C=(0,1,2,2,...).

The types mean respectively singleton; pair; even order>=4; odd order>=3 with even Hamilton-cycle parity; and odd order>=3 with odd Hamilton-cycle parity. All coordinates from4 onward vanish. In every case normalization retains H as an odd scalar. The actual C3 has H=3; its normalized C type is not mistaken for a literal profile with H=1.

Write u_i=q2/H, let c_i be1 for O or C and0 otherwise, and let v_i be0 for S and1 otherwise. For a specified physical block query, define exact integer run-word coefficients:

- K0: every child occurs once, and each port once.
- K_i: child i occurs twice; every other child once.
- L_i: child i occurs three times; every other child once.
- K_ij: children i and j occur twice; every other child once.

Every run transition follows the quotient/port arcs and adjacent identical child symbols are forbidden. Both port orders are included whenever directed; no old-name endpoint condition is discarded. Then the COMPLETE normalized response is

    F_norm = K0 + sum_i u_i K_i + 2 sum_i c_i L_i
                     + 4 sum_(i<j) v_i v_j K_ij mod8.  (1)

The actual response is (product_i H(A_i))F_norm modulo8, with the same odd scalar for Q and Q^op.

Here is the complete truncation argument. Every repeated-child cover factor is even; three repeated children therefore contribute zero modulo8. A run multiplicity>=4 contributes zero directly. L_i is even by swapping its two nonempty return segments, whose other symbols occur once. If i occurs three times and only one other child j repeats (twice or three times), the same swap is still fixed-point-free: each return word contains at least two distinct non-i symbols, since a one-symbol return is a directed two-cycle and adjacent equal symbols are forbidden. Equal swapped return words would require at least two other repeated symbols, but only j can repeat. Thus every such coefficient is even and its even×even profile weight kills it modulo8. The only two-repeated-child term left is K_ij, whose weight is4v_i v_j. Finally the q3 term's extra4 in type O or E multiplies even L_i and vanishes. This proves(1) for each response, not merely for contrasts.

This is the same small-model reduction developed independently in `../coarse_mod8_affine/PROOF.md`, building on the inherited return-segment involution in `../parity_profile_transfer/PROOF.md`. Its exact physical coefficient extraction is

    K_i=(F_TT2(i)-F_single)/2,
    L_i=(F_TT3(i)-3F_TT2(i)+2F_single)/6,
    K_ij=(F_TT2(i),TT2(j)-F_TT2(i)-F_TT2(j)+F_single)/4.

All other children in these formulas are actual singletons. The L_i numerator is divisible by12. These are offline functional-reconstruction formulas, not a query count on an unknown original tournament.

## 3. Complete finite classification for five prime quotient names

The complete prime classification in `../prime_five_even_profiles/PROOF.md` has264 labeled prime quotients and three representatives, in outgoing-row notation:

    Q10=(20,1,2,7,14), orbit size120;
    Q41=(18,8,3,5,14), orbit size120;
    Q76=(24,17,3,6,12), orbit size24.

For each representative the verifier retains all1,024 block queries, both orientations, and all21 integer run channels modulo16. Their complete256-coordinate Boolean coefficient banks of bidegree at most(2,2) are also retained. Every full query is reconstructed from those banks, including every cubic/quartic term and zero.

Substituting every one of5^5=3,125 normalized type assignments into(1) gives the exact classification:

- Q10: no modulo-eight-blind assignment.
- Q41: no modulo-eight-blind assignment.
- Q76: exactly one blind assignment, EEEEE.

The receipts retain the complete type list, a separating query when one exists, the full-bank contrast hash and the number of nonzero queries. The full values are determined losslessly by the retained run channels and formula(1). No sampling, isomorphism assumption beyond the verified orbit enumeration, or floating arithmetic is used.

For Q76 with EEEEE, the separately reviewed eight-profile modulo-sixteen local certificate in `../prime_five_even_profiles/PROOF.md` applies to every actual child profile. Its fixed two-query battery cannot vanish simultaneously. This is precisely the remaining branch. Together these facts prove the local theorem for every nonempty child assignment.

A type S is an actual singleton, not an absent child. The finite classification includes all singleton-position patterns explicitly. The verifier also finds a modulo-four separator for all31 nonempty singleton-position masks at every representative, independently of the other nontrivial type choices. This is a statement about baseline quotient converse; it does not apply an exposed-anchor theorem directly to arbitrary crossing descendant families.

## 4. Scope and reproduction

The complete256-query Boolean basis is a fixed finite battery for the known five-child partition; the local EEEEE battery lies in it. If complete source profiles are known, the tables select a separating query, with the EEEEE branch using its extra modulo-sixteen bit information or its fixed two-query fallback.

The result is an unbounded-child-order consequence of a finite normalized profile domain. It is not a census of all large named tournaments, arbitrary prime quotient orders, exact response collisions, or unrestricted reconstruction. Modular blindness of the EEEEE branch at eight is not an exact collision and is resolved here at sixteen.

Direct physical subset counts on actual mixed canonical children check the physical-mask expansion and normalization scalar, including actual C3 interiors. They are checks on the implementation, not replacements for complete domain coverage. Every generated file and every check is replayable with standard Python and the installed C++ compiler, with checks retained under optimization.
