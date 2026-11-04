# Five-child exact classification and its first overlapping-module obstruction

This is a small producer census, not a nine-quotient or size-29 separation theorem. All graphs have five named old vertices. `five_child_records.json` contains one least binary edge code per isomorphism class and every simultaneous family of its original nonsingleton modules. The bit order is (01,02,03,04,12,13,14,23,24,34); a set bit means the increasing arc. Profile coordinate h counts ordered spanning covers by h nonempty directed paths. Isomorphism transports both source and selected family together.

## Complete classification

There are 12 isomorphism classes, of which 11 are nontransitive; their labeled multiplicities sum to 1024 and 904 respectively. Complete nontransitive source-family enumeration has 7,408 outputs after multiplicities. Full enumeration including the transitive source has 130,288 outputs, represented compactly by the representative families rather than duplicating labeled records.

| Code | Sorted degrees | Labeled | q1,q2,q3,q4,q5 | Possible target class codes |
|---:|---|---:|---|---|
|0|0,1,2,3,4|120|1,30,150,240,120|all 12|
|2|1,1,1,3,4|40|3,42,162,240,120|2,5,11,40|
|4|1,1,2,2,4|120|5,54,174,240,120|4,9|
|5|0,2,2,2,4|40|3,42,162,240,120|2,5,11,40|
|8|1,1,2,3,3|120|9,66,186,240,120|8,40|
|9|0,2,2,3,3|120|5,54,174,240,120|4,9|
|10|1,1,2,3,3|120|9,66,186,240,120|10|
|11|0,1,3,3,3|40|3,42,162,240,120|2,5,11,40|
|12|1,2,2,2,3|120|11,78,198,240,120|12|
|40|1,2,2,2,3|40|15,78,198,240,120|40|
|41|1,2,2,2,3|120|13,78,198,240,120|41|
|76|2,2,2,2,2|24|15,90,210,240,120|76|

The exact rows, all module masks, and family-by-family profiles are in the JSON. No classification is inferred from degree sequences alone: codes 8/10 and 12/40/41 demonstrate their ambiguity.

Seven of the eleven nontransitive classes retain their complete profile for every family. Codes 2,5,11 can increase from (3,42,162,240,120) to (15,78,198,240,120); code 8 can increase from (9,66,186,240,120) to (15,78,198,240,120). These are all changes, and every coordinate is nondecreasing.

## Explicit minimal obstruction

Take source code 8 with rows [16,1,3,7,14]. The set {1,2,3} is a transitive triple inside a cyclic substitution C3[TT3,1,1]. Both A={1,2} (mask6) and B={2,3} (mask12) are original modules. Simultaneously reverse their internal edges, producing code152, rows [16,5,9,3,14], isomorphic to code40, the cyclic substitution C3[C3,1,1]. The profile changes from (9,66,186,240,120) to (15,78,198,240,120). Neither operation alone changes the profile; the overlapping family does. This cannot occur in a nontransitive child of size at most four, as the preceding size28 classification establishes.

For code2, rows [4,1,2,7,15], toggling original modules masks15 and24 gives code693, also in class40, and changes (3,42,162,240,120) to (15,78,198,240,120). Codes5 and11 have the same transition, with example families [15,29] and [6,29] respectively.

These failures do not invalidate a collision-based reduction: the imported hereditary equality theorem forces corresponding complete child profiles to match under an assumed full-response collision. They do invalidate an unconditional claim that all five-child original families preserve profile. To obtain a family-uniform root battery from a quotient-only battery, use the separately established monotonicity plus paired converse-complement queries, not a false profile-invariance assertion.

## Modular profiles

Every q1 is odd. Divide by q1 in the indicated ring. Modulo eight there are only two normalized profiles:

- (1,6,6,0,0): codes 0,2,4,5,9,11,41,76.
- (1,2,2,0,0): codes 8,10,12,40.

Modulo sixteen there are four profiles:

- (1,14,6,0,8): codes 0,2,4,5,9,11.
- (1,2,10,0,8): codes 8,10,40.
- (1,10,2,0,8): code12.
- (1,6,14,0,8): codes41,76.

Fourth visits vanish modulo sixteen, but fifth visits have weight8 and must be retained in a mod16 word recurrence. This is an important difference from a three-visit truncation. Modulo eight both fourth and fifth visits vanish; the second type is exactly the C3 run profile, while the first is obtained by multiplying the C3 second- and third-visit weights by3. In profile-changing families q1 also changes, so normalized comparisons must retain separate physical odd prefactors.

For two four-children, the preceding complete profiles (3,18,36,24) and (5,22,36,24) both normalize to (1,6,4,0) modulo8. Modulo16 they normalize to (1,6,12,8) and (1,14,4,8), respectively; fourth visits do not vanish modulo16. Each is preserved by its own source-relative original families. These are local inputs only; no prime-nine separation follows from this local classification alone.

## Complete named polynomial controls and verification

`named_controls.json` retains all 1,024 Boolean values and all 1,024 multilinear coefficients for each canonical representative and its literal converse, 24 vectors and 24,576 coefficients total. Variable order is a0,...,a4,c0,...,c4, with actual ports x->y and a_v=[v->x], c_v=[v->y]. Even converse-isomorphic controls are retained on their literal physical names. This includes every zero coefficient and all cubic and quartic coordinates. There are 2,184 nonzero cubic and 2,160 nonzero quartic coordinates.

`prepare.py` computes profiles by ordered permutations and cut masks and responses by Hamilton subset dynamic programming plus Boolean Mobius inversion. `verify.py` independently recomputes profiles via ordered subset convolution of exact induced Hamilton counts, and obtains all named coefficients by directly enumerating seven-name permutations and expanding their signed Boolean edge factors. It checks all representative families, complete profiles, mod8/mod16 normalization, and all named values by zeta reconstruction. Run both `python3 verify.py` and `python3 -O verify.py`. The verifier uses explicit exceptions rather than assertions.
