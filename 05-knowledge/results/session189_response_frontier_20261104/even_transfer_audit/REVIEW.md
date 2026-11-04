# Independent review: all-even-child modulo-sixteen ambient transfer

**APPROVED for the stated all-order transfer theorem, not for universal local nonvanishing.** The proof has been checked independently and its complete retained kernel and physical data have been recounted by independently written algorithms. The frozen producer proof is SHA-256 `05921da21eebcfc97cdf462fd247615cc0fe2cb45899ab9dab868589f3f1fde3`.

## 1. Exact premise and conclusion

Let Q be any finite tournament on at least one child label. Every actual source child A_i has even order at least four, and its corresponding target B_i has exactly the same complete ordered-cover profile on the same physical child set. Different child positions need not have the same orders or profiles. The target quotient is Q converse. The actual ports satisfy x to y and are individually uniform on each child, with no independent relabeling of their bits or the physical children.

For every such query, the augmented core contrasts obey Delta q_2 = 2 Delta H and Delta q_h = 0 for all h>=3 modulo sixteen. Hence a nonzero local Hamilton residue survives every fixed ambient root substitution when both ports inherit precisely that root's exterior incidences. The multiplier is c_1+2c_2=H(R[TT2]), an odd integer. An odd factor preserves nonvanishing but need not preserve the actual residue.

This is not a claim that a nonzero local residue always exists. It does not extend the domain to arbitrary module pairs, odd or smaller children, nonuniform internal ports, unmatched child profiles, arbitrary exterior changes, or all-pairs orbit reconstruction. The proof correctly distinguishes this theorem from a false generic mod-sixteen transfer principle.

## 2. Independent all-order proof review

The ordered path-slot action is handled with its actual stabilizers. Because all labels except the distinguished i occur once, repeated identical path words can only be singleton-i words. If d paths have that form, the stabilizer is exactly S_d, not a free action. At least one other path exists. The orbit size h!/d! is even whenever at least two non-singleton-i slots exist. The only potentially odd higher-cover orbit therefore has h=d+1 odd and one remaining path.

That remaining path contains every other label and r=k-h+1 copies of i. For r=0 or 1 its count is an odd tournament Hamilton count, including the one-vertex complement case. For r>=3, the r-1 intermediate return segments have distinct nonempty disjoint supports. Their permutation is free and preserves the joins, so their count is even. The exceptional r=2 occurs at h=k-1; for k in {2,3,5}, this is either h=1 or even. Thus the higher-cover parity table is exactly the three asserted entries: (h,k)=(3,2),(3,3),(5,5). The scalar k=3,5 parities are separately even. No return-segment permutation is incorrectly applied through arbitrary cuts.

The inherited even-child normalized domain has geometric baseline (1,14,4,8,0,...) and independent formal perturbations 8e at q_2, 8f at q_3, and 8g at q_5. Two perturbations vanish modulo sixteen; one perturbation forces every other child to occur only once. This produces exactly the displayed three derivative kernels. Their higher-cover parities are independent of every arc, so they cancel between the two quotient orientations. Their scalar k=3,5 parities vanish too. The only possible correction is 8 eta in Delta H, which disappears after multiplication by two. The common odd product of child Hamilton counts is preserved throughout.

The general geometric full-profile identity is supplied by the separately reviewed Seidel theorem, not inferred from finite quotient lists. The proof controls the geometric tail by profile valuations. Independently, every actual ordered q_h is divisible by h!, so h>=6 is zero modulo sixteen. Finally the complete profile identity, rather than scalar nonvanishing alone, justifies the ambient substitution. Substituting the actual transitive pair proves the multiplier's oddness. These steps establish every higher coordinate and every ambient without a quotient-order cutoff.

## 3. Independent exact kernel recount

`kernel_audit.py` uses an exact suffix recurrence whose state consists of the unused unique labels, remaining copies of i, and the preceding label. Its cut polynomial distinguishes an arc continuation from an unrestricted cut. Coefficients are packed into noninteracting integer fields, with an explicit factorial bound preventing carries. This is independent of the producer's forward mixed-radix state recurrence.

It enumerates every labeled tournament of orders two through five, every distinguished label, and k=2,3,5. There are 16,212 complete kernel records, 117,940 possible-length coordinates, or 145,908 coordinates after explicit zero padding through h=9. Every exact coordinate agrees with the producer. Complete independent data, not just parity histograms, are retained in `kernel_recount.json`. Higher h outside the word length are identically zero by nonempty slots.

## 4. Independent physical and named checks

`producer_audit.py` imports no producer implementation. It reconstructs every physical child, quotient arc, exterior arc, port name and mask. It independently checks all 14 distinct complete actual child profiles. The separate `literal.cpp` uses backward endpoint subset recurrence, adding either an allowed predecessor arc or a new path cut. All 944 producer physical sides and all 4,720 exact q_1 through q_5 coordinates agree, including all 400 local queries, eight cyclic ambient controls and the complete 64-query outside-domain bank. The 18! binomial(17,4) bound is checked before fixed-width exact arithmetic.

The independent `physical_audit.py` additionally computes every ordered-cover coordinate, not merely the first five, on 40 local actual queries and three cyclic ambient lifts. It uses exact packed cut polynomials on literal vertices, including nontransitive children, unequal orders four and six, and profile-preserving converses. It verifies 1,012 complete exact coordinates through augmented order sixteen. Twenty-two of the local residues are nonzero; all three selected ambient controls preserve a nonzero residue with the correct odd multiplier. All actual row lists, child sets, port/exterior labels and masks are retained.

`named_audit.py` retains the complete Boolean Mobius banks for all twenty full block-query cubes present in the producer controls: nineteen all-even two-child families and the complete three-child hostile family. These comprise 3,680 coefficient coordinates over 200 ordered-cover channels, including all 960 cubic and 340 quartic positions and every zero. Every cube value is reconstructed from those coefficients, and the full bidegree bound is checked. Physical bit-expansion masks are retained. These are complete block-uniform observers; no nonuniform within-child bank is claimed.

The hostile directed-triangle example with child orders (2,2,3) has the asserted nonzero third contrast modulo sixteen. It is outside the theorem and confirms that the hypothesis has not been silently erased.

## 5. Reproduction, pins, and finite scope

Place the producer at sibling `even_transfer/` (the original local name `session189-even-transfer/` is also supported). Run each of `kernel_audit.py`, `physical_audit.py`, `producer_audit.py`, and `named_audit.py` under normal Python and Python -O, saving the corresponding receipts. Each normal/optimized pair must be byte-identical. The producer cross-audit requires the independent kernel recount generated by the first program. All existing retained data are compared exactly; checks use explicit exceptions and remain active under optimization. Only standard Python and the already installed C++17 compiler are used.

The publication allowlist, checksum manifest, and dependency pins freeze the sources and full evidence. These finite computations corroborate the independent all-order proof; they are not the basis for extrapolating beyond order five or beyond the displayed physical sizes. The inherited S188 profile-domain theorem and the reviewed Seidel geometric identity remain explicit dependencies.
