# Every prime five-name quotient converse is visible in every fixed ambient

## 1. Physical setting and scope

Let Q be any prime five-name quotient and let its five named children be arbitrary nonempty tournaments, with arbitrary unequal orders and complete profiles. Compare T=Q[A_i] and U=Q^op[A_i], retaining the child names and profiles. Substitute T or U for the root of any fixed ambient tournament R.

The physical two ports x->y use the selected quotient block masks inside the core. Both ports imitate the original root at every exterior vertex: v->root gives v->x,v->y, and root->v gives x->v,y->v. The port-augmented core is therefore a common module. All exterior names and arcs remain fixed.

**Theorem.** Some child-uniform query with these exterior incidences distinguishes the actual full named responses. The proof below covers every child-type assignment and every total core parity. It combines separately proved low-modulus transfers with a new finite complete ordered-cover certificate. It does not infer an ambient theorem merely from the local modulo-sixteen classification.

## 2. A general modulo-four lift and the singleton branch

For any two tournaments of the same positive order N, the ordered-cover residues q2 and q3 modulo four depend only on N. Indeed

    qk/k! = S(N,k) modulo2,

because each unordered set partition is weighted by a product of odd induced Hamilton counts. For k=2,3 the factorial is two modulo four; for k>=4 it is divisible by four. Hence any same-order pair with Hamilton contrast two modulo four has complete ordered-cover contrast(2,0,0,...) modulo four. Substitution into a fixed rooted exterior multiplies its first coordinate by b1=H(R), an odd number, so the ambient contrast is still two modulo four.

The complete local type bank in `../prime_five_all_profiles/PROOF.md` gives a modulo-four baseline separator for all31 nonempty singleton-position masks at every prime representative. Other nontrivial child types cannot change these residues: q2/H is always two modulo four and each three-run coefficient is even. Thus every source with a singleton child is covered by this direct finite baseline certificate and the preceding lift. No laminarity assumption about a later arbitrary family is used here.

The representative Q41=(18,8,3,5,14) is also covered without singleton children: query(a,c)=(0,1) has normalized contrast two modulo four for every nontrivial type assignment. The actual odd normalization scalar preserves two, so it survives every ambient regardless of total core parity.

## 3. Even-total nontrivial children

It remains to consider Q10=(20,1,2,7,14) and Q76=(24,17,3,6,12), with all children nontrivial. The normalized types P,E are even-order, and O,C are odd-order. Their total parity is therefore determined exactly by the number of O/C types.

For even total core order, the local modulo-eight classification gives a nonzero gap except at Q76 with EEEEE. The port-augmented core is even, so the inherited even-order transfer in `../mod8_transfer/PROOF.md` preserves nonzero modulo-eight gaps in every fixed ambient.

The remaining EEEEE branch is covered by the separate full-cover modulo-sixteen theorem `../prime_five_even_ambient/PROOF.md`: one of its two fixed queries has complete augmented-profile contrast(8,0,0,...) modulo sixteen. This branch is not passed through a nonexistent general modulo-sixteen transfer principle.

## 4. Odd-total nontrivial children: complete finite cover certificate

There are exactly512 odd-total assignments from the four normalized nontrivial types P,E,O,C on five children. For each of Q10 and Q76, the present verifier computes the first three ordered-cover response coordinates for every block query in the256-coordinate Boolean query basis. It finds, for EVERY one of these512 assignments, a physical block query with normalized contrast

    (Delta q1,Delta q2,Delta q3)=(4,0,0) modulo8.         (1)

The complete cut-aware kernels, every type assignment and selected query, and all relevant response values are retained. The normalized type profiles are actual complete residues, not approximations. Every coordinate qk for k>=4 is divisible by eight. Multiplication by the common odd normalization factor preserves four and zero, so (1) gives the complete actual profile contrast(4,0,0,...) modulo eight. Common-module substitution then gives ambient Hamilton contrast4H(R)=4 modulo eight, for every fixed exterior.

The calculation is different from Hamilton-only affine reconstruction. A global ordered cover is a run word with cuts between its nonempty path slots. An uncut transition must follow an arc and cannot repeat a child symbol. A cut permits any transition, including equal child symbols. Ports appear once each, and each child's runs across all slots form one ordered cover. Thus coefficients retain every port order and placement, every empty-old path slot allowed by a singleton port, and all cubic/quartic interactions.

Modulo eight, only multiplicity1,2,3 is needed in a child and at most two children can have multiplicity above one. Each such profile factor is even. The51 possible multiplicity shapes are: all ones; one specified child at2 or3; or two specified children each at2 or3. All other products vanish modulo eight. The verifier retains both orientations' complete153-entry arrays (51 shapes times3 cover coordinates) at every basis query. This exact truncation does not apply a return-segment involution across path cuts.

## 5. Exhaustion and original-module families

Sections2–4 exhaust all nonempty child assignments and both augmented-core parities. The complete264-quotient classification and physical relabeling convention are inherited from the separately audited local packets. The theorem is a fixed-ambient quotient-converse result, not unrestricted named reconstruction.

A separate source-relative family consequence follows by the same monotonicity argument as in `../prime_five_even_ambient/PROOF.md`, Section4. Restrict operations to original modules inside this prime core and include the core itself after parity reduction. Every proper original module lies in a child. Hereditary original-module monotonicity makes each target child profile coordinatewise no smaller; no laminarity of descendants is assumed. Pair a certified baseline query with its core converse-complement, leaving the lifted exterior incidences fixed. Complete augmented-profile duality gives opposite baseline ambient contrasts, so choose the negative source-minus-baseline query. Its uniformity on each child and nonnegative run substitution prevent any target profile increase from cancelling it.

This argument does not apply the exposed-anchor theorem directly to crossing families. It explicitly uses baseline separation, prime-child containment, full-profile duality and hereditary monotonicity. It is source-relative detection in a fixed exterior, not pairwise injectivity among arbitrary family outputs. Operations involving exterior names are outside scope.

## 6. Reproduction

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. The source uses standard Python and the installed C++17 compiler, and every check remains active under optimization.

The complete odd-type certificate is complemented by direct physical subset counts on unequal canonical children inside a cyclic ambient. These counts verify the actual name masks, odd normalization factor and root-incidence lift. They do not substitute for full finite-domain coverage. The local classification and this ambient certificate require separate independent audits before being combined into a final claim.
