# A finite-state construction target and its limits

Status: the arithmetic induction and ordinal-product branch are proved. No construction of the required state bank or two offset gadgets is asserted. This is a sufficient structural proof target, not a proof of the spectrum.

A decorated tournament has one of finitely many state labels. A proposed bank is an explicitly defined class of such objects. Suppose:

1. For each odd h from 23 to 67 there is a bank object T with H(T)=h.
2. There are three fixed, algorithmically specified bounded gadget operations U_0,U_1,U_2, taking any bank object T with H(T)>=23 to another bank object, and
   H(U_j(T))=3H(T)+2j, j=0,1,2.

Then every odd N>=23 is realized in the bank. For N>=69, take the unique j in {0,1,2} with N-2j divisible by 3. Then h=(N-2j)/3 is odd, >=23, and smaller than N; induct and apply U_j. The base is 23,25,...,67. Together with existing verified witnesses for odd values below 23 other than 7,21, and the established nonattainment of 7,21, this proves the entire spectrum. The ternary induction requires no overlap hypothesis. Once explicit gadgets and a local finite-state invariant are provided, their verification is an algebraic preservation problem; the target does not assume arbitrary spectrum coverage. The state family, witnesses, and offset constructions remain unresolved and may fail in the proposed form.

The count identity for U_0 is unconditionally available: U_0(T)=T→C_3 has 3H(T) paths by the ordinal-product formula. Its membership in any proposed bank is a separate issue. No assertion is made that either offset can be achieved by a single vertex or by a bounded gadget.

Minimal obstruction to naive ear intervals: exhaustive direct permutation counting for all labelled strong tournaments at n=3,4,5 finds every nonconstant one-vertex ear image solid at n=3,4, while order-five lex-upper-pair code 8 has H=9 and image {15,17,19,23,25,27,29,33,37,41}. Its holes are 21,31,35,39. The n=5 obstruction and image were already documented in THM-4118; the small-order minimality check adds context, not novelty. `check_minimal.py` prints the first nonsolid strong parent.

The idea that a fixed parent always has single-vertex ears at both offsets is false: for T=C_3→C_3, H=9 and the nonconstant ear image is {9,15,31,41,45,61,75,99}; 3H+2=29 is absent. For code 8 the value 3H+4=31 is absent. These examples lie below the bank threshold and therefore do not refute the finite-state target. They show why universal local assertions and unlabelled interval overlap cannot substitute for an invariant bank and valid gadgets.
