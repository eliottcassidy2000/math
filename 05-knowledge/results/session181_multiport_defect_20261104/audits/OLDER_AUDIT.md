# Independent audit: root-module defect-slot packet

2026-11-04. **PASS after one explicit scope repair.**

Reviewed `../older/PROOF.md` and its independent-control implementation. The positive slot formula, inclusion-exclusion version, arbitrary-module substitution, transitive-clone exponential formula, characteristic recurrence, and strongness preservation/reflection are sound. In particular, an ordered path cover is assigned in the increasing order of the occupied slots, so no missing factorial or unintended extra slot permutation occurs.

The arbitrary-module target bound uses a triangle through each quotient vertex, Hamilton-path insertion monotonicity, and `H(C3[B,K1,K1])=H(B)+P_B(2)`. Its initial wording omitted the nontrivial quotient hypothesis. A singleton quotient is strong by a common convention and admits arbitrarily large transitive modules at H=1. The author repaired the statement to require quotient order q>=3; the current proof is valid with that condition.

Source attribution was checked against pinned main `9f5bdfaea3a2debeaac8d076c0367d5006aca71a`: THM-4145 supplies the zero/one/two-defect rooted-pair formula and strongness argument, THM-3121 supplies arbitrary substitution and its triangle specialization, and THM-4162 supplies the rooted-pair capacity tensor/coset filter. The local Session175 clone-obstruction proof explicitly already supplies the strong fixed-template polylogarithmic bound. The new packet accurately labels its specialization and branch-tree application and does not claim a new general substitution principle or global completeness.

`older_independent.py` imports no author or repository code. It uses backward/successor memoization for literal Hamilton counts and explicit permutation/cut enumeration for the factor path-cover data. It verifies:

- 1,136 rooted arbitrary-module substitutions, including every rooted labelled quotient of orders2 through4 and four structurally different modules.
- 64 triangle substitutions with three independently varying arbitrary modules, checking every module-size lower bound.
- 14 characteristic-recurrence equations using literal child counts rather than slot-formula values.

Normal and optimized execution pass. The author's full 21,300-case arbitrary-module control was separately replayed and matched its frozen receipt. These finite controls support the reviewed all-order proofs; they are not substituted for them. No native guard, all-depth two-digit closure, or global spectrum theorem follows from this packet.
