# Independent review: vertex-transitive extreme children

PASS. The frozen proof and verifier require no mathematical corrections. This review binds the final wording-corrected proof 23aaf2f3ff68bf1f57bed21528a19e1d517d18afb02e1e6c2973ad541618c5f7 and verifier cd7385d116621654233f4f523ef58ac2820fa46327efaffce13543c0b6ab5264.

## Analytical review

1. Regularity passes to every module because its exterior outdegree is constant. Genuinely overlapping modules have a uniform transitive three-part quotient, so the two nonempty parts of either module form an ordered join; their mean internal degrees differ by half the module order. Such a join cannot be regular. This establishes strongness without assuming vertex transitivity.
2. An automorphism taking one name of a strong module to another maps that module to an equal-sized intersecting strong module, hence to itself. Its restriction supplies the required internal automorphism. Singleton and whole modules are included.
3. The run decomposition for a balanced-module converse is exact. An interior physical mark occupies one prescribed run slot, whose slot permutation contributes no additional multiplicity. Exterior marks belong to the unchanged skeleton. This preserves complete one-marked E/S and unmarked q profiles, but does not imply full two-port invariance.
4. Original modules in the vertex-transitive child are laminar. Increasing-inclusion operations keep subsequent selected modules homogeneous. Earlier descendant operations preserve those modules' balanced marked profiles, so the induction is legitimate even if an intermediate tournament is not known to remain vertex-transitive. Ancestors perform whole converse and also preserve these balanced profiles.
5. Prime-root proper-module containment is used only inside the original substitution. Scalar equality supplies arbitrary middle-child q restoration, while Theorem4 supplies physical marked-profile restoration at the two extremes. The hub-selected quartic and hub-unselected uniform anchor therefore apply exactly. The scalar (0,0) query is already included, and m+20 is a valid safe union bound. No full arbitrary-child orbit injectivity or arbitrary strong outer lift follows or is claimed.

The author corrected one boundary sentence before the final freeze: endpoint balance does imply regularity and laminar modules, but it has not been shown here to imply endpoint balance of every submodule. The final proof retains vertex transitivity for that stronger step. No theorem or numerical data changed.

## Independent computational evidence

The audit imports no producer Python code. It independently checks all 27 regular tournaments among all 1,099 labelled tournaments through order5, their induced regularity and module laminarity, and the retained literal automorphism witnesses.

A distinct C++ enumerator goes through every vertex permutation and sums all valid cuts by binomial coefficients. It independently reconstructs 34 complete q/E/S profiles, including all16 original-module XOR families of C3[C3,C3,C3], and all eight arbitrary-ambient balanced-module controls. It checks256 marked-profile restorations in128 root families.

A separate predecessor subset DP checks all2,860 actual physical two-port counts in the retained root-family banks. All128 families are detected:16 scalar-strict,56 scalar-equal/hub-selected,56 scalar-equal/hub-unselected. Physical XOR rows, complete original module lists, selected branches, coefficients, residues and battery bounds are independently reconstructed.

Normal and optimized Python modes pass with identical receipts. The receipts contain source pins and numerical summaries, not a second copy of the complete producer data.
