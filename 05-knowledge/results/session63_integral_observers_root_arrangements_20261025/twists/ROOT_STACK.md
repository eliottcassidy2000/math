# The quotient stack retains the simultaneous branch-root construction

**Status: PROVED from the local quotient presentations below.** The root-stack construction is classical, not newly invented here. Characteristic two is excluded. X and G are as in `RESULT.md`; no F1 model is asserted.

For each branch line D_i=(L_i=0), its square-root stack parametrizes a line bundle M_i with a section s_i and an isomorphism M_i²≅f*O(D_i) carrying s_i² to the pulled-back defining section. Let R be the fiber product over P² of these seven root stacks.

## Proposition

There is a natural equivalence of stacks over P²

    R ≅ [X/G],       G=(μ₂)^7/μ₂.

Its coarse space is P², its geometric stabilizer on exactly r branch lines is (μ₂)^r, and its restriction to the branch complement U is U itself. Thus it retains the additive branch equations and inertia that the coarse quotient forgets.

**Proof.** Cover P² by U_j={L_j≠0}. A root of a nowhere-zero section is canonically trivial as a root object: its root section trivializes its line bundle. On U_j trivialize O(1) by L_j. The six other roots are roots of the regular functions L_i/L_j, and the product root stack has presentation

    [ Spec_(U_j)(O_(U_j)[r_i:i≠j]/(r_i²−L_i/L_j)) / μ₂^6 ].

This affine scheme is π⁻¹(U_j)⊂X with r_i=x_i/x_j. The group is the projective sign group normalized at j. On overlaps the normalization changes r_i↦r_i/r_k, with the corresponding relabeling of signs. These are exactly the transitions inherited from both the root-line objects and the projective root coordinates, so the local equivalences glue. Stabilizers independently flip the zero root coordinates, proving the remaining assertions. QED.

At a triple point the atlas is the singular node w²=a u²+b v². Hence the stack is **not claimed smooth** merely because each divisor is smooth: three lines on a surface are not a simple normal crossing divisor. Roots of a line bundle WITHOUT a section are a different construction; retaining sections is essential for trivial stabilizer on U.

Over F_q its groupoid mass equals #P²(F_q), as proved in `RESULT.md` Theorem A. Over a rational point on r branch lines there are2^r local root forms, each with automorphism group of order2^r, giving mass1. Numerical equality does not identify the stack with its coarse space. On U, freeness does give [V/G]≅U. No unspecified stack/orbifold Euler convention is identified with the formal q=1 value: the named Euler values in `RESULT.md` are those of the coarse varieties.

This stack and the inherited Session59 additive blueprint are constructive carriers. Neither supplies an F1 descent, a preferred q=1 point functor, a Frobenius-lift system, or an RH theorem.

## Classical sources and scope

- [Stacks Project, Tag04V8](https://stacks.math.columbia.edu/tag/04V8) records roots of effective divisors and credits Cadman and Abramovich–Graber–Vistoli.
- Charles Cadman, *Using Stacks to Impose Tangency Conditions on Curves*, Definition2.2.4 and Remark2.2.5 define simultaneous roots and iteration; Examples2.4.1–2.4.2 give the affine quotient presentation and triviality for a nonvanishing section. [Verified author-hosted primary text](https://www.charlescadman.com/raw/pdf/stacks.pdf). These precise definitions, not the paper's stable-map applications, are used here.
- [Stacks Project, Tag04UV](https://stacks.math.columbia.edu/tag/04UV) gives the quotient-stack torsor presentation. The specific equivalence is proved above.
