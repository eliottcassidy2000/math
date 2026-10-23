# Independent equality-rigidity audit

Status: PASS; major errors0. Scope: EQUALITY_RIGIDITY.md and equality_checks.py/json only, building on the previously frozen homogeneous-completion audit. No initial audit files altered. No Git mutation.

## Detailed proof audit

1. The only extra dimension beyond247 lies in degree0, because nonzero grades occur in equally sized opposite pairs. Thus the claimed profile is forced, without a hidden assumption restricting support to five grades.
2. Under adjoint action of s, L0 contains the adjoint submodule s and has a one-dimensional quotient. Complete reducibility splits that quotient as a trivial module, yielding h with [s,h]=0. This establishes the Lie direct sum s⊕Ch; it does not assume reductivity merely from codimension1.
3. Invariance gives B(s,h)=B([s,s],h)=B(s,[s,h])=0. Homogeneity makes B|L0 nondegenerate, so both its adjoint block and one-dimensional h block are nondegenerate. Thus h is not isotropic. The Schur scalar action and all signs ±a,±2a follow from Jacobi and skew-adjointness.
4. If a=0, h is central on every grade. Its orthogonal complement is a graded ideal with a nondegenerate restricted form. It contains all of g: positive degrees are orthogonal to degree0 h and s is orthogonal to h. Its dimension247 contradicts the earlier audited theorem. This validates the crucial a≠0 step.
5. The radical's degree0 component is an s-submodule. Complete reducibility and the nonisomorphic adjoint/trivial factors force it to lie in Ch, since including s contradicts solvability. The s-projection of [V+,V−] is nonzero by invariance, using B(s,h)=0. Irreducibility therefore excludes R±1. Nonzero a excludes R0. The retained positive bracket paired with z− excludes R−2. Finally [z+,z−] centralizes s and hence is in Ch; its h pairing is2aB(z+,z−), nonzero. Thus R2=0. Every radical component has been eliminated.
6. For semisimple L, simple ideals are graded because the grading derivation is inner. Projections of s into any simple ideal lie in its degree0 part and, when nonzero, give an adjoint submodule. Two nonzero projections would give two independent adjoint copies in L0, impossible. Therefore s is contained wholly in one simple ideal I. All other ideals commute with s, so the nontrivial irreducible V± lie in I. Their positive bracket gives z+ in I. Distinct simple ideals are orthogonal for any invariant symmetric form, and the one-dimensional L−2 splits across graded ideals in only one summand, forcing z− into I. Then [z+,z−] gives h. Thus I=L.

These arguments also prove the stated general equality⇒simple assertion for any finite-dimensional complex simple s and nontrivial irreducible symplectic V. They do not prove general equality existence or marked embedding uniqueness.

## Identification and diagonal probe

Classical-family dimensions cannot equal248: type A would require249 to be a square; types B,C,D require discriminant1985 to be a square. Neither is. The five exceptional dimensions leave E8 uniquely. The classification is a cited standard input, not derived from numerical scans.

Verified primary author survey: Garibaldi, https://arxiv.org/pdf/1605.01721 , introduction and §2. The old maths.ed.ac.uk mirror failed to open; the accessible arXiv replacement was sent to the author. This is a link-quality issue, not a mathematical blocker.

The diagonal probe is correctly limited to fixed homogeneous metric and retained positive bracket. Opposite exponents sum to0 and w2=2w1, giving wj=jw1 on every present grade. Hence each transported bracket exponent is0. The executable's bounded integral search is a finite control, not proof of the abstract theorem; the linear equations themselves hold over the rationals or reals. Independent controls include rational exponents and analytic rank-discriminant exclusions.

## Verdict

The new equality rigidity result is supported: every248-dimensional completion satisfying the precise graded homogeneous metric hypotheses is simple E8. Consequently a genuinely nonisomorphic contraction must lose at least one hypothesis. No conclusion about all ungraded completions, marked embeddings, normalization of B, or real forms follows.
