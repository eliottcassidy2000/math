# Independent curve / genus-two / Fano–Richelot audit

**Verdict: PASS for the geometric packet.** Audited 2026-10-25. The incorporated
clarifying guard is that a transverse **base** line pulls back to `2H`, not
`H`; it does not change any formula in the original packet. This separates
its genus-17 curve from the ordinary genus-5 hyperplane section of the
six-root K3. The audit's arithmetic sidecar below is explicitly additional
and does not enlarge the producer's original C-only claim.

## 1. Mathematical gate

- Independent branch valuations give a connected root curve of degree
  `2^(n-1)`. Riemann–Hurwitz yields genus `1+(n-4)2^(n-3)`, hence 17 and 49.
  The same result follows from adjunction with class `2H`:
  `2g-2=(2H).((n-6)H+2H)`.
- Characters really are **even subsets** of the n labels, without quotienting
  by complements. Their index-two quotient curves have branch set S and
  genus `|S|/2-1`. For n=6 this gives 15 elliptic and one genus-two factor;
  for n=7 it gives 35 elliptic and seven genus-two factors. The pullback
  differential decomposition proves an isogeny, not a polarized product
  isomorphism or simplicity of the factors.
- The complement quotient occurs instead in the different model
  `Jac(D)[2]=even subsets of six labels / complements`. The packet does not
  confuse these two spaces. Its pair-intersection form, matching subgroups,
  four-element cuff kernel, and seam reconstruction are exact.
- Residue characters on the root surfaces are the full six-element support
  at n=6 and the seven complements of single labels at n=7. Coordinate
  numerator signs cancel the odd ambient sign in the latter case. This is
  common character indexing, not restriction of a two-form to a one-form.
- The full-root pants lift is a four-sheeted sphere-with-six-holes over each
  base pair of pants. The lifted gluing graphs have `(V,E,b1)=(32,48,17)`
  and `(80,128,49)`. Parallel edges count and were retained.
- The marked Fano pairings give maximal isotropic two-torsion subgroups and
  hence degree-four (2,2)-isogenies. The Borel21 formulas are statements about
  labels. The explicit order-three Möbius failure prevents the invalid
  inference that this is a symmetry of the fixed geometric curve.
- All nine original magic line restrictions have distinct roots on the
  displayed rational slice. Infinity is correctly retained via homogeneous
  forms. Every one of the seven displayed Richelot determinants is nonzero,
  every target is a squarefree sextic, and the split-target hostile control
  correctly shows why a matching alone does not guarantee a genus-two target.

## 2. Reproduction and independent controls

The producer's checker was run normally and with `python3 -O`; both outputs
match the frozen `normal.json` exactly, with **392,767 active checks**.

`check_curve_independent.py` imports no producer code. It rebuilds the seven
splittings from homogeneous binary forms and computes the brackets as
half their homogeneous Poisson brackets. It verifies:

- all seven determinants by a separate cofactor formula;
- target squarefreeness via each bracket discriminant and every pairwise
  resultant, rather than only the source checker's derivative gcd;
- all seven target polynomials against the recorded artifact;
- the signed polynomial identity stated below;
- 144 normalized source/target Fp point-count comparisons, and 41 Fp²
  comparisons for p<=31. Its prime filter conservatively requires the entire
  seven-point configuration to have good reduction as well as the individual
  isogeny. This explains why other less restrictive filters can count more
  cases.

Normal and optimized independent runs pass with identical output. All audit
checks are explicit exceptions, not optimization-sensitive assertions.
The bounded checks supplement the cited theorem and proofs.

## 3. Arithmetic normalization sidecar — convention matters

The original geometric target `v²=product B_i` is correct over C. It is not
in general the correct Q-twist: our checks produce **45** finite-prime
counterexamples to treating that raw equation as the Q-isogenous target.

With coefficient rows ordered `(1,u,u²)` and brackets
`B_i=F_j'F_k-F_jF_k'`, the independent exact identity is

    sum_i F_i(x) B_i(z) = -Delta (x-z)^2.

Put `G_i=B_i/Delta`. The Bruin–Doerksen correspondence is defined over Q
when the individual factors are in Q[u]. Its normalized target is
`w²=product G_i`; equivalently, after `v=Delta² w`,

    v² = Delta * B_1 B_2 B_3.

The equal-leading-coefficient setup in that source causes no obstruction
for arbitrary nonzero leading coefficients `a_i`: put `A=a_1 a_2 a_3` and
replace `F_i` by `(A/a_i)F_i`. The source rescales by the square `A²`, and
the resulting normalized target rescales over Q to the displayed formula.
The homogeneous version or the exact correspondence identity handles a
linear affine factor corresponding to infinity. No new square-point
existence claim follows from this normalization.

**Citation-convention warning:** Castryck–Decru–Smith's printed determinant
uses columns in reverse order `(u²,u,1)`. Do not mix that determinant sign
with Bruin–Doerksen's displayed normalization. The original C-only packet
is unaffected; any Q implementation should keep one internally checked
convention throughout.

## 4. Primary-source checks

- [Auffarth–Lucchini Arteche–Rojas](https://arxiv.org/html/1905.12690v2),
  Theorem 2.1 and Remark 1.2, support the classical Humbert–Edge isogeny
  decomposition after shifting the type by one. The packet independently
  proves the even-character specialization it uses.
- [Castryck–Folsom–Hubrechts–Sutherland](https://math.mit.edu/~drew/proclms104.pdf),
  §6, Lemma 4 and the paragraph after Lemma 6, explicitly give pair classes
  and the Weil pairing. These support the six-label symplectic model.
- [Bruin–Doerksen](https://arxiv.org/html/0902.3480v2), §4, equations 4.2–4.3
  and Lemma 10, give the determinant gate, normalized bracket target,
  squarefreeness and base-field correspondence for individually rational
  factors. No geometric-to-arithmetic twist inference was silently made.
- [Castryck–Decru–Smith](https://nutmic2019.imj-prg.fr/confpapers/HashRichelot.pdf),
  §3.2, Definition 1 and Proposition 1, explicitly include the degree-five
  case and distinguish genus-two from split elliptic-product codomains.

No blocker was found in the packet's mathematical claims. Its limits on
geometric symmetries, polarized decompositions, and arithmetic square points
are all necessary and correctly maintained.

## 5. Final frozen source revision

The parent-directed narrow clarification was applied after the original audit.
Only the theorem prose and provenance/manifest changed; source code and both
392,767-check results are unchanged. No arithmetic scope was added to the source.

- Original `THEOREM.md`: `b81f159316edebb9e44ee0eb292db03fc87523c66fafe9450c0e7a3ff9255809`.
- Final `THEOREM.md`: `20030d4fadd02733b765866833c1bdf89f5c5e262b965cedc337c7d89054d271`.
- `check.py`: `184642280d799faff0c6c73f0ea2ba3ebbc452d15b8a78c5cceecfcfadbddcbe`.
- `normal.json` and `optimized.json`: `6fdb6febc394149edafb03b074053f54752288fafe8afdcc6f103337e001a687`.

The final paragraph was re-read against `f^*O(1)=O(2)` and adjunction; both
general-H sectional genera (5,17) and base-line-2H genera (17,49) agree.
**Final gate remains PASS.**
