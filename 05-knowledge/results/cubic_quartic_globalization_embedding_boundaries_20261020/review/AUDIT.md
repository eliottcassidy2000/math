# Independent audit: cubic/quartic globalization and embedding boundaries

20 October 2026.

## Verdict

**PASS within the stated hypotheses.** No substantive mathematical defect or
scope contraction was found. The relative hostile, primitive-cubic lattice,
whole-algebra globalization criterion, and inequivalent embeddings are valid.
The optional characteristic-sequence enumeration is also justified by the
primary text, including its ordering convention.

The accepted consequence is a failure-boundary packet and a conditional
globalization result. It does not exclude the actual Keller conductor-six
sector, raise a universal bound to eight, or prove the planar Jacobian
conjecture. The inherited conductor-six alternative and moving quadratic
exclusion were supplied as already audited premises; they were not re-proved
in this review.

All reviewer writes are confined to this review directory. No producer file,
shared checkout, or remote repository was modified by the reviewer.

## 1. Saturation and the entire conductor

The monic hypothesis is essential and is used correctly. For
`D=C[x,f]`, the normalization candidate `B=C[x,y]` is finite free of rank
`n=deg_y f`, with the indicated powers of `y` as basis. Slice primitivity,
together with fixed degree `n`, makes the specialized characteristic
polynomial of multiplication by `g` its minimal polynomial of degree `n`.
Therefore the specialized hypersurface presentation is exactly the
polynomial image, rather than a quotient with an unobserved nilpotent kernel.

Generic primitivity follows from the hypotheses: the common field cannot
have degree greater than one while a slice has full minimal degree `n`.
Equivalently, a generic relation of smaller `g`-degree would specialize to
such a relation away from its finitely many bad coefficients. The standard
integrality/normality argument places the monic minimal polynomial in
`C[x,u][v]`. Both power-basis modules in the proof are consequently free
over the same two-dimensional base.

The saturation deduction is valid one linear factor at a time. If
`(x-a)b` is in `T`, injectivity of specialization puts it in `(x-a)T`;
cancellation occurs inside the ambient domain. Factoring a nonzero
polynomial in `C[x]` proves the asserted full saturation.

The trace-dual proof does not require a PID. For a monic polynomial `M` of
degree `n`, the trace of `q(theta)/M'(theta)` is the coefficient of degree
`n-1` in the remainder of `q`. Hence the power-basis pairing is
anti-triangular with unit anti-diagonal over the base itself. These perfect
pairings give the two dual formulas and double duality.

In particular, the polynomiality step is not circular. The inclusion
`B† ⊂ T†` implies `H_v/f_y` belongs to `T`, by applying it to `1/f_y`.
Writing that polynomial as `c`, the equivalences are

    bB ⊂ T  iff  bT† ⊂ B†  iff  (b/c)T ⊂ B  iff  b/c ∈ B.

The last step uses both `1∈T` and `T⊂B`. Thus the full conductor is `cB`,
not merely a contained principal ideal or a divisor with the right degree.
No Smith normal form over `C[x,u]` is needed or used.

## 2. The relative hostile survives specialization

The pair `P=y^4`, `Q=xy^3+y^2+y` is immersed on every slice, including
`x=0`. Coprime degrees prove primitivity away from zero; the displayed
rational reconstruction proves it at zero. These facts satisfy every
hypothesis of the saturation lemma.

For each nonzero `a`, the standard power basis over the cubic coordinate
has leading residues `0,4,8` modulo three. It gives exactly the degree
semigroup `<3,4>`, not just an inclusion of that numerical semigroup.
Distinct leading residues prohibit hidden leading cancellation. Codimension
is consequently three.

The revised direct collision is exact for every nonzero parameter: choose
`r^2=-1/a`. Then `r` and `-r` are distinct and both map to
`(1/a^2,-1/a)`. At zero the translation in the note gives
`C[t^2,t(t^2+1/4)]`, whose two conductor roots produce a genuine node.
Thus no injective slice is hidden in the globalization counterexample.

The relation and derivative identity are correct. Applying the same
trace-dual formula to each specialized pair is legitimate, since its
characteristic polynomial remains minimal there. Accordingly the global
conductor specializes to the actual fiber conductor in this example.
The degree falls from six to two. This is not an appeal to a general
unproved theorem saying arbitrary conductors commute with base change.

The nonexistence of a monic cubic is intrinsic to the entire algebra:
the generic degree-at-most-three space is the span of `1` and `Q/x`.
A constant term cannot remove either coefficient `1/x` of `y^2` and `y`.
The displayed integral cubic filtration is correct. The affine-source
change argument uses only a constant nonzero linear coefficient, precisely
the automorphisms of `C[x,y]` that fix `x`.

The Jacobian `-4y^6` explicitly prevents interpreting this family as a Keller
counterexample.

## 3. Primitive lattices and whole-algebra generation

The cubic lattice proof is valid for an arbitrary saturated `R`-subalgebra
in the stated generic sector; finite generation is not silently assumed.
Clearing a denominator supplies a cubic in `S`. Saturation allows division
by the gcd of its nonconstant coefficients. Evaluation at `y=0` separates
the constant coefficient, and the Bezout identity for the primitive
coefficient triple forces its scalar multiplier into `R`.

Therefore every zero of the primitive leading coefficient produces a
nonconstant element of degree at most two, rather than the zero
specialization. If linear elements are forbidden, it is exactly quadratic.
The cited one-line theorem has the correct role in the Keller corollary:
a linear element of the polynomial-image slice algebra makes the original
map injective on that source line. It does not need a selected generator
itself to be linear.

The quartic construction in Section 3.2 also preserves membership in `S`.
After the cubic is monic, subtraction removes the cubic coefficient
without leaving the generic algebra; denominator clearing and saturation
then give a primitive quartic coefficient triple. A zero of its leading
coefficient would again give a prohibited degree-at-most-two element.

The final equality is equality of entire algebras. The monic cubic/quartic
subalgebra is saturated by the all-fiber primitivity lemma. Its generic
semigroup has the same three gaps as `S_K`; inclusion and equal finite
codimension force generic equality. Clearing a base denominator followed
by saturation then gives equality over `R`. The subsequent source
translation and generator subtractions have polynomial base coefficients.

The universal quantifier is correctly over every complex fiber, not just
the generic fiber or a finite test set. The sufficient condition is not
assumed for all actual Keller sources.

Optional strengthening, not a correction: an entire presentation by monic
degrees three and four itself forces every specialized image to have
semigroup `<3,4>`. Thus the no-low-degree-fiber condition is also necessary
for that exact kind of entire presentation. The note's decision not to
claim necessity is a harmless underclaim.

## 4. Prescribed lifts versus equivalence of embedded images

The generation identity `z^2-64k^3-9=8h` holds exactly and proves equality
of the coordinate rings. The conductor of the standard presentation is
`(y^3+1)^2`; its support is the three distinct preimages of the single
ordinary triple point. Translating its implicit equation by `U=-1`
gives tangent cone `V^3+W^3`, with three distinct tangent lines.

For a prescribed ambient lift, each coordinate differs from the explicit
map `T` by a multiple of the irreducible equation. At the singular point
both the equation and its gradient vanish, so the local derivative is
forced and its determinant is four. This uses the entire prescribed
coordinate functions, not a claim that the selected pair's determinant
governs all possible generating pairs.

The period identity has the correct sign. Differentiating `U dV` gives
`dU∧dV`. Therefore a determinant-`lambda` automorphism makes
`Phi*(U dV)-lambda U dV` closed, hence polynomial-exact. Pullback of its
primitive has equal endpoint values at a collision. The computed periods
force `lambda=-28/5`, contradicting four.

The abstract-automorphism argument genuinely closes the gap from
prescribed lifts to ambient equivalence of images. Any algebra
automorphism extends to an affine automorphism of the normalization and
preserves the conductor ideal. Preservation of its three-point support
forces translation zero and scaling `alpha^3=1`. Every such scaling
has the stated ambient linear lift in the standard presentation.
Precomposition therefore reduces any unparametrized equivalence to the
already impossible prescribed lift.

As an independent check, composing the second parametrization with an
arbitrary such scaling forces local determinant `4alpha`, while the
period ratio is `(-28/5)alpha`; they still disagree for every nonzero
`alpha`.

One type-clarity suggestion was accepted and verified in the final producer
wording: the generating pairs define closed embeddings of `Spec A`; the
displayed maps from the normalization parameter `y` are not themselves
closed embeddings of `A^1`, since three points collide. This resolved
exposition clarification was not a flaw in the argument.

## 5. Primary-source and characteristic-sequence audit

The downloaded Sathaye–Stenerson text, Section 3.1, printed pages 8–9,
explicitly provides:

- The three characteristic-sequence axioms.
- The conductor formula.
- For a nonprincipal sequence, initial entries at most `c+1` and later
  entries at most `c-1`.
- Every primitive semigroup generator occurs in the sequence.

For `<3,4>`, the conductor is six, the initial entries belong to
`{3,4,6,7}`, and later positive entries belong to `{3,4}`. Both primitive
generators three and four must appear.

The ordering convention is justified, rather than assumed from a
different convention in the source. Swapping the first two entries
preserves their gcd and lcm. The first relation holds in either order;
the first later strict inequality is the same bound by that lcm.
All later initial semigroups and gcds remain unchanged. The initial
conductor contribution `1-delta0-delta1+lcm(delta0,delta1)` is symmetric.
Thus either order gives a characteristic sequence, and the larger
initial entry may be placed first.

After ordering, the nondividing initial pairs are exactly
`(4,3),(6,4),(7,3),(7,4),(7,6)`. Each pair involving seven already has gcd
one, so it cannot be followed by another entry under the strict
gcd-drop axiom and cannot acquire the missing primitive generators.
The pair `(4,3)` is complete. The pair `(6,4)` forces the later entry
three and then ends. Hence the two patterns in the note are exhaustive
as a numerical statement about characteristic sequences.

This does not classify all embeddings or transport a parameterized
classification over `C(x)`. The distinction is maintained in the note.

The downloaded Gurycz original, Section 7, explicitly distinguishes
equivalence of images from extension of a specified isomorphism.
The attribution and its use to motivate Section 4.3 are accurate.

A fresh web search returned indexed original Abhyankar–Sathaye page 1062
at the stated mirror. It confirms the distinct-prime hypothesis of
Theorem 1.3, the composite coprime exponent hypothesis of Theorem 1.4,
and the reduction to nonprincipal initial degree pairs. Direct mirror
fetch returned HTTP 429 and the AMS route was unavailable; no fresh full
PDF download is claimed. None of that paper's nonuniqueness constructions
is assumed rational or immersed. The new witness is proved directly.

The Gwoździewicz one-line theorem is inherited from the cited completed
low-delta source audit. It is used only within its planar complex
constant-nonzero-Jacobian hypotheses.

## 6. Independent finite controls and reproducibility

The final producer checker passed normally and with `-O`, with
**74 checks** and identical output. Both replay files equal the frozen
producer receipt. After the accepted final wording clarification and
review-status promotion, the refreshed seven-file producer hash manifest
was verified in full. The review summary pins the final producer hashes.
No mathematical or publication-bookkeeping issue remains in this audit.

The independent checker imports no producer code and uses explicit
exceptions. It passes **56 checks** normally and with `-O`, with identical
output. Besides separate elimination and period calculations, it checks
the polynomial trace-dual Gram determinant over the two-dimensional base,
conductor containment by a separate linear solve, direct all-nonzero-
parameter collisions, and recursive gcd-descent enumeration of the
characteristic patterns.

For seven rational parameter fibers, it also computes the actual finite
algebra modulo the proposed conductor and then the largest ideal of the
normalization quotient contained in that algebra. The latter ideal is
zero in each case, independently testing conductor maximality rather than
only its containment or degree. It performs the same maximality test
for both presentations of the triple-point algebra.

These are finite exact controls. The all-fiber, saturation, and universal
generation statements are established by the proofs above, not by the
sampled rational parameters.

From the packet directory:

    python3 review/independent_checks.py > /tmp/cq-independent.json
    python3 -O review/independent_checks.py > /tmp/cq-independent-O.json
    cmp /tmp/cq-independent.json /tmp/cq-independent-O.json
    python3 verify_cubic_quartic_boundaries.py
    python3 -O verify_cubic_quartic_boundaries.py
    sha256sum -c SHA256SUMS.txt

The review summary records the reviewed producer hashes and exact counts.
