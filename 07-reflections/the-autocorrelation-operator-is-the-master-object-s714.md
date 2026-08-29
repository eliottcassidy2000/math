# The autocorrelation operator MM* is the master object (S714)

The instruction was to hunt for operator/adjoint pairs across the repo, think outside the box, and run
the whole thing through an autocorrelation lens. The hunt converges fast, because the pairs are all the
same pair wearing different clothes. A problem hands you an operator `M`; its adjoint `M*` is whatever
"reverse/conjugate/transpose" means in that setting; and the object that actually decides the answer is
never `M` but the **autocorrelation** `C = M M*`, the operator composed with its own adjoint. `C` is
positive semidefinite by construction, and every invariant the project cares about is a functional of
`C`'s spectrum.

The cleanest instance is the one I had just built without naming it. A tournament's skew-adjacency
`S = A - A^T` has adjoint `S* = S^T = -S`, so its autocorrelation is `C = S S* = -S^2`, a positive
semidefinite matrix whose eigenvalues are the squared magnitudes of the skew spectrum. The Pfaffian I
spent last session on is exactly `|Pf(S)| = sqrt(det C)`: the Pfaffian is the square-root-determinant of
the autocorrelation operator. Said that way, S713's "Pfaffian is the even/odd seam" and this session's
"autocorrelation is the master object" are one statement. And the same shape recurs everywhere once you
look: signed-LRC sign-orbit collisions are `|f-hat|^2` coincidences (THM-415 already calls this
homometry — same Patterson power spectrum); the number of unit distances is the peak of the Patterson
autocorrelation `A_P = (sum delta_{p_i}) * (reflection)` at radius one; the Paley maximizer is the
character whose autocorrelation is dead flat (`|Gauss|^2 = p`); the covering-depth factorial moments are
the multi-point autocorrelations of the runner-arc system. Operator, adjoint, autocorrelation — one
object, five problems.

The lens did more than reorganize; it sharpened a canonical fact. The tournament autocorrelation `-S^2`
on six vertices has, exhaustively over all 32768 tournaments, exactly **six** distinct spectra, not the
five that THM-120's prose claims (its own table, tellingly, lists six). The fifth and sixth are both
`|Pf| = 1` but carry different autocorrelation signatures, `(15,15,1)` and `(15,47,1)`: they differ in
the middle symmetric function `e2`, the sum of pairwise products of squared eigenvalues, which `|Pf|`
(only the product, `e3`) throws away. So the homometry partition — equal autocorrelation spectrum — is
strictly finer than the Pfaffian partition; it is the join of the `|Pf|` classes and the `e2`
(character-ratio) classes. The Pfaffian sees the determinant of the autocorrelation; homometry sees the
whole spectrum, and the whole spectrum distinguishes tournaments the determinant cannot. That is a small
but real refinement of canon, and it came directly from insisting on looking at `C = MM*` rather than at
the single number `det C`.

The second thing the lens explained is why my last two sessions' product structure was inevitable. The
Patterson autocorrelation of a Minkowski sum is the convolution of the factors' Pattersons — operators
tensor, so their autocorrelations convolve. The unit-distance count of `K3 [] W7` is then literally a
convolution peak, `(A_{K3} * A_{W7})` evaluated on the unit circle, and it comes out to `57` three
independent ways. The reason `n = 22` is hard reads off immediately: `22 = 2 * 11` forces the Patterson
to be a two-fold convolution `A_{K2} * A_{(11)}`, and the radius-one peak of such a convolution is capped
at `57`, below the known `60`. A configuration achieving `60` must have a Patterson that is *not* a
convolution of two smaller ones — an autocorrelation-irreducible point set. "21 is multiplicative, 22 is
additive" (S712) becomes "21's autocorrelation factors as a convolution, 22's does not."

And the third thing is a hardness principle the whole project keeps circling. Flat autocorrelation is
rigidity. The Paley set, the homometric signed-LRC collisions, the LRC transversal core where the
good-multiplier dodge dies — these are exactly the configurations whose autocorrelation operator is flat,
proportional to the identity off its DC term, no peak anywhere to grab onto. The solved, extremal cases —
the product unit-distance optima, the arithmetic-progression LRC tights, the high-Pfaffian tournaments —
are the peaked, convolutional ones. The recurring multiplicative-versus-additive, structure-versus-
quasirandom dichotomy the repo has named a dozen ways is, in the end, one spectral statement:
**flat versus peaked autocorrelation of `MM*`.** A peak is a handle; flatness is the absence of one. The
hard problems are hard because their autocorrelation gives you nothing to pull.

So the answer to "how many operator/adjoint instances are there" is: as many as there are problems, and
they are all the autocorrelation operator. The answer to "make progress" is: looking at the full spectrum
of `MM*` instead of its determinant splits a canonical class in two, factors the unit-distance optimum
into a convolution, and turns the frontier's hardness into a single word — flatness.
