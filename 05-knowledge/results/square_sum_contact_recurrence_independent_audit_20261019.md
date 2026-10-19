# Independent adversarial audit: logarithmic elliptic contact bound

**Verdict: PASS for the stated signed, even-index upper bounds.** The global
claim holds for every admissible primitive anchor and for either mixed-sign
port. No numerical integration, equidistribution rate, or unproved
Diophantine approximation estimate is needed. This audit does not claim a
matching lower bound, a positive-index-only bound, or an output-height bound.

Audited the recurrence-bound note and its two mathematical predecessors,
then inspected the companion checker. This audit made no changes outside
the recurrence staging directory. One minor exact-arithmetic qualification
was reported and incorporated: a local width in the rational-budget formula
must itself be rational, or be replaced by a positive rational lower bound.

## Proof-critical checks

1. **Both components and the period.** The cubic has three distinct real
   roots, and `P=(0,k)` lies on its bounded, nonidentity oval because
   `-1<0<min(K,1)`. Hence `2P` acts within each component. The inherited
   non-torsion proof applies to both port choices, so `2P` gives an irrational
   circle rotation. With the specified differential `dX/(2Y)`, a full turn
   traverses both signs of `Y`; its length is exactly
   `Omega=2 K(sqrt(1-f0/D))/sqrt(D)`. The potentially dangerous factor of two
   is correct. Independently, the bounded oval has the same period: the
   substitution `X=-1+(e2+1) sin(theta)^2` gives the identical complete
   integral. Translation therefore uses this period even when the anchor
   is on the other component.

2. **Strict circular recurrence, including wrap-around.** For
   `N=floor(Omega/ell)`, the `N+1` distinct prefix points yield a gap at most
   `Omega/(N+1)<ell`. Its oriented endpoint difference has the form `jA`
   with `0<|j|<=N`, including the last-to-first gap. Reversing `j` reaches
   the other orientation. Translating that short displacement by `Q` enters
   the open arc attached to `Q`, never an endpoint. The proof genuinely uses
   signed indices; it says nothing comparable about positive indices only
   or an arbitrary remote target arc. Passing from `A=2P` to `P` correctly
   multiplies the index budget by two.

3. **Differential, charts, and strict distinctness.** Direct differentiation
   of the birational map gives `dX/(2Y)=-du/z`. In fact the numerator of
   `z dX/du+2Y` is a multiple of the quartic equation. Setting `X=0` in
   the forward map and quartic gives
   `4u^2(2k^2 u+1)/k^2=0`, so its nonzero finite exception has `|u|<1`.
   Setting `X=K` gives `-4u^2(k^4-1)/k^2=0`, impossible here.
   The proposed open arcs have finite `|u|=t>1`, positive `z`, and no
   missing inverse chart. Their images are simple arcs because the inverse
   recovers `u`. They cannot silently cross infinity. The derivative of
   `f` vanishes only at `1+sqrt(2)`, where `f=1`, outside admissibility.
   Thus each interior point changes spread strictly. Correct primitive
   normalization, including division by two for odd/odd parameters, then
   gives a genuinely new unordered primitive petal.

4. **Arc lengths and the global constant.** The two angular intervals both
   have width `pi/24`. The radicand lies between `1/2` and `3/2`; the longer
   anchor-to-boundary arc therefore has length at least
   `pi/(24 sqrt(6))`. The spread-coordinate differential and the formulas
   for `L_-` and `L_+` are correct on both swapped-leg branches and both
   ports. Their combined length is at most `sqrt(2) pi/24`, safely less
   than `Omega/2`. Choosing the longer arc proves the global bound; asking
   for a particular direction correctly uses that direction's own length.
   A finite forbidden set only shortens an anchor-attached arc to its nearest
   forbidden spread. It does not invalidate this argument.

5. **Elementary period estimate and rational stopping.** Splitting the
   tangent-form complete integral at `1` and `1/r` gives the stated two
   `asinh(1)` tails and `log(1/r)` middle. All constants in formulas
   (14), (15), and (18) agree with this calculation. Since
   `log(D/f0)<=J log(2)<J` and `sqrt(D)>7/5`, the rational quantity
   `W=(5/7)(4+J)` strictly exceeds the period. The spread density is strictly
   greater than `1/5`, so a rational width `Delta` gives `ell>Delta/5`.
   Floors preserve the desired upper bound even when the coarse rational
   ratio is an integer. The global width `1/4` is available on one side
   of every anchor and gives exactly
   `2 floor(100(4+J)/7)`. Finite exact enumeration of both signs of all even
   indices to that budget is consequently certified to succeed.

6. **Cusp terms and scope.** With complementary modulus
   `sqrt(f0/D)`, the complete-integral expansion gives
   `Omega=log(32/f0)/sqrt(2)+O(f0 log(1/f0))` for either port. Both explicit
   arc constants are correct. Substituting `v=sqrt((1+s)/2)` for the plus
   port and `v=sqrt((1-s)/2)` for the minus port gives the stated `artanh`
   expressions directly. Uniform smoothness of the integrand on
   `0<=s<=1/2` near `f0=0` justifies `L_+=C_epsilon+O(f0)` and
   `L_-=f0/4+O(f0^2)`. The coefficients `sqrt(2)/C_epsilon` and
   `4sqrt(2)` in the two directional budgets are correct. These describe
   upper budgets, not first-success asymptotics. The predecessor's cusp
   divergence theorem is specifically for its plus-port search; the new
   upper theorem independently covers both ports.

## Independent controls and checker review

The separate audit control script retains **49,622 explicit checks** under
both ordinary and optimized Python. Outputs are byte-identical. They include
symbolic differential/chart checks; **7,624 exact rational-circle cases**,
with **251 wrap-gap cases**; and separately diagnostic 70-digit period,
oval, arc, bound, and cusp evaluations. Rational circle tests use prefixes
shorter than the rotation's prime order, which is enough to test the finite
distinct-prefix combinatorics. They do not pretend to certify irrationality.

Independent limiting constants:

- `C_+=0.15400230965624083193415914881489055...`
- `C_-=0.11740357026709195070581580989334068...`
- `sqrt(2)/C_+=9.18306722496473331048251963847045...`
- `sqrt(2)/C_-=12.04574579082879087767822589573117...`

The main companion checker correctly separates `exact_symbolic`,
`exact_rational`, and `numerical` checks. Exact witness equations and the
global stopping budget use rational/integer arithmetic. Every mpmath check
is classified as numerical. Its chosen 32-even-step witness cap is explicitly
diagnostic, and selected successful witnesses are not offered as a proof of
the universal theorem. Numerical floors are not used as certified search
stopping criteria. No additional implementation blocker was found.

Sources checked directly: [NIST DLMF 19.12.1 and 19.12.3](https://dlmf.nist.gov/19.12)
for the complete-integral expansion, and
[Elkies's real elliptic-group discussion](https://people.math.harvard.edu/~elkies/M223.24/index.html)
for the circle/two-component structure. The new quantitative recurrence
argument is proved directly rather than imported from those sources.
