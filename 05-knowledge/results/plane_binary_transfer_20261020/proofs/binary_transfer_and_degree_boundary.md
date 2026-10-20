# Binary transfer: robust cell rounding and the exact degree boundary

**Status:** PROVED below as standalone analytic/algebraic statements. Applying the transfer to an arbitrary coloring of an algebraically closed carrier is CONDITIONAL on the supplied blueprint's bounded projected-field and nonzero pair-correlation conclusions. Those source conclusions and the claimed Lean developments were not rebuilt or independently validated here.

## 1. The faithful object

For real `a,b,c ∈ [-1,1]` put

    Q(a,b,c) = 1 + ab + ac + bc.

Then

    2Q = (1+a)(1+b)(1+c) + (1-a)(1-b)(1-c) >= 0.                 (B1)

For Boolean signs, `Q/4` is exactly the indicator of monochromaticity. This is a special cancellation at arity three, not a general principle that a pattern is controlled by pair correlations.

Let `F : R² × Ω -> [-1,1]` be jointly measurable on a probability space. Let `T_i=(v_i0,v_i1,v_i2)` be finitely many templates with distinct vertices. A placement is `y=(p,u)`, with `p∈R²` and `u∈SO(2)`; write `yT_i=(p+uv_i0,p+uv_i1,p+uv_i2)`. Let `Π` be any probability measure on the disjoint union of these template-placement spaces. Assume the evaluations in question are measurable. In the application `Π` can use bounded translation windows and normalized angular measure, with assigned template weights.

Define the projected quadratic defect

    D_F(Π) = (1/4) ∫ E_ω Q(F(yv_i0,ω), F(yv_i1,ω), F(yv_i2,ω)) dΠ(i,y). (B2)

No stationarity is needed for the rounding theorem.

## 2. Loss-one cell rounding

**Theorem B1.** There exists a measurable two-coloring `c:R²->{-1,1}` whose `Π`-averaged monochromatic defect is at most `D_F(Π)`.

**Proof.** Put `δ=min_{i,j≠k}|v_ij-v_ik|>0`. Partition the plane into half-open squares of side `h` with `sqrt(2)h<δ`. Let `(U_z)_{z∈Z²}` be independent uniform variables on `[0,1]`, independent of `ω`. At a point in cell `z`, set

    c_(ω,U)(x) = +1 iff U_z <= (1+F(x,ω))/2.                   (B3)

Every placed triple occupies three distinct cells, for every translation and rotation: two points in one cell have distance less than `sqrt(2)h`. Conditional on `ω`, the three rounded signs are therefore independent with respective means `a,b,c`. Consequently

    P_U(monochromatic | ω) = [Π_j(1+a_j)+Π_j(1-a_j)]/8 = Q/4. (B4)

The coloring is measurable because it is defined on a countable Borel partition using measurable values of `F`. Tonelli gives expected `Π`-defect exactly `D_F(Π)`. Some realization has defect at most that expectation. ∎

**Simultaneous-coordinate hostile.** Take five fixed distinct points, ten separate atomic objectives (one for each three-element subset), and `F≡0`. Every original objective is `1/4`. Every two-coloring has three points in one color, so at least one output objective is `1>1/4`. Thus coordinatewise simultaneous improvement is genuinely false in this general framework, even though every chosen weighted objective admits the conclusion of B1.

The same proof covers any countable family with a common positive lower bound on its pairwise within-template distances. An uncountable family is also allowed when the combined parameterization and placement measure are measurable and this separation bound holds. The theorem does not assert simultaneous coordinate-wise improvement of several nonzero defects; a realization that optimizes one weighted objective need not optimize all objectives.

**The coefficient one is asymptotically sharp in this abstract formulation.** Take `2m` distinct fixed points, use every three-element subset with equal weight at its fixed placement, and set `F≡0`. Then `D_F=1/4`. Every coloring has defect at least

    2 binom(m,3) / binom(2m,3) -> 1/4,

because a balanced split minimizes `binom(r,3)+binom(2m-r,3)`. Thus no universal coefficient strictly below one can replace the coefficient in Theorem B1. This sharpness is for the stated arbitrary-placement-measure formulation, not a sharpness claim for Haar-uniform plane placements.

## 3. Application to pair-correlation transfer

Suppose an invariant sign process `s(z)` on a discrete carrier has correlation `ρ(w)=E[s(0)s(w)]` and produces a bounded measurable field with

    E[F(x)F(y)] = ρ(y-x)                                      (B5)

at the distinct-point evaluations used, interpreted by continuous extension and rotation invariance as in blueprint §6. For each template let

    δ_i = [1+ρ(v_i1-v_i0)+ρ(v_i2-v_i0)+ρ(v_i2-v_i1)]/4.        (B6)

For an actual invariant sign process this is its exact monochromatic probability. If the three correlations of each rotated and translated template agree with those in (B6), B2 equals the weighted average of the `δ_i`. Theorem B1 therefore transfers every chosen finite weighted defect objective without increasing it. This extends the blueprint's zero-defect transfer to nonzero defects and arbitrary distinct three-point templates, including collinear ones; equilateral geometry is not used in B1--B6.

**Conditional boundary.** The blueprint supplies carriers for the algebraic numbers and the algebraic closure of `Q(x)` for one real parameter `x`. A finite family with several algebraically independent shape parameters does not automatically lie in one of those carriers. No extension of the source template property to higher transcendence degree is assumed here. B5 is not a consequence of positive definiteness alone. It is the substantive rigidity/projection input. Nor does the conclusion give a coloring with all individual defect inequalities at once.

For exact avoidance of countably many templates one can avoid the common minimum-separation condition entirely. Define `c(x)=sign_+(F(x))`, where zero is assigned `+1`. If this threshold gives one color to three arguments, all have the same weak sign and

    1_mono(sign_+(a),sign_+(b),sign_+(c)) <= Q(a,b,c).         (B7)

Thus `E Q=0` forces zero monochromatic defect. Countably many template constraints can be intersected after Tonelli, using countably many bounded placement windows. This is the general zero-defect mechanism behind §6.

The deterministic pointwise coefficient in B7 is sharp: at `a=b=c=0`, its left side and `Q` are both one. Relative to `Q/4`, this is the sharp factor four for that deterministic threshold bound. Cell rounding removes the factor by adding independent local threshold data at a scale below the geometric separation.

## 4. Boolean degree: where pair transfer stops

For `k>=2`, write `M_k(x)=1` when `x_1=...=x_k`, and zero otherwise. Its unique multilinear representative on `{±1}^k` is

    M_k(x) = 2^(1-k) Σ_(S⊆[k], |S| even) Π_(i∈S) x_i.       (B8)

This follows by expanding `2^(-k)[Π(1+x_i)+Π(1-x_i)]`. Hence its Fourier degree is `k` for even `k` and `k-1` for odd `k`. In particular, it is quadratic precisely for `k=2,3`.

**Theorem B2 (quadratic certificate obstruction).** For every `k>=4`, no polynomial of degree at most two can be nonnegative on the Boolean cube, vanish at all nonmonochromatic Boolean vertices, and be strictly positive at at least one monochromatic vertex.

**Proof.** Multilinearize the polynomial using `x_i²=1`; this cannot increase its degree. Let its values at all-plus and all-minus be `A,B>=0`, with `A+B>0`. Since all other values vanish, its Fourier coefficient at any four-element subset is `(A+B)/2^k>0`. Degree at most two requires this coefficient to be zero. ∎

Nonnegativity on the whole continuous cube is stronger and is therefore ruled out too. If nonnegativity is dropped, cancellation at the two exceptional vertices does not rescue a degree-two certificate: for `A=-B≠0`, every three-element Fourier coefficient is nonzero. The nonnegative form is the relevant one for expectation-zero transfer.

**Exact hostile.** Let `X` be uniform on the eight sign vectors of even parity `ΠX_i=+1`, and `Y` uniform on the eight vectors of odd parity `ΠY_i=-1`, for four coordinates. Both laws have uniform marginals on every proper subset, so all moments of orders one, two and three agree. But

    E M_4(X)=1/4,       E M_4(Y)=0.                           (B9)

No method retaining only those proper-subset marginals can recover four-point monochromaticity. This is a distributional obstruction, not a claim that the two laws are stationary geometric processes.

## 5. New hypotheses and decisive tests

- **H-B1, OPEN (geometric fourth-order transfer).** Identify a class of four-point templates and invariant processes for which projection preserves the specific fourth moment in B8, not all fourth moments. Cheapest hostile: the parity laws B9; any proposed hypothesis must explicitly exclude their mechanism.
- **H-B2, OPEN (simultaneous defect region).** Characterize the attainable deterministic vectors of measurable placement defects and their closed convex hull (the closure of the randomized defect region). Atomic placement measures show that the deterministic region need not be convex. B1 places every nonnegative weighted support functional below the corresponding original objective; it does not identify a single point below every original coordinate. Study two finite templates with competing threshold objectives before a universal claim.
- **H-B3, OPEN (separation-free robust transfer).** Extend B1 to a prescribed probability measure on triples of distinct points without a uniform separation lower bound. Shrinking cells gives collision error at most the `Π`-mass of triples with a pair closer than the cell diameter, which tends to zero. Thus an `ε`-loss follows immediately for every `ε>0`; exact attainment requires an additional compactness argument and is not claimed.

**Proved H-B3 survivor.** With finite total placement mass and no repeated vertices almost surely, for every `ε>0` a measurable coloring has defect at most `D_F+ε`: use a sufficiently fine partition, bound the exceptional collision event by one, and use monotone convergence of the small-distance exceptional probabilities. No uniform geometric scale is needed for this approximate statement.

## 6. Connection contract

- Source: blueprint §§4,6, bounded projected sign field and pair coefficient identity.
- Target: stochastic rounding / NAE-3 constraint optimization on geometric placements.
- Map: one independent uniform threshold per fine square, equation B3.
- Preserved: the expectation of every separated triple's monochromatic indicator, then a specified weighted average objective.
- Destroyed: individual sample correlations, all simultaneous nonzero coordinate bounds, microscopic continuity, and fourth-order data.
- Sidecars: minimum distance, the placement probability measure, and independent cell addresses.
- Hostiles: coincident points, competing objectives, and B9's fourth-parity ambiguity.

This is a genuine arity-sensitive transfer, not a tournament encoding: the native relation is a symmetric ternary NAE constraint. An orientation would discard the relevant Boolean degree.
