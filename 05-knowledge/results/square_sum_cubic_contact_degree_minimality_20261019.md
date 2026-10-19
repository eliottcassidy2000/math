# Cubic degree is minimal for nontrivial global rational AP contact maps

**Status: PROVED scoped, with exact certificate replay.** The degree-two exclusion combines elementary divisor arguments with three explicitly supplied rational polynomial-combination certificates. No bounded census, numerical approximation, unproved Mordell–Weil basis, or modular lifting assumption enters the theorem.

## 1. The precise minimality theorem

Let

\[
 f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad
 k_\sigma(t)=\frac{t^2+2\sigma t-1}{t^2+1},
 \quad k_\sigma(t)^2=1+\sigma f(t).
\]

The degree of a reduced rational function `h=P/Q` is `max(deg P,deg Q)`, its degree as a map of the projective line. A **global rational contact map** is a rational function `h in Q(t)` for which, at prescribed signs `epsilon,eta in {±1}`, the identity

\[
 y(t)^2=1+\varepsilon f(t)+\eta f(h(t))                 \tag{1}
\]

holds for some `y in Q(t)`.

**Theorem 1.** Every such map satisfies the following classification at low degree.

1. If `deg h=0`, its output has zero spread: `h in {0,1,-1,infinity}`.
2. If `deg h=1`, then `f(h)=-epsilon eta f(t)` and `y²=1`. Whenever input and output spreads are both positive, they are equal and correspond to the same unordered primitive Pythagorean petal.
3. There is **no** degree-two global rational contact map.

Consequently every global rational contact map that advances between distinct positive-spread petals has degree at least three. This lower bound is sharp: the maps

\[
 g_+(t)=-\frac{2t^3+t^2-2t+1}{t(t^2-4t+1)},\qquad
 g_-(t)=-\frac{t^2-4t+1}{t^3-2t^2+t+2}                \tag{2}
\]

have degree three, give respectively `++` and `--` contacts, and advance between distinct admissible primitive petals in the fixed cusp chart established in [the cubic contact-chain theorem](square_sum_cubic_contact_chains_20261019.md).

The word “global” matters. Equation (1) is an identity in `Q(t)`, not merely a rational-square condition at one input or along a specially chosen set of inputs. Adaptive choices, algebraic correspondences, and parameterizations by another curve are outside this minimality statement.

## 2. Complete coverage of the degree-two case

Write `x=k_epsilon(t)` and `z=k_eta(h(t))`. Then

\[
 y^2=x^2+z^2-1,\qquad (y-z)(y+z)=\varepsilon f(t).     \tag{3}
\]

The poles of `x` are the simple conjugate pair `±i`. Write `h=P/Q` in coprime form with maximum degree two. The degree-four polynomial `P²+Q²` is its contact pole polynomial.

The residue argument in [the low-degree obstruction](square_sum_low_degree_contact_obstruction_20261019.md), Section 5, first proves that this quartic must split over `Q` into two quadratic factors, with multiplicity allowed, and every root belongs to `Q(i)`. An irreducible quartic would force `y=alpha x+beta z+c`; solving this relation makes a degree-four function `z` fractional-linear in the degree-two function `x`, a contradiction.

There are now exactly three exhaustive configurations.

### A. One repeated quadratic factor

Section 7 of the low-degree obstruction excludes this case. If its roots are disjoint from `±i`, one of `y±z` has four zeros counted with multiplicity but at most two poles, impossible. If they coincide with `±i`, the same factor is constant; then `z` is affine in `f`, contradicting the branch-value sets of `k_eta o h` and of `f`.

### B. Two distinct quadratic factors, one equal to `t²+1`

The divisor of `f` forces, up to exchanging the factors,

\[
 y-z=c\frac{q(t)}{t^2+1},\qquad
 y+z=\frac{4\varepsilon}{c}\frac{t(t^2-1)}{q(t)(t^2+1)}.
\]

The [shared-pole exclusion](square_sum_quadratic_shared_pole_exclusion_20261019.md) expands the necessary Gaussian quadratic form for `z`. Three exact two-by-two determinants reduce a possible rational solution either to a forbidden coincident pair, a zero normalization constant, or a rational root of `2b⁴+4b+1`. The rational-root theorem excludes the latter. This proof is elementary and handles both signs directly.

### C. Two distinct quadratic factors, neither equal to `t²+1`

The two new pole pairs must carry opposite residue signs. Each of `y−z,y+z` has simple poles at `±i` and at one new pair, and simple zeros at the other new pair. Their remaining zeros partition the four simple real zeros `{0,1,-1,infinity}` of `f`, with exactly two assigned to each factor.

There are six assignments. Complementary assignments are equivalent by exchanging `y−z` and `y+z` and replacing `h` by `-1/h`, since `k_eta(-1/h)=-k_eta(h)`. Hence there are exactly three representatives containing zero:

\[
 \{0,1\},\qquad \{0,-1\},\qquad \{0,\infty\}.         \tag{4}
\]

Replacing `h` by `eta h` reduces the output sign to `eta=+1`, because `k_eta(h)=k_+(eta h)`. Replacing the input variable by `-t` reduces `epsilon=-1` to `epsilon=+1`; it interchanges the first two representatives and preserves the third. All three are treated separately, so no unproved symmetry merges inequivalent patterns.

The [adjacent-pole exclusion](square_sum_quadratic_disjoint_adjacent_exclusion_20261019.md) supplies necessary coefficient ideals for the first two cases and explicit certificates proving respectively

\[
 m^2+4=0,\qquad
 (m-2)^3(m+2)(m^2+20m+52)=0.                          \tag{5}
\]

The first has no real solutions. In the second, rationality forces `m=±2`: `m=2` makes the new pole pairs coincide, while the coefficient equations at `m=-2` force an input pole pair to be shared. Both contradict this configuration.

The [opposite-pole exclusion](square_sum_quadratic_opposite_pole_exclusion_20261019.md) treats the third case and supplies an explicit certificate forcing

\[
 m^2+16=0,                                           \tag{6}
\]

which again has no real solutions. Its normalization variable `m` differs from the adjacent cases; only the stated contradiction within each separately derived system is used.

Configurations A, B, and C exhaust all degree-two maps. Each is impossible, proving Part 3. Parts 1 and 2 are proved in Sections 3–4 of the low-degree obstruction. The rational square identities and primitive hypotenuse transport for (2) are proved in the cubic contact-chain note, establishing sharpness. QED.

## 3. What the computational certificates prove

The three ideal-membership certificates contain 349, 204, and 402 nodes. A node is either one of the displayed necessary input equations or a sum of polynomial multiples of earlier nodes, with exact rational coefficients. Thus each node belongs to the input ideal by induction. Their final nodes are precisely the three polynomials in (5)–(6).

The verifiers reconstruct the input equations from the mathematical models and check every polynomial addition and multiplication. They do not assume that the discovered Groebner bases are correct, that a numerical root search is exhaustive, or that a modular calculation lifts. The certificates prove ideal membership directly; the final real/rational contradictions are elementary.

For reproducibility the mathematical packet consists of:

- the low-degree obstruction note and its elementary formula checks;
- the elementary shared-pole note and determinant checker;
- the two adjacent systems, their two polynomial-combination certificates, and their independent checker;
- the opposite system, its polynomial-combination certificate, and its independent checker;
- this coverage and minimality theorem;
- the cubic construction supplying an actual sharp example.

The publication verifiers and input-equation generator use only the Python standard library. Every check remains active under `python -O`, and normal/optimized controls agree byte-for-byte. The scripts that originally discovered the ideal-membership relations are useful for regeneration but are not prerequisites for verifying the proof.

## 4. Consequence for the research direction

The cubic maps are not merely smaller formulas found by a lucky search: degree three is the first possible global rational degree for a genuinely new positive contact. This is a precise algebraic obstruction behind the repeated height expansion in fixed-map constructions. It does not by itself rule out compact chains from adaptive maps or correspondences, and it does not identify square-sum or graceful labeling problems with a Collatz dynamics.


## Repository replay

From the repository root run:

- `python3 05-knowledge/experiments/verify_low_degree_contact_identities_20261019.py`
- `python3 05-knowledge/experiments/verify_quadratic_shared_pole_exclusion_20261019.py`
- `python3 05-knowledge/experiments/verify_quadratic_contact_certificates_20261019.py`

All use only the Python standard library. The three input certificates and three published output receipts are in `05-knowledge/results`; the unified verifier reads that directory. No exploratory CAS scripts are needed.
