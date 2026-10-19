# A finite rational contact catalogue necessarily has exponential bit growth

**2026-10-19. PROVED using the companion's complete degree-at-most-two classification.** The height argument below is elementary. The cubic two-map chain has the optimal exponential base among all finite catalogues of global rational advancing contact maps, even with adaptive nonperiodic switching.

## 1. Precise scope

Let f(t)=4t(t²-1)/(t²+1)² be the signed Pythagorean spread. A global rational contact map h∈Q(t) has an identity

    y(t)²=1+epsilon f(t)+eta f(h(t)),

for a rational function y and fixed signs epsilon,eta∈{±1}. Actual orbit points must additionally represent admissible positive primitive petals and must use valid ports. Those guards are retained, not inferred from the rational identity alone.

The companion [minimal-degree theorem](square_sum_cubic_contact_degree_minimality_20261019.md) proves:

- constant maps give only the zero-spread boundary;
- degree-one maps preserve the unordered petal whenever both endpoints are admissible positive-spread petals;
- degree-two maps do not exist;
- useful degree-three maps do exist, explicitly.

Thus every step that advances to a different admissible petal has rational-map degree at least three. A degree-one change of Euclidean parameter cannot be used as artificial padding in the chain length.

Fix any finite catalogue H={h_1,...,h_s} of these advancing maps. A chain may choose its next map arbitrarily, depending on the complete history. Neither periodicity nor a fixed choice sequence is assumed. Require its petals to be pairwise distinct.

**Theorem 1.** There are constants c>0 and k_0, depending only on the catalogue, such that every such finite chain with k≥k_0 transitions has final primitive hypotenuse bit length at least c·3^k. Any common-center integer realization of the chain has at least this order of label bits. In particular no fixed finite catalogue yields arbitrarily long distinct chains of polynomial bit length.

For a fixed infinite chain, the same exponential lower bound holds eventually. The [alternating cubic construction](square_sum_cubic_contact_chains_20261019.md) attains Theta(3^k) for every fixed initialized orbit. Hence base three is optimal within this catalogue model.

## 2. An elementary two-sided rational height bound

For t=a/b in reduced form put h(t)=log max(|a|,|b|), using (1,0) for infinity. For a fixed degree-d rational map F over Q, there is a constant C_F such that

    |h(F(t))-d h(t)|≤C_F                                 (1)

at every projective rational point.

Here is a self-contained proof. Represent F by coprime homogeneous integer forms P(X,Y),Q(X,Y) of degree d, with no common projective zero. Coefficient bounds give

    max(|P(a,b)|,|Q(a,b)|) ≤ A max(|a|,|b|)^d.

On the compact real set max(|x|,|y|)=1, the continuous function max(|P(x,y)|,|Q(x,y)|) has a positive minimum c, because P,Q have no common zero. Homogeneity gives the matching lower bound c max(|a|,|b|)^d.

A nonzero integer resultant multiple R has homogeneous Bezout identities

    U P+V Q=R X^(2d-1),
    U'P+V'Q=R Y^(2d-1).

They follow by clearing denominators in the inverse Sylvester matrix. If gcd(a,b)=1, any common divisor of P(a,b),Q(a,b) consequently divides R. Reducing the output fraction therefore changes the logarithmic height by at most log|R|. Taking logarithms in the preceding bounds proves (1). The constants are computable from the fixed forms; no statistical cancellation assumption is used.

For a finite catalogue, choose C≥0 bounding every C_F. Since each advancing degree is at least three,

    h(t_(i+1)) ≥ 3h(t_i)-C.                              (2)

## 3. Finite-height escape works even under adaptive switching

Let B=C/2+1. There are only finitely many rational projective points with h(t)≤B: their reduced integer numerators and denominators are bounded by exp(B). Let M be any upper bound for their number.

Pairwise distinct petals imply pairwise distinct parameters. Thus among the first M+1 points of any sufficiently long chain, some t_j has h(t_j)>B, with j≤M. Iterating (2) from that point gives

    h(t_(j+r)) ≥ C/2+3^r(h(t_j)-C/2) > C/2+3^r.          (3)

In particular the height can never re-enter the bounded region. For k≥M,

    h(t_k)>3^(k-M).

This is an exponential lower bound with constants depending only on the fixed catalogue, independent of the input seed and of the adaptive map choices. The threshold M can be extremely large; the theorem is an asymptotic obstruction, not a practical bound for a small search.

For a rational Euclidean parameter a/b, primitive Pythagorean normalization gives

    C(t)=(a²+b²)/delta,  delta∈{1,2},
    log C(t)=2h(t)+O(1).

Hence (3) gives the asserted primitive hypotenuse bit growth. A common-center integer scaling has H divisible by each primitive hypotenuse and center 2H², so its largest label has logarithmic size at least 2log C(t_k)+O(1). The same lower bound therefore survives integer normalization and lcm sharing.

This proves Theorem 1. In particular adaptive switching among finitely many maps cannot defeat the lower bound by arranging repeated lucky cancellations: the resultant bound controls all such cancellations uniformly.

## 4. Periodic catalogues have an exact arithmetic growth rate

If the choice sequence has period r, let F be its r-step composition and D its degree. Multiplication of rational-map degrees gives D≥3^r. Equation (1) implies that

    lambda=lim_(n→infinity) D^(-n) h(F^n(t_0))

exists: consecutive normalized heights differ by at most C_F D^(-n-1), a summable geometric sequence. Moreover

    h(F^n(t_0))=lambda D^n+O(1).                          (4)

If lambda were zero, (4) would bound all rational heights. The period subsequence would then take only finitely many values, hence be eventually periodic, contradicting distinct petals. Therefore lambda>0. Each intermediate phase is a fixed rational composition, so (1) transfers (4) to it as well.

Thus the arithmetic exponential rate of a periodic scheme is D^(1/r), at least three. The two cubic maps have D=9 and attain the lower rate exactly. Their special Gaussian factorization additionally identifies every old and new hypotenuse-prime contribution, rather than merely bounding total height.

A fixed birational change of the rational parameter is a degree-one Möbius transformation. It preserves map degree and changes rational height by only O(1), so it cannot remove this obstruction. More generally a rational semiconjugacy π∘F=G∘π with nonconstant π forces deg F=deg G by multiplication of degrees. A finite rational cover alone therefore cannot lower the dynamical degree.

## 5. The remaining polynomial-bit problem now needs a different mechanism

This theorem excludes more than a particular cubic orbit or a single repeated section. It excludes every fixed finite global rational catalogue, including nonperiodic adaptive switching, as a polynomial-bit construction for growing distinct chains.

It does not rule out:

- an unbounded family of maps whose coefficients depend on the evolving arithmetic data;
- algebraic correspondences with rational branches selected only at the particular input, rather than one global rational function;
- integer-specific contact operations that do not extend to a rational identity for all parameters;
- constructions building many vertices in a shallow branching tree instead of following one long orbit.

The last option is already realized by the polynomial-bit balanced-tree theorem. The first three remain potential long-chain routes. None automatically gives consecutive-integer coverage, which is a separate packing constraint.

The new research distinction is therefore exact: polynomial bit growth cannot come from improving the fixed rational atlas alone. A successful long-chain proof must leave that model or abandon pairwise distinct advancement. No lower bound for every conceivable contact-chain construction is asserted.
