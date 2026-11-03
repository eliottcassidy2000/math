# An integral carry scheme and the Pythagorean two-valued correspondence

**Status:** explicit proved models and consequences, extending the S164 bookkeeping. Gaussian factorization, the rational-circle group, and quotient two-valued laws are classical. The points here are the exact integral model of the carry, its rational trivialization, and a rational-selector obstruction with the absolute-leg quotient correctly retained. These are not constructions of the missing strong ternary branches.

## 1. The carry is the multiplication table of a concrete integral algebra

Fix a finite set S of distinct rational primes congruent to 1 modulo 4. Choose Gaussian primes `pi_p` with `pi_p conjugate(pi_p)=p`. For an integer vector `v∈Z^S`, let

    alpha(v)=∏ pi_p^max(v_p,0) conjugate(pi_p)^max(−v_p,0),
    g(v,w)=∏ p^((|v_p|+|w_p|−|v_p+w_p|)/2).

S164 proves `alpha(v)alpha(w)=g(v,w)alpha(v+w)`. Introduce the finitely presented ring

    A_S = Z[X_p,Y_p : p∈S] / (X_pY_p−p : p∈S).             (1)

For v define its normal monomial `e_v=∏ X_p^max(v_p,0) Y_p^max(−v_p,0)`.

**Theorem.** The e_v form a Z-basis of A_S, and

    e_v e_w = g(v,w)e_(v+w).                                (2)

The assignment `X_p→pi_p, Y_p→conjugate(pi_p)` is a ring map to Z[i] taking (2) to Gaussian multiplication.

**Proof.** Replacing every paired `X_pY_p` by p spans the displayed monomials. Embed the resulting normal forms into the Laurent polynomial ring `Q[T_p^(±1)]` by `X_p→T_p, Y_p→p T_p^(−1)`. The images are distinct Laurent monomials with nonzero integer coefficients, so are Z-linearly independent. Pairing opposite powers in a product gives exactly the exponent of p in g. The Gaussian specialization respects every defining relation. ∎

For nonempty S, the specialization to Z[i] is not injective and is not an equivalence between the formal torus and the actual Gaussian arithmetic. For example, if `pi_5=2+i`, the specialization kills `X_5+Y_5−4`. Keeping this distinction avoids claiming that an entire positive-dimensional scheme has been identified with a discrete collection of triples.

Let `P=∏_(p∈S)p`. Over `Z[1/P]`, each X_p becomes a unit with inverse `Y_p/p`, so

    A_S[1/P] = Z[1/P][X_p^(±1):p∈S].                       (3)

Thus the localized scheme is a split torus. The integral scheme in (1) is not being asserted to be a group scheme over Z: inverse prime factors are needed for the torus multiplication in these coordinates. The carry is explicitly a coboundary after localization. Set `c(v)=∏p^max(−v_p,0)`. Then

    g(v,w)=c(v)c(w)/c(v+w),
    f_v=e_v/c(v),             f_v f_w=f_(v+w).              (4)

This is a genuine distinction: the positive integral normal form incurs a carry, while the localized Laurent coordinates multiply without one. The cocycle should not be advertised as a nontrivial rational group extension.

For a particular p∈S, reduction modulo p gives

    A_S ⊗ F_p ≅ F_p[X_p,Y_p]/(X_pY_p)
                  ⊗ F_p[X_q^(±1):q∈S, q≠p].              (5)

It has two reduced components, meeting along `X_p=Y_p=0`. This is the precise algebraic degeneration behind local cancellation. From (2),

    e_v e_w = 0 mod p   iff   v_p w_p < 0.                 (6)

Indeed opposite nonzero signs contribute at least one p; otherwise the structure constant is nonzero modulo p and the surviving normal monomial is nonzero. Under Gaussian specialization a nonzero product of normalized factors also acquires a coordinate-common p exactly in this case, by the exact gcd theorem in S164. This distinguishes a raw product killed by a carry from its renormalized primitive output; normalization is not a ring homomorphism modulo p.

## 2. The absolute primitive triple quotient is an explicit correspondence

The rational norm-one torus is

    T: x²+y²=1,
    (x,y)*(u,v)=(xu−yv,xv+yu).

A primitive triple `(a,b,c)`, with the odd leg listed first, gives the rational point `(a/c,b/c)` once signs are chosen. Forgetting both leg signs is the quotient by the four transformations `(x,y)→(±x,±y)`. Its invariant coordinate is `U=x²`; its invariant ring over Q is exactly `Q[U]`. Equivalently, first quotient the central two-torsion `(x,y)→(−x,−y)` and then quotient inversion `(x,y)→(x,−y)`.

For two absolute input classes U,V, the two possible product classes are the roots Z of

    Z² − 2(1−U−V+2UV) Z + (U+V−1)² = 0.                  (7)

**Proof.** Put `U=x²,V=u²`. The two relative sign choices give

    Z_±=(xu ± yv)²
       =1−U−V+2UV ± 2xyuv.

Their sum is twice the displayed center and their product is `(U+V−1)²`. This proves (7), including repeated roots. The discriminant is

    16UV(1−U)(1−V).                                      (8)

It is nonsquare in `Q(U,V)`, since its valuation at the irreducible divisor U=0 is one. Consequently (7) is irreducible over that field, and there is **no rational function of U,V selecting one branch generically**. The same applies to any rational formula intended to select a branch on every pair of positive rational primitive triples: their U-values are infinite, their Cartesian square is Zariski dense, and any proposed formula would force the polynomial identity (7). Pointwise square-root conventions or a supplied relative sign do select a branch, but they are extra structure rather than rational descent of the group law.

For actual input triples, U and 1−U (and V and 1−V) are rational squares; the odd-leg ordering additionally selects the rational-circle points whose primitive first numerator is odd. Consequently both roots in (7) are rational squares for every admissible pair, even though the polynomial does not split over the generic function field Q(U,V). An arbitrary rational number U is not automatically a triple class.

This is a finite algebraic correspondence, not a group law on the affine quotient line. Associativity is understood as multiset composition inherited from the group, exactly as in S164's `[v]⋆[w]=[[v+w],[v−w]]`. Degenerate outputs must be retained for cancellation. No tournament is called a Kummer variety here.

For `(3,4,5)` and `(5,12,13)`, the two roots of (7) are

    (33/65)²,  (63/65)².

Thus the two rational points over this one quotient pair encode precisely the two output triples `(33,56,65)` and `(63,16,65)` from S164. The relative bit is the choice of point in this degree-two correspondence.

## 3. Where the literal tournament map sits

For r,s≥1, S164's endpoint-reversal tournament E(r), with r ordered interior vertices, has one increasing Hamilton path and 2^r nonincreasing paths. A nonincreasing path is addressed by a subset of the interior vertices. Attaching Gaussian prime-power marks to the r slots turns that address into a choice of pi or its conjugate in each slot. Identifying the two endpoint terminals of E(r) and E(s), and ordering the first interior block before the second, literally gives the strong tournament E(r+s). Its nonincreasing path addresses concatenate.

This path concatenation can be evaluated in the normal monomial algebra (1), then specialized to Gaussian multiplication. Distinct repeated-prime slots remain distinct vertices; the relation `X_pY_p=p` acts only on their arithmetic decoding. Primitive normalization removes the scalar in (2). Passing to absolute triple classes then lands in (7), where the retained relative sign chooses one branch. Hence the scheme and correspondence describe actual maps of this marked strong-tournament path family.

The maps have different information losses: concatenation is bijective on ordered address pairs; integral reduction collects prime carries; primitive normalization discards scales; and the absolute quotient discards global signs. None of these maps provides the three physical successor states required for the full ternary tournament program. In particular, the two components in (5), the two roots in (7), and three affine targets `3h,3h+2,3h+4` are not interchangeable counts of states.

## Sources and controls

- [S164 literal endpoint gluing and Gaussian signs](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session164_emergent_20261103/triple_gluing/NOTE.md): path gluing, sign-interface theorem, exact gcd/carry and two-valued lattice law.
- [Amnon Yekutieli, *Pythagorean Triples, Complex Numbers, Abelian Groups and Prime Numbers*](https://arxiv.org/abs/2101.12166): the prior rational-circle group and Gaussian organization of triple fibers.
- [V. M. Buchstaber and V. Dragovic, *Two-Valued Groups, Kummer Varieties, and Integrable Billiards*](https://amj.math.stonybrook.edu/pdf-Springer-final/018-0085.pdf): prior quotient and multivalued-group framework. The explicit elementary correspondence above is proved directly and is not offered as a new general theory.

The accompanying checker verifies the integral normal-form multiplication, coboundary, prime-local cancellation, Gaussian evaluation, and rational correspondence on finite grids. Those checks corroborate the formulas; the ring basis, generic nonsplitting, and physical source maps are established by the proofs.

**Verification receipt.** Both normal and optimized Python runs pass all 729,799 active checks and produce byte-identical JSON. The geometric controls use 15,625 ordered valuation-vector pairs on three marked primes, 19,683 triples for reassociation, and the explicit norm-65 example.
