# Weighted Wronskians, Frobenius jet budgets, and prime-power cellular templates

**Session:** 2026-10-19. **Status:** elementary proofs and exact finite controls; scoped reusable instruments. No claim of a new resolution of a major conjecture, or of priority over the classical theories cited below. The supplied Barreto text is an inspiration and an object of analysis, not an assumed theorem dependency.

## 0. Results in brief

1. A weighted-support Wronskian bound gives a universal sufficient certificate for injectivity of arbitrary matched-point Hermite evaluation. It retains the coefficient-space rank, the vanished specializations, and each node's jet depth. For weight `(1,k)` the Wronskian degree bill is `ell*D-(k+1)*ell*(ell-1)/2`, sharper than the unweighted bill.
2. The Frobenius graph ideal `(p,z^p-w)` is stable under every ordinary `z` derivative. Mixed and Hasse derivatives have an exact loss budget: one unit for each `w` derivative and at most `floor(a/p)` additional units for Hasse `z` order `a`. A sorted assignment formula turns these losses into a lower bound for the valuation of the same determinant minor.
3. Every scalar additive cellular automaton modulo `p^m` has a fixed finite prime-power-time template after time `p^(m-1)`. The template depth `m-1` and the erasure time for different integral lifts of the same mod-p rule are both universally sharp.
4. For the Pascal rule `1+X`, the surviving kernel has exactly `p^(m-1)+1` terms at those times. Seed content, overlaps, torus aliasing, noncommuting matrix rules, and nonlinear Rule 30 are explicit boundaries.
5. A bounded weighted signed-prime functional alone cannot control an unweighted race sign. An exact Abel transform and a constructive sign-sequence hostile identify the missing local information; this is not a counterexample involving a Dirichlet character.

## 1. Inheritance and source status

### Closest mechanisms

- `THM-4711-legendre-wronskian-and-binomial-cube-lucas-carry-divisibility.md`: a finite-degree Wronskian, an actual carry companion, and a retained characteristic-p integration constant. We use its methodological boundary, not its hypergeometric conclusions.
- `THM-3882-rational-dual-one-place-wronskian-projection-criterion.md`: the exact base-divisor and ramification ledger for a two-dimensional linear series. Our higher-rank polynomial Wronskian is the same degree-versus-vanishing mechanism on an affine line; it does not imply Jacobian-conjecture entry.
- `THM-2172-frobenius-collapse-of-mahler-and-twisted-cyclic-packets.md`: whole-layer Frobenius in characteristic p. Section 4 retains the extra p-adic precision coordinate rather than silently transporting a mod-p identity to mod-p^m.
- `THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy.md`: carried reset thickness is an integral coordinate. Here a distinct scalar polynomial carrier has an explicit precision-depth law.
- `prime_power_bias_20261014d.md`: the full residue-power transform retains powers and lower cutoffs; it is not a single race scalar.

The canonical hostiles are `1,z^p` for characteristic-p ordinary Wronskians; `D_w(z^p-w)=-1` for differential direction; Rule 30 at time two for nonlinear dilation; and a noncommuting matrix lift for power-lifting failure. The corrected near miss is the repo's warning that a zero derivative in characteristic p need not mean a constant. The least-used sidecar here is the **per-column differential loss budget**, rather than another bare rank statistic.

### Primary checks

- Alin Bostan and Philippe Dumas, *Wronskians and linear independence*, https://arxiv.org/abs/1301.6598; author PDF https://mathexp.eu/dumas/publications/Bostan_Dumas_2010_WLI.pdf. This supplies a classical characteristic-zero Wronskian/Vandermonde baseline and explicitly warns about positive characteristic. All polynomial estimates below are proved directly.
- Satoshi Takahashi, *Self-similarity of linear cellular automata*, JCSS 44 (1992), 114–140, https://doi.org/10.1016/0022-0000(92)90007-6. The publisher record establishes relevant classical prime-power CA literature. We do not attribute our precise finite-time formula to a theorem not inspected in that paper, nor claim priority for it.
- The primary arXiv records https://arxiv.org/abs/2610.15483 and https://arxiv.org/abs/2610.21099 resolve to the cited Banks and Kerr preprints. They confirm those bibliographic records and their claimed subject, not peer review or an independent validation of every headline. The supplied Barreto manuscript's public version was not located in this scoped search.

The source attachment prints `p sum_T floor(r/p) | D`. The row-factor argument mathematically gives `p^(sum_T floor(r/p)) | D`; this exponent is required for the logarithmic determinant estimate. This note uses the exponent form and does not silently treat the printed expression as correct notation.

### Portfolio and concept board

Anchor: characteristic-zero rank plus integral determinant valuation. Niche: prime-power scalar cellular automata. Wildcard: what a weighted prime statistic forgets.

The live concepts are weighted support, linear-series ramification, differential filtration, p-adic precision, spatial aliasing, and the counting cutoff. Each connection below specifies its carrier and loss rather than identifying similarly named objects.

## 2. Weighted matched-point Hermite interpolation

Let K have characteristic zero. Fix integers `D>=0`, `k>=1`, put `E=floor(D/k)`, and define

`V_(D,k) = span_K { z^r w^s : r,s>=0, r+k*s<=D }`.

Take M points `(alpha_i,beta_i)` such that the alpha coordinates are pairwise distinct and the beta coordinates are pairwise distinct. At node i observe the ordinary z derivatives of orders `0,...,m_i-1`, where `m_i>=0`.

### Proposition 2.1. Weighted coefficient-space Wronskian budget

If a nonzero `F in V_(D,k)` has exactly e identically zero specializations `F(z,beta_i)`, then `0<=e<=E`. Factor

`F(z,w) = H(w) G(z,w),  H(w)=product_(F(z,beta_i)=0)(w-beta_i)`.

The weighted degree of G is at most `D-k*e`. Let ell be the dimension of the K-span of G's coefficient polynomials in w. Then

`1<=ell<=floor((D-k*e)/k)+1=E-e+1`.

Its nonzero Wronskian W satisfies

`deg W <= B(D,k,e,ell)`

where

`B(D,k,e,ell) = ell*(D-k*e) - (k+1)*ell*(ell-1)/2`.       (2.1)

**Proof.** Weighted degree is additive on nonzero polynomials, so division by H costs exactly `k*e`. Write `G=sum_s f_s(z)w^s`. The bound on support gives `deg f_s <= D-k*e-k*s`. Echelonize the span of the f_s by polynomial degree and write its degrees in decreasing order `delta_1>...>delta_ell`. Modulo polynomials of degree at most `D-k*e-k*j`, only the first j coefficient polynomials can contribute. Thus

`delta_(j+1) <= D-k*e-k*j` for `0<=j<ell`.

A monic degree-echelon basis has Wronskian leading coefficient, up to its ordering sign, `product_(a<b)(delta_b-delta_a)`, which is nonzero in characteristic zero. Hence

`deg W=sum_j delta_j-ell*(ell-1)/2`,

and summing the degree bounds proves (2.1). The bound on this degree budget itself is attained by the coefficient space spanned by `z^(D-k*e),z^(D-k*e-k),...,z^(D-k*e-k*(ell-1))`. This does **not** assert that the resulting interpolation criterion is optimal. ∎

### Proposition 2.2. Nonuniform jet certificate, with discarded nodes retained

For each ell define

`a_i(ell)=max(m_i-ell+1,0)`

and sort these M numbers in increasing order. If `M>E` and, for every

`0<=e<=E,  1<=ell<=E-e+1`,

one has

`sum_(j=1)^(M-e) a_(j)(ell) > B(D,k,e,ell)`,             (2.2)

then the matched-point Hermite map on `V_(D,k)` is injective.

**Proof.** Suppose F lies in its kernel. Apply Proposition 2.1, using the actual e and ell. At every nonzero specialization, G specializes to a nonzero element of its coefficient span with a zero of order at least m_i at alpha_i. Complete that element to a basis. Its Wronskian column is divisible by `(z-alpha_i)^max(m_i-ell+1,0)`. Changing basis changes W only by a nonzero constant. Distinct alpha coordinates therefore give

`sum_(i: F(z,beta_i)!=0) a_i(ell) <= deg W <= B(D,k,e,ell)`.

Deleting any e summands leaves a sum at least that of the smallest M-e summands, contradicting (2.2). ∎

### Uniform version

If every `m_i=m>E`, it suffices that

`M > max_(e,ell) { e + B(D,k,e,ell)/(m-ell+1) }`.        (2.3)

All denominators are positive. For `D=99,k=10,m=20`, (2.3) certifies **M=46** matched nodes for all choices of distinct coordinates, using 550 monomials. The corresponding unweighted maximal-degree Wronskian estimate certifies M=82. The source's simple choice uses roughly `2x` nodes; reducing their number by a constant does not improve its leading logarithmic prime-sum coefficient.

### Failure boundaries

- Pairwise distinct alpha coordinates cannot be dropped: with every alpha_i=0, the polynomial `z^m` has all requested jets zero whenever `m<=D`, however many beta values are used.
- Pairwise distinct beta coordinates cannot be dropped: with every beta_i=b, `w-b` has all z jets zero at every node.
- The characteristic-zero hypothesis can be replaced by `char(K)>D`, with ordinary jets interpreted in the range in which factorials are units. In unrestricted characteristic p, `1,z^p` have zero ordinary Wronskian although they are linearly independent over K.
- Distinct leading degrees modulo p suffice for the Vandermonde leading term to survive. Colliding degrees modulo p do not force the Wronskian to vanish: `1,z^p+z` has nonzero constant Wronskian. The lower coefficients remain relevant.
- The criterion is sufficient, not necessary, and it proves linear evaluation rank, not algebraic independence, global injectivity of a polynomial map, or a Jacobian-conjecture case.

### Typed Jacobian connection

Source: a bounded polynomial coefficient space. Target: a Wronskian polynomial. Map: choose a basis and take its derivative determinant. Preserved: linear independence and quantified total ramification when characteristic zero and the degree ledger apply. Lost: the original multivariable map, nonlinear fibres, branch entry, and global behavior at infinity. Sidecar: complete coefficient support, field characteristic, base factors, and coordinate matching. Decisive test: the inseparable pair `1,z^p`, and the repo's `THM-3882` base-divisor factor `W(hf,hg)=h^2 W(f,g)`.

## 3. Frobenius graph ideals and exact derivative loss

Let p be prime and `Q=z^p-w` in `Z[z,w]`. Let R be a commutative ring that is finite free over Z, with a fixed integral basis. Suppose the evaluation nodes obey

`alpha_i^p-beta_i in pR`.                                  (3.1)

No field or inertness assumption is needed for the following local statement. In the supplied quadratic setting, inert Frobenius is one mechanism producing (3.1).

Write `D_z^[a] D_w^[b]` for integral Hasse derivatives, equivalently the coefficient of `t^a u^b` after substituting `(z+t,w+u)`.

### Proposition 3.1. Jet-divisibility budgets

For all nonnegative r,s,h,a,b, evaluate `P=z^r Q^h w^s` at a node satisfying (3.1). Then

`partial_z^a partial_w^b P(alpha,beta) in p^max(h-b,0) R`,  (3.2)

and

`D_z^[a] D_w^[b] P(alpha,beta)
   in p^max(h-floor(a/p)-b,0) R`.                         (3.3)

**Proof of (3.2).** Every positive ordinary z derivative of Q has coefficient divisible by p, because its coefficient contains the factor p from `(p)_j`; derivatives above p are zero. The undifferentiated Q has p-divisible value. Thus every Q factor costs at least one p under arbitrarily many z derivatives. A positive w derivative can remove a Q factor, but total w order b can remove at most b such factors. Integral product-rule coefficients do not decrease divisibility.

**Proof of (3.3).** At a node,

`Q(alpha+t,beta+u)=p*rho + sum_(j=1)^(p-1) binom(p,j)alpha^(p-j)t^j + t^p-u`.

Every term except `t^p` and `-u` is p-divisible. To form the coefficient of `t^a u^b` in the hth power, at most `floor(a/p)+b` factors can use these exceptional terms. All other factors cost p. The outer factors `z^r w^s` and the Hasse product formula cannot increase the available a or b budget. ∎

This is the precise failure of silently normalizing the source's derivatives:

`partial_z^p Q=p!`, but `D_z^[p]Q=1`.

It also distinguishes direction: `partial_w Q=-1`.

The bounds are locally sharp whenever `h<p`: at `(alpha,beta)=(0,-p)`, the coefficient with `a=t*p`, `b=u`, `t+u<=h` equals

`(-1)^u multinomial(h;t,u,h-t-u) p^(h-t-u)`.

Its multinomial coefficient is a p-adic unit. This proves Hasse sharpness. For the corresponding ordinary derivative, `v_p((t*p)! u!)=t` because `t,u<=h<p`, so its valuation is exactly `h-u`, proving sharpness of (3.2) as well. For larger h, additional factorial carries can improve the bound; no universal equality is asserted there.

### Proposition 3.2. Same-minor tropical valuation certificate

Let `T={(r,s):r+k*s<=D}` with `1<=k<=p`. Use its monomials as rows and any collection of node/jet/integral-coefficient observations as columns. Set `h_(r,s)=floor(r/p)` and replace each row polynomial by

`z^(r-p*h) (z^p-w)^h w^s`.                               (3.4)

All terms remain in T, and every non-leading term has smaller z exponent. Thus (3.4) is an integral unitriangular row transformation, and it preserves **each square minor using all rows**, with the same selected columns.

Choose such a square minor D_0, with N rows and N columns. Assign a column defect

- `c_j=b_j` for ordinary mixed jets of orders `(a_j,b_j)`;
- `c_j=floor(a_j/p)+b_j` for Hasse mixed jets.

Sort the N row h-values and N column c-values in increasing order. Then

`v_p(D_0) >= sum_(i=1)^N max(h_(i)-c_(i),0)`.             (3.5)

Zero determinant has infinite valuation by convention. In the source's ordinary z-only case, all defects are zero, recovering `v_p(D_0)>=sum_T floor(r/p)`.

**Proof.** Proposition 3.1 makes every transformed matrix entry divisible by `p^max(h_i-c_j,0)`. Every determinant term therefore has valuation at least the corresponding assignment sum. The minimum assignment is the similarly sorted one. Indeed, for `h<=H` and `c<=C`,

`(h-c)_+ +(H-C)_+ <= (h-C)_+ +(H-c)_+`.

Repeatedly uncrossing an assignment proves (3.5). ∎

**Additional factorial ledger.** Ordinary derivatives equal `a!b!` times Hasse derivatives. Factoring these integers from columns gives the independent bound

`v_p(D_0) >= sum_j v_p(a_j! b_j!)
             + sum_i max(h_(i)-c_(i),0)`,                (3.6)

where now `c_j=floor(a_j/p)+b_j`. One may use the maximum of (3.5) in its ordinary version and (3.6). Normalization changes the arithmetic matrix and cannot be ignored.

### What this transfers and what remains open

The map sends a weighted monomial basis to an integral basis adapted to a Frobenius graph ideal. It preserves the selected determinant exactly and exposes divisibility. It forgets the size and nonvanishing of that determinant; Section 2 and an Archimedean height estimate are separate obligations. Prime-power congruences stronger than (3.1) cannot be inferred from a single mod-p Frobenius relation. Nor may one add the bounds for alternative row transformations at the same p unless a single simultaneous transformation or filtration proves that addition.

A useful next question is to optimize nonuniform jet acquisition jointly with (2.2) and (3.5), for a concrete arithmetic sampling family. The finite search can optimize a stated acquisition cost, but no globally optimal determinant architecture is claimed here.

## 4. Fixed prime-power templates for scalar additive cellular automata

Work in the commutative Laurent polynomial ring

`Z[X_1^±1,...,X_d^±1]`.

For a scalar local rule A, a finitely supported state S evolves by `S_t=A^t S`. This is the algebraic carrier of a translation-invariant additive cellular automaton. The finite-torus version is its quotient; identities survive quotienting, while support counts can change through aliasing.

### Lemma 4.1. Integral power lifting

In any commutative ring, if `U-V in p^e R` with `e>=1`, then

`U^(p^n)-V^(p^n) in p^(e+n)R` for every `n>=0`.           (4.1)

**Proof.** Write `U=V+p^e C`. In the pth-power binomial expansion, the terms with `1<=j<p` have p-valuation at least `1+e*j>=e+1`, and the last term has at least `e*p>=e+1`. Iterate. This includes p=2 without an extra assumption. ∎

### Theorem 4.2. Universal finite template

For every `m>=1` and `n>=m-1`, coefficientwise modulo `p^m`,

`A(X)^(p^n)
 = A(X^(p^(n-m+1)))^(p^(m-1))`.                         (4.2)

Here every spatial variable is dilated by the displayed exponent. Equivalently, once precision m is fixed, the complete rule kernel at every later prime-power time is a spatial dilation of the single row at time `p^(m-1)`.

**Proof.** Freshman's dream modulo p gives `A(X)^p=A(X^p) mod p`. Applying Lemma 4.1 yields

`A(X)^(p^r)=A(X^p)^(p^(r-1)) mod p^r` for `r>=1`.

For `n>=m`, apply this at r=n, then at r=n-1 after substituting `X^p`, and continue down to r=m. Every congruence is valid modulo `p^m`, and the final exponent is `p^(m-1)`. The case n=m-1 is an identity. ∎

### Theorem 4.3. Erasure of integral rule-lift data, with sharp time

If `A=B mod p^e`, `e>=1`, then

`A^(p^n)=B^(p^n) mod p^m` whenever `n>=max(m-e,0)`.       (4.3)

Thus the eventual template at precision m depends only on A modulo p. In particular two different lifts of the same finite-field local rule synchronize at all prime-power times `p^n` with `n>=m-1`.

This cutoff is uniformly sharp when `m>e`: take `A=1+p^e X`, `B=1`. At `n=m-e-1`, the X coefficient of `A^(p^n)-1` is exactly `p^(m-1)` and is nonzero modulo `p^m`.

**Proof.** Apply Lemma 4.1 in the Laurent polynomial ring. The coefficient calculation proves sharpness. ∎

This erases the *lift of the rule*, not the state: S is still present in `A^(p^n)S`. It does not erase initial-condition information, prove mixing, or make the evolution invertible.

### First correction layer

Define the integral polynomial

`C_p(A)=(A(X)^p-A(X^p))/p`.

For every n>=1,

`A(X)^(p^n)=A(X^(p^n))
              +p*C_p(A)(X^(p^(n-1))) mod p^2`.          (4.4)

It obeys the exact product law

`C_p(AB)=A(X^p)C_p(B)+B(X^p)C_p(A)+p C_p(A)C_p(B)`.

Formula (4.4) follows either from (4.2) at m=2 or directly from the defining identity. It identifies a concrete carry polynomial, rather than treating the mod-p zero error as zero at higher precision.

### Proposition 4.4. Pascal support and optimal template depth

For `A=1+X`, `n>=m-1`, the nonzero coefficients of `(1+X)^(p^n)` modulo `p^m` occur **exactly** at

`j = b*p^(n-m+1),  0<=b<=p^(m-1)`.                      (4.5)

There are `p^(m-1)+1` of them. At every interior index `0<j<p^n`,

`v_p binom(p^n,j)=n-v_p(j)`.                             (4.6)

The value at index j in (4.5) equals `binom(p^(m-1),b) mod p^m` by Theorem 4.2.

**Proof of (4.6).** Use

`binom(p^n,j)=(p^n/j) binom(p^n-1,j-1)`.

Modulo p, `(1+X)^(p^n-1)=(1+X^(p^n))/(1+X)` has every coefficient through degree `p^n-1` equal to `(-1)^j`, hence each is a unit. This proves (4.6); the two endpoints are units separately. Inequality `n-v_p(j)<m` is exactly (4.5). ∎

For m>=2, no universal template depth r<m-1 can replace m-1 in (4.2): a depth-r dilation would have support in multiples of `p^(n-r)`, but the surviving index `p^(n-m+1)` is not such a multiple. The coefficient at that index has valuation exactly m-1.

**Concrete rows.** Modulo 8, all rows at times `2^n`, n>=2, are dilations of

`(1+X)^4 = 1+4X+6X^2+4X^3+X^4 mod 8`.

Modulo 9, all rows at times `3^n`, n>=1, are dilations of

`(1+X)^3 = 1+3X+3X^2+X^3 mod 9`.

The intermediate terms are the precision sidecar erased by reduction modulo p.

### Seed content and separated replicas

Suppose S is nonzero modulo `p^m`, and let s<m be the minimum p-valuation of its nonzero coefficients. Put `q=m-s`. Multiplication by S converts a kernel identity modulo `p^q` into a state identity modulo `p^m`. For the Pascal rule, if `n>=q-1` and the spacing `p^(n-q+1)` is larger than the diameter of the one-dimensional seed support, the evolved state is a disjoint collection of exactly `p^(q-1)+1` **nonzero scaled seed packets**.

Each packet can have fewer occupied cells than S: a nonunit scalar can kill the higher-valuation entries of the seed. The term “replica” must retain this scaling and truncation. If packets overlap, their cells can cancel; on a periodic lattice, quotient aliasing can produce additional cancellation.

### Sharp boundaries: scalar, commuting, and linear

1. **Rule 30 is not covered.** With the single occupied cell at position 0 and the convention `new(x)=old(x-1) XOR (old(x) OR old(x+1))`, time-one support is `{-1,0,1}` and time-two support is `{-2,-1,2}`. Doubling time-one support gives `{-2,0,2}`, which differs. The nonlinear `center*right` term prevents a scalar convolution-power representation.
2. **Noncommuting lifts need not synchronize.** Modulo 4 take

`V=[[1,0],[0,0]], C=[[0,1],[0,0]], U=V+2C`.

Then `U=V mod2`, but `U^2=U`, `V^2=V`, so their difference stays `2C`, not zero modulo4. Commutativity is a real hypothesis in Lemma 4.1, not a cosmetic choice of notation.
3. **A monomial mod-p rule has a stronger collapse.** If `A=aX^d+pB`, then `A^(p^n)=a^(p^n)X^(d*p^n) mod p^(n+1)`. This is again Lemma 4.1. At fixed precision m and n>=m-1, the kernel is a single scaled shift. If a is a unit the shift is nonzero; if a is not a unit, nilpotence must be tracked separately.
4. **A nonlinear polynomial local rule is not a polynomial multiplier.** Formal use of the word polynomial does not identify these two categories.

### Typed transfer from the source

Source: `z^p-w in pR`, with ordinary differentiation preserving the congruence filtration. Target: scalar spatial substitution `A(X)^p-A(X^p) in pR`. Map: Frobenius error followed by the integral power-lifting lemma. Preserved: exact coefficient divisibility and the whole kernel. Lost by mod-p reduction: the correction polynomial and precision depth; restored by (4.2) or (4.4). Cheapest tests: `(1+X)^4 mod8`, the matrix example above, and the Rule30 singleton. This is a literal ring operation transfer, not an analogy between visual fractals.

## 5. Prime-race transfer: an exact formula and a no-go boundary

For a real Dirichlet character chi, define

`S_chi(x)=sum_(p<=x) chi(p)*log(p)/p`,

`C_chi(x)=sum_(p<=x) chi(p)`.

The exact partial-summation identity is

`C_chi(x) = x/log(x) * S_chi(x)
             - integral_2^x S_chi(t)*(log(t)-1)/(log(t))^2 dt`. (5.1)

The jump at 2 is included by using the left limit zero before the first prime. This weighted-to-unweighted map needs the whole lower-cutoff function, just as the repo's prime-power Möbius inverse does. A one-sided bound at a single cutoff is not enough; a one-sided bound for the full function still does not in general control the sign because (5.1) subtracts an integral.

### Proposition 5.1. Bounded weighted sums can coexist with arbitrarily negative count excursions

For every positive nonincreasing sequence `w_n -> 0`, there is a sign sequence `epsilon_n in {-1,1}` whose weighted partial sums are bounded, but whose unweighted partial sums are unbounded below.

**Proof.** Start with alternating signs, whose weighted partial sums are bounded by the alternating-series criterion and whose unweighted sums are bounded. For each j choose a later disjoint block of even length `2L_j`, with `L_j->infinity`, starting at an even baseline boundary and sufficiently late that `4L_j*w_(start)<=2^-j`. Replace the signs in its first half by all -1 and in its second half by all +1. Choose L_j even if necessary. Each completed replacement has zero total count, so the unweighted count resets to its bounded baseline after the block; halfway it is the baseline minus L_j. The change to every weighted partial sum caused by block j has absolute value at most `4L_j*w_(start)<=2^-j`. Summing these perturbations proves boundedness. ∎

For prime weights `w_n=log(p_n)/p_n`, discard the first finitely many terms to make the weights decreasing and apply the proposition. These signs are **not asserted to be Dirichlet-character values**. The example proves only that a weighted inequality alone, stripped of the character's multiplicative and analytic structure, cannot force the sign of a prime-count race.

The supplied prime-sum bound, if established with its hypotheses, controls an analytic functional relevant to real zeros. It neither removes the square correction in `prime_power_bias_20261014d.md` nor determines the next race crossover or a permanent winning class. A genuine new race consequence needs an explicit analytic remainder or local cutoff control, not just a new name for the weighted scalar.

## 6. Validation and open questions

The companion `verify_wronskian_frobenius.py` is pure standard-library integer/rational arithmetic. It checks weighted degree budgets, exact uniform threshold examples, modular full-rank witnesses under the sufficient certificate, Hasse and ordinary derivative valuations, exhaustive assignment-vs-sorting equality, prime-power kernel identities, lift erasure, binomial support and its failure at smaller depth, seed-content controls, noncommuting failure, and the Rule30 hostile. Its finite universes and counts are printed explicitly. Modular full row rank certifies rational full row rank; failure of a modular test would not by itself prove rational dependence.

The mathematical claims are the proved formulas above, not extrapolations from the finite tests. The script is a reproducibility and hostile-control companion.

Concrete next questions, left OPEN here:

- Jointly optimize per-node jet depths, weighted support, and the sorted valuation loss in a specified arithmetic family. What is the optimal tradeoff when derivative acquisition has a fixed total cost?
- The scalar depth question is closed in the companion `SCALAR_TEMPLATE_DEPTH_CLASSIFICATION.md`: every normalized nonmonomial has full depth, and a Frobenius-compatible lift modulo p^m exists exactly for support in one residue class modulo p^(m-1). Classification of finer prime-to-p spacing and all coefficient patterns remains separate.
- For matrix-valued CA rules, characterize a commutative coefficient subalgebra or a controlled commutator filtration sufficient to recover a precision-depth law. The explicit 2x2 example rules out an unrestricted extension.
- For a specific prime-race application, identify which local analytic estimate supplies the information lost by passing to the weighted sum.

None of these questions is promoted merely because its low-dimensional tests pass.

## 7. Integral divided-power duality: the exact common module

A parallel session establishes a source-splitting decomposition of the integer-valued function lattice on unlabelled tournaments. On each primitive string, its divided vertex-deletion operator acts as multiplication by `x^[m]` in the integral divided-power algebra `Gamma_Z(x)`, where

`x^[a] x^[b]=binom(a+b,a)x^[a+b]`.

This yields an exact bridge, independent of visual or terminological analogy. Pair `Gamma_Z(x)` with `Z[X]` by `<x^[j],X^k>=delta_(j,k)`. The adjoint of multiplication by `x^[m]` is the Hasse derivative, because

`<x^[m]x^[j], X^(j+m)> = binom(j+m,m)
 = <x^[j], D_X^[m] X^(j+m)>`.

Consequently the adjoint of multiplication by `m! x^[m]` is the ordinary derivative `partial_X^m`. On every finite integral graded piece, matrix transposition preserves all determinantal divisors and Smith factors. The factor `m!` is therefore a genuine lattice-index change, exactly the factor retained in (3.6). The valuation of the divided-power multiplication entry `binom(j+m,m)` is the base-p carry count for adding j and m; it is not automatically the floor-loss bound (3.3).

The target/source bases and the integral pairing are the preserved sidecars. A Frobenius graph ideal, arithmetic evaluation nodes, and the weighted support filtration are additional structures in this note, not structures supplied by the tournament transpose. In particular the bridge does not transfer a prime-sum estimate to tournaments. The theorem on tournament primitives is recorded in the companion session's `PADDING_FROBENIUS_PROOF.md`, §5; the abstract divided-power/Hasse duality above is proved here directly.
