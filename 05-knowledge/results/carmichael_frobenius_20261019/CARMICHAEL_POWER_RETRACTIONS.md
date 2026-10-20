# Carmichael power maps: exact retraction clocks, arithmetic jets, and hidden sheets

**Status: PROVED within the stated elementary scopes; FINITE-EXACT companion controls.** This packet does not establish a new Carmichael-counting estimate, a Jacobian conjecture result, or a new priority claim. Teichmüller representatives and simple-root lifting are classical. The contribution is a self-contained synthesis and an exact iteration formula for monomial thickenings, including coefficient torsion.

## 0. Inheritance, sources, and the question

The supplied Howe–Huffman–Kohlmeyer–Kruer paper, *The exponent of growth of the Carmichael counting function* (2026-10-08), uses in Proposition 5.2 the chain

    p−1 | L | n−1, with n a squarefree product of at least two selected primes.

This is the only input needed here. Its analytic inputs and its claimed counting theorem are not dependencies of this packet. The chain says that every residue is fixed by the nth-power operation on Z/n. What survives if one retains infinitesimals, higher arithmetic precision, or extension-field points?

Closest previous mechanisms:

- `05-knowledge/results/wronskian_frobenius_20261019/unique-lift-addendum/UNIQUE_LIFTS_AND_EXACT_COUNTS.md`: uniqueness and multiplicativity of compatible scalar coefficient lifts, with a nonadditive carry.
- The same packet's `SCALAR_TEMPLATE_DEPTH_CLASSIFICATION.md` and `FINITE_LOGARITHM_BRIDGE.md`: coefficient-level precision, support obstructions, and the loss incurred by scalar evaluation.
- The same packet's `WRONSKIAN_FROBENIUS_TRANSFER.md`: ordinary versus Hasse jet budgets.
- `01-canon/theorems/THM-3548-planar-keller-conductance-shadow-gates.md` and the current planar-Jacobian guardrails: finite observables do not replace a polynomial or geometric condition.
- Incoming `THM-4718-galois-group-of-the-problem-of-apollonius-is-the-parity-cube.md`: its full-group conclusion is explicitly under audit and is not used. The present example independently shows how rational-point data can miss geometric sheets.

Classical attribution: William Stein's [Hensel and Teichmüller exposition](https://wstein.org/books/bsd/bsd/node10.html) identifies the unique root-of-unity lift of a nonzero residue. The [Stacks Project's formally étale lifting definition](https://stacks.math.columbia.edu/tag/00UP) supplies the standard geometric language. The proofs below are given directly; no originality claim is made for Hensel lifting, Teichmüller sections, LTE, or the binomial valuation identities.

Related primary literature already studies these broader objects: Bae–Hu–Sha, [*On the Carmichael rings, Carmichael ideals and Carmichael polynomials*](https://arxiv.org/abs/1809.05432), and Burnette, [*On power maps over weakly periodic rings*](https://arxiv.org/abs/2207.14283). The former studies the identity power map with exponent equal to ring cardinality; the latter studies power-map periodicity and finite-ring examples. These references are context and priority guardrails, not unproved inputs to Theorem 2. No exhaustive novelty audit is claimed.

Live concept board: group exponents; power-map contraction; coefficient precision; Hasse jets; finite-field extensions; polynomial-versus-point evaluation. The hostile controls are nonsquarefree moduli, a nonzero mixed binomial coefficient, a surviving pth jet, nonunit inputs, and extension fields on which a pointwise identity fails.

## 1. The elementary nilpotent contraction

### Theorem 1: a canonical multiplicative section

Let A be a commutative unital ring, I an ideal with I^k=0, and e≥2 an integer. Suppose

1. e·1_A belongs to I;
2. b^e=b for every b in B=A/I.

Write F(x)=x^e. Then F^(k−1) is independent of the choice of lift of a residue in B, takes values in Fix(F), and induces a bijection

    τ : B → Fix(F) ⊂ A

inverse to reduction. It preserves 0, 1, and multiplication. It need not preserve addition. The construction commutes with homomorphisms of such rings preserving the ideals and exponent.

**Proof.** If x−y belongs to I^j, j≥1, expand (y+(x−y))^e. The linear term belongs to I^(j+1), because e belongs to I. Every higher term belongs to I^(2j), hence to I^(j+1). Thus

    x=y mod I^j  ⇒  F(x)=F(y) mod I^(j+1).               (1)

The assumption on B gives F(x)−x in I. Applying (1) repeatedly gives F^k(x)=F^(k−1)(x). Applying (1) to two lifts proves independence. Two fixed lifts congruent modulo I must be congruent modulo every I^j and therefore equal. Multiplicativity follows from F(xy)=F(x)F(y), and functoriality follows from the power formula. ∎

This is a set/monoid retraction, not an assertion that the quotient map splits as a ring map. The polynomial defining the fixed locus is X^e−X, whose derivative is eX^(e−1)−1, a unit because e is nilpotent. Thus its fixed locus is finite étale of degree e over A, even though the power map itself need not be étale. A-valued fixed points and geometric points must still be counted separately.

### Corollary 1.1: every Carmichael residue has a canonical lift

Let n be a squarefree composite integer with p−1 dividing n−1 for every p|n. Such n is odd. For every k≥1,

    τ_(n,k)(a) = a^(n^(k−1)) mod n^k                      (2)

is the unique solution t of t^n=t modulo n^k lifting a modulo n. Consequently there are exactly n such solutions. Formula (2) is independent of the integer representative a, multiplicative, and compatible under reduction from k+1 to k.

At each p|n its component is the classical Teichmüller lift, including 0: a nonzero component satisfies t^(p−1)=1. Indeed those p simple roots satisfy t^n=t, and uniqueness identifies them with the fixed lifts. This also explains why the limiting lift is independent of the chosen admissible power exponent.

For k=2 the formula is particularly simple: a↦a^n mod n². It is never additive for k≥2: an additive multiplicative section has τ(1)=1, so nτ(1)=τ(0)=0, whereas n is nonzero in Z/n^k.

**Boundary.** The residue identity holds for n>1 exactly when n is prime or squarefree with every p−1|n−1. If p²|n, the residue p violates it because v_p(p^n−p)=1. For a squarefree modulus, the assertion follows in both directions from the cyclic groups F_p^×. Nothing here turns a composite-exponent polynomial into a Frobenius ring endomorphism.

## 2. Exact arithmetic timing, including zero-divisors

Fix p|n and work modulo p^k. Let t be the Teichmüller lift of the residue of x.

- If p does not divide x, put s=min(k,v_p(x−t)). Then the exact first hitting time is k−s. More precisely, before saturation,

      v_p(x^(n^r)−t) = s+r.                              (3)

  To prove this, write x=t(1+u) with v_p(u)=s and use the odd-prime lifting-the-exponent identity, with v_p(n)=1. The value s=k means the point is already fixed.
- If p divides x, put s=min(k,v_p(x)). Then t=0 and the exact first hitting time is

      min{r≥0 : n^r s≥k}.                                (4)

  This follows directly from v_p(x^(n^r))=n^r v_p(x), with truncation at k. In particular x=0 has hitting time zero.

The global hitting time modulo n^k is the maximum of the local times (3) and (4). The uniform bound k−1 is sharp for every k≥2: x=1+n has t=1 at every p and valuation exactly r+1 after r iterations. Zero-residue components typically disappear much faster; it is incorrect to assign their valuations the unit rule s+r.

The odd-prime qualification is real. For exponent 2 modulo 2^k, extra 2-adic cancellation changes the sharp unit timing. The general upper bound in Theorem 1 still holds, but the stated equality witness is for odd Carmichael n.

## 3. The exact clock on arithmetic monomial thickenings

The next theorem retains arithmetic precision and all surviving jet coefficients simultaneously. The polynomial quotient is a concrete input, not merely an evaluated function.

### Definition: coefficient-precision monomial rings

Fix an odd prime p and finitely many variables T_1,…,T_d. For every exponent vector γ in N^d choose an integer w_γ≥0 such that:

- w_0=a≥1;
- γ≤δ coordinatewise implies w_γ≥w_δ;
- only finitely many w_γ are positive.

Let A_w be the quotient of Z[T_1,…,T_d] by the monomial relations

    p^(w_γ) T^γ = 0  for every γ.                        (5)

The case w_γ=0 means T^γ=0. Equivalently, an element has unique coefficients at T^γ modulo p^(w_γ). This includes flat truncated polynomial rings and nonflat thickenings such as pT=0. Reduction of the constant coefficient gives A_w→F_p, with nilpotent kernel.

Define the finite axis maximum

    H = max { h+w_(p^h e_i) : w_(p^h e_i)>0 },             (6)

where e_i is the ith coordinate vector; set H=0 if every variable is zero. Let e≥2 satisfy p|e and e≡1 mod p−1, and put s=v_p(e).

### Theorem 2: the axis-precision formula

For F(x)=x^e on A_w, the least r such that F^r sends every element to the constant Teichmüller lift of its residue is exactly

    R = ceil( max{a−1,H} / s ).                          (7)

Thus all mixed monomial relations are irrelevant to this worst-case clock once the constant precision and p-power axis precisions are fixed.

**Upper bound.** Let E=e^r and v=v_p(E)=sr. Choose any integral polynomial representative P and write P^E=Σ c_γ T^γ. The integral Euler-derivative identity

    T_i ∂_i(P^E) = E P^(E−1) T_i ∂_i P

gives, for every γ_i>0,

    v_p(c_γ) ≥ v−v_p(γ_i).                               (8)

For a surviving nonconstant γ choose i with h=v_p(γ_i) minimal among its positive coordinates. The axis exponent p^h e_i is coordinatewise at most γ, so w_(p^h e_i)≥w_γ. Hence H≥h+w_γ. If v≥H, (8) kills c_γ modulo p^(w_γ). Every nonconstant coefficient is therefore zero.

The remaining constant is P(0)^E modulo p^a. If its residue is nonzero, write it as its Teichmüller lift times a principal unit; raising to E kills the principal unit when v≥a−1. If its residue is zero, E≥p^v≥a when v≥a−1, so its E-th power is zero. This proves the bound, including a=1 and r=0 when the ring has no nonzero variables.

**Sharpness of the constant constraint.** For p odd, the element 1+p has valuation v_p((1+p)^(e^r)−1)=1+sr. If sr<a−1, it has not reached 1 modulo p^a.

**Sharpness of the axis constraint.** We use the exact elementary identity

    v_p(binomial(E,j)) = v_p(E)−v_p(j),
    1≤j≤p^(v_p(E)).                                    (9)

Indeed j·binomial(E,j)=E·binomial(E−1,j−1), and the latter binomial coefficient is a p-adic unit: in its product formula all numerator factors E−i, 1≤i<j, have the same valuation as i.

If sr<H, choose i,h with h+w_(p^h e_i)>sr and let j=p^min(h,sr). In (1+T_i)^E the coefficient of T_i^j has valuation max(0,sr−h) by (9). If h≤sr this is smaller than w_(p^h e_i); if h>sr it is zero and the coefficient precision at T_i^j is at least w_(p^h e_i)>0. In either case the coefficient survives. Distinct monomials cannot cancel in (5). This proves (7). ∎

### Corollary 2.1: flat one-variable and multivariable clocks

For A=(Z/p^a)[T]/(T^b), b≥2, one has

    R = ceil( (a+floor(log_p(b−1))) / v_p(e) ).           (10)

For an ideal J generated by monomials with unit coefficients over Z/p^a, with finite quotient, suppose T_i^(b_i) is the least pure power of variable i in J. Ignore variables with b_i=1. If at least one variable survives, the clock is

    R = ceil( (a+max_i floor(log_p(b_i−1))) / v_p(e) ).    (11)

The proof uses only the longest surviving axis at a p-power degree. A large total-degree nilpotence index need not imply a slow power-map clock.

### Corollary 2.2: the Carmichael arithmetic/jet clock

Let p_min be the least prime divisor of a Carmichael n. On

    A = (Z/n^a)[T]/(T^b),  a≥1, b≥2,

the exact worst-case iteration count is

    R = a+floor(log_(p_min)(b−1)).                        (12)

The image is the n-element set of constant CRT Teichmüller lifts. The witness 1+T is sharp at the p_min component. If b=1, the clock is instead a−1. At a=1 the terminal retraction is the additive constant-term projection; at a≥2 it is not additive, by the characteristic obstruction in Corollary 1.1.

For comparison, the ideal (n,T) has nilpotence index a+b−1, so Theorem 1 gives the weaker bound a+b−2. Formula (12) replaces the linear jet-length contribution by a logarithmic one.

### Nonflat example

For p=3, e=561, and

    A=(Z/27)[T]/(T^4,3T),

the coefficient precisions are (3,1,1,1). Here H=max(1,1+1)=2 and a−1=2, so R=2. The flat ring (Z/27)[T]/(T^4) instead has R=4. Thus the ideal's thickness alone loses the coefficient-torsion coordinate; the axis precision formula retains it.

**Scope boundary.** The exact theorem is for commutative monomial coefficient-precision quotients, with p odd. Nonmonomial relations or a noncommuting coefficient ring require a separate cancellation analysis. A quotient of a stated flat truncated ring inherits its upper bound, but need not inherit sharpness. This is the next precise extension question, not a proved general classification.

## 4. Why first-order tangent information is insufficient

Put R=Z/n. For every a,b in R,

    (a+bε)^n=a in R[ε]/(ε²).                             (13)

So the power operation is the identity on R-points but kills every tangent vector. The two polynomials X and X^n induce the same function on R, while their derivatives are 1 and 0. Differentiation does not descend through the quotient from polynomials to evaluated functions: the null polynomial X^n−X has derivative −1.

The first local Hasse jet restores missing information. For p|n write m=n/p; p does not divide m. In characteristic p,

    (a+T)^n = (a^p+T^p)^m.

For a in F_p^×, all positive coefficients below degree p vanish, and

    [T^p](a+T)^n = m mod p.                              (14)

Indeed m≡1 mod p−1, so a^(n−p)=1. At a=0 the first term instead has degree n. The distinction between unit and zero base points matters here too.

The operation on F_p[T]/(T^b) is constant-term projection for b≤p. At b=p+1 it fails additivity:

    (1+T)^n = 1+mT^p,
    1^n+T^n = 1.                                        (15)

The second equality uses n>p, which is automatic for composite squarefree n. For exponent e=p, Frobenius remains additive and this hostile is unavailable. Ordinary derivatives in characteristic p miss the surviving pth jet; Hasse coefficients see it immediately.

This supplies the exact loss ledger: reduced-point evaluation preserves the residue identity and multiplicative exponent condition; it discards the differential, the first surviving Hasse degree, and the algebraic degree of the power map. The necessary sidecars are nilpotent test algebras, coefficient precision, and extension-field points.

## 5. Hidden separable sheets and a transverse fixed locus

For each p|n, the algebraic map X↦X^n over F_p factors as the pth Frobenius followed by X↦X^m, m=n/p. Its total degree is n, inseparable degree p, and separable degree m. The whole map is **not purely inseparable** when n is composite.

On G_m, its kernel group scheme decomposes

    μ_n ≅ μ_p × μ_m.

The μ_p part is infinitesimal in characteristic p; μ_m is étale. The latter has only the identity as an F_p-point because m≡1 mod p−1. Thus both infinitesimal and extension-field information can be invisible in the base-field point table.

For q=p^f, cyclicity of F_q^× gives exact tests:

    #Fix(x↦x^n on F_q) = 1+gcd(n−1,q−1);                 (16)
    the map is identity ⇔ q−1 | n−1;                     (17)
    the map is a permutation ⇔ gcd(n,q−1)=1.             (18)

Identity and bijectivity are different. The first extension on which noninjectivity appears has degree

    min_{prime ℓ | n/p} ord_ℓ(p).                        (19)

This is greater than one for a Carmichael n, since p−1|n−1 precludes any ℓ|n from also dividing p−1.

For n=561:

- over F_9 the map is identity;
- over F_27 it is a permutation with exactly 3 fixed points;
- over F_121 it has 41 image points and 41 fixed points;
- the first noninjective extension degrees at p=3,11,17 are respectively 5,2,2.

Meanwhile the fixed polynomial X^n−X has derivative −1 over every characteristic p dividing n. It has n distinct geometric roots, even though only p roots are F_p-rational. Its fixed scheme is étale; the power map is inseparable. These are different statements. For example, X^561−X has 561 geometric roots in characteristic 3 but only 3 F_3-rational roots.

The planar-Jacobian relevance is a warning about transport, not a theorem about characteristic zero: a pointwise permutation or identity over a finite field does not certify a polynomial inverse, an invertible Jacobian, or a geometric degree. Likewise a finite-field point census alone cannot establish a full generic Galois group.

## 6. The composite carry and the old scalar-lift bridge

At precision n² let

    c_n(a,b) = ((a+b)^n−a^n−b^n)/n mod n.                 (20)

The numerator is divisible by n for integer inputs. The result depends only on a,b modulo n: shifting an input by n changes its nth power by a multiple of n². It measures the additive failure of the canonical section. It is symmetric, normalized at zero, and satisfies the cocycle identity

    c_n(a,b)+c_n(a+b,c)=c_n(b,c)+c_n(a,b+c).              (21)

For p|n, the exact local normalization is

    c_n(a,b) mod p = (n/p)^(-1) c_p(a mod p,b mod p),      (22)

where

    c_p(X,Y)=Σ_(j=1)^(p−1) [binomial(p,j)/p] X^j Y^(p−j)

is the classical integral pth-power carry polynomial reduced modulo p.

**Proof of (22).** Both a^n and a^p modulo p² are the unique Teichmüller lift of a modulo p, so their differences are multiples of p². Divide the three-term difference by p, then note that division by n is multiplication by (n/p)^(-1) after reducing modulo p. ∎

A crucial boundary is that the formal expression in (20) is generally not an integral polynomial. For n=561,

    binomial(561,3) mod 561 = 187,

so its coefficient divided by 561 is not integral. Integer-valued divisibility after evaluation is weaker than coefficientwise divisibility. One can represent the resulting function over Z/n by CRT-combining the local polynomials in (22), but that does not make the displayed rational-coefficient polynomial integral.

This is the direct bridge to the previous scalar Frobenius-lift packet. Both constructions are unique multiplicative lifts with additive carries. The earlier equation A(X)^p=A(X^p) is a coefficient-polynomial identity with a spatial-substitution sidecar and a support obstruction. The present equation is t^n=t in specified finite rings. Scalar evaluation relates their carries; it does not remove the earlier support obstruction or transport a full polynomial identity.

## 7. Verification, proof boundaries, and next question

Run from this directory:

    python3 -B verify_power_retractions.py > /tmp/power-retractions.json
    python3 -B -O verify_power_retractions.py > /tmp/power-retractions-optimized.json
    cmp verification.json /tmp/power-retractions.json
    cmp verification.optimized.json /tmp/power-retractions-optimized.json
    python3 -B verify_axis_precision.py > /tmp/axis-precision.json
    python3 -B -O verify_axis_precision.py > /tmp/axis-precision-optimized.json
    cmp axis-verification.json /tmp/axis-precision.json
    cmp axis-verification.optimized.json /tmp/axis-precision-optimized.json

The first verifier uses independent Hensel-digit/CRT and modular-power paths, exact truncated-polynomial convolution, binomial valuations, local hitting times including nonunits, all inputs in Z/561², all 561² carry pairs, and actual finite-field models. The second verifier enumerates coefficient-torsion monomial rings using their normal forms and compares measured first hitting times with (7). Neither uses removable Python assertions.

These computations support the proofs, not an extrapolation to untested classes. The next precise question is which nonmonomial mixed-characteristic quotients retain an axis-determined clock, and which relations introduce earlier cancellation. The monomial result already shows that total nilpotence index is too coarse: one must retain where p-adic coefficient precision is available.
