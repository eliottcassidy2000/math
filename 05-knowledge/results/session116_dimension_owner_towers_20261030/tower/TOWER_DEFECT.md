# Doubling defects, complete sedenion singularity, and arithmetic bad primes

**Status: PROVED, with finite-exact reproducibility checks.** This is a scoped synthesis/refinement, not a claim to have discovered sedenion zero divisors or settled a major open problem. The full real 16-dimensional left-multiplication singularity classification is proved below. The explicit integer family is complete only in its stated subspace. No complete higher-stage zero-divisor classification is claimed.

## Inheritance, convention, and corrections

Use real Cayley–Dickson algebras A_0=R, A_(n+1)=A_n⊕A_n with

    (a,b)(c,d)=(ac−bar(d)b, da+b bar(c)),
    bar(a,b)=(bar(a),−b),       N(a,b)=N(a)+N(b).

N is the **squared** Euclidean norm; < , > is its real inner product. Our associator is [u,v,w]=(uv)w−u(vw). The order and parentheses are essential. This convention is explicitly used in [Geometry of the Octonions, sedenions](https://sites.science.oregonstate.edu/coursewikis/GO/book/go/sedenions.html); Baez's [Cayley–Dickson discussion](https://math.ucr.edu/home/baez/octonions/node5.html) uses a different but isomorphic convention, so its coordinate formulas must not be copied without transport.

The inherited `session112_counting_interfaces_20261030/norms/NORM_SIGNS.md` proves the rectangle criterion and a five-rectangle 0=1 certificate prohibiting Euclidean norm composition for any normalized signed-XOR table at dimension≥16. Its coherent eight-dimensional gauges preserve composition, but cannot fix the dimension16 obstruction. The present note quantifies the canonical doubling failure rather than contradicting that result. Session115 odd-eighth Boolean gauges are not themselves dimension16 algebras; no numerical identification is imposed.

Historical `cayley_dickson_tower_s136.out` incorrectly says power-associativity disappears at32. In fact, at every stage xx*=x*x=N(x)1, x+x*=2Re(x)1, and

    x²−2Re(x)x+N(x)1=0.

The real span of 1,x is an associative subalgebra (its multiplication is the quotient by this quadratic relation). Hence all powers of one element associate at every stage. Every nonzero x has the two-sided inverse x*/N(x); this does NOT imply cancellation in a nonassociative algebra. Terminology sometimes reserves “invertible” for invertibility of the multiplication operators; here “two-sided inverse” means only the two displayed product equations.

For the classical full zero-divisor cone see [Reggiani, The geometry of sedenion zero divisors, arXiv:2411.18881](https://arxiv.org/abs/2411.18881). Its distinction matters: normalized zero-divisor **elements** form a Stiefel manifold; normalized zero-product **pairs** form the G2 space. We do not claim these known classifications as new.

## 1. Uniform recursive defect and the first failure

Put D_n(u,v)=N(uv)−N(u)N(v). Expansion, valid at every stage, gives

    D_(n+1)((a,b),(c,d))
      =D_n(a,c)+D_n(bar d,b)+D_n(d,a)+D_n(b,bar c)
       +2(<da,b bar c>−<ac,bar d b>).                    (1)

Thus there are four inherited defects and one cross-term defect. Below dimension16 the inherited terms vanish. For a,b,c,d∈O, the composition-algebra adjoint identities

    <uv,w>=<v,bar u w>=<u,w bar v>

turn the last difference into <(da)c−d(ac),b>. Consequently

    D_4((a,b),(c,d)) = 2<[d,a,c],b>.                   (2)

This explains precisely why nonassociativity of the lower stage first becomes a norm defect in its double. For H the associator is zero, so doubling to O still composes the norm. For O the associator can be nonzero. Formula (2) is NOT asserted unchanged at32 and above; use (1), which retains the inherited defects.

## 2. Complete singular-value classification at dimension16

For x=(a,b)∈S, write A=Im(a), B=Im(b), and

    N=N(a)+N(b),   Δ=N(A)N(B)−<A,B>² ≥0.

**Theorem.** For real left multiplication L_x:S→S, the eigenvalues of L_x*L_x, counted with multiplicity, are

    N                   (8 times),
    N+2 sqrt(Δ)         (4 times),
    N−2 sqrt(Δ)         (4 times).                     (3)

When Δ=0 these coincide; the multiset statement remains valid. In particular,

    det L_x = N⁴ (N²−4Δ)².                           (4)

For nonzero x, the following are equivalent:

* x has a nontrivial right annihilator, i.e. xy=0 for some y≠0;
* L_x is singular;
* a,b are pure imaginary, have equal nonzero norms, and are orthogonal.

In that case dim_R ker L_x=4.

**Proof.** Any two octonions a,b lie in a quaternion subalgebra H (if the generated associative subalgebra is smaller, extend it). Choose a unit ℓ perpendicular to H. The usual octonion quaternion splitting O=H⊕Hℓ has the displayed Cayley–Dickson convention. Both H⊕H and Hℓ⊕Hℓ are invariant under L_x. On H⊕H, this is octonion multiplication, so the Gram operator equals N I_8.

Write c=uℓ,d=vℓ in the second block, u,v∈H. The product there is

    (u,v) ↦ (ua+v bar b, v bar a−ub).                 (5)

Quaternion associativity and adjoints give its Gram operator

    (u,v) ↦ (Nu+v[b,a], Nv+u[a,b]),                  (6)

where [a,b]=ab−ba is the commutator, not the three-input associator. Let K=[a,b]=[A,B]. It is imaginary and N(K)=4Δ. The off-diagonal part T(u,v)=(−vK,uK) is self-adjoint, T²=4Δ I_8, and tr(T)=0. Thus its eigenvalues are ±2sqrt(Δ), each four times if Δ>0, and all zero if Δ=0. This proves (3).

The squared determinant is N^8(N²−4Δ)^4. Both det(L_x) and Q=N⁴(N²−4Δ)² are real polynomials in the16 coordinates. Since (det−Q)(det+Q)=0 identically in an integral domain, one factor vanishes identically; x=1 fixes det=Q. This proves (4) without a sign ambiguity.

Finally

    2sqrt(Δ) ≤ 2sqrt(N(A)N(B)) ≤ N(A)+N(B) ≤ N.

Equality throughout holds precisely when <A,B>=0, N(A)=N(B), Re(a)=Re(b)=0. If x≠0 the common imaginary norm is positive, and (3) gives exactly four zero eigenvalues. This proves completeness.

**All annihilators, not just their dimension.** In the singular case let r=N(a)=N(b)>0 and H=span_R{1,a,b,ab}. Then

    ker L_(a,b) = { (c, b(ca)/r) : c∈H^perp }.       (7)

Indeed put c=uℓ; formula (5) says v=−u ab/r, and the second equation follows from ab=−ba, a²=b²=−r. Translating back gives b(ca)=−(uab)ℓ. The H⊕H component is zero because its Gram operator is 2r I. Parentheses in b(ca) may not be dropped.

## 3. Integer matrices: the exact bad-prime set

In the standard integral CD basis, integer x gives an integer16×16 multiplication matrix. Formula (4) is an identity over Z, so for every prime p, including2,

    L_x mod p is singular  iff  p divides N(N²−4Δ).   (8)

This is a genuine arithmetic transfer: source is the real block/spectral calculation; map is its integer determinant polynomial; preserved predicate is singularity after reduction; the lost structure is positive definiteness. Consequently the real cone criterion “pure, equal, orthogonal” must NOT be substituted for (8) over finite fields.

For example x=(1,i) has N=2, Δ=0 and is nonsingular over R, but L_x is singular modulo2. Higher-stage embedding preserves exact zero products, but formula (4) is a dimension16 result, not an all-dimensional determinant formula.

## 4. Complete zero-product relation in a fixed arithmetic subspace

Now fix the standard H⊂O=H⊕Hℓ, and the8-dimensional linear subspace

    W={ (a,sℓ) : a,s∈H } ⊂ S.

This W is an **eight-dimensional subalgebra**, but is NOT the canonical octonion subalgebra and does not compose the Euclidean norm. Thus dimension eight alone never guarantees the octonion laws: the multiplication tensor matters. For x=(a,sℓ), y=(c,tℓ), direct multiplication gives

    xy=(ac−bar s t, (t bar a+sc)ℓ).                 (9)

For nonzero x,y∈W, equation xy=0 holds exactly when

    a,c are nonzero pure imaginary quaternions,
    <a,c>=0,  N(s)=N(a),  t=sac/N(a).               (10)

**Proof.** If any of a,s,c,t is zero, (9) and quaternion division force one entire factor to vanish. Otherwise the first equation yields t=sac/N(s). Substitution in the second yields ac bar a=−N(s)c. Taking quaternion norms gives N(a)=N(s); multiplying on the right by a then gives ac=−ca. Writing scalar and vector parts shows that two nonzero anticommuting quaternions are precisely two orthogonal imaginary quaternions: if a's scalar part were nonzero, the vector equation would make c's vector part a scalar multiple of a's and the scalar equation would force c=0. Therefore a,c have zero scalar parts and orthogonal vector parts. All steps reverse, proving necessity and sufficiency.

For fixed eligible x∈W, its annihilator **inside W** has real dimension2 (c lies in the plane perpendicular to a in Im H). The full annihilator in S has dimension4 by §2; confusing these dimensions would lose half the solutions.

**Integral construction without divisibility assumptions.** Choose nonzero a∈Z³⊂Im H and s∈Z⁴ with N(s)=N(a)=r; choose nonzero c∈Z³ with a·c=0. Then

    x=(a,sℓ),    y=(r c,(sac)ℓ)                    (11)

are nonzero integral sedenions with xy=0. Such s always exists, e.g. s=a. This is an unlimited family indexed by an ordinary orthogonal lattice and a four-square shell. Formula (10), not (11), is the complete rational classification in W; (11) only clears the denominator uniformly.

## 5. Explicit failed cancellation and stable higher-stage witnesses

With quaternion ij=k, take

    x=(i,ℓ),  y=(j,kℓ).

Both have N=2; xy=0, so their norm defect is−4. Yet x^−1=−x/2 satisfies x^−1 x=x x^−1=1. The associator gives exactly

    [x^−1,x,y]=(x^−1 x)y−x^−1(xy)=y.               (12)

Thus the missing cancellation is measured by an explicit nonzero associator, not by loss of the conjugation formula. The embeddings z↦(z,0) preserve multiplication and Euclidean norm at every doubling; (11) and (12) therefore persist at every dimension2^n≥16. This does not classify all new zero divisors at later stages.

## Checks and failure controls

Run from this directory:

    python3 check_tower.py
    python3 -O check_tower.py
    python3 check_spectrum.py
    python3 -O check_spectrum.py

`check_tower.py` verifies4096 exact cross-term coefficients; all signed two-term input pairs at dimensions4,8,16; the entire9^4 signed-basis-plus-zero quaternion factor bank;2064 integral cone pairs; and fixed witnesses through128. It checks degenerate factors as well as nonzero ones, exact quaternion reduction, inverse equations, and the cancellation associator. `check_spectrum.py` independently constructs16D multiplication matrices, checks the Gram cubic and traces, computes determinants by exact Bareiss elimination, and compares six-prime ranks by independent row reduction. The finite tests are corroboration; the displayed arguments prove the real/rational families. Both scripts use explicit exception gates, never Python `assert`.

## 6. Exact shell counts and a genuine modulo-eight obstruction

For m≥1, let r_k(m) count ordered signed integer representations as a sum of k squares. The complete cone criterion in §2, restricted to W, gives an exact bijection

    {integer zero-divisor elements x∈W : N(x)=2m}
       ↔ {a∈Z³:N(a)=m} × {s∈Z⁴:N(s)=m}.

Hence the count is r_3(m)r_4(m). By the classical Jacobi four-square theorem and Legendre three-square theorem, this is

    8 r_3(m) sum_(d|m, 4∤d) d,

and is positive exactly when m is NOT of the form 4^j(8k+7). These two classical number-theory theorems are external dependencies; the new transfer here is the explicit shell bijection, not a new proof of those theorems. Primary references matching these exact classical statements are [Doležálek, introduction, printed page2](https://arxiv.org/pdf/2409.16970) for Jacobi and four-square existence, and [Danilkin–Chevalier, Three Squares Theorem](https://isa-afp.org/entries/Three_Squares.html) for Legendre.

For instance the norm-squared2 shell has48 elements, whereas the norm-squared14 shell has none **inside W**. This does not say there are no integer zero divisors of that norm in the entire16-dimensional algebra. The modulo-eight obstruction is a true representation-theoretic arithmetic predicate; it is not an identification with the Session115 odd-eighth phase probabilities.

`check_integer_shells.py` and its optimized run independently enumerate ordered signed three- and four-square representations through m=128, verify both classical formulas in that finite range and record all shell counts in `integer_shell_results.json` (258 gates). The final matrix checker additionally includes12 non-coordinate Householder-frame zero divisors and tests the full kernel parametrization (7); final totals are90 cases and109,028 gates. The main checker remains102,990 gates. All normal/optimized outputs agree.
