# From matched null boundaries to a determinant singlet: Gaussian nullcones, color gluing, and a planar obstruction

2026-11-26. Read against Venkatesa Chandrasekaran, [*An algebraic quantization of gravity from the null cone*, arXiv:2611.34924v1](https://arxiv.org/html/2611.34924v1). **Proved here:** the factorial-weighted rank-one nullcone theorem (using the inherited one-variable constant-term theorem), the exact radial image of the polynomial Weyl algebra on determinant singlets, classification of a marked classical and quantum boundary shear, and the explicit resummation in the three-color model. The replay gives independent exact finite checks. This is a finite-mode bosonic color model, not QCD dynamics, gravity quantization, a mass-gap result, or a proof of planar Jacobian or Dixmier.

## 1. What actually transfers from the preprint

- [§3.6](https://arxiv.org/html/2611.34924v1#S3.SS6), (3.60)–(3.66): the overlap matching uses actual frames and a triple-overlap consistency relation, not equality of aggregate charges alone. With opposite orientations, the common shift becomes gauge when the two charges match; the *relative* boost and angular matching retain information. This suggests testing a diagonal invariant projection, not independently zeroing each edge state.
- [§4 and Appendix C](https://arxiv.org/html/2611.34924v1#S4), especially (4.1) and (4.29)–(4.30): deformation insertions have an ordering dictated by their side of the cut. For prepared states, absolute gluing follows from square summability by Cauchy–Schwarz. That is not itself a trace bound or a reason to remove the ordering/contact data. The finite-dimensional causal-diamond claim is a separate §5 argument under the paper's hypotheses. The null vectors in §5.4 are in a gluing Gram kernel.
- [§4.4](https://arxiv.org/html/2611.34924v1#S4.SS4), (4.33)–(4.34): for the displayed Gaussian two-graviton example the full sum converges absolutely but the cubic term of the expansion before removing the cutoff produces a harmonic sum. §5.4 and this zero/limit issue have no automatic connection to the Gaussian moment nullcone.

The repo's [NC2 Wick-face proof, THM-2022](../../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md) retains a tied balanced face before reducing modulo a prime; the [typed color work](../session77_ternary_color_20261027/COLOR_FUSION.md) retains actual invariant contractions beyond triality; the [rank-two compensator firewall](../../../01-canon/theorems/THM-4408-rank-two-poisson-hamiltonian-primitive-nondescent-compensator-firewall.md) separates an inverse on a larger carrier from descent to a specified polynomial subring. The construction below combines these three ideas in an exactly solvable setting.

## 2. A whole higher-dimensional Gaussian sector has the rank-one nullcone

Let Z be a d×d matrix of independent normalized circular complex Gaussians and write u=det Z. Left SU(d) acts on rows (color), with d distinguishable bosonic columns. For the compact color-singlet projection, invariant holomorphic polynomials are exactly C[u]: SU(d)-invariance of a holomorphic polynomial implies SL(d,C)-invariance by differentiating and complexifying; on the dense invertible set the latter group acts transitively on each fixed-determinant fiber. Evaluation on diag(u,1,...,1) gives the polynomial in u. This does not identify quark fields with commuting variables.

The usual **orbit-closure nullcone** on configurations is det Z=0. Indeed a singular matrix of image dimension r<d can be sent to zero by scaling the image with t^(d−r) and a complementary subspace with t^(−r), in a determinant-one one-parameter subgroup. A nonzero determinant cannot tend to zero along an orbit. In contrast a **moment-null** element is a polynomial whose positive Gaussian moments vanish. For example u is moment-null on account of the global phase of the Gaussian, even though it is a nonzero SU(d)-invariant observable. The global phase charge is not color.

By Gaussian column elimination, the successive squared perpendicular lengths are independent gamma variables of shapes d,d−1,...,1. Conditional on earlier columns, rotation invariance gives the next gamma law independent of the earlier values. Their product is |u|². The global phase then gives

    E[u^a conj(u)^b] = 0              for a ≠ b,
    h_n = E[|u|^(2n)] = ∏_(j=0)^(d−1) (n+j)!/j!.

**Theorem (factorial-weighted rank-one transfer).** For any linear moment functional that kills a≠b and at a=b has weight w(a) equal to a nonempty product of factorials (c_j a+h_j)!/h_j!, where c_j are positive integers and h_j nonnegative integers, the moment-nullcone in C[u,v] is exactly the union of the positive and negative strict one-sided charge loci, plus zero. Each charge is exponent of u minus exponent of v. In particular, this holds on C[det Z, conj(det Z)] for every d, despite the ambient Gaussian having 2d² real coordinates. For the determinant sector the eventual vanishing for a fixed polynomial multiplier holds even when the multiplier is **not** a color singlet: global phase sees charges d(a−b), and a polynomial multiplier has only finitely many opposing phase charges.

Proof. The forward direction for a one-sided support follows by adding charges. For the converse, fix an exactly nonzero complex monomial support. If all its moments vanish for some coefficient tuple in the coefficient torus, Noetherianity and a maximal ideal over Q give an algebraic coefficient tuple with the same support and vanishing. Use THM-2022's lowest balanced-height face: minimize the u exponent per letter on balanced convex combinations, and choose its supporting linear functional. There is no collision of charges among the distinct monomials of this face because the supporting equality determines the u exponent from the charge. The inherited [one-variable compound-root/constant-term theorem](../../../01-canon/theorems/THM-2111-effective-compound-root-bound-for-one-variable-constant-terms.md) gives a nonzero face constant term Q at some order m0, all of whose monomials have the same balanced exponent A0.

Choose a good prime p>m0 for the algebraic coefficients and Q. In order p m0, all balanced channels have exponent A≥p A0. Divide the entire moment by the common integer w(p A0) *before* reduction. Every factorial ratio is an integer since all its arguments increase with A. The multinomial carries kill allocations not all divisible by p. For a dilated allocation outside the face, A≥p(A0+1); each factorial interval of length at least c_j p contains a multiple of p, regardless of the offset h_j. The ratios on the face are all one. Lucas and Frobenius thus give the residue Q^p, nonzero, a contradiction. This supplies the load-bearing extension of the original lowest-face proof: several factorials are harmless because they have **one common monotone radial coordinate**. For independent Gaussian coordinates without this constraint the proof does not apply.

An exact color/gluing control: if Π is Haar projection to the singlet, then Π z_ia Π = Π ∂_jb Π = 0 (fundamental and dual are nontrivial). But on u^n,

    Π z_ia ∂_jb Π = δ_ij δ_ab N/d,
    Π ∂_jb z_ia Π = δ_ij δ_ab(1+N/d),

where N u^n=n u^n. Indeed the color average is its trace divided by d, and the remaining flavor infinitesimal action on det^n is n times the flavor trace. This gives an explicit nonzero intermediate colored contribution that would disappear if the two pieces were projected individually. The algebraic fact here is compact SU(d) color projection, not an identification of SU(d) with the Lorentzian surface symmetry.

## 3. The exact reduced operator algebra and its planar chart

Take the holomorphic Gaussian inner product, U multiplication by u, D=det(∂/∂Z)=U† on polynomials, and N=(1/d)Σ z_ia∂_ia. Let G(n)=n(n+1)···(n+d−1). The orthogonal basis u^n and the weight ratio h_n/h_(n−1) give

    D u^n = G(n) u^(n−1),
    [N,U]=U,  [N,D]=−D,  UD=G(N),  DU=G(N+1).

**Radial image theorem.** All SU(d)-equivariant polynomial differential operators on the d×d variables, restricted to the singlet, have image exactly the algebra generated by U,D,N. Split such an operator into homogeneous total-degree shifts; this preserves equivariance since SU(d) commutes with scaling. A component of shift ds sends u^n to a(n) u^(n+s) for all sufficiently large n, and any other shift acts as zero. The coefficient a(n) is polynomial in n: evaluate the differential operator at Z=I and expand the bounded jet of det^n there using binomial coefficients. For a raising or zero shift this is U^s a(N). A lowering shift by j is the adjoint of a raising shift by j, so its coefficient is divisible by h_n/h_(n−j) as a polynomial in n. That is exactly D^j times a polynomial of N. The same descriptions agree at the few small n where the lowering result is zero. Linear independence of different shifts and of polynomials evaluated on all sufficiently large n gives the normal form.

There is **no** singlet polynomial differential operator P with [P,U]=1 when d≥2. Only its degree −d component contributes to shift zero. If it lowers with coefficient a(n), it would need a(0)=0 and a(n+1)−a(n)=1, giving a(n)=n, whereas that coefficient must be divisible by G(n). The rational correction D[(N+1)···(N+d−1)]^(-1) *is* the canonical derivative on the singlet, as an operator on polynomials, but does not belong to the radial polynomial algebra. This is a specified descent failure, not a Weyl embedding obstruction in general.

Giving U,D weight d and N weight 2 in the normal form yields a commutative associated graded algebra

    R = C[u,v,s]/(uv=s^d) = C[x,y]^(μ_d),
    u=x^d, s=xy, v=y^d, with {y,x}=1/d.

The leading commutators give {s,u}=u, {s,v}=−v and {v,u}=d s^(d−1), consistent with the displayed plane cover. On u≠0 the conjugate coordinate is s/u. This is a canonical chart, not a global polynomial coordinate.

**Matched-shear theorem (classical and quantum).** Among Poisson endomorphisms of R fixing u, the only possibilities send s to s+f(u) with f∈u C[u], and all are automorphisms. The exact analogue holds for endomorphisms of the radial algebra fixing U: N goes to N+f(U), again with f∈U C[U], and all are automorphisms.

- Classically, localize at u. The bracket with u forces ∂_s(s_new)=1. Regularity in R makes the added Laurent term polynomial in u. The image of v would have to be (s+f(u))^d/u; reduction modulo (u,s) forces f(0)=0. Conversely if f(u)=u h(u), this is the descent of the **planar Keller shear** (x,y)↦(x,y+x^(d−1)h(x^d)), with explicit triangular inverse and Jacobian determinant one.
- Quantum mechanically, the centralizer of U on the polynomial module consists of multiplication by polynomials of u, determined by the value on 1. Thus N_new=N+f(U). From U D_new=G(N_new), the lowering coefficient of D_new on u^n is G(n+f(0)). It must be divisible by G(n); since the polynomials are monic of equal degree, comparing their next coefficient gives f(0)=0. Conversely choose a polynomial g with g′(u)=−f(u)/u. The finite-on-each-operator conjugation exp(ad g(U)) gives this shear: commutation with multiplication reduces ambient differential order until zero. D_new is determined uniquely by the product relation.

This is a quantified, marked extension of the planar-Jacobian theme: a shift allowed in a local canonical chart need not descend across its singular boundary. It makes no statement about arbitrary Keller pairs. The comparison to the paper's relative boost is at the level of retaining a conjugate coordinate and respecting the boundary, not a physical identification.

## 4. Two exact checks of summation before expansion

In the preprint's example set t=1/s. Its series is (t²/2) Li₂(e^(-2it)); the cubic cutoff divergence is replaced in the resummed result by **i t³ log t**. Indeed for 0<t<π the derivative of Li₂(e^(-2it)) is 2i log(1−e^(-2it)), with argument π/2−t. Integrating from zero gives

    (π²/12)t² + i t³(log(2t)−1) − (π/2)t³ + t⁴/2 + O(t⁵),       t→0+.

In the three-color bosonic model take the balanced real observable P=u+conj(u). Condition on the first two columns. The variance of the remaining conditional complex Gaussian determinant is the product of Gamma(3) and Gamma(2) variables. Hence, for z=t²≥0, its actual characteristic function is

    Φ(z)=E[exp(i t P)]= (1/2) ∫_0^∞ e^(-x) x²/(1+zx)² dx.

It is finite, and has the formal series Σ(-1)^n h_n z^n/(n!)². That series has radius zero, even though every moment exists. The exact Borel transform is (1−2ξ)/(1+ξ)^4, since h_n/(n!)³=(n+1)²(n+2)/2. Its Laplace integral on the positive ray equals Φ(z): (1−2zx)/(1+zx)^4 is the second x derivative of x²/[2(1+zx)²], and integration by parts twice gives the equality. The two examples exhibit different analytic singularities but the same warning: fixed-moment algebra or absolute convergence of a glued amplitude does not justify termwise parameter expansion.

Replay: **python3 replay.py** (SymPy) checks determinant differential identities directly, two simultaneously surviving Wick face channels with exact prime residues, the classical Jacobian-one lifts, quantum constant-shift rejection, and the Borel formula coefficients. The arguments above, not the finite checks, supply the general statements.
