# Unit-ball volume, algebraic norm, and the rank-two Thue slice

**Status:** PROVED statements with proofs below; FINITE-EXACT verification in `verification.json`; numerical illustrations explicitly marked NUMERICAL. No global conjecture is claimed.

**Session:** 127, 2026-10-31. **Inherited source:** read-only repository snapshot `522535a7b6f519fd0c1e4304a2c5b1c909fce4c9`. The startup packet reported a cached one-commit lag; this note makes no claim about that unseen commit. GSMA means the user's “Generative Source Morphism Architechture”; no pre-existing axioms or categorical definition were recovered in this lane. The polynomials below are a specified reconstruction, not a claim to have authenticated the user's old `normcollapse` definition.

## 1. The useful conclusion

The formula `pi^(n/2)/((n/2)!)` is the volume of an **n-dimensional Euclidean unit ball**, with half-integer factorial interpreted through Gamma. It enters a norm-form problem only after specifying the real vector space, its metric, its measure, and its lattice.

For degree `k=3,4,5`, these are three different objects:

1. the full algebraic-integer lattice has real rank `k`;
2. the candidate set `p-q theta` in a binary Thue equation has rank **2**;
3. the polynomial norm has homogeneous degree `k`.

Thus `V_k` belongs to the full-lattice Euclidean bound, while the binary metric section uses `V_2=pi`. Neither is the volume of the algebraic norm sublevel. The full norm sublevel has infinite volume; the binary norm sublevel has finite but polynomial-dependent area and unbounded cusps. At the tribonacci input `(103,56)`, norm `-1` coexists with Euclidean length about `185.23`. This is a concrete example of anisotropic cancellation, not a small Euclidean vector.

## 2. Inheritance and concept board

Closest proved mechanism: discriminant/covolume plus Minkowski's convex-body theorem, not a numerical coincidence between sequence terms. Canonical hostile: a norm-one unit with an expanding real embedding. Corrected near miss: the algebraic trace `Tr(theta^2)=3` is **not** the positive Euclidean squared length `sum |sigma(theta)|^2`. Least-used sidecar: the rank-two Gram matrix and the shape of a binary norm cusp.

The relevant repository warning is the connection-contract/type audit in `00-navigation/RESEARCH-PROTOCOL.md`, and the scalar/operator warning under the 2026-10-28 Lucas correction in `01-canon/MISTAKES.md`. The `META-PATTERNS.md` type/direction card applies: scalar, dimension, rank, and polynomial degree are not interchangeable. No theorem-ID promotion is attempted here.

| Object | Exact invariant | Operation | Lost coordinate / hostile |
|---|---|---|---|
| Euclidean balls | `V_n`, radius, metric | dimension shift by 2 | arithmetic and lattice absent |
| Number-field order | degree, discriminant, signature | Minkowski embedding | none if the full embedding is retained |
| Binary candidate lattice | rank 2, Gram matrix | section of full embedding | arbitrary order elements usually leave section |
| Algebraic norm | product of embeddings | multiply conjugates | individual sizes, arguments, signs except product |
| Norm sublevel | homogeneous star body | scalar norm threshold | convexity/compactness can fail |
| GSMA architecture | declared source/map/predicate | typed composition | unspecified morphism cannot carry an asserted theorem |

## 3. Ball volume and a recurrence obstruction

### VOL-1 (PROVED): ball formula and exact dimension shift

For every integer `n>=0`, with the zero-dimensional ball given volume 1,

`V_n = pi^(n/2)/Gamma(n/2+1)`,

`V_(n+2) = (2 pi)/(n+2) V_n`.

Proof: for `n>=1`, integrate `exp(-||x||^2)` in Cartesian coordinates and radially. The case `n=0` is the separate convention `V_0=1`, which agrees with the displayed Gamma formula. The Cartesian integral is `pi^(n/2)`. The radial one is

`n V_n integral_0^infinity exp(-r^2) r^(n-1) dr = V_n Gamma(n/2+1)`.

The recurrence follows from Gamma's functional equation. The standard Gamma conventions are recorded in [DLMF 5.4](https://dlmf.nist.gov/5.4) and [DLMF 5.5](https://dlmf.nist.gov/5.5).

In particular,

| `n` | `V_n` |
|---|---|
| 2 | `pi` |
| 3 | `4 pi/3` |
| 4 | `pi^2/2` |
| 5 | `8 pi^2/15` |
| 6 | `pi^3/6` |

The integer-dimensional volume is largest at `n=5`: each parity's maximum is at 5 or 6 by the two-step recurrence, and `V_6/V_5=5 pi/16<1`. This is a geometric fact only; it does not rank the difficulty or number of solutions of cubic/quartic/quintic norm equations.

### VOL-2 (PROVED): no direct tribonacci/Lucas recurrence identity

The sequence `(V_n)_(n>=0)` satisfies no eventual nonzero constant-coefficient linear recurrence. Indeed, each parity subsequence has successive ratio tending to zero, so its ordinary generating series, and therefore `sum V_n z^n`, is entire. An eventual constant-coefficient recurrence would make that ordinary generating function rational. An entire rational function is polynomial, contradicting `V_n>0` at every index. This also rules out direct identification with a nonzero tribonacci or Lucas sequence over the complete infinite index range. It does not rule out a specifically declared transform or selected finite identity.

## 4. The full lattice: the legitimate `V_k` bridge

Let `K` be a number field of degree `d=r+2s`, with real embeddings `sigma_i` and one representative `tau_j` of each complex-conjugate pair. Use the **scaled** real embedding

`Phi(alpha)=(sigma_1(alpha),...,sigma_r(alpha), sqrt(2) Re tau_1(alpha), sqrt(2) Im tau_1(alpha),...)`.

In these scaled coordinates the norm polynomial is

`Ncal(x,u,v)=product_i x_i product_j ((u_j^2+v_j^2)/2)`.

In particular, omitting the factors of 2 would change the norm. Complex-block multiplication is transported as `w star w_prime = w w_prime/sqrt(2)`, not ordinary multiplication of the scaled complex coordinates.

Set

`Q(alpha)=||Phi(alpha)||^2=sum_i sigma_i(alpha)^2+2 sum_j |tau_j(alpha)|^2`.

### VOL-3 (PROVED): norm/metric inequality and covolume

For all `alpha` in `K`,

`|N_(K/Q)(alpha)| <= (Q(alpha)/d)^(d/2)`.

Equality holds exactly when all `d` embedding absolute values are equal (with each complex absolute value repeated twice). Apply AM-GM to their squares. The constant is sharp, including on every order, since nonzero rational integers attain equality.

For any full-rank order `O` with a Z-basis `b`,

`covol(Phi(O))=sqrt(|disc(b)|)`.

Proof: the determinant of the matrix of all complex embeddings has squared determinant equal to the basis discriminant. Replacing each conjugate pair by its real/imaginary rows changes the determinant's absolute value by `1/2`; subsequently scaling both rows by `sqrt(2)` multiplies it by `2`. These factors cancel. This also proves `det(G_b)=|disc(b)|`, where `G_b` is the positive Euclidean Gram matrix. **It does not identify `G_b` with the algebraic trace matrix**, which can be indefinite.

Consequently the coefficient ellipsoid `Q(sum a_i b_i)<=R^2` has coefficient-space volume

`V_d R^d/sqrt(|disc(b)|)`.

If `I` is an integral ideal of the maximal order `O_K`, then `covol(Phi(I))=sqrt(|D_K|) N(I)`. The unscaled embedding has the different factor `2^(-s)`. The standard lattice and ideal facts are [Milne, Algebraic Number Theory, Proposition 4.26](https://www.jmilne.org/math/CourseNotes/ANT.pdf), printed pp. 79-80. Order discriminant and field discriminant may differ by an index square; no maximal-order claim is needed for the power bases used below.

### VOL-4 (PROVED using the cited convex-body theorem): a sufficient bound, not an iff

For every full lattice `L=Phi(O)`, a nonzero `alpha in O` exists with

`Q(alpha) <= (2^d covol(L)/V_d)^(2/d)`

and hence

`|N(alpha)| <= C_d sqrt(|disc(O)|)`, where `C_d=2^d/(d^(d/2) V_d)`.

Use the closed Euclidean ball and Minkowski's compact convex symmetric-body theorem. If using its strict-volume formulation, take radii decreasing to the displayed threshold: a fixed slightly larger closed ball contains only finitely many lattice points, so one nonzero point occurs along a subsequence and lies in the limiting closed ball. This justifies the non-strict bound. See [Milne, Theorem 4.19](https://www.jmilne.org/math/CourseNotes/ANT.pdf), printed p. 76. The corresponding ideal bound multiplies the right side by `N(I)`.

| `d` | `C_d` |
|---|---|
| 3 | `2/(sqrt(3) pi)` |
| 4 | `2/pi^2` |
| 5 | `12/(5 sqrt(5) pi^2)` |

This guarantees an unspecified element of the **full** lattice. It does not guarantee a point of the binary slice, a specified sign, primitivity, specified prime support, or an element in a given recurrence orbit. The usual signature-adapted Minkowski body gives the often sharper constant `(4/pi)^s d!/d^d`, not the Euclidean-ball constant; see [Milne, Proposition 4.27](https://www.jmilne.org/math/CourseNotes/ANT.pdf), printed p. 80. A volume threshold is a sufficient existence test, never an exact description of norm collapse.

## 5. The reconstructed degree-3/4/5 family and its exact binary Gram matrix

Define explicitly

`f_k(X)=X^k-sum_(j=0)^(k-1) X^j`,

`K_k=Q(theta_k)`, where `f_k(theta_k)=0`,

`F_k(p,q)=p^k-sum_(j=0)^(k-1) p^j q^(k-j)`.

The independent determinant and resultant computations in `verify_volume_bridge.py` certify

`N_(K_k/Q)(p-q theta_k)=F_k(p,q)`.

Irreducibility is certified by reduction modulo the primes shown; real-root counts are certified by exact Sturm calculations. All discriminants below are **power-order/polynomial** discriminants.

| `k` | irreducible reduction prime | signature `(r,s)` | `disc(Z[theta_k])` |
|---|---|---|---|
| 3 | 3 | `(1,1)` | `-44` |
| 4 | 2 | `(2,1)` | `-563` |
| 5 | 5 | `(1,2)` | `9584` |

These are separate fields, not a tower: a unital Q-algebra homomorphism between fields is injective, and the source degree must divide the target degree. No two distinct degrees here divide one another. Hence there is no unital Q-algebra map between distinct members of this list. Cross-order GSMA transport must name a different category/map, or put the fields in a common target.

### VOL-5 (PROVED): rank-two section formula

Write the roots of `f_k` as `theta_1,...,theta_k`, including conjugates, and set

`S_k=sum_i |theta_i|^2`.

Newton's identities give `Tr(theta_k)=1` and `Tr(theta_k^2)=3`. But

`S_k = 3+4 sum_(complex-pair representatives j) (Im theta_j)^2 > 3`.

The strict inequality follows from the nonreal roots certified in the table. The correct binary quadratic form is

`Q_k(p,q)=||Phi(p-q theta_k)||^2 = k p^2-2pq+S_k q^2`

`= k(p-q/k)^2+(S_k-1/k)q^2`.

Its Gram matrix is `[[k,-1],[-1,S_k]]`, with determinant `k S_k-1>0`. The lattice `Phi(Z+Z theta_k)` has rank 2 and induced two-dimensional covolume `sqrt(k S_k-1)`; it has no finite full-dimensional covolume in `R^k`.

Therefore:

* an intrinsic radius-R disk in its real plane has area `pi R^2`;
* its coefficient ellipse in `(p,q)` coordinates has area `pi R^2/sqrt(k S_k-1)`;
* `|F_k(p,q)| <= (Q_k(p,q)/k)^(k/2)`, still only in the displayed direction;
* the norm's degree remains `k`, while the parameter lattice's rank is 2.

### The 305 predecessor as a hostile control

Let `tau` be the positive real root of `f_3`. Then `S_3=tau^2+2/tau`: the other two roots are a conjugate pair with squared modulus `1/tau`.

At the verified arithmetic input `(p,q)=(103,56)`,

`F_3(103,56)=-1`,

`Q_3(103,56)=20291+3136(tau^2+2/tau)`.

For an exact lower bound, `9/5<tau<2`, so `S_3>106/25` and `Q_3>839691/25>180^2`. Thus a norm-one point already has Euclidean length greater than 180.

**NUMERICAL illustration, not an exact certificate:**

`Q_3(103,56)=34310.02949575866...`,

`||Phi(103-56tau)||=185.22966688886...`,

`103-56tau=-0.0000582919930234...`.

The other two conjugate factors compensate for the tiny real factor. The appearance of 305 as the next continued-fraction coefficient is handled in the [arithmetic lane](../tribonacci/REPORT.md); ball volume alone neither predicts nor explains that coefficient.

## 6. Norm bodies: full space versus binary slice

### VOL-6 (PROVED): full norm sublevels have infinite volume

On `K tensor_Q R = R^r x C^s`, extend the norm as

`P(x,z)=product_i x_i product_j |z_j|^2`.

If `r+s>=2`, the set `{|P|<=1}` is unbounded and has infinite Lebesgue volume. Choose two distinct archimedean factors of real dimensions `w_1,w_2 in {1,2}`. Scaling the first factor by `t` and the second by `t^(-w_1/w_2)` preserves both `P` and volume. Start with a small box away from every zero factor inside `0<|P|<1`. Applying a sufficiently separated geometric sequence of `t` gives disjoint positive-volume boxes inside the same sublevel. This proves infinite volume, not merely unboundedness. The hypothesis holds for all degrees `d>=3`.

There are also concrete lattice hostiles in the reconstructed family. Every `f_k` has a real root `rho_k in (1,2)` and constant term `-1`. Hence `theta_k` is an algebraic-integer unit and

`N(theta_k)=(-1)^(k+1)`, `N(theta_k^(2m))=1`,

while `||Phi(theta_k^(2m))|| >= rho_k^(2m) -> infinity`.

These are points of the full power-order lattice. They generally leave the binary slice; the unit orbit must **not** be offered as infinitely many binary Thue solutions. Full norm bodies cannot be turned into bounded Euclidean balls by an invertible linear map.

### VOL-7 (PROVED): the binary norm sublevel has a different finite area

Let `f` be a real monic polynomial of degree `d>2` whose real roots, if any, are simple, and `F(x,y)=y^d f(x/y)` the associated homogeneous form. Then

`area{(x,y) in R^2 : |F(x,y)|<=1} = A_f := integral_R |f(t)|^(-2/d) dt < infinity`.

Proof: in `y>0`, put `x=ty`; the Jacobian is `y`. The radial range is `0<y<=|f(t)|^(-1/d)`, giving half the displayed integral. The half-plane `y<0` contributes the same amount by homogeneity and absolute values. The line `y=0` has area zero. At a simple real root, the integrand is bounded by a constant times `|t-r|^(-2/d)`, which is integrable because `d>2`. At infinity it is `O(|t|^-2)`. This proves the claim.

If `f` has a real root, this finite-area set is unbounded along its real root line. Every reconstructed `f_k` has such a line. It is not convex: it contains a neighborhood of the origin, and an unbounded convex planar set containing a disk has infinite area (the convex hull of that disk and arbitrarily distant points already has unbounded area). The usual convex-body theorem therefore does not directly apply to this norm sublevel.

For a threshold `M>0`, coefficient-space area scales as `M^(2/d) A_f`; induced embedded-plane area scales as `sqrt(d S-(Tr theta)^2) M^(2/d) A_f`. This correct scalar depends on the entire polynomial/root geometry, not degree alone. No numerical integral or binary-solution census is claimed here. In particular, finite area of an unbounded cusp is not by itself a finite-lattice-point or effective-height theorem.

## 7. GSMA-compatible map contract

Within each specified order, the following is a commuting **typed data diagram**:

`(p,q) in Z^2 -> alpha=p-q theta in Z[theta] -> Phi(alpha) in R^k`

and

`F_k(p,q) = N(alpha) = Ncal(Phi(alpha))`.

Here `Ncal` includes the factor `1/2` for each scaled complex block, as specified in Section 4.

| Map / derived operation | Preserved predicate | Lost data | Necessary sidecar / decisive test |
|---|---|---|---|
| `(p,q)->p-q theta` | exact addition, binary membership, norm equation | none when basis retained | named basis and rank-two image |
| `alpha->Phi(alpha)` | addition, injectivity, multiplication with transported complex-block product `w w_prime/sqrt(2)` | none | signature, embedding order, scaled metric |
| `Phi(alpha)->N(alpha)` | multiplicative norm, signed product | individual sizes, complex phases, source coordinates | full conjugate tuple or logs/arguments and coordinates |
| `Phi(alpha)->Q(alpha)` | Euclidean size | norm sign, product, arithmetic support | norm and prime valuations; `(103,56)` is the hostile |
| ball/lattice -> volume ratio | a sufficient full-lattice existence condition | labels, shape-specific lattice occupancy, binary slice, prime support | lattice basis, Gram matrix, actual candidate search |
| binary norm body -> `A_f` | exact real area and homogeneity | integral points, primitivity, local solvability | coefficient lattice plus Diophantine/local arguments |
| cross-order `k=3,4,5` | only a declared common data predicate | no field homomorphism exists between distinct degrees | common ambient target or explicitly non-algebraic map |

This is a concrete way to merge the geometric formula into a Generative Source Morphism Architecture without asserting that a scalar volume is a recurrence, a norm, a solution count, or a proof of Thue–Mahler finiteness. A complete old GSMA definition, if recovered elsewhere, must be checked against this contract before inheriting additional claims.

## 8. Verification and scope

Run `python3 verify_volume_bridge.py`. The exact universe is `k in {3,4,5}`, symbolic `p,q`, the specified polynomial family, and the finite unit-power controls recorded in JSON. Independent paths compare polynomial resultants, multiplication-matrix determinants, trace-Gram discriminants, modular irreducibility, and Sturm root counts. Positive control: `q=0` attains the AM-GM norm bound. Hostiles: the norm-one `(103,56)` vector, full unit powers, incorrect replacement `S_k=3`, and attempting to treat a rank-two lattice as a full rank-k lattice.

No assertion is made about all Thue–Mahler solutions, all continued-fraction partial quotients, Lucas-47 classifications, maximal orders, class numbers, or any global open conjecture. The missing historical normcollapse definition does not obstruct the explicitly conditional reconstruction; it does prevent claiming that this is the authenticated old object.
