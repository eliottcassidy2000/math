# Session127: a generative source architecture for norm collapse

**Status:** proved uniform implications; cited global cubic classification;
exact bounded quartic/quintic controls. No claim of historical priority or
complete quartic/quintic Thue–Mahler solution lists.

GSMA is the user's **Generative Source Morphism Architecture**. This session
instantiates it with explicit sources, maps, predicates, and failure tests.
Only the name was recovered from the user; the construction below is proposed
here rather than attributed to an unseen prior formalism.

## 1. One source family, several genuinely different observations

For each order k, use

    f_k(X)=X^k−X^(k−1)−...−X−1,
    alpha_k = its positive real root,
    F_k(p,q)=p^k−p^(k−1)q−...−q^k.

At k=3,4,5 the irreducible polynomials define number fields. The exact source
is the coefficient vector of alpha_k^(-m) in its named power basis. Applying
the high-coordinate projection tests whether this vector belongs to the
rank-two additive lattice Z+Z alpha_k. This is the **sparse-unit gate**.
The lattice is not a subring, and the gate is not equivalent to norm one.

A passing vector supplies a signed linear element p−q alpha_k. Its norm,
its rational approximation p/q, the reciprocal rational-function input q/p,
and its Euclidean embedding are separate observations. Retaining sign,
basis, order, embedding, and predecessor is essential.

The [GSMA packet](gsma/INHERITANCE_AND_KERNEL.md) proves the maps and their
exact domains. All orders share the identity

    alpha_k^(-k)=2−alpha_k,
    (p−q)F_k(p,q)=q^(k+1)−(2q−p)p^k.

This transports an algebraic identity across the family; it does not give
field embeddings between degrees 3,4,5. Such embeddings do not exist, since
none of these distinct degrees divides another.

## 2. The exceptional305 is fully explained and globally isolated

For tau=alpha_3,

    tau^(-16)=56tau−103,
    F_3(103,56)=−1,
    tau−103/56=tau^(-16)/56.

The preceding convergent is46/25, so the next complete quotient is exactly

    T=(tau^16−25)/56,       305<T<306.

This proves the coefficient305 without numerical rounding. The identity
itself is classical and credited. A published Delone–Nagell cubic Thue
bound, together with four explicit solutions, gives the complete result

    m>0 and tau^(-m) in Z+Ztau   iff   m in {3,16}.

See the [exact certificate and cited global proof](tribonacci/REPORT.md).
The old floating-point continued-fraction output first fails at a_12:
14 replaces17. The coefficient305 and its preceding convergent survive;
the remaining old floating tail is not inherited as evidence.

## 3. A uniform theorem that strengthens the isolated example

The [all-order approximation theorem](uniform_norm/UNIFORM_NORM_TO_CF.md)
proves, for **every k>=3 and all positive integers p,q**,

    |F_k(p,q)|=1  =>  p/q is a continued-fraction convergent of alpha_k.

This is not a finite-search inference. Positivity and unit norm force
1<p/q<=2. Normalize by

    h_k(x)=f_k(x)/x^k=1−sum_(j=1)^k x^(-j).

On [1,2], h'_k>1/2. For 1<p/q<2, the mean-value theorem gives

    |alpha_k−p/q| < 2/p^k <= 2/(q+1)^3 < 1/(2q^2).

Legendre's classical criterion finishes the proof. The endpoint2/1 is
handled separately. A general residual m=F_k(p,q) satisfies the useful
sufficient gate

    1<p/q<2 and 4|m|q^2<=p^k  =>  p/q is a convergent.

The converse fails:11/6 is a cubic convergent with norm−7. Prime support
alone also fails: F_3(3,2)=−11 is supported on S={11}, but3/2 is not a
convergent. Even a large continued-fraction coefficient need not come from
a unit norm; an exact next-digit503 hostile is supplied.

For an actual convergent with predecessor denominator s and next complete
quotient T, the exact common formula is

    T+s/q = q^(k−2) C_k(p/q)/|F_k(p,q)|,
    C_k(x)=f_k(x)/(x−alpha_k)>0.

The exponent k−2 compares the degree-k norm residual with the determinant-
controlled q^(-2) continued-fraction error. It gives scales q,q^2,q^3 at
orders3,4,5 when the other factors are retained. It does not assert an
infinite family of fixed-norm solutions.

## 4. Thue–Mahler is an exact denominator condition

The recurrence generating function

    U_k(z)=z^(k−1)/(1−z−...−z^k)

obeys, for primitive integers p,q with p!=0,

    U_k(q/p)=p q^(k−1)/F_k(p,q),    gcd(F_k(p,q),pq)=1.

Thus integrality is exactly F_k=+/-1. For a fixed finite prime set S,
S-integrality is exactly the condition that F_k has only prime factors
in S. This is a genuine uniform Thue–Mahler connection, not a similarity
of terminology. Classical finiteness applies to each of the three fixed
irreducible forms and each fixed S; it does not supply the missing explicit
quartic/quintic solution lists.

Evaluation domain is a necessary sidecar. At the exceptional reciprocal
input56/103, U_3=-323008. That point lies outside the power series' radius
1/tau: the negative integer is a value of the rational continuation, not
a convergent sum of positive recurrence coefficients.

## 5. A precise candidate Lucas47 triad

A recovered exact junction is

    L_8=47=2^7−3^4=103−56.

The [Lucas packet](lucas47/README.md) treats this as a candidate interpretation
of the requested triad, not an authenticated recovery of its earlier name.
Its source objects are a quadratic-unit trace, a power difference, and the
coefficient difference of a cubic unit. These maps share an output but
forget different information.

The unimodular shear r=p−q is a genuine norm-preserving coordinate change:

    F_3(r+q,q)=r^3+2r^2q−2q^3,
    47^3+2*47^2*56−2*56^3=−1,
    56^4−9*103^3=−47.

Using the cited complete cubic classification, the intersection is global:
for q>0, n>=0, p−q=L_n and |F_3(p,q)|=1, the only triples are

    (n,p,q)=(1,2,1), (8,103,56).

This is stronger than the accompanying index100 census. It is not a
field isomorphism or a new proof of the cited cubic classification.
A hostile retained from prior work is Norm(phi^8−1)=−45, not−47.

## 6. Where the ball-volume formula really enters

[Geometry and dimension](volume/THEORY.md) distinguishes three quantities:
field degree k, binary parameter rank2, and norm-polynomial degree k.

    V_n=pi^(n/2)/Gamma(n/2+1),
    V_(n+2)=2pi V_n/(n+2).

The full Minkowski lattice has dimension k. Its Euclidean ball and covolume
provide a sufficient small-norm existence bound. The binary candidate
lattice is a plane: its metric sections use V_2=pi and their own Gram
matrix. Neither statement identifies norm one with Euclidean smallness.
The point103−56tau has norm−1 but length greater than180 in the
specified sqrt(2)-weighted Minkowski embedding. The metric convention matters.

The full real norm sublevel has infinite volume. In contrast the binary
homogeneous norm sublevel has finite area

    A_f = integral_R |f(t)|^(-2/k) dt,

because k>2 and the real roots are simple. It still has unbounded cusps.
The degree-versus-rank gap underlying this integrability threshold is
consistent with the k−2 exponent in the approximation formula, but area
alone is not an integral-point finiteness or height theorem.

The polynomial discriminants are−44,−563,9584. These are order-basis
invariants; maximal-order field discriminants cannot be substituted without
an index computation. Likewise V_5 being the largest integer-dimensional
unit-ball volume does not explain or classify quintic norm solutions.

## 7. What closes, and what remains open here

- Complete: cubic inverse-power collapse and the stated Lucas intersection,
  using an explicitly cited classical theorem.
- Uniform proved: denominator/S-integrality gate, positive-unit-to-convergent
  theorem, exact next-quotient identity, and the typed geometric bounds.
- Finite exact: inverse powers through500 give only m=4 at order4 and m=5
  at order5 (apart from m=0); no global exclusion follows from that census.
- Separate: all binary norm units, the cyclic alpha-power source, arbitrary
  S-units, and full lattice elements are not interchangeable solution sets.

The session's anchor was the norm-collapse mechanism, its niche was exact
continued-fraction certification, and its wildcard was the geometric
rank/degree comparison. Each successful connection now has a map, retained
hypotheses, a hostile example, and a reproducible test. No general tournament
spectrum theorem, new GSMA axiomatization, or complete higher-order Thue–Mahler
classification is asserted.

[Verification and reproducible full-session runner](VERIFICATION.md).
