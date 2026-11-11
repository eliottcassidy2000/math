# Session 267: Gaussian null packets, uniform escape, and the planar inverse

## Status and inheritance

**PROVED HERE, conditional only on inherited NC2:** uniform-family charge orientation; explicit radical of the finite-support moment ideal; finite moment certificates; polynomial marked partition functions. **PROVED HERE without NC2:** the rank-one covariance sampling lower bound and formal inverse packet implication. **OPEN:** existence of bounded null packets for planar-origin nilpotent Keller reductions. No claim that JC(2) is settled, that high degree alone excludes counterexamples, or that a topological anomaly has been classified.

Pinned source: `1432c0ac0d987f1f7d8219a04c5326b2826e552d` in `eliottcassidy2000/math`. Exact connector reads recovered:

- `01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md`: inherited NC2 and GMC(2), including scalar radial-factor mechanism and higher-dimensional carry obstruction.
- `04-computation/lean/TournamentH7/TournamentH7/GMC2Main.lean`: exposes unconditional `GMC2.nc2` and `GMC2.gmc2`. This session uses the repository theorem; it did not rebuild or re-audit it and makes no external-consensus claim.
- `01-canon/theorems/THM-2801-sharp-special-image-boundary-and-beta-shift-witness.md`: inherited SIC(1), explicit SIC(2) failure, and the surviving xi-linear cubic and Hessian restrictions.
- `01-canon/theorems/THM-2811-special-image-linear-intertwiner-rigidity-and-jc-degree-wall.md`: affine intertwiner rank obstruction, bidegree wall, fixed versus growing stencil distinction.

Closest proved mechanism: strict charge escape. Canonical hostile: THM-2801's bidegree-(4,4) null tower with a persistent coordinate observer. Corrected near miss: compressing all two-pair contractions to one pair by a linear substitution is impossible. Least-used sidecar selected here: the *parameter-wide orientation* of a null family and the fixed complexity of its marked exponential generating function.

Live board: (1) charge geometry; (2) moment ideals; (3) formal inverse residues; (4) finite Wick packets; (5) covariance sampling rank; (6) recursive interfaces. The positive connection is between (1), (3), (4), and (6); (5) supplies its exact failure boundary.

## 1. Three different meanings of dimension

For a circular complex Gaussian Z and W=conjugate(Z), use the algebraic functional

    L(Z^a W^b) = a! if a=b, and 0 otherwise.

This is two real Gaussians, or one complex Gaussian pair. In this convention NC2 says that a nonzero P with L(P^m)=0 for all m>=1 has either only positive charges a-b or only negative charges. GMC(2) says L(QP^m)=0 eventually for every fixed Q.

The special-image operator in n pairs is instead

    E_n(xi^alpha z^beta) = partial_z^alpha z^beta.

Its scalar constant-term contraction uses n complex Gaussian pairs, hence 2n real Gaussians. E_n is polynomial-valued and scalar Wick expectation is not generally the same predicate. They agree on the balanced bihomogeneous slice. Consequently ordinary GMC(2) naturally reaches SIC(1), not unrestricted SIC(2). Inherited THM-2801 already refutes unrestricted SIC(2).

The cubic homogeneous reduction relevant to JC uses f=xi dot H(z), with xi-degree one and z-degree three, usually in a dimension larger than the original plane. Calling the original map planar does not turn that reduction into two real Gaussian variables.

## 2. Uniform escape and marked partition functions

### Proposition 267-G1 (uniform observer bound)

Assume inherited NC2. Let P be nonzero and L(P^m)=0 for all positive m. Choose its charge sign epsilon in {+1,-1}, and let delta=min epsilon(a-b) over its support. Thus delta>=1. For any Q, let B=max(0,max[-epsilon(a-b)] over supp Q). Then

    L(Q P^m)=0 whenever m delta > B.

Proof: every product monomial has signed charge at least m delta-B>0 and is killed by L. In particular deg Q<=E implies vanishing for m>E, independent of deg P. This strengthens an existence threshold to an explicit uniform one; it is a direct corollary of the inherited classification, not a new proof of NC2. It is sharp: P=Z and Q=W^E give L(QP^E)=E!.

Define the formal marked partition function

    A_Q(t)=sum_{m>=0} L(QP^m)t^m/m!.

It is a polynomial of degree at most floor(B/delta). There are no convergence assumptions: this is a formal series becoming a polynomial.

### Proposition 267-G2 (irreducible parameter-family orientation)

Let X be an irreducible complex affine variety and P in O(X)[Z,W] have finite Gaussian support. Suppose for every x in X(C), L(P_x^m)=0 for every m>=1. Then globally either all nonpositive-charge coefficients vanish or all nonnegative-charge coefficients vanish. The identically zero family belongs to both alternatives.

Proof: split coefficients into positive, negative, and neutral charges. By NC2, X is the union of the two closed loci where all nonpositive coefficients vanish or where all nonnegative coefficients vanish. Irreducibility forces one locus to be all X. Regular functions vanishing at every complex point are zero in the reduced coordinate ring. This proof is entirely algebraic and avoids selecting an orientation separately at each parameter value.

Equivalently, over any integral finitely generated complex algebra R, the identities L(P^m)=0 in R imply a common sign. If Q in R[Z,W] has fixed finite Gaussian support, G1 gives one threshold valid on every fibre, including special fibres where coefficients disappear. A reducible parameter space may require different signs on different components; finite-component bounds still work by taking their maximum. Nonreduced parameter rings require a separate nilpotent-coefficient analysis and are not claimed here.

This closes a serious quantifier gap: pointwise eventual vanishing alone need not supply one polynomial inverse coefficient bound. Here finite support plus irreducibility does.

### Corollary (polarized null spaces)

If a finite-dimensional linear subspace V of C[Z,W] consists entirely of null polynomials, then V lies entirely in the positive-charge subspace or entirely in the negative-charge subspace. Apply G2 to the generic linear combination of a basis. Conversely either containment gives nullness.

Individual null basis vectors do not suffice. Z and W are separately null; their span contains Z+W with nonzero second moment. The correct hypothesis is nullness of every linear combination, equivalently all polarized moment identities.

## 3. The explicit moment-ideal radical and finite certificates

Fix a finite monomial support S, allow zero coefficients, and write P_c=sum_{s in S} c_s Z^{a_s}W^{b_s}. Set M_m(c)=L(P_c^m) and I=(M_1,M_2,...) in C[c_s]. Let C_0, C_+, C_- be its neutral, positive and negative coefficient variables.

### Proposition 267-G3

    radical(I) = ( C_0, c_p c_n : p positive, n negative ).

Proof: NC2 identifies V(I) as the union of the two coordinate subspaces defined by (C_0,C_-) and (C_0,C_+). Its vanishing ideal is their intersection, which is exactly the displayed squarefree monomial ideal. Apply the Nullstellensatz.

Consequences:

1. There exists N(S) such that vanishing of M_1,...,M_N already forces all higher moments to vanish and fixes a charge side. No effective bound is asserted.
2. For each neutral coefficient c_0 and each opposite-sign product c_p c_n, some positive power belongs to a finitely generated moment ideal. These are explicit *targets* for algebraic certificates, not merely a search for more moments.
3. On the exact-support coefficient torus, when S is not one-sided, there exists a finite certificate

       1=sum_{m<=N} A_m(c) M_m(c)

   with Laurent-polynomial A_m. Clearing denominators yields a monomial multiple on the left. Rational certificates exist because the moment equations have rational coefficients and finite linear ideal-membership systems descend from C to Q.

Proof of finite cutoff: Noetherianity gives finitely many generators from the moment family; choose the largest index. The same reasoning applies to the localization. This is an existence theorem, not a claim that a chosen prefix already generates the ideal.

There is no degree-independent first-return bound. For P=Z^r+W^s, the first nonzero moment occurs at L=(r+s)/gcd(r,s), with value

    binom(L,s/gcd(r,s)) (rs/gcd(r,s))!.

This is a useful hostile for any fixed-depth recursive or tournament certificate. A certificate budget must retain support width or another growth-controlling invariant.

## 4. Why low-rank covariance sampling does not automatically solve the bridge

Let E_C denote scalar complex Wick contraction with covariance C, characterized by E_C exp(xi dot u+z dot v)=exp(u^T C v). For a rank-one C=a b^T, it is realized on one circular pair by xi=a W, z=b Z. For polynomial-valued E_C at an external point y, use z=y+bZ instead.

For bidegree (m,m) in two pairs, the rank-one contraction of xi_1^k xi_2^{m-k} z_1^l z_2^{m-l} is

    m! a_1^k a_2^{m-k} b_1^l b_2^{m-l}.

The identity-covariance target, divided by m!, is the (m+1)-by-(m+1) diagonal matrix

    D_{kl}=delta_{kl}/binom(m,k).

### Proposition 267-G4 (sampling-rank lower bound)

Any expression of E_I on the whole bidegree-(m,m) space as a linear combination of J rank-one covariance evaluations needs J>=m+1. Every evaluation contributes a rank-one matrix, whereas D has rank m+1. This allows arbitrary complex nodes and weights, even depending on m.

Thus fixed finite covariance quadrature cannot work on all degrees. This result does NOT exclude compression on the much thinner span of a single tower {f^m,z_i f^m}; that is precisely the loophole the conjecture below uses.

A second obstruction is more basic: pure nullness at C=I does not descend to each rank-one sample. For f=xi_1 z_2, E_I(f^m)=0 for all m>=1, while C=ones(2,2) gives E_C(f^m)=m!. This hostile is already xi-linear and nilpotent-linear. Interpolation cancellation is not a collection of individually null packets.

For cubic H, the unshifted rank-one pullback f=(a dot H(b))W Z^3 has automatic one-sided charge. That observation is useless for the polynomial-valued operator: the external shift y is essential and destroys this automatic charge argument. One must preserve output polynomials, not merely their constant terms.

## 5. A precise conditional chain to JC(2)

### Formal inverse interface (standard identity, proof included)

Let H:C^n->C^n have order at least two and nilpotent Jacobian, f=xi dot H(z), and F_t(z)=z-tH(z). Write G_t for its formal inverse. Then

    sum_{m>=0} E_n(g(z) f^m) t^m/m! = g(G_t(y)),

where the output coordinates of E_n are renamed y. In particular E_n(f^m)=0 for m>=1 and the marked tower g=z_i is the inverse-coordinate series.

Formal residue proof: expansion of the left side is

    sum_{alpha>=0} t^{|alpha|}/alpha! partial_y^alpha(g(y) H(y)^alpha).

Expand the formal Cauchy kernel for z-y-tH(z) to identify this with the formal residue of g(z)/prod_i(z_i-y_i-tH_i(z)). Formal substitution w=F_t(z) evaluates it as g(G_t(y))/det(I-t JH(G_t(y))). Nilpotence makes this determinant identically one. All operations are t-adic, so no analytic integration is used. Without nilpotence the determinant factor must remain.

Let K_pl be a specified class of homogeneous cubic nilpotent-Jacobian maps containing a stable reduction of every planar Keller map, with polynomial invertibility reflecting back through the reduction. The standard dimension-raising cubic homogeneous reduction is a **CITED reduction input**, not a theorem newly proved in this session. The subclass must carry its planar provenance; replacing it by arbitrary two-pair polynomials is invalid.

### Conjecture FGNP_pl (finite Gaussian null packets, OPEN)

For each H in K_pl and each coordinate i, there exist finitely many P_ij,Q_ij in C[y_1,...,y_n][Z,W], independent of m, such that

    L(P_ij^m)=0 in C[y] for every m>=1,
    E_n(z_i f^m)(y)=sum_j L(Q_ij P_ij^m) for every m>=1.

One may also permit fixed polynomial prefactors in y, absorbed into Q. An independently certified finite prefix is harmless: the identity need only hold for all m beyond an explicit fixed index. A number of packets or Gaussian degree that grows with m is NOT this conjecture.

### Proposition 267-G5 (proved implication)

Inherited NC2 + FGNP_pl + the stated planar stable-reduction input imply JC(2).

Proof: G2 gives each packet a parameter-wide orientation. G1 supplies a uniform threshold for every marked packet, and finite summation supplies a threshold for each inverse coordinate. The formal inverse interface then makes G_t polynomial in both y and t. Substitute t=1 into the two formal inverse identities to obtain a polynomial inverse for z-H(z). Reflect invertibility through the specified stable reduction to the original planar map.

The weaker conjecture that the marked generating functions merely admit a finite sum of polynomial Gaussian-null marked partition functions is equivalent to the needed packet implication but may permit more flexible identities. The coefficient-level FGNP form is preferable for exact tests.

There is no implication FGNP_pl from GMC(2) alone supplied here. FGNP_pl is the new missing construction, not a renamed proved theorem. At bare existence level this conjecture is **equivalent to eventual marked vanishing**, not a strictly stronger intermediate principle. Indeed any finitely supported sequence a_m(y) can be encoded by P=Z and Q=sum_{m>=1} a_m(y) W^m/m!: then L(QP^m)=a_m(y). Thus a polynomial inverse supplies packets after the fact. This is a useful sharp warning: bare packet existence alone only repackages the inverse problem. The genuinely additional research proposal is the local, forward-data construction in Section 6, which must not use the inverse or its eventual truncation as an oracle. The new proved advances here are the uniform-family theorem, explicit moment-ideal radical and covariance sampling obstruction, not a claim that this reformulation by itself advances JC.

## 6. Recursive finite-interface strengthening

An ambitious but testable sufficient strengthening is to encode each planar-origin reduction by a finite rooted decomposition. At a leaf, produce finite null packets. At each composition node, prove an exact rule taking the child marked sequences to a finite sum of sequences of the form

    sum_{r=0}^R a_r(y) L(Q_r P_r^m),

with P_r null over C[y], fixed finite support, and no m-dependent increase in R or observer degrees. This closure rule implies FGNP_pl by structural induction, even when the source map has arbitrarily high degree.

This is a genuine route by which recursive structure could rule out high-degree counterexamples: arbitrary depth is allowed but each *finite* map must admit a finite lawful packet derivation. A uniform absolute bound independent of map degree is unnecessary. The exact algebraic closure rule is currently OPEN, as is the assertion that every planar Keller map has such a derivation.

For a tournament analogue, use packet obligations as vertices and let a signed charge budget be the observable. Strictly ordered budgets give an acyclic orientation; equal budgets remain ties, not arbitrarily oriented tournament edges. This ordering certifies termination only after the contraction identities and signs are proved. Pairwise dominance alone cannot establish a missing packet identity or cure opposite-charge cancellation.

The appropriate topological sidecar is still the actual algebraic inverse/nonproperness data of the source map. A surface sketch or recursive combinatorial classification cannot imply global polynomial invertibility without the exact map from that data into the packet rule. The conjecture localizes the gap in a concrete algebraic construction rather than postulating an unspecified exotic topology.

## 7. First research experiments and boundary controls

1. Start with triangular homogeneous H=(p(z_2),0). Then f=xi_1 p(z_2); the first coordinate tower is a_1=p(y_2), a_m=0 for m>1, represented explicitly by P=Z, Q=p(y_2)W. The second marked tower is zero. Record how this forward-data construction transforms under a legal coordinate gauge.
2. Study tame compositions after stable cubic reduction. An inverse exists, but the challenge is a bounded compositional *construction* of null packets from the forward description, not using the already-computed inverse as an oracle.
3. Compute moment-ideal certificates for sparse candidate packet supports. G3 specifies exactly which neutral/opposite-sign monomials the ideal must detect.
4. For a proposed Gaussian compression, test the covariance hostile xi_1 z_2 and the THM-2801 persistent-observer base. The first detects unjustified individual-nullness transfer; the second proves that any construction extending to every special-image null tower is false.
5. Record output parameters before taking Gaussian expectations. Specializing y=0 can erase the entire inverse obstruction.

Executable companion: `python3 gaussian/gaussian_bridge_exact.py` (from the packet root). It gives 6,655 counted exact checks (including the shear extension below) of primitive first-return formulas, the common observer bound on explicit parameter families, the covariance rank lower-bound matrix, and the nullness-transfer hostile. These controls support the new formulas; they neither prove FGNP_pl nor re-audit the inherited theorem.

## 8. Connection card

Source: NC2 one-sided nullcone. Target: termination of the formal inverse tower of a planar-origin nilpotent cubic reduction. Proposed map: a finite marked null-packet realization over the polynomial output-parameter ring. Preserved predicate: every marked coefficient identity at every power, not only scalar moments. Lost information: arbitrary full bidegree spaces and arbitrary special-image null bases are deliberately not transported. Required sidecars: planar provenance, output parameters, fixed support, individually null packets, stable-reduction reflection. Cheapest decisive test: construct a packet for a nontrivial composed family and verify the all-m identity symbolically; bounded numerical moment agreement alone is not decisive.

## 9. A constructive nonempty sector: shear words and their exact observer recursion

This section supplies a forward-data algorithm, not an encoding that first computes an unknown inverse truncation. It applies to a specified finite word of triangular shears. It does **not** assert that every Keller map admits such a word.

Take

    A_t(x,y)=(x+t p(y),y),
    B_t(x,y)=(x,y+t q(x)),
    F_t=A_t composed with B_t,

where p,q are arbitrary complex polynomials, also allowed to have coefficients in a reduced irreducible parameter ring. The inverse is obtained by reversing the two given elementary operations, so its coordinates are

    U_t=u-t p(v),
    V_t=v-t q(u-t p(v)).

Construct a single positive-charge packet P=Z and its two marked observers directly from the *forward* shear polynomials:

    Q_1(W;u,v)=u-W p(v),
    Q_2(W;u,v)=v-W q(u-W p(v)).

No moment solving, unknown inverse coefficients, or unproved truncation enters this construction. Since L(Z^m W^k)=m! delta_mk,

    sum_{m>=0} L(Q_1 P^m)t^m/m! = U_t,
    sum_{m>=0} L(Q_2 P^m)t^m/m! = V_t.

For m>=1 the exact coefficients are

    L(Q_1 Z^m) = -p(v) if m=1, and 0 otherwise;
    L(Q_2 Z^m) = -m (-p(v))^(m-1) q^(m-1)(u).

The second formula is zero whenever m>deg(q)+1. Its threshold is uniform in all coefficients of p and q, and strikingly independent of deg(p) for this ordering. If p is nonzero and q has degree d>=1, the last coefficient is generically nonzero at m=d+1. Special parameter values can decrease the degree, never increase the bound. The m=0 values are u and v.

### Exact local closure rule

For a given shear word, maintain two observers R_1(W;u,v),R_2(W;u,v), initially u,v. Traverse the word in the order that applies its elementary inverses to the output point. The inverse shear subtracting t p from coordinate i updates

    R_i <- R_i - W p(R_j),     R_j unchanged.

Thus every update stays inside a single common null base P=Z with polynomial observers. If the old W-degrees are D_i,D_j and deg p=d, the new degrees satisfy

    D_i_new <= max(D_i,1+d D_j),     D_j_new=D_j.

This is a proved recursive closure with a computable finite budget at every node; it requires neither a uniform degree bound over all words nor a global topological assumption. For nonzero p, the formula is a valid upper bound even when cancellations lower the actual degree. A zero shear is skipped. Initial budgets are D_1=D_2=0. This realizes the finite-interface proposal on all specified triangular shear words.

### Operator interpretation and precise scope

The map F_t has determinant one, but its displacement h_t(z)=z-F_t(z) normally depends nonlinearly on t. The formal residue identity extends verbatim to

    E_2(z_i exp(xi dot h_t(z))) = (F_t inverse)_i.

The exponential is t-adically defined because h_0=0. Therefore the packet identities above exactly match the **coefficient of each t power** in this contraction. They do not assert an identity E_2(z_i f^m)=L(Q_i P^m) for one t-independent xi-linear f: the higher-t terms in h_t mix its contraction powers. This distinction preserves the dimension/reduction boundary. A stable cubic reduction would need an additional transport rule before this shear-word construction becomes an FGNP_pl instance.

Reproduction from the published packet root: `python3 gaussian/gaussian_bridge_exact.py`, or `python3 -O gaussian/gaussian_bridge_exact.py`. All checks use explicit exceptions rather than Python assertions.
