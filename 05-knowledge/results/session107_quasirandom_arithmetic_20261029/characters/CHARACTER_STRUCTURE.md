# Squareclass colors inside the prime Legendre tournament

Session 107, 2026-10-29. **PROVED deductions from CITED classical inputs; exact finite guards; independent proof audit in `../audit_structure/AUDIT.md`.** No novelty or major-open-problem claim.

## 1. Independent audit of the inherited quasirandomness proof

Audited P1 of `05-knowledge/results/tournament_open_problem_routes_20261029/PROOFS.md` at commit `96b41ec57202ea9cc89dc301abd7e8053fbce179` ([immutable source](https://github.com/eliottcassidy2000/math/blob/96b41ec57202ea9cc89dc301abd7e8053fbce179/05-knowledge/results/tournament_open_problem_routes_20261029/PROOFS.md)). Its deduction is valid.

Let P(X) be primes at most X congruent to 3 modulo 4, N its size, and K(p,q)=(p/q), with diagonal zero. Quadratic reciprocity gives K transpose = -K. Heath-Brown's quadratic large sieve, restricted in both indices to P(X), yields B_X:=||K|| <= C_epsilon X^(1/2+epsilon). The coefficient vector can be arbitrary complex, so this really is an operator bound, not merely an estimate for the all-ones vector. Since N ~ X/(2 log X), B_X/N tends to zero for epsilon<1/2. The diagonal convention W=(J+K)/2 assigns 1/2 on the diagonal; collisions therefore must be removed, exactly as in P1.

In a telescoping replacement of one edge factor, fixing all other vertices leaves separate functions of the two endpoints. This separation is valid for simple directed graphs/tournaments, not arbitrary higher-order observables. It gives labelled k-pattern density 2^(-binom(k,2))+O_k(B_X/N+1/N). Multiplication by k!/|Aut(F)| is correct only for the unlabelled isomorphism class. It is not a law about the likelihood of one specified labelled pattern.

**Source verification:** the original author's Oxford PDF returned HTTP 403 in this audit. The exact cited quadratic large-sieve inequality was independently checked in Soundararajan–Young, *The second moment of quadratic twists of modular L-functions*, JEMS 12 (2010), §2.3, Theorem 2.4, printed p.1101 ([publisher PDF](https://ems.press/content/serial-article-files/31734)), which attributes it to Heath-Brown. We use the theorem unconditionally; the later GRH results in that paper are not used. Original reference: Heath-Brown, *A mean value estimate for real character sums*, Acta Arith.72 (1995),235–275, [DOI](https://doi.org/10.4064/aa-72-3-235-275).

## 2. Multi-fiber counting, with an explicit deterministic bound

**Theorem A.** Let S_1,...,S_k be nonempty subsets of P(X), n_i=|S_i|, and independently choose x_i uniformly from S_i. For any labelled tournament F on [k], put m=binom(k,2). Then

    |E product_(i->j in F) W(x_i,x_j) - 2^(-m)|
       <= (B_X/2) sum_(i<j) (n_i n_j)^(-1/2).                 (A)

The same estimate plus

    D = sum_(i<j) |S_i intersect S_j|/(n_i n_j)              (B)

bounds the absolute error for the probability of distinct vertices inducing F (without renormalizing). Conditional on distinctness, the bound is at most (right side of (A)+D)/(1-D), provided D<1.

**Proof.** Telescope edge factors to 1/2. Fix all other vertices for an edge ij. The remaining product is c f(x_i)g(x_j), where c,f,g lie in [0,1]. Extend f,g by zero off their respective fibers. Their Euclidean norms are at most sqrt(n_i),sqrt(n_j). Thus the normalized bilinear contribution of K/2 is at most B_X/(2 sqrt(n_i n_j)). Sum over edges. The union bound gives (B), and conditioning follows by dividing by the noncollision probability. QED.

Consequently any fixed number of fibers with n_i>=X^(1/2+eta) have independent fair limiting moving-vertex edges. Fibers may overlap, be arithmetically selected, or be adversarial. This assertion does not assume random vertex selection beyond the stated uniform sampling. More generally it suffices that the sum in (A) tends to zero after multiplication by B_X, and D tends to zero.

## 3. Exact named-port dependence is a squareclass code

For fixed nonzero integers a_1,...,a_r, exclude primes dividing their product and write

    a_j = (-1)^s_j 2^e_j product_(ell odd prime) ell^v_(j,ell) times a square,

with exponents in F_2. On prime q=3 mod4, chi_q(-1)=-1. Encode chi_q(a_j)=(-1)^y_j. Let V be the binary matrix whose columns are 2 and the distinct odd primes present, with rows their squarefree exponent vectors; let s=(s_j). Then

    y = s + V z.                                           (C)

**Theorem B.** As q runs through P(X), the sign vectors are asymptotically uniform over the affine image s+im(V). Every attainable vector has relative density 2^(-d), d=rank_F2(V); unattainable vectors never occur away from the excluded divisors.

**Proof of exact constraint.** Multiplicativity gives (C); square factors contribute +1. For every u in ker(V transpose), product_j chi_q(a_j)^u_j=(-1)^(u dot s). These relations also characterize the affine image by elementary linear algebra.

**Proof of sufficiency and density.** Set Q=8 product ell over the distinct odd prime columns. Among reduced residue classes modulo Q with q=3 mod4, chi_q(2) takes each sign equally: q=3 or7 mod8. For an odd prime ell, reciprocity gives

    chi_q(ell)=(-1)^((ell-1)/2) chi_ell(q).

Exactly (ell-1)/2 nonzero classes modulo ell have each Legendre sign. CRT independently combines these choices, proving uniformity of the latent sign vector z among the reduced classes with q=3 mod4. (An unused column for 2 may harmlessly be included.) The prime number theorem in each fixed reduced progression transfers this class counting to prime densities. Each vector in im(V) has exactly 2^(number_of_columns-d) preimages. QED.

The analytic input here is PNT in **fixed** arithmetic progressions, not merely Dirichlet infinitude; a precise statement is in [Tao's exposition, equation (3)](https://terrytao.wordpress.com/2009/09/24/the-prime-number-theorem-in-arithmetic-progressions-and-dueling-conspiracies/). No uniformity in growing Q is asserted.

**Corollary (colored rooted quasirandomness).** Fix the bases a_j and prescribe an attainable color to each of k moving prime vertices. By Theorem B every prescribed color fiber has size asymptotic to 2^(-d)N. Theorem A then says that their moving-vertex tournament is asymptotically uniform on labelled tournaments, conditionally on these colors and on distinctness. For colors c_1,...,c_k and labelled F, the unconditional probability is

    2^(-dk-binom(k,2)) + o(1).

The fiber-density error is only o(1) as stated: do not replace it by P1's power-saving error without an appropriate quantitative fixed-modulus PNT input. For fixed roots that themselves are distinct primes3 mod4, color signs are precisely incidences from root to moving vertex. Root-root incidences remain fixed, not random. Thus fixed named-root extensions ARE obtainable with an additional CRT/PNT argument; they do not follow from the spectral estimate alone. Roots growing with X or an unbounded number of ports require new estimates.

## 4. Cheap hostile examples and the requested numbers

* Bases (2,3,6): rank2, not3; chi_q(6)=chi_q(2)chi_q(3). Exactly four of eight patterns occur.
* Bases (4,2,5,3,7,13): rank5; the 4-port is identically +1, while the other five signs are asymptotically independent. There are 32, not64, possible patterns.
* Bases (2,3,7,5,2): rank4; first and last ports agree exactly. Cyclic placement does not change this rank or create a return-map identity.
* Bases (-1,2,-2): the first sign is -1 and the last is the negative of the second. This is an affine, not linear, sign code.
* Distinct primes (1093,3511,13,139,19) as **external bases** give rank5 and32 asymptotic colors. As tournament vertices only3511,139,19 are3 mod4:1093 and13 are not in this tournament's carrier. An external color is not a new vertex.

This gives an actual common mechanism (squareclasses, independent local residue choices, and spectral mixing) rather than a numerological merger. No implication about Monster/Baby Monster actions, Sylvester/Fermat primality, or Collatz convergence follows.

## 5. Multiplicative-order bridge: what a quadratic character retains

For odd prime p and p not dividing a, let d=ord_p(a). Then

    chi_p(a)=+1 iff d divides (p-1)/2,
    chi_p(a)=-1 iff v_2(d)=v_2(p-1).                         (D)

Indeed Euler's criterion identifies +1 with a^((p-1)/2)=1; the latter holds exactly when d divides (p-1)/2. Since d divides p-1, its failure is exactly equality of the2-adic valuations. In particular for p3 mod4, chi_p(a)=+1 iff d is odd. Character signs retain only this single2-primary gate, losing the odd part of the order and almost all lifting information.

Exact checks: ord_1093(2)=364 and ord_3511(2)=1755; both satisfy2^(p-1)=1 mod p^2 (the stated Wieferich condition), but chi_1093(2)=-1 and chi_3511(2)=+1. The common lifting condition does not force a common quadratic sign. These are verified finite examples, not a claim that there are no other Wieferich primes.

An affine3n±1 dynamics introduces addition. Multiplicativity of chi does not turn chi(3n±1) into a function only of chi(n). For example modulo7, n=1 and2 are both quadratic residues, but chi_7(3*1+1)=+1 whereas chi_7(3*2+1)=0; n=1 and4 give +1 and-1 with no zero exception. Thus even the sign quotient is not a deterministic Collatz state. Exact intercepts and admissible division valuations remain necessary sidecars.

## 6. Transfer ledger and verification

Source: integer bases and prime-indexed character matrix. Map: squareclass sign colors and normalized tournament kernel. Preserved: exact multiplicative sign relations; averaged finite patterns within large attainable fibers. Destroyed: full orders, prime-power lifting, additive intercepts, named-prime identities, rare/growing-root fibers. Sidecars: squareclass matrix, modulus Q, root scope, numerical height X, fiber cardinalities, and automorphism groups for unlabelled counts.

Run `python3 verify_characters.py` in this directory. It uses integer arithmetic only, verifies order gates and examples, exhaustively checks finite CRT color counts for manageable base lists, checks squareclass constraints on prime samples, and checks exact labelled/unlabelled three-vertex normalization. Prime-prefix frequencies are diagnostics, not proofs of convergence or an effective large-sieve constant.

### Inheritance and sparse conditioning hostile

The earlier `root27_machine_bridges_20261026/anchor/LEGENDRE_SIDECAR.md` at the audited commit was also read directly. Its CRT/Dirichlet argument establishes extension existence, and its change-of-enumeration example correctly shows why existence is weaker than prime-prefix densities. Theorem B upgrades fixed-root existence to density with PNT in fixed progressions; Theorem A then combines that density with P1's spectral input. This is the specific added analytic operation.

The read local `session100_selective_collatz_targets_20261029/bijections/VALUATION_CONJUGACY_AND_INTEGER_SIDECAR.md` proves valuation-word conjugacies while explicitly retaining positive-integral-source conditions. The additive hostile above is a weaker, elementary instance of the same loss-of-information warning; it does not subsume that conjugacy theorem.

The parent session supplied an exact Sylvester hostile: later Sylvester terms are1 modulo each earlier term. Hence their prime terms3 mod4, whenever present, have a fixed reciprocity orientation and form a transitive subtournament, not a quasirandom one. No infinitude assertion is required. The cardinality hypothesis of Theorem A is essential and is not supplied for this sparse recursive subsequence. In particular1807=13*139 has chi_p(13)=chi_p(139) for earlier prime terms p=3,7,43. Reciprocity makes chi_13(p)=-chi_139(p), because the two factors differ modulo4. The script checks these six incoming/outgoing identities exactly. The composite product constraint is real, but changing the direction of a named port needs its reciprocity sign sidecar.
