# Session52 independent MRB / CM / mod-15 audit

Status: FINAL PASS. The producer applied both scope repairs and integrated the audited Fourier and Bolza additions. The final 14-file producer manifest is 5cf13f7a96cee50502f81f0e9e74624a6454f5702a76aa72996d5a34075d6102. The main mathematics, numerical certificates, finite universes, and source-dependent implications pass. No publication was performed.

## Verdict and useful result

The final producer bundle was copied using only its 14-entry SHA-256 manifest plus the manifest itself, then replayed in both normal and optimized Python. Both gave **157,050 active checks**, identical output and result bytes, and unchanged manifest hashes. The earlier 12-entry snapshot is retained under `pre-amendment-manifest-copy/` as explicitly superseded audit history. All eight producer negative-control executions rejected their deliberately false claims.

An independently written checker adds **3,725 active checks**. It computes Hasse invariants rather than Legendre-sum point counts, checks the entire claimed all-j census, proves the finite root-of-unity projection in an exact cyclotomic quotient, and evaluates the MRB constants with a different Euler split and explicit finite differences. Normal and optimized executions have identical results and logs. Twelve additional semantic-mutation executions reject false Ogg membership, the wrong curve, an inclusion-exclusion sign error, a density substitution, omission of j=0, and conflation of a quadratic order with the maximal order.

Both originally requested wording/domain repairs are now present and verified in the final producer text:

1. Theorem 2 now says “each nonzero one-sign residue subseries diverges when separated from its opposite-parity partner.” Zero-weight residue subseries converge to zero, including the nonunit residues of the unit mask and every residue when w=0.
2. Corollary 3 now specifies “For positive odd d.” The nth-root profile here is defined on positive arguments.

Neither repair changes a substantive formula, certificate, or conclusion. The final candidate input bytes are pinned in `actual-manifest-copy/`; the parent owns integration. All final additions have been reviewed: the cyclotomic reduction computes the stated projector, the twist is verified at every claimed prime, and the bounded inherited excerpt is literally present in the source with the recorded SHA-256. The final `RESEARCH.md` SHA-256 is c2ff3f80d7108289ced1469b845aab14238d6e9754bb9ddff5d202e161e1ba39.

## Important recovered prior-work connection

The inherited canon already contains the same **geometric** CM curve, with a rational-twist distinction that must not be erased. In `THM-4012-weighted-leading-face-good-elliptic-factor-observer.md`, Sections 5.1–5.2, the Bolza curve

    C_B: y^2=x(x^4+1)

has quotient

    v^2=(u+2)(u^2-2).

After shifting the coordinate, this is

    E_minus: y^2=x^3-4x^2+2x.

That theorem computes c4=160, Delta=512, j=8000 and constructs a degree-two endomorphism alpha with alpha^2=[-2]. It also proves the geometric isogeny Jac(C_B) ~ E_minus^2, with carefully limited use in its planar-Jacobian argument.

The present E_plus: y^2=x^3+4x^2+2x is the quadratic twist of E_minus by -1. The explicit complex isomorphism is

    x_minus=-x_plus, y_minus=i y_plus.

Substitution proves it directly. At an odd prime p, changing x to -x in the Legendre sum gives a_p(E_minus)=(-1/p) a_p(E_plus); therefore the zero-trace and supersingularity predicates agree. The independent checker verifies the corresponding Hasse-coefficient relation for all 167 odd primes below 1000.

This is a real inherited connection through a concrete elliptic isomorphism class and its endomorphism algebra, not through the number 15. The final producer synthesis now credits this prior result. Its face-stability assumptions and Jacobian-conjecture conclusion do not transfer to the MRB constant. The present bridge still starts from A(8), not from C or D.

Source: `session50-publication/clean-repo/01-canon/theorems/THM-4012-weighted-leading-face-good-elliptic-factor-observer.md`, especially equations (33), (37)–(46). Current inherited checkout: `f50a8c2f3a307869ec51650851dbc43eaedfc176`.

## Mathematical proof audit

### The actual MRB functional

For f(n)=exp(log(n)/n)-1, f decreases to zero for n>=3 and f(n)=log(n)/n+O(log(n)^2/n^2). Thus C=sum (-1)^n f(n) converges conditionally. The unregularized series with summands (-1)^n n^(1/n) has no ordinary limit. Its even partial sums have limit C, but that is a different statement.

For periodic real a_n with mean mu, write a_n=mu+b_n. The partial sums of b_n are bounded. Summation by parts makes sum b_n f(n) convergent. The difference f(n)-log(n)/n is absolutely summable, while the integral-test error for sum log(n)/n converges. Consequently

    sum_{n<=N} a_n f(n)=(mu/2)(log N)^2+c_a+o(1).

This proves, and slightly sharpens, the producer's bounded-error statement. With an unsigned unit mask the coefficient is 4/15, not 8/15; the factor 1/2 comes from integrating log(x)/x.

### Odd-modulus projection and legal regrouping

Translating a 15-periodic w by 15 reverses (-1)^n w(n). The signed mask therefore has zero mean on 30 terms. For the unit mask the producer excludes every proper divisor of 30 as a period; since 30 is a period this proves the fundamental period exactly.

For each r in 1..15, the sequence f(15k+r) decreases eventually to zero. Its alternating sum converges. Pairing adjacent terms in that subseries gives differences f(30j+r)-f(30j+r+15)=O(log(j)/j^2); therefore the paired series is absolutely convergent. Taking a finite number of residue limits is lawful. No general rearrangement theorem for conditionally convergent series is invoked.

The identity u(n)=1-1_(3|n)-1_(5|n)+1_(15|n) holds for every n. Because all divisors used are odd, (-1)^(dk)=(-1)^k. Applying the identity first at finite N, retaining floor(N/d), and then taking four convergent limits proves both C_u=M(1)-M(3)-M(5)+M(15) and its complementary formula.

The rational generating function R(z) in the packet is the generating function of the signed **mask**, not of f(n). Its coefficients follow by a geometric series on each residue. Both rational forms and all signs are correct.

### Final integrated fixed-root-of-unity formula

A genuine full-profile Fourier formula is available without any problematic boundary interchange. Define

    G(z)=sum_{n>=1} f(n) z^n,  |z|<1,
    zeta=exp(2 pi i/15),
    c15(j)=sum_{r in (Z/15Z)^*} zeta^(-jr).

Then

    C_u=(1/15) sum_{j=0}^{14} c15(j) G(-zeta^j).

Here the values on the unit circle mean their ordinary convergent series. Every -zeta^j differs from 1 because 15 is odd. Geometric partial sums are bounded by 2/|1+zeta^j|, so Dirichlet's test gives convergence of each series. The Fourier projection holds at every finite cutoff; the outer sum has only fifteen terms, so taking its limit is legitimate. Abel's theorem also identifies these boundary values with radial limits if that formulation is desired.

The exact coefficient vector for j=0..14 is

    [8,1,1,-2,1,-4,-2,1,1,-2,-4,1,-2,1,1].

It follows from c15(j)=sum_(d|gcd(15,j)) d*mu(15/d). Fourier orthogonality proves the identity for all n; the audit additionally checks the projection for n=1..120 by integer polynomial reduction modulo Phi_15. The j=0 contribution is (8/15)C, but the remaining fourteen modes do not vanish. Thus the density rule would wrongly delete real information. This formula has now been integrated as Corollary 4 of the final producer packet, with its exact finite certificate.

### CM bridge and quantifiers

The AGM calculation gives K(k')/K(k)=sqrt(2) at k=sqrt(2)-1, with lambda=k^2=3-2sqrt(2). The polynomial lambda^2-6lambda+1=0 reduces the standard j formula to 8000. The normalization from A(8) is correct. The lattice Z+Z sqrt(-2) is a module over the maximal quadratic order Z[sqrt(-2)], so its elliptic curve has CM discriminant -8. The rational E has c4=160 and Delta=512; hence every odd prime is a good prime.

Deuring's criterion applies at p>=5. If needed, formulate it after extending the number field so that the geometric CM endomorphisms are defined; supersingularity is unchanged by that base-field extension. There are no odd ramified primes in Q(sqrt(-2)). Therefore the non-splitting criterion is exactly inertness, or (-8/p)=(-2/p)=-1, namely p=5 or 7 modulo 8. The p=3 ordinary calculation is a separate finite check.

The points at p=17 and p=37 respectively have counts 24 and 38, giving traces -6 and 0. They disprove both implications between Ogg membership and supersingularity of this one curve. No MRB-defect identity follows.

### Ogg counts and mod 15

All eight unit classes contain both an Ogg prime and a non-Ogg witness. This proves the impossibility of an exact mod-15 membership function on primes, and also proves that a sound mod-15 rejection rule retaining every Ogg prime >5 retains every prime >5. The same residue classes contain both CM reduction types. The 32 mod-120 classes are an exact CRT product, not a probabilistic assertion about independence of primes.

The total/rational supersingular counts are correctly distinguished. For p>3 the reduced-form calculation uses the nonmaximal order of discriminant -4p when p=3 mod 4. The source explicitly gives the combined class-number formula as well as its piecewise equivalent. No class number at p=2 or p=3 is certified by that p>3 routine; those Ogg-list entries remain cited.

Aricheta's full Atkin–Lehner level 15 quotient and Fricke-only level 15 quotient are different. Their cited rationality lists are {2,7} and {2}, respectively, for p not dividing 15. These are level-structure rationality predicates, not residue-class projections.

### Interval proof and entire-perturbation boundary

The Cauchy disk |z-x|=x/2 stays in the principal-log domain. For x>=M>=4, |Log(z)| is bounded by log(3x/2)+pi/6, and its quotient by x decreases. Cauchy's derivative estimate, the finite-difference unit-cube identity, and first-term-plus-integral summation yield the printed Euler remainder, including its factor 2 and 1/d scaling. The finite Euler sum uses exact integer weights; interval transcendental operations account for rounding. The alternate audit uses M=512,N=128 and explicit successive differences instead of those weights. Its four resulting M(d) intervals lie strictly inside the producer intervals. The interval library is not a proof kernel; this limitation is already disclosed.

The positivity comparison D<q^4 rules out the unweighted positive even-theta tail and only that proposed identity. The entire perturbation h_K=e^(-z) product_(j=1)^K(z-j) preserves the listed branch/prefix/asymptotic data and changes the alternating functional by the stated nonzero exact amount. It need not remain an nth-root profile or obey stronger uniform half-plane bounds. The producer retains these essential restrictions.

## Primary-source verification

All URLs below were retrieved during this independent audit on 2026-10-24.

- [Crandall, Section 7.5](https://www.marvinrayburns.com/UniversalTOC25.pdf), printed page 28: the MRB definition uses (-1)^k(k^(1/k)-1). The derivative representation is historical context, not used as a proof of modularity.
- [DLMF 19.8.5](https://dlmf.nist.gov/19.8#E5) and [DLMF 23.15](https://dlmf.nist.gov/23.15): verify the AGM, modulus, nome and lambda conventions. Both moduli here lie strictly between 0 and 1.
- [Mula–Murru–Pintore, Theorem 4.1](https://link.springer.com/article/10.1007/s10231-024-01528-x): states the p>=5 CM good-reduction criterion. The article's preceding prose has an apparent reversed discriminant-valuation phrase; the present work correctly uses Delta=512 to establish good reduction and does not import that typo. Its Theorem 4.18 provides the Hasse-invariant criterion used by the independent finite check.
- [Delfs–Galbraith](https://arxiv.org/pdf/1310.7789), Sections 1–2, especially Proposition 2.5: supplies total and prime-field supersingular counts, the p>3 trace-zero setting, and the order-class-number expression. This resolves the maximal/nonmaximal-order issue without assuming h(-4p)=h(-p).
- [Aricheta](https://arxiv.org/html/1809.07421v2), Introduction and Appendix A: verifies the Ogg equivalence and Monster prime support, defines full versus partial Atkin–Lehner notation, and tabulates the two distinct level-15 lists. The prime-to-level restriction is explicit in the setup.
- [Elkies–Ono–Yang, Section 3](https://people.math.wisc.edu/~tonghaiyang/deuring.pdf): independently corroborates the CM reduction mechanism. It is not needed to extend the stated p>=5 range.

- [Florit–Smith, Section 4.14](https://arxiv.org/html/2101.00917v2#S4.SS14), independently retrieved on the final pass: corroborates the Bolza normal form and its isogeny to a square of a j=8000 elliptic curve. The characteristic-zero application is within its stated characteristic restriction; this is not used to extend any small-characteristic claim.

## Reproduction and limits

From this audit directory (commands preserve the pinned file closure):

    python3 verify_audit_manifest.py
    python3 actual-manifest-copy/verify_manifest.py
    python3 actual-manifest-copy/check.py > final.normal.out
    cp actual-manifest-copy/results.json final.results.normal.json
    python3 -O actual-manifest-copy/check.py > final.optimized.out
    cmp actual-manifest-copy/results.json final.results.normal.json
    cmp final.normal.out final.optimized.out
    python3 actual-manifest-copy/test_rejections.py > final.rejections.normal.out
    python3 -O actual-manifest-copy/test_rejections.py > final.rejections.optimized.out
    cmp final.rejections.normal.out final.rejections.optimized.out
    python3 actual-manifest-copy/verify_manifest.py
    python3 independent_check.py > independent.normal.out
    cp independent_results.json independent_results.normal.json
    python3 -O independent_check.py > independent.optimized.out
    cmp independent_results.json independent_results.normal.json
    cmp independent.normal.out independent.optimized.out
    python3 semantic_mutations.py > semantic_mutations.out
    python3 verify_audit_manifest.py

`verify_audit_manifest.py` verifies every outer payload hash, byte length, and complete recursive file-set equality. The outer manifest includes the final nested producer manifest and all fourteen of its payloads. Final replay logs are `final.normal.out`, `final.optimized.out`, and the `final.rejections.*.out` pair; they match the producer’s pinned outputs. The final normal/-O independent and semantic checks were rerun against the amended candidate. `independent_results.json` stores Hasse and interval data. `semantic_mutations.json` records the twelve rejection reasons. Direct ordinary partial sums in the independent checker are numerical consistency checks, not additional high-precision interval certificates.

No finite computation proves the all-prime CM theorem, the Ogg classification, or analytic convergence; the proofs and imports above do that work. No result determines the arithmetic nature of C or identifies its defect D with a modular object. No git state or producer artifact was changed by the auditor.
