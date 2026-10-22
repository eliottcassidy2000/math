# Independent hostile audit: additive prime-power closure

Verdict: **PASS, proof complete with the explicit conventions below.** The candidate proves, for every odd prime p and r >= 1,

s(n + phi(p^r)) = (2/p) s(n) modulo p^r for every n >= 2r - 1.

This is exactly the currently open statement HYP-9358(d′), not the already obtained multiplicative dilation congruence. The fuller proof in `session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md` supplies the explanatory details identified below. No repository files were edited and no Git action was taken by this audit.

## Proof audit

1. **Coefficient extension.** An unramified quadratic extension R of Z/p^r is sufficient because 2(p-1) divides p^2-1. Its constants inject, and an element is a unit exactly when its residue is nonzero. Taking the quadratic extension uniformly avoids the false assertion that Z/p^r contains every 2(p-1)st root when epsilon=1. It generally does not. Distinct Teichmuller roots have unit differences.
2. **Teichmuller scalar.** The lift omega of 2 has omega^(p-1)=1 and omega^((p-1)/2)=epsilon. The ratio 2/omega is a principal unit; choose u with 2=omega(1+pu). The nonuniqueness of u modulo p^(r-1) has no effect. For r=1 any u works.
3. **Exact root support.** The exponent difference in t0(n+p-1)/t0(n) is (p-1)n+(p-1)(p-2)/2, giving epsilon because p-2 is odd. The length L=2(p-1) Fourier transform is invertible as L is a unit. The shift identity forces support in alpha^(p-1)=epsilon: the only other eigenvalue is -epsilon and the difference is the unit 2. There is no need for all allowed coefficients to be nonzero. Their sum is t0(0)=1.
4. **Additive line.** Fix a root theta with theta^(p-1)=epsilon. Each supported alpha/theta is a (p-1)st root and belongs to the base ring. Hence all support is inside G=theta(Z/p^r), an additive cyclic group of order p^r. Multiplication and inversion may introduce any point of this group, not just its unit points. The group is closed under addition; the Fourier node set by itself is not.
5. **Augmentation-unit inverse.** In R[G] modulo p, Frobenius to the p^r power gives A^(p^r)=sum c_alpha^(p^r)[p^r alpha]=(sum c_alpha)^(p^r)[0]=[0]. Coefficients need not be fixed by that Frobenius: summing first gives 1 anyway. Thus A^(p^r)=1+pV in the unreduced group ring. Raising to p^(r-1) annihilates the error by the ordinary principal-unit induction, yielding A^(p^(2r-1))=1. This is a sufficient exponent, not an asserted exact order. It also works at r=1. Consequently A inverse is finite and has support in G.
6. **Hurwitz carrier.** Define E_beta(n)=beta^n with 0^0=1 and delta_d supported at d. The identities E_beta*E_gamma=E_(beta+gamma) and delta_d*delta_e=binom(d+e,d)delta_(d+e) hold over every commutative ring. The group-ring map is unital, since E_0=delta_0. Injectivity is never required: a homomorphic image of an actual inverse remains an inverse. This avoids any unjustified assumption that the kernel of an image reduction has square zero.
7. **Integral Newton expansion.** For f_j(n)=binom(binomial(n,2),j), finite differences at 0 are integers because all its nonnegative integer values are integers. Since its rational degree is at most 2j, the expansion in binom(n,d) stops at d=2j. This remains valid for j>=p and d>=p; there is no illegal division by j!, d!, or p in the coefficient ring. Explicitly, binom(n,d)alpha^n=alpha^d(delta_d*E_alpha)(n). The factor alpha^d must be included in writing out the lift.
8. **Filtration and inverse.** The span of p^j delta_d*E_beta for j<r,d<=2j is a ring. Product grades and degrees add; terms of grade >=r vanish. The correction U=T-A is represented entirely in positive grades. Every power (A inverse U)^k has grade at least k, hence the geometric inverse terminates at k=r-1. Its degree never exceeds twice its total grade. This is a representation argument, not a uniqueness argument for grades. Nonunique representations cannot invalidate closure or nilpotence.
9. **Identification and signs.** A tournament has uniquely ordered strong components. Taking its first nonempty strong component gives T_n=sum_(k=1)^n binom(n,k)s(k)T_(n-k), n>0. With S_0=0 this is S*T=T-e, so S=e-T inverse. The constant term is 0 and the sign is exactly as in the candidate. The inverse is unique because T_0=1.
10. **Nonunit tail.** Write beta=theta b. If p divides b, a generator has valuation at least j+n-d >= n-j >= r at n>=2r-1. Here n>d automatically as d<=2r-2, eliminating any 0^0 issue for beta=0. The same applies at n+P. This includes the identity's point mass at 0.
11. **Unit tail.** If b is a unit then b^P=1 and theta^P=epsilon, hence beta^P=epsilon exactly in R. For j=0 necessarily d=0. For j>=1, Vandermonde and l*binom(P,l)=P*binom(P-1,l-1) give v_p(binomial(P,l))>=r-1-v_p(l). For l>P the term is zero and harmless. Since l<=2j<p^j, v_p(l)<=j-1. Multiplication by p^j therefore kills every difference modulo p^r. The inequality p^j>2j holds for all odd p and j>=1, including p=3,j=1. The proof would not work unchanged at p=2.
12. **Descent and boundaries.** The base-ring injection proves an ordinary congruence of integers, not just an extension-ring equality. At r=1 only grade 0 survives and the zero exponent disappears at n>=1. Do not infer validity at n=0: p=5 gives s(4)=24=4 modulo5, whereas s(0)=0. Nor claim sharpness for each pair: p=3,r=3 also works at n=3,4 below the stated sufficient bound 5.

## Source reconciliation

Read the current `math/05-knowledge/hypotheses/HYP-9358-...md` and `math/05-knowledge/results/tournament_pseudoprimes_audit_20261021.md`, especially Lemma B and its distinction between the prime-square case and the higher-power conjecture. Incoming commit 3596c5d0b1990b05ae960e980a863956b47cc6b5 already proves the r=2 additive case and records the exact general threshold proposed here. The new proof extends that finite derivative correction using integral Newton coefficients and a finite exact additive group carrier. The previous multiplicative theorem was not itself a closure of (d′). Other open parts of HYP-9358 remain open, and connected-graph prior-art cautions remain intact. No external novelty claim has been checked or made.

## Independent executable controls

`check.py` is fresh independent Python code, importing no producer or prior-audit code. It uses explicit runtime checks rather than assert, so optimization does not suppress checking.

- A direct modular strong-component recurrence verifies 706 additive shifts: p=3,r=1..5; p=5,r=1..4; p=7,r=1..3; p=11,13,17,19,r=1..2. Every tested index in the claimed range passes. All lower-index failures are retained separately.
- Finite differences produce integral Newton coefficients for j=0..12; polynomial degree and values beyond the interpolation interval are checked exactly.
- Exact rings (Z/p^r)[sqrt(D)] with nonsquare D verify Fourier expansion and support, augmentation 1, the group-algebra inverse exponent p^(2r-1), and the complete Newton-to-Hurwitz reconstruction of T. Cases: (p,r)=(3,1),(3,2),(3,3),(5,1),(5,2),(7,1),(7,2), covering both epsilon signs and genuinely nilpotent/non-field coefficient rings.
- Normal and `python -O` output files agree byte-for-byte. These are controls, not substitutes for the proof.

Reproduce:

    python session25-additive-review/check.py > session25-additive-review/check.out
    python -O session25-additive-review/check.py > session25-additive-review/check.optimized.out
    cmp session25-additive-review/check.out session25-additive-review/check.optimized.out

All hashes, including pinned source proof/candidate and incoming canonical sources, are in SHA256SUMS.json.

## Final frozen-file extension: §§7–8

**Expanded verdict: PASS for the entire final file**, SHA256 `ad3c9b993debf5cb2403d414eb532df9a19a8a260d8ddc1705f41545664451c6`. The producer explicitly confirmed this final freeze. This supersedes the earlier core-only source hash.

- Any integer base a coprime to p works, including negative a: Teichmuller decomposition and finite binomial expansion are unchanged; the nonnegative exponent n(n-1)/2 creates no ambiguity. The scalar multiplier becomes (a/p). The a=3 semicomplete interpretation excludes p=3.
- The shift identity is exact, including d=0 and n<d: Delta(delta_d*E_beta)=delta_(d-1)*E_beta+beta delta_d*E_beta. Degrees do not increase. Thus F is stable under shifts as claimed.
- For arbitrary f in F with f(0) a unit, a chosen grade-zero lift A_0 has unit augmentation modulo p. Its reduction is a unit in the cyclic p-group algebra (normalize the augmentation to 1 if desired), hence A_0 itself is invertible. The remaining terms have positive grade and the same terminating geometric inverse works. No unique filtration representation is needed.
- The closure claim is only for the stated Hurwitz operations, with scalars multiplying delta_0. It does not extend to pointwise products or to adjoining arbitrary sequences such as the all-ones sequence. Every resulting expression remains in the same fixed F and inherits the same epsilon, despite closure under Hurwitz multiplication; this is consistent because exponential frequencies add within G.
- The connected-graph recurrence gives Delta b=(Delta T)*T inverse. Since its index is one less than the vertex count, the threshold is n>=2r, exactly as stated. This is sufficient and not claimed minimal or stronger than every known literature bound.
- The CRT representative t has lower bound 2r-1 and upper bound phi(p^r)+2r-2<=p^r: the last step uses 2r-2<=p^(r-1), including r=1. Since p^r divides n, t<=n and n-t is a nonnegative multiple of phi. Repeated congruences multiply D(t) by a unit, preserving vanishing.
- For even n the rotation subgroup of order 2^a, where 2^a exactly divides n, acts freely on tournaments: a nontrivial stabilizer in a cyclic 2-group contains its unique involution, the half-turn, which swaps two endpoints of an edge and cannot preserve a tournament. Strongness is invariant. Thus 2^a divides s(n); c(n)=0. This justifies the automatic 2-part and the empty list for pure powers of two without invoking an odd-prime theorem at p=2.
- For n=p^r the representative is p^(r-1), using p^(r-1)>=2r-1 (equality at p=3,r=2; r=1 also works). For general m, the same bound holds for mp^(r-1), and the two multiplicative indices differ by m phi(p^r), producing epsilon^m. These are predecessor-level relations, not the invalid replacement of every prime-power test by D(p) at full precision.
- D is explicitly the aggregate s-c, and no individual Hamiltonian or metagraph statistic is implicated. The criterion has no squarefree/cube-free restriction. Other unrelated open conjectures remain open.

Additional independent `corollaries.py` controls pass: 714 arbitrary-unit-base shifts (including negative bases and base1), 96 connected-graph shifts, and the full CRT criterion at every composite from4 through200 (153 cases). Both normal and optimized runs agree byte-for-byte. These controls are separate from the producer's larger tests.
