# Cooper weight-one prime dichotomy: recovery, audit, and characteristic-three consequence

Status: RECOVERED REPOSITORY PROOF + independently reproduced FINITE-EXACT tests. The characteristic-three corollary below is PROVED from the cited Lucas theorem; no novelty claim is made. No repository files edited.

## Exact object and quantifiers

Let s7(n)=sum_{k=0}^n binom(n,k)^2 binom(n+k,k) binom(2k,n), and let G(t)=sum c(n)t^n be the unique normalized formal square root of F(t)=sum s7(n)t^n, with c(0)=1. Then c begins 1,2,22,336,6006. It is not the sequence of numerical square roots of individual s7(n).

The claim is:

For a prime p, [for every integer n>=0, c(pn)=c(n) modulo p²] iff [p=7 or p mod7 is 1,2,4].

The negative direction means there exists n for each excluded prime, not failure for every n. In particular n=0 always succeeds. Prime 2 is INCLUDED because 2 is a square modulo7. Prime3 is excluded: c(3)-c(1)=334 is 1 modulo3. Prime7 is a separate ramified case, not a nonzero quadratic residue.

Integrality is genuine over Z. Write h=(eta(7tau)/eta(tau))^4, t=h/(1+13h+49h²), and theta=sum_{x,y}q^(x²+xy+2y²). The known modular parametrization G(t(q))=theta(q), t(q)=q+q²Z[[q]], implies q(t) in t+t²Z[[t]] and G=theta(q(t)) in Z[[t]]. Merely taking the square root of an arbitrary integral series would not establish this, particularly at 2.

## Inheritance and source status

The exact full iff is already Theorem D, sections 5.1–5.8, of 05-knowledge/results/cooper_s7_weight_one_root_frobenius_20261018.md (copied here as inherited-proof.md). It identifies itself as a proof of Chan–Cooper–Sica Conjecture5.4 together with O'Brien Conjecture9.3. This session must not call it a new theorem. THM-4669 is a different half-period comparison with Apéry numbers. THM-4686 is a different CM harmonic result. Neither proves the present all-n statement by itself. MISTAKE-629 warns against mistaking classical equivalent guises for novelty; its named affected claims do not retract Theorem D.

Primary-source check on 2026-10-24:

- Cooper, arXiv:2302.00757v2, section11, equations11.2–11.3: explicitly distinguishes the square-root sequence from s7, attributes the split-prime assertion to CCS2010 and p7 to O'Brien, and reports the square-root conjecture open as of that version. https://arxiv.org/html/2302.00757v2#S11
- O'Brien's 2016 thesis, chapters6 and9, gives the quadratic recurrence and Conjectures9.2/9.3. The printed phrase “p=0” means the zero residue class modulo7 in the surrounding classification, hence prime7, not a prime zero. https://mro.massey.ac.nz/server/api/core/bitstreams/16e1c23a-dd52-4f08-93bc-6fb17ed7d79b/content
- Beukers–Tsai–Ye, arXiv:2408.16616v2, Theorem1.3: sufficient weight-one Lucas congruences, including the permitted primes here. This is modulo p, NOT the full p² theorem. https://arxiv.org/html/2408.16616v2#S1.Thm3
- Beukers, arXiv:2403.03301v3: Proposition4.2 provides the weight-one polynomial mechanism; Remark1.14(2) still treats the s7 excellent-lift congruence as experimental beyond his main theorem. Thus the inherited proof needs its double-cover adaptation, not blind citation of that main theorem. https://arxiv.org/html/2403.03301v3

Searches found no later primary publication proving this precise p² iff. That is a bounded negative search, NOT proof of current literature openness. The CCS author-hosted link reached via OEIS failed to open; its exact attribution is independently verified in Cooper's own survey. The repository proof and externally published theorem status must remain distinct.

## Audit of inherited mechanism

Put s²=Q=(1+t)(1-27t), s(0)=1, H=G/s, C_p(sum a_nt^n)=sum a_{pn}t^n. The inherited weight-one Frobenius theorem gives G=G_p G(t^sigma) modulo p² for split p; at inert p it instead gives G=s H_p G(t^sigma). The Fricke involution on X0(7) supplies exactly this two-sheet distinction. The excellent lift t^sigma=t(q(t)^p) must NOT be silently replaced by t^p modulo p².

For split p, write t^sigma=t^p(1+p epsilon). Cartier extraction yields
C_p G = G + p C_p(epsilon G_p) (t d/dt)G modulo p².
The inherited Eisenstein/Fricke argument proves epsilon G_p modulo p is a polynomial of degree<p with zero constant coefficient, killing the error term. The degree bound is floor(p/3)+floor((p+1)/3)<p. All split odd primes begin at11. At2 the separate lift identity is t(q²)=t²+2t³ modulo4. At7 use f=theta³s, U7f=f modulo49, and t(q^7)=t^7+14t^8+28t^9 modulo49. These are not covered by the split-prime Eisenstein division argument.

The inert impossibility is particularly robust: if C_pG=G modulo p, the inert factorization forces C_p(sH_p)=1. But
C_p(sH_p)=s P/Q, P=C_p(H_p Q^((p+1)/2)) in F_p[t].
It follows s=Q/P is rational, impossible because Q has distinct simple roots for inert p>=5, and Q=1+t for p3. This establishes failure somewhere without needing the experimentally observed first witness n=1 for every inert p.

I checked the algebraic reductions, the p2/p7 separation, and the Beukers polynomial statement against the primary source. I independently recomputed all three exceptional lift identities to q^150 and the weight-three Eisenstein identity to q^149. The inherited claimed Sturm bounds are24,32,140 and2 respectively. The first three modular forms used to clear denominators have weights12,12,30 and levels14,21,49. These checks support the stated Sturm proof, but are not a fresh full audit of every modular transformation, pole divisor, or integral-model assertion of the inherited LemmaT. Its original PARI checker was recovered from git, but PARI/GP is absent here and that complete checker was not rerun. Treat the full iff as an inherited proof with a documented partial independent audit, not as newly independently certified in every dependency.

## PROVED characteristic-three corollary

This uses only the established Lucas property for s7, not the stronger inherited p² theorem. Since F_3=1+t modulo3, Lucas gives F(t)=(1+t)F(t³). In F3[[t]], normalized square roots are unique, so G(t)=s(t)G(t³), s²=1+t. Let H=G/s. Then H=G(t³), hence C3H=G. Moreover s=(1+t)²/s(t³), so C3s=1/s; consequently C3G=H.

Thus C3²G=G, and for every n>=0:
- c(9n)=c(n) modulo3;
- c(27n)=c(3n)=d(n) modulo3, where H=sum d(n)t^n;
- even powers of3 return G and odd powers select H.

This is a two-state swap, NOT multiplication by -1 on G. At n1, c3=336=0 modulo3 whereas -c1=1; c27 is also0. Modulo9 the simple swap already fails: c3-d1=336-15=321=6 modulo9.

Iterating G=sG(t³) once gives G=(1+t)²G(t^9). Therefore
G(t)=product_{j>=0}(1+t^(9^j))² in F3[[t]], and G(t)^4=(1+t)^(-1).

Precisely, c(n) is zero modulo3 unless every base9 digit of n belongs to {0,1,2}. On that support it equals (-1)^(number of base9 digits equal to1). The product proof avoids ambiguity about fractional powers in characteristic3. Among 0<=n<9^r, exactly3^r coefficients are nonzero: a rigorously defined digit-support self-similarity of density tending to0. It is not a Collatz statement and supplies no orbit coverage or convergence claim.

Targeted searches of the materialized theorem/results files for c(9n), base9-digit and Cartier-swap formulations found no exact prior statement. A full git-object grep was stopped after it stalled fetching omitted objects; that search was incomplete. Web searches for A279619 modulo3 and Cooper square-root Cartier found no primary match. The corollary is elementary from known Lucas data, so no originality claim is justified.

## Reproduction and boundaries

Run python3 check.py and python3 -O check.py; outputs match byte-for-byte. Exact integer recurrence through index8000; independent binomial-convolution square verification through100; all46 primes<200 and n=0..40 (1886 pairs); all included primes pass p², every excluded prime fails modp already at n1 in this finite universe. The digit law is checked through8000, and c9n=c_n through n888. A deliberately corrupted c22 is rejected in both optimization modes, without reliance on assert.

Run python3 check_modular.py and python3 -O check_modular.py; outputs match byte-for-byte. All arithmetic is exact modular integer arithmetic. These finite tests alone do not prove any infinite quantifier. Positive, inert, ramified, small-prime, wrong-normalization, signed-law, and corruption controls are included.

