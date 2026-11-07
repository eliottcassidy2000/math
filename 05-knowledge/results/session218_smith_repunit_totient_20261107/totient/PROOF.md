# Session 218: a certified exponent-to-nontotient bridge

Status: PROVED below, with finite exact regression. These are elementary derivations, not claims of historical novelty. Parent integrates/publishes; this scratch packet does not modify the source repository.

## Inheritance and scope

Read-only source: session216-stage at parent-specified 88ff4b0afd38315a4ad96db6164b80b55488e649. AGENTS, research protocol, active guards and core-papers register consulted. Closest mechanism is THM-4746-camion-pseudoprimes..., section 5: simultaneous Mersenne/Wagstaff primes give a triangle-class Korselt certificate for 4^q-1. Hostile controls: 341 fails Camion while 21 passes; the old squarefree-only Mersenne scan omitted 63 and 4095 (MISTAKE-651). THM-4884-definability-lenses... distinguishes actual strong-class multiplicity from numerical grammar-profile multiplicity. MISTAKE-559 forbids treating a small-number coincidence as Chan–Pak coverage. Least-used sidecar here: prime support of a totient preimage, rather than the exponent alone.

Live concepts: Camion residual classes; simultaneous Mersenne/Wagstaff primes; totient fibers; coprimality/support-dependent scaling; strong tournament path spectrum. Anchor is the first two, niche is the exponent-to-totient bridge, wildcard is honest comparison with Chan–Pak.

## A. Prime-index classification of the fiber at 2q

For an odd prime q>3,

  phi^{-1}(2q) = {2q+1, 2(2q+1)} if 2q+1 is prime, and is empty otherwise.

Proof. Any prime p dividing a preimage n has p-1 dividing 2q, so p belongs to {2,3,q+1,2q+1}; q+1 is even greater than 2. Since v_2(phi(n))=1, at most one odd prime divides n, and v_2(n)<=2. If the odd prime is 3, its contribution is 2*3^(a-1); obtaining 2q forces q=3, excluded. With no odd prime, phi(n) is a power of 2. Thus the odd prime must be 2q+1. Its exponent is 1, and the 2-part is either 1 or 2. These two possibilities work. At q=3 the fiber is {7,9,14,18}; at q=2 it is {5,8,10,12}.

## B. Double-prime Camion exponents are nontotients, with sharp exceptions

Suppose q>=3 is odd and a=2^q-1 and b=(2^q+1)/3 are prime. Then q is prime. If q>5, 2q+1 is composite and consequently 2q is a nontotient.

Proof (factor-exhaustion, simplified by the parent). Assume ell=2q+1 prime. Fermat gives ell | 2^(2q)-1=3ab. Since ell>3, primality of a,b forces ell=a or ell=b. For q>5, 2^q>6q+2 (holds at q=6 and propagates by doubling); hence both a,b exceed ell, a contradiction. Apply A.

Optional Euler perspective: if q=3 modulo 4 then ell=7 modulo 8 and Euler's criterion forces ell|a; if q=1 modulo 4 then ell=3 modulo 8 and it forces ell|b. The equality cases give q=3 and q=5, respectively. The simpler Fermat proof already establishes the needed exclusion.

The exceptions are sharp: q=3 gives k=6 and q=5 gives k=10, both totients, while 4^3-1=63 and 4^5-1=1023 pass Camion (THM-4746; q=3 uses its prime-power correction). q=2 yields 15 at exponent 4 outside this odd-q formulation. Exponent 12 gives 4095 by the corrected census, not by this sufficient double-prime family.

Map ledger: source is a certified double-prime pair (a,b) with a=3b-2; map to exponent k=2q and separately to n=2^k-1=3ab. Target predicates are phi^{-1}(k)=empty (q>5) and Camion(n), respectively. The first proof uses the extra Sophie-Germain-prime factor-exhaustion obstruction; the second uses THM-4746's zero residue 3. A bare list of exponents loses primality certification and is not enough for either converse. Sidecar: q plus exact primality certificates for a,b. Cheapest hostile: q=5 (Camion exponent but a totient). The reverse implication is false: q=37 gives the nontotient 74, but M_37=223*616318177 is composite. This refutes recovery of the double-prime certificate, without deciding its Camion test.

This proves a bridge between the two appearances of 14,26,34,38,..., not an identification of their full sequences. The full even nontotient sequence also includes 50,62,68,74,...; the corrected Camion Mersenne exponent list through 64 is 4,6,10,12,14,26,34,38,62.

## C. Support-aware scaling of totient fibers

For positive integers t,n,

  phi(tn) = t phi(n) product_{p|t, p not|n} (1-1/p).

Proof: cancel prime-power factors in the Euler product. Thus n -> tn preserves the scaling phi(n) -> t phi(n) exactly when every prime of t already divides n. More generally partition a fiber phi^{-1}(m) by support intersection S={p:p|t and p|n}; each stratum injects into the fiber at t*m*product_{p|t,p notin S}(1-1/p). A common coprime multiplier injects a fiber into phi^{-1}(phi(t)m); a support-contained multiplier injects it into phi^{-1}(tm). These are injections, not surjections: new target preimages can appear.

Hostile: phi(3)=phi(4)=2, but phi(6)=2 and phi(8)=4, so scaling an entire fiber by 2 does not preserve equality. Nontotient status also does not survive doubling: 14 is a nontotient, but 28=phi(29); 34 and 68 are both nontotients. These are target-value scalings, distinct from multiplying preimages.

For an odd n, phi(2n)=phi(n); for v_2(n)=1, phi(n/2)=phi(n). A singleton totient preimage, if one exists, must therefore have 4|n. This familiar elementary necessary condition neither proves nor disproves Carmichael's totient conjecture: each nonempty fiber has at least two elements. Nontotients have empty fibers and are not counterexamples. Korselt Carmichael numbers instead concern squarefree composite n with p-1|n-1, equivalently lambda(n)|n-1, and are not this conjecture.

## D. Chan–Pak boundary

The legitimate common language is fibers of integer-valued maps, and multiplication of source objects with a sidecar. Ordered sums of tournaments give H(A->B)=H(A)H(B) without a coprimality restriction; phi(ab)=phi(a)phi(b) only when gcd(a,b)=1. Thus there is no automatic fiber-preserving functor. Division k->k/2 sends this certified nontotient subfamily to primes q, but q=7 is a forbidden tournament H value and q=13 is attained (THM-4884 gives a singleton strong class). This is an immediate hostile to identifying gaps after halving. Doubling likewise cannot identify the full nontotient gap set with twice {7,21}. Chan–Pak completeness remains OPEN. Strong-class Carmichael analogues already fail at h=3,5,11,13 (THM-4884), so a singleton analogy is not a transfer theorem.

## Reproduction and limits

Run python3 session218-scaling/verify.py. Exact inverse fibers are independently checked by two complete methods: candidate-support enumeration and a totient sieve through 2m^2. This latter cutoff is rigorous: phi(n)^2/n is a product of p^(a-2)(p-1)^2, each odd-prime factor >=1, the 2-factor >=1/2, hence phi(n)>=sqrt(n/2). No numerical library primality probable-prime claim is used beyond 64-bit inputs in this finite verifier. No conjecture solved, no all-base extension claimed, no infinitude of double-prime indices assumed.
