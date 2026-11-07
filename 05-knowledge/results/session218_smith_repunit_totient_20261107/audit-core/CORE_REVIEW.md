# Frozen checkpoint 1 independent audit

PASS. Reviewed `session218-scaling/PROOF.md` (SHA256 5d2ceb3dfa9db0531439d46eb80b948d3a2baa09a94f8d05b520a6a6603e65b4), `verify.py`, and `verify.out`. Re-ran the producer verifier; output matches byte for byte. Independent standard-library verifier `verify_core.py` passes; receipt `core.out`.

The complete inverse-fiber classification for prime q>3 is correct. The Fermat factor-exhaustion proof is valid and simpler than the optional Euler proof: if ell=2q+1 were prime it would divide 3MW, forcing ell=M or W, contradicted by 2^q>6q+2 for q>5. Primality of M itself forces q prime. Sharp q=3 and q=5 exceptions are included and correct.

The support-aware scaling identity, support-stratum injections, and failure of uniform scaling are correct. The singleton-fiber necessary condition 4|n is valid. Neither zero-sized totient fibers nor Korselt Carmichael numbers settle Carmichael's totient conjecture. The complete finite verifier cutoff n<=2m² is rigorous from the Euler-product bound phi(n)²/n>=1/2. The recursive support enumeration correctly includes powers of2 and n=1.

Independent checks: all prime q from5 through97, with complete inverse fibers of2q by a totient sieve to80000; simultaneous prime pairs for odd q<=31; support-aware scaling for1<=n<=100 and1<=t<=50. Positive and hostile controls include q=3,5 and phi(3)=phi(4) with phi(6)!=phi(8).

Scope: the exact Camion census and tournament statements are inherited cited claims, not freshly re-proved by this audit. The theorem certifies a common subfamily, not equality of the full Camion exponent and even-nontotient sequences. Publication must retain that distinction and the source provenance correction.
