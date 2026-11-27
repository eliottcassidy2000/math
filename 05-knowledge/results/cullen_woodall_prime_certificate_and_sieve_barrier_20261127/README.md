# A prime-certificate machine and a limit on finite sieves for Cullen and Woodall

This follows the [two-sign digit controller](../cullen_woodall_signed_collatz_clocks_20261127/README.md). Write \(N_\pm(n)=n2^n\pm1\). Infinitude of primes in either family is not established here. The distinction below is between an exact per-term certificate, which we can make symmetric in the two signs, and proving that the certificate machine returns new prime hits at unbounded indices.

## 1. Both signs have an exact, one-seed test

Put \(e=n+\nu_2(n)\) and \(k=n/2^{\nu_2(n)}\), so \(N_\pm=k2^e\pm1\) with \(k\) odd. In both cases a prime certificate is a half-exponent equality. For comparison with published searches see [Keller, New Cullen primes](https://cir.nii.ac.jp/crid/1362262944614601984).

**Cullen (the Proth test).** A seed \(a\) certifies primality if \(a^{(N_+-1)/2}\equiv-1\pmod{N_+}\). Any prime factor \(q\) must then have \(2^e\mid q-1\), hence \(q\ge2^e+1>\sqrt{N_+}\). Conversely if \(N_+\) is prime, a multiplicative generator is a seed. Here \(k<2^e\).

**Woodall (norm-one analogue).** Take \(n\ge2\) and a seed \(x\) modulo \(N_-\), put \(D=x^2-1\), and work in the rank-two ring with \(u^2=D\). Form \(z=x+u\), whose norm is one. The test is that the scalar coordinate of \(z^{(N_-+1)/2}\) is \(-1\bmod N_-\). Equivalently use a doubled trace \(T\) initialized by \(T(z)=2x\); multiplication and squaring can be performed by pairs of residues. No factorization of \(k\) is required.

To prove the Woodall criterion, reduce modulo any prime divisor \(q\) of \(N_-\). If \(D=0\bmod q\), then \(x=\pm1\bmod q\) and the exponent \(k2^{e-1}\) is even, so the output scalar would be \(+1\), not \(-1\). Otherwise the quadratic algebra over \(\mathbb F_q\) is either split or a field. A norm-one element with scalar coordinate \(-1\) must be \(-1\) itself, since the square of its \(u\) coefficient times \(D\) is zero. The multiplicative order of \(z\) therefore contains \(2^e\). The norm-one group has order \(q-1\) or \(q+1\) according to the two cases. Thus \(q\ge2^e-1>\sqrt{N_-}\), using \(k\le2^e-2\) for \(n\ge2\). This forces primality.

Conversely if \(N_-\) is prime, take a generator of the norm-one group of its quadratic field. It is not a scalar since \(N_-+1>2\). Write it as \(x+yv\), so \(y\ne0\) and \(x^2-1=y^2v^2\). Changing the quadratic generator to \(u=yv\) gives exactly the seed \(x\) above, and the half exponent returns \(-1\). In particular one may scan just \(x\), not independent discriminants and coefficients.

### Turing-machine interpretation

For each \(n\) and sign, interleave trial factors with trials of the appropriate seed. A factor is a composite certificate; a successful seed is a prime certificate. Every input finishes; Woodall \(n=1\) is the special value 1. To look for the next prime beyond \(B\), scan \(n>B\) and halt on the first prime certificate. For each sign, the infinitude conjecture is precisely totality of that next-hit program on all \(B\), not a consequence of termination on any finite set of \(B\)'s. The replay checks exact certificates for Cullen indices 1, 141 and 4713, and Woodall indices 2, 3, 6, 30, 75, 81, 115 and 123.

## 2. Any finite sieve state can be followed by an arbitrarily long composite block

**Factor-planting theorem.** For any progression \(a\bmod D\), any length \(H\), and any size cutoff, one can construct a larger start \(s\equiv a\bmod D\) so that **both** signs are composite at all \(s,s+D,\ldots,s+(H-1)D\). In fact there is a progression of such starts, and each of the \(2H\) terms can be assigned its own distinct prime divisor, larger than any prescribed bound.

Take \(2H\) distinct primes \(p_{i,\epsilon}\) in one interval \((Y,2Y)\), with \(Y>D\) and beyond the desired bound. The prime number theorem supplies enough in that interval once \(Y\) is large. No selected prime divides another's \(p-1\): a multiple would have to be the prime itself by the interval, but \(p-1\) is even. Set \(L=\operatorname{lcm}(D,\{p_{i,\epsilon}-1\})\). Require \(s\equiv a\bmod L\) and for each chosen prime prescribe
\[
s\equiv-iD-\epsilon\,2^{-(a+iD)}\pmod{p_{i,\epsilon}}.
\]
All the moduli \(L,p_{i,\epsilon}\) are pairwise coprime, so CRT gives a unique simultaneous class. The clock ensures \(2^{s+iD}\equiv2^{a+iD}\) at its assigned prime, proving the divisibility. Choose \(s\) larger than all the chosen primes so that every factor is proper.

For a finite prime sieve, take \(D\) a period of all the sieve tests (for a test at \(p^r\), \(p^r(p-1)\) suffices). Any class the sieve accepts can be used for \(a\), so all members of the planted block still pass the original sieve. This gives a rigorous limitation: a periodic local state cannot guarantee a prime nearby **along its own progression**, no matter how many fixed prime-power tests the state records. For example, in phase \(0\bmod30\) the replay plants distinct primes \(61,67,71,73,79,83\) for three successive steps in that progression; all six terms survive the old 3 and 5 tests.

There is also no fixed-cover proof of eventual compositeness. To avoid a given finite set of primes, take \(n\) divisible by all of them; both signs are automatically \(\pm1\) modulo each odd prime, and both are odd. This yields infinitely many new **prime divisors** in either sequence by iteration, not necessarily infinitely many prime **terms**.

## 3. The proof target left by both constructions

The finite-state side controls phase. The half-exponent tests control primality of a specific term by forcing each hypothetical factor to lie above the square root. What is missing is a construction that forces a successful seed at unbounded indices, or a distribution theorem strong enough to force such hits. Allowing longer modular clocks alone is not that argument: every finite stage can be matched by the composite constructions above. This is the same bounded-witness distinction that matters in the Pollock and Chan–Pak work.

Run `node replay.js` for the witness checks, the total per-term machine on small inputs, and three modular factor-planting replays. It verifies planted factors by modular exponentiation without expanding the enormous terms.
