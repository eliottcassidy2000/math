# Collatz-style prime certificates and a complete inverse for the odd-run exit

This follows the [Cullen/Woodall certificate machines](../cullen_woodall_prime_certificate_and_sieve_barrier_20261127/README.md) and the [shared midpoint calculation](../cullen_woodall_signed_collatz_clocks_20261127/README.md). The useful link is stronger than an analogy: odd-index Cullen and Woodall terms can realize **every finite legal parity word after a prescribed exit**, even when both terms are forced composite. This also marks the limit of parity-only arguments.

## 1. The certificates are countdown machines

Write \(n=k2^{v_2(n)}\), \(e=n+v_2(n)\), \(N_\pm=k2^e\pm1\). The integer \(e\) is exactly the length of the initial synchronized signed odd run. It also controls the depth of the fast certificate:

* Cullen starts at \(R=a^k\pmod{N_+}\) and repeats \(R\gets R^2\) exactly \(e-1\) times. It accepts a prime certificate at \(R=-1\).
* Woodall starts with \(R\) equal to the trace of \(z^k\pmod{N_-}\), for the norm-one seed \(z=x+\sqrt{x^2-1}\) of the preceding proof. It repeats \(R\gets R^2-2\) exactly \(e-1\) times, accepting at \(R=-2\). The identity is \(\operatorname{Tr}(z^2)=\operatorname{Tr}(z)^2-2\) for norm one.

This is one countdown register with two different update lanes. The Collatz odd run consumes the same valuation budget one factor of two at a time, while multiplying the center by \(3/2\). The certificate updates square *residues* and never choose the next step by an integer parity test. The endpoint implications are different: the half-power witness excludes small prime factors, whereas Collatz needs an actual return below a previous integer value. Neither state machine implicitly proves that return.

## 2. Prescribe the whole exit and any finite continuation modulo powers of two

Use signed halved maps \(T_\pm(x)=x/2\) on even \(x\), and \((3x\pm1)/2\) on odd \(x\). For odd positive \(n\), start at \(W_n=n2^n-1\) in \(T_+\) or \(C_n=n2^n+1\) in \(T_-\). After exactly \(n\) consecutive odd steps, the first even values are
\[
n3^n-1,\qquad n3^n+1.
\]
For either sign \(\epsilon\in\{-1,1\}\) on the center, let \(s\ge1\) be the desired number of following halvings. Define \(X=(n3^n+\epsilon)/2^s\); forcing \(X\) odd ensures the number of halvings is exactly \(s\). The subsequent map is \(T_{-\epsilon}\).

**Post-exit inverse.** For each \(s\ge1\), \(k\ge1\), and odd residue \(r\bmod2^k\), there is exactly one odd class \(n\bmod 2^{s+k}\) such that the selected exit has \(s\) halvings and \(X\equiv r\bmod2^k\). The test is
\[
n3^n\equiv 2^s r-\epsilon\pmod {2^{s+k}}.
\]
Indeed \(n\mapsto n3^n\) permutes the odd residues at every precision. From a lifted odd \(n\) at precision \(2^j\), \(j\ge1\), the candidates \(n,n+2^j\) produce outputs differing by \(2^j\bmod2^{j+1}\), since \(3^{2^j}\equiv1\bmod2^{j+1}\). Thus exactly one next bit works, beginning from the unique odd class mod 2.

For either signed map a parity word of length \(k\) determines one starting residue mod \(2^k\). To see this, once the first \(j\) inputs have the same parities, translating the starting value by \(2^j\) translates the \(j\)-th iterate by \(3^{b}\), where \(b\) is the number of odd inputs in that prefix. This flips the next input parity. So any word beginning in 1 is realized by an odd residue, and composing this inverse with the exit inverse gives an explicit index for that word. A quantitative caution: the index can be chosen below \(2^{s+k}\), but the actual exit grows like \(n3^n/2^s\), not like the least residue. For example on the Woodall side with \(s=1\), all such exits are \(1\bmod3\); covering every odd 2-adic class is not covering all positive odd numbers.

## 3. The same words occur on provably composite inputs

Fix the index class from the inverse. Keep its 2-power modulus and also fix the exponent mod 12 by choosing any compatible common phase. The new moduli 5 and 7 are coprime to this combined clock, and the orders of 2 at those primes divide 12. The two-sign CRT rule therefore allows \(5\mid C_n\) and \(7\mid W_n\) simultaneously, without changing the selected exit or its word. The indices are unbounded, so these factors can be made proper. Replacing 5 and 7 by new primes as in the [factor-planting theorem](../cullen_woodall_prime_certificate_and_sieve_barrier_20261127/README.md) gives analogues beyond any fixed factor cutoff if desired.

For instance, choose an odd Woodall exit after one halving, with six following input parities \(111010\). The inverse gives \(n=53\bmod128\); \(n=2997\) realizes the word and has the planted factors on both sides. Taking \(k\) arbitrarily large and choosing the all-odd word shows that one cannot guarantee a drop in a fixed number of steps even **after** removing the initial forced odd run. This is an obstruction to uniform finite-horizon arguments, not to eventual descent.

## 4. What would close Collatz

The universal positive-integer conjecture is equivalent to eventual strict descent below the starting value for each \(x>1\): descent then iterates by well-ordering. A usable quantitative proof target follows directly from the exact multiplier. If \(N>1\) never goes below itself, and \(o_j\) odd steps occur among its first \(j\) steps, then
\[
\frac{T_+^j(N)}{N}=2^{-j}\!\prod_{\substack{i<j\\T_+^i(N)\ {\rm odd}}}\left(3+\frac{1}{T_+^i(N)}\right)
\le 2^{-j}(3+1/N)^{o_j}.
\]
Thus it would suffice to force a prefix with \(o_j/j<\log 2/\log(3+1/N)\) for every \(N>1\) (or prove another bound on the same product). The exact constructions above explain why merely listing admissible short words, or deciding which special starting numbers are prime, supplies no such control. This is consistent with the gap between finite local dynamics and global size estimates in [Tao's almost-bounded theorem](https://arxiv.org/abs/1909.03562); there is no claim of a proof of Collatz or of a reduction from Cullen/Woodall infinitude.

The JavaScript replay checks 8184 inverse cases, an actual composite pair with the prescribed exit, and five prime countdown certificates.
