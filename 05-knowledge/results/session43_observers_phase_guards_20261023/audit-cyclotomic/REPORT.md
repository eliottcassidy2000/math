# Independent audit of Session43 cyclotomic path-cover bounds

Status: PASS. No mathematical correction required. The audited producer packet is `session43-cyclotomic/`, frozen manifest SHA256 `8cc181b333bf0ebad59efd20c65493932c5446ff54fd35a9dec5fe2619fdb43c`. No producer files or Git state were modified.

## Proof audit

1. THM-4756 explicitly contains both the factorial path-cover expansion and the all-k ordered-block principal specialization. Its retained cycle convention is directed cycles modulo rotation, not reversal. The producer's falling-factorial conversion to the exact odd-cycle Eulerian expansion is valid, including d>=1 and n−d even.
2. Legendre's formula gives termwise pi-valuation n−s_p(m). Taking the maximum digit sum in1..n gives a universal lower bound; cancellation can only increase divisibility.
3. For odd p, multiplying by the correct root-of-unity phase makes the evaluated polynomial real. Conjugation fixes the residue field and negates the uniformizer to leading order, so every nonzero valuation is even. This proves the rounded bound. The zero case is separately handled.
4. At p=3 the real-phase coefficient is a rational integer because the other two residue-class coefficient sums coincide. This rationality does not extend automatically to p>=5, and the producer correctly says so.
5. For odd n, the empty cycle collection has valuation n−v_2(n+1). Every nonempty collection has valuation n+r−v_2(d+1), with r>=1 and d+1<=n−1. The bound B(n) follows exactly.
6. In the nontransitive extremizer, all directed cycles contain the same unique reversed arc. Hence no collection contains two disjoint cycles. Its unique longest q-cycle supplies a unique term of valuation B(n); every shorter cycle gives d+1 strictly between2^a and2^(a+1), while the empty term is higher by the case choice. There is no same-valuation term available to cancel it. This proves sharpness for every odd order, not merely the tested orders.
7. Transitive attainment occurs at odd n with n+1 a power of2 or three times a power of2, with the stated small-order convention. Away from powers of2, the condition is v_2(n+1)=a−1 in the open dyadic interval; its only solution is n+1=3·2^(a−1). This is a binary valuation classification, not a Collatz theorem.
8. The Gaussian radial expression is a finite factorial-moment transform. It does not assert vanishing of all powers of one fixed polynomial. No Gaussian moment conjecture or Jacobian consequence is established, as correctly stated by the producer.

## Classical tangent input checked against primary references

NIST DLMF24.15.E4 gives the tangent/Bernoulli formula, and24.10.E1 gives the von Staudt–Clausen denominator statement:

- https://dlmf.nist.gov/24.15#E4
- https://dlmf.nist.gov/24.10#E1

Together they give v_2(T_(2j−1))=2j−1−v_2(2j). The identification A_(2j−1)(−1)=(-1)^(j−1)T_(2j−1) can be checked without an extra citation: the Stirling expansion yields the exponential generating function

Σ_(d>=0) A_d(x)t^d/d! = (x−1)/(x−exp((x−1)t)).

At x=−1 this is2/(1+exp(−2t))=1+tanh(t). Comparing tanh with tan gives the required sign. Thus the proof's complete tangent-valuation chain is supported.

## Independent computational architecture

The auditor does NOT import or invoke the producer checker.

- Computes F_T by direct enumeration of all vertex permutations, rather than subset dynamic programming.
- Enumerates actual directed cycles and disjoint collections directly to reconstruct the Eulerian expansion.
- Evaluates the cyclotomic norm as the determinant of the exact integral multiplication matrix, using fraction-free Bareiss elimination. Because p is totally ramified with residue degree1 in the pth cyclotomic field, the rational p-valuation of this norm equals the pi-valuation. This is independent of the producer's repeated division algorithm.
- Builds tangent numbers from the differential equation tan'=1+tan² and binomial convolution, rather than the Eulerian recurrence.

Results:

-1108 tournaments counted by direct permutations, including every labelled tournament of orders1..5 and larger deterministic witnesses through9.
-1103 exact cycle-expansion comparisons.
-4441 independent cyclotomic norm checks at primes2,3,5,7 and additional larger-witness checks at11.
-200 odd orders3..401 verify the explicit sharp-family scalar and transitive-attainment classification using tangent ODE data.
-Actual extremizer tournament witnesses through9 attain B(n), not merely its lower bound.
-Normal and optimized outputs agree byte-for-byte. No assertions or floating literals are used.

## Provenance and remaining limits

All8 files in the frozen producer manifest match their checksums. All11 source pins match repository revision8b71b746f7a9d90378c702e01608207227b1f83e. `source_audit.json` records this check.

This audit does not claim global literature novelty, optimal odd-prime valuations, a universal odd part of the gcd, a Hamiltonian-count spectrum theorem, native port realizability, or equivalence with the8/9 observer coupling. The all-size new conclusion audited here is the exact universal binary valuation and its explicit attaining family on the forward-step polynomial carrier.

## Publication portability replay

The checkers also run from an unrelated working directory. The affine audit resolves its independent basic-audit sibling relative to its own file, accepting both workspace and published directory names. The basic audit has no runtime packet dependency. The producer affine checker was repinned after the same layout-only repair; no mathematics or numerical results changed. `portability.json` records normal and optimized executions of both producer packets and both independent audit packets in a mock published directory, all launched from `/tmp`, plus workspace-layout independent-audit replays.
