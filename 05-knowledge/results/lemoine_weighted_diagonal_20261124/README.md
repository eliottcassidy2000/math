# Lemoine: the missing diagonal and the two-anchor obstruction

2026-11-24. Continuation of HYP-2218/S642 and HYP-2219/S643.
Three separate outputs: an unconditional obstruction to short anchor banks
even on **represented** inputs, an unconditional obstruction to a
one-Fourier-sum proof transfer, and a square-root-scale exceptional set
obtained from the same Dirichlet-GRH input used by the linked development.
None asserts Lemoine for every odd integer.

## 1. What Glazer's proof transfers, and what it does not

Read the pinned [README](https://github.com/elliotglazer/ternary-goldbach/blob/d16af8a0ef27610cbd25c29189d19e463c43bdc8/README.md),
[NumericAssembly](https://github.com/elliotglazer/ternary-goldbach/blob/d16af8a0ef27610cbd25c29189d19e463c43bdc8/proof/NumericAssembly.lean),
[ProfileRoute](https://github.com/elliotglazer/ternary-goldbach/blob/d16af8a0ef27610cbd25c29189d19e463c43bdc8/proof/ProfileRoute.lean),
[CompletedFiniteRows](https://github.com/elliotglazer/ternary-goldbach/blob/d16af8a0ef27610cbd25c29189d19e463c43bdc8/proof/CompletedFiniteRows.lean)
and [CompletedGoldbach](https://github.com/elliotglazer/ternary-goldbach/blob/d16af8a0ef27610cbd25c29189d19e463c43bdc8/proof/CompletedGoldbach.lean).
The published route isolates a finite binary prefix, a prime-gap bridge
through 2·10^19, 223 analytic intervals, and a uniform tail. The final
declaration substitutes the imported RH/Dirichlet-GRH theorems into the
conditional assembly. This session uses the mathematical interface, not a
fresh replay of those Lean objects.

S642/S643 already give E=p+q, N=p+2q, and swap N↔3E−N. Think of a
fixed-N slice as a **different** symmetry problem: its two endpoints satisfy

    3p − N = −2(3q − N).

This is a dilation, not the fixed-E reflection. Iterating inside the
prime set produces paths with doubled distance from N/3; the sole possible
fixed point is the prime diagonal p=q=N/3. The number 21 has three
representations (p,q)=(17,2),(11,5),(7,7); the last is the old E=14
diagonal, the second survives after removing both the diagonal and the
prime-2 boundary. In contrast 7 has only (3,2). Neither is a Lemoine hole.

For a finite prime weight w and exponential sum S, the ternary integral
with two auxiliary frequencies is the coefficient of
S(α) S(α+β) S(α−β). Its β=0 **slice** sums all ternary representations;
its β-average selects equality of the last two primes. With log weights
that average selects log² q on the repeated prime, instead of log q.
That difference is irrelevant for positivity, but slice and average
are not interchangeable. Helfgott's third prime provides an L∞
minor-arc saving alongside two L² factors; Lemoine has just two
factors. Glazer's finite gap step likewise leaves two free primes
after choosing a nearby prime. Choosing one prime in Lemoine leaves
a single prescribed integer, which a prime-gap bound cannot make prime.

## 2. Unconditional: represented targets with neither small endpoint

Write log_j for j-fold iterated natural logarithm. For every sufficiently
large X there is an odd N≤X **with** a Lemoine representation, but every
representation p+2q=N has

    min(p,q) ≫ log X · log_2 X · log_4 X / log_3 X.

This is a transfer of the covering form of the
[Ford–Green–Konyagin–Maynard–Tao prime-gap theorem](https://arxiv.org/pdf/1412.5029)
and [Linnik's least-prime theorem](https://arxiv.org/abs/0906.2749). It does not use a conjectural Lemoine
witness. A consecutive gap by itself would obstruct small q but say
nothing about whether the target is represented, or about small p.

**Proof.** At scale y the covering theorem covers [1,U] by one class
a_r mod r for every prime r≤y, with U≫y log y log_3 y/log_2 y. Choose
b∈{−1,0} so 2m+b avoids a_2 mod 2. The classes
c_r=(a_r−b)/2 mod r at **odd primes only** now cover every
m∈[1,L], L=floor(U/2). Use a truncated U of the asserted scale.
Impose N≡c_r mod r and N odd. If N>2L+2y and 2q≤L, then
N−2q>y and has a divisor r≤y. If p≤L, then (N−p)/2>y and has
one too (since r is odd). The entire covered block does both jobs.

To ensure N is represented, for each r pick a residue b_r for q avoiding
both 0 and c_r/2. This is possible even at r=3. Choose a fresh prime
s>2L+2y+2 of size O(L+y), force q≡s−1 (mod s), and include odd parity.
The progression is reduced. By Linnik there is a prime Q in it with
Q≤(2s∏r)^O(1)=exp(O(y)), and Q≥s−1. For p take the reduced progression
p≡c_r−2Q mod r with odd parity. Linnik supplies a prime P≤exp(O(y)).
Now N=P+2Q is odd, satisfies every prescribed residue, and is above
2L+2y. Take y a sufficiently small fixed multiple of log X.
All constants are independent of X. This rules out uniform O(log X)
anchor cutoffs even when either prime may be used as the anchor.

**GRH refinement.** For any fixed ε>0 the same cover produces
≫ε X^(1−ε)/log X represented hard targets up to X, with a constant
depending on ε in the lower bound on min(p,q). Choose y small enough
that the CRT modulus M=2∏r is ≤X^min(ε,1/4) and Q≤X^(1/4). Keep Q
fixed and put P in [X/2,3X/4] in its prescribed reduced class.
Dirichlet GRH gives weighted count ≫X/φ(M) with error
O(√X log²(XM)), sufficient here. The targets form a sparse but
certifiably represented family, not a lower bound on exceptions.

## 3. Unconditional: remove one target without appreciably moving S

For any odd N≤2X, let S(α)=∑_{p≤X} log p·e(αp). There exist
weights w_p∈{0,2} on the **same primes** such that the coefficient
of e(Nα) in S_w(α) S_w(2α) vanishes, while

    sup_α |S_w(α) − S(α)| ≪ √X (log X)^(3/2).

This also holds simultaneously over every interval of prime indices.
If 3∤N the error improves to O(√X log X). This is a local
countermodel at one selected N, not a claim about multiplicative
Euler products or GRH for the altered weights.

**Proof.** Let the graph on primes have an edge q→N−2q when the
endpoint is also in the list. Remove the possible loop at N/3
(weight zero). The graph has indegree and outdegree at most one
and has no nonloop directed cycle by the dilation identity. Thus
its remaining components are simple paths. Independently choose
one of each path's two color classes to keep and double. A fixed
frequency/interval discrepancy is a sum of independent signed
complex numbers with total squared magnitude at most

    (maximum component order) · ∑_{p≤X}(log p)^2 ≪ X log² X.

The maximum is O(log X) by the doubling identity and the bounds
1≤|3p−N|≤3X. If 3∤N, residues along any chain advance by N mod 3;
six vertices would force two occurrences of the prime 3. So the
maximum is then five. Hoeffding on the real and imaginary parts,
a union bound over O(X²) intervals and a mesh of O(X³) frequencies,
and the O(X²) derivative bound yield the assertions; the removed
loop costs O(log X). Character sums on progressions are not
assumed: AP marginals follow directly by additive orthogonality,
and fixed smooth weights by partial summation.

The same perturbation changes **every** ternary Fourier coefficient
by O(X^(3/2) log^(5/2) X), by telescoping S_w³−S³, taking the
difference in sup norm and bounding the other two factors in L².
Thus any ternary coefficient with a positive main term of order X²
stays positive. On arcs of total length X^(-1+o(1)), it changes the
mixed major coefficient by at most X^(1/2+o(1)); this follows by
expanding the difference and bounding each unperturbed sum by O(X).
In any bulk range where the true binary coefficient is ≫X, the
mixed minor **signed coefficient** consequently moves by order X.
Such coefficients exist, e.g. pigeonhole pairs with both primes
in [X/3,X/2] using the prime number theorem.
This isolates the information needed: a bound on S alone at
square-root accuracy, even uniform in frequency, does not replace
control of the correlation of S(α) and S(2α) at the prescribed phase.
It does not preclude using the exact multiplicative structure of Λ.

## 4. Analytic positive transfer and the exact open estimate

[Details and proof](GRH-EXCEPTIONAL.md): with the RH/Dirichlet-GRH
interface read in Glazer, for every ε>0 at most Oε(X^(1/2+ε))
odd targets up to X lack a representation, even when both primes
are required to lie in fixed bulk windows at each dyadic scale.
The Lemoine singular series is **exactly the ordinary Goldbach
singular series at 2N**; the geometric factor differs. This is a
weighted adaptation of classical GRH binary methods, not a claim
to priority for the exponent.

The next proof test should be a one-sided bound on the mixed
minor coefficient Re∫_minor S_f(α) S_g(2α)e(−Nα)dα, uniformly in N,
strong enough that it cannot cancel the positive major lower bound.
The square-root exceptional-set result is a second-moment statement:
even a vastly better finite certificate does not close that gap.
One concrete structural question is whether multiplicative identities
for the *same* prime weight at α and 2α reject the independent
chain-color passes in §3 without discarding the signed coefficient.
The Chan–Pak analogy is precise in scope: two projections can
preserve their coarse readouts while simultaneously suppressing
one exact target. The missing object is a compatible **joint** witness,
not another check of the marginal readouts.

## 5. Finite replay and the same-pair bridge

    g++ -O3 -std=c++17 replay.cpp -o replay
    ./replay 100000000

Deterministic Eratosthenes and ascending scans, no floating point.
Among 49,999,997 odd targets 7 through 99,999,999 there are zero
misses. The odd-prime bridge fails only at 7; deleting both q=2
and the diagonal leaves 7, 9, 15. Maximum first doubler q is 1231
(first at N=45,521,269, p=45,518,807); maximum first undoubled p
is 2833 (first at N=87,170,987, q=43,584,077). The largest
min(p,q) required, allowing either endpoint to be small, is 829,
first at 74,475,571. The replay is a reproducible same-pair
census, not a record claim for ordinary Lemoine verification.
No tournament fingerprint is reported: the exact readout here uses
a graph with in/out-degree ≤1, and filling missing orientations
would dominate the fingerprint with arbitrary choices.
