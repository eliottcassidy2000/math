# Weighted-to-binary contextual simulation at 3, 5, 7, 11

**Status: PROVED using the independently audited finite prime-state coverage; exact positive recipes supplied for the same-session five-vertex cospectral pair.** This preserves modular context responses, not real stability, graph order, or zeta.

## 1. The general simulation theorem

Let p be one of3,5,7,11. Let A be any nonempty rational weighted tournament: Aii=0, Aij+Aji=1, with every entry's denominator coprime to p. Nonnegative real weights are allowed but positivity is not needed for the algebraic statement. Weighted path covers use products of the edge weights inside each path; define p_A(c) and F_A as in the binary case.

> **Theorem.** There is an actual binary tournament T of substitution width at most3 such that, for **every** binary homogeneous one-port context Q,
>
>     H(Q[T]) = H(Q[A]) mod p.                              (W1)
>
> Reduction of the rational right side is canonical because all denominators are p-units. Each T can be supplied by a finite positive constructor expression.

**Proof.** Expand A as a formal multilinear mixture of binary orientations. For each binary tournament R on the same vertex set give it coefficient

    w_R=product_(i<j) [Aij if i->j in R, else Aji].

The coefficients sum to1 because Aij+Aji=1. For any directed path cover its distinct used edges have expectation equal to their product, while unused edges sum out. Therefore

    F_A(x)=sum_R w_R F_R(x).                              (W2)

This identity is polynomial and remains valid over the localization at p. Since every R has the same vertex count n, parity gives

    F_A(0)=0,       F_A(-x)=(-1)^n F_A(x).

After reduction, F_A is in one of the two parity function spaces from `PROOF.md`, Theorem2. Its finite-exact all-state coverage at p=3,5,7,11 supplies an actual width-three T with the same function on F_p. At the zero function a parity tag is irrelevant; source padding can match it if desired.

The homogeneous substitution kernel extends to weighted input profiles without change: select the maximal input runs, multiply their internal weights, and count the exterior quotient words. Its coefficients retain the exact divisors d_c. For c<p these are units modulo p; for c>=p they vanish. Thus the same triangular evaluation proof as Theorem2 shows that equality of F-functions is exactly equality of every weighted/binary homogeneous-context response. This proves(W1). QED.

This is an actual modular simulation theorem, not a statement that the weighted model itself becomes an ordered/probabilistic object over F_p. At other primes or general composite moduli, arbitrary weighted-to-binary contextual simulation is not proved here.

## 2. The five-vertex cospectral pair becomes two explicit binary modules

Use A,B from the same-session `../riemann_five_eleven/WEIGHTED_ZETA_BOUNDARY.md`. Their entries have denominators supported at2 and3, so reduction at5 and11 is permitted. Their determinant zetas are identical over the rationals, but H(B)-H(A)=delta=1/3456.

A literal weighted permutation-and-cut enumeration, independently of that determinant argument, gives

    p_A=(15/2,30,30,10,1),
    p_B=(25921/3456,30,30,10,1).

Equivalently,

    F_A(x)=x^5+5x^3+(3/2)x,
    F_B(x)=x^5+5x^3+(5185/3456)x = F_A(x)+delta*x.        (W3)

### Modulo5: the H5/H11 join and a singleton

The first function is zero on F5, while the second is x. Hence

    T_A=W5->W11       (order9, H55),
    T_B=K1           (order1, H1)

are actual width-three simulations of A and B modulo5 in **every** homogeneous context. The equality of the original weighted zetas does not transport to these binary representatives, and is not claimed to do so.

### Modulo11: explicit positive words

Words are chronological from K1, with P,C,D exactly as in `PROOF.md`.

| input | positive word | order | actual integer H | F(1),...,F(5) mod11 |
|---|---|---:|---:|---|
|A|`PPDDCDDCDPCDP`|29|100671287595046269|2,9,3,8,1|
|B|`PPCDCCPCDDDD`|27|13291517610922779|8,10,10,10,9|

Full integer polynomial replay checks these functions at every residue of F11, not only their scalar H. Their difference is6x as a function, because delta=6 mod11. The large integer H values illustrate exactly why modular simulation does not imply exact integer realization at the small weighted input values.

## 3. Input/context duality with the Session138 one-edge reflector

The pair's difference is particularly rigid: (W3) says only p(1) changes. The one-run coefficient of any context Q is H(Q), because collapsing its input to a singleton gives precisely Q. Consequently, over the rationals,

    H(Q[B])-H(Q[A])=delta*H(Q) for every Q.                (W4)

This is stronger than separation by the identity observer, but it does not assert that every modular observer separates the pair: observers with H(Q)=0 modulo p can be blind.

In the inherited Session138 reflector construction, the two actual contexts R_plus,R_minus have coefficient-row difference `(2,0,0,...)`. The weighted cospectral pair has profile-column difference `(delta,0,0,...)`. They are dual first-channel isolators in the exact substitution pairing

    H(Q[T])=sum_c a_c(Q) p_T(c).

Their mixed rectangle is therefore exactly

    [H(R_plus[B])-H(R_minus[B])]
     -[H(R_plus[A])-H(R_minus[A])] = 2*delta.              (W5)

At p=5 and11, the actual binary simulator modules give the same mixed identity modulo p. This is an explicit map between old-edge context cancellation and the weighted five-vertex zeta-loss sidecar; it is not a syntax-only analogy.

### The spectral loss survives EVERY homogeneous context

In the weighted witness B=U A U^T, where U is rational orthogonal and U1=1. Orthogonality also gives U^T1=1, hence both row- and column-ones identities. Let Q be any binary homogeneous one-port context and order the substituted module first. The cross blocks between the module and every exterior vertex are either the zero vector or the all-ones vector. Therefore

    Q[B] = diag(U,I) Q[A] diag(U^T,I).                       (W6)

The conjugating matrix fixes the global all-ones vector. Thus, for every Q, the two weighted substitutions have the same full adjacency and skew characteristic polynomials, determinant zeta, and all-one marked resolvent (equivalently all moments `1^T Q[A]^k 1`). Their Hamiltonian-path sums nevertheless differ by the positive rational amount `delta H(Q)` in(W4). Since the determinants u and v also agree, their real stability deficits Delta differ by the same amount.

This is **contextual cospectrality**: even allowing all homogeneous binary exterior contexts does not turn these spectral observations into an exact Hamiltonian-path or stability-deficit reconstruction. The loss can be amplified by any attained H(Q), not merely preserved under ordinal singleton padding. In the reflector rectangle(W5), the exact rational value is `2 delta=1/1728`.

The claim remains about failure of exact reconstruction, not failure of spectral sufficient lower bounds or impossibility of a proof using spectral data together with a non-spectral sidecar. The minimum-order-five weighted blindness boundary is unchanged. The actual binary simulators in Section2 preserve the context-H data modulo a good prime; they are NOT asserted to inherit the weighted contextual cospectrality.

## 4. Verification and scope

`verify_weighted_simulation.py` is standard-library-only. It enumerates all120 permutations and16 cut patterns for each of the two rational5-vertex matrices, divides by c! to obtain unordered cover profiles, checks(W3), and checks the modulo5 simulations. It then finds the explicit modulo11 words in the actual P/C/D quotient and independently replays their **full integer polynomials**. It additionally checks the exact block orthogonal conjugacy and first-channel response law for all75 labelled port-zero contexts of orders1..4 plus both five-vertex reflector contexts, using independent rational Held–Karp path counts on the resulting graphs of order at most9. The mixed rectangle is exactly1/1728. Normal and optimized outputs agree.

    python3 -B verify_weighted_simulation.py > weighted.normal.json
    python3 -B -O verify_weighted_simulation.py > weighted.optimized.json
    cmp weighted.normal.json weighted.optimized.json

Connection contract: rational weighted A -> its parity F-function modulo a good prime -> an actual positive binary expression. Preserved: every homogeneous-context Hamiltonian response modulo that prime. Lost: real weights, order, exact H, determinant/zeta, V/tau, Delta, absolute values, and inequalities. Therefore this theorem does not transfer the unresolved real stability range, prove RH, or solve the integer Hamiltonian-path spectrum. The finite prime-state realization hypothesis is explicitly limited to3,5,7,11.
