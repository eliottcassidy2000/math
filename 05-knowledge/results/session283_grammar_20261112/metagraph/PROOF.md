# Separator heat flow: an exact finite response grammar

**PROVED, independently reviewed in this packet.** This is the expectation law for
uniform labelled arc-reversal dynamics, not a closed Markov process on the
compressed separator/rank statistics or an A000568 enumeration formula.

## Inheritance and coordinates

The homogeneous-pair path forest and its maximal transitive-module quotient are
inherited from `weighted_core_transport_20261020/README.md` and its cited core
classification. Write b(T) for its number of components. The sibling
[local-rewrite proof](../grammar/PROOF.md) supplies its weighted-core separator
polynomial and exact one-step rank kernel. Generic hypercube/Krawtchouk diagonalization
is classical; repo inspiration includes `07-reflections/krawtchouk-coordinates-of-tournament-space.md`
and THM-3417's exchangeable shell inversion. The old reflection is inspiration only: its identification of H with a linear Hamming coordinate is not used (and is not a valid identity). Our cube is the full labelled arc cube, not a fixed-path tiling quotient. We apply it to this specific faithful
response, not claim a new family of orthogonal polynomials.

For each unordered pair uv, D_T(uv) is the number of exterior vertices w that
see u and v differently. Set N=n-2, M=binom(n,2), and

    P_T(z)=sum_{uv} z^{D_T(uv)}=sum_d e_d z^d.

A pair is homogeneous iff D=0. Its graph is the inherited path forest, so
b(T)=n-e_0. This includes n=2; we state dynamics below for n>=3.

## 1. Exact invariant expectation space

Let T' be obtained by reversing one of the M arcs, uniformly. Then

    E[P_T'(z)|T] = P_T(z)
        + (2/M)[N(z-1)P_T(z)+(1-z^2)(d/dz)P_T(z)].       (1)

For a fixed uv, reversing uv changes no separator bit.
For each exterior w, exactly two arc flips toggle its bit, uw and vw. All
other flips leave the bit vector unchanged. A pair at separation d therefore
has transitions to d+1 with probability 2(N-d)/M and to d-1 with 2d/M.
Summing over the actual pairs proves (1), without independence assumptions
between different pairs. Linearity then proves that the expected polynomial
at every time t is L^t P_T, where L is the explicitly displayed linear map.

Thus the n-1 coefficients are a closed expectation grammar. They need not
form a Markov state for the actual random histogram. In particular, the
sibling same-histogram hostile pair has different rank transition laws.

## 2. Exact all-time spectral rank decoder

Define K_s(d)=[u^s](1-u)^d(1+u)^(N-d), and Q_s=sum_d e_d K_s(d).
The separator bit walk is the lazy coordinate-flip walk on {0,1}^N with
probability 2/M for each coordinate. A size-s character has eigenvalue
lambda_s=1-4s/M. Fourier inversion at the all-zero vector gives

    E[b(T_t)|T_0=T]
       = n - 2^(-N) sum_{s=0}^N Q_s lambda_s^t.      (2)

No full tournament enumeration is required once P_T is known. All eigenvalues
except lambda_0 have absolute value<1 for n>=3, and Q_0=M. Hence

    lim_t E[b(T_t)] = n-M/2^(n-2).

This stationary law is for uniform LABELLED tournaments, not uniform
isomorphism classes. A000568 requires stabilizer weights and cannot be read
from this limit.

## 3. Poissonization and a universal macroscopic fracture law

Let every arc flip independently at rate1. Put

    a(s)=(1+exp(-4s))/2,  c(s)=(1-exp(-4s))/2.

For each fixed pair, its N separator bits flip independently at rate2,
because their two-edge sets are disjoint. Therefore

    E[b(T(s))]=n-a(s)^N P_T(c(s)/a(s)).              (3)

Equation(3) also follows by Poissonizing (2); it is not an assumption of
independence between pairs.

Consider any sequence of weighted transitive blowups of pair-free cores with
at most B vertices, for fixed B, and total size n→infinity. The sibling
separator formula yields e_0=n-b and e_1<=n+binom(B,2)B. Since sum e_d=M,
for s=tau/n with fixed tau>=0, q=c/a=O(1/n) and

    sum_{d>=1} e_d q^d <= e_1 q+M q^2/(1-q)=O(1).

Also a(s)^N→exp(-2tau). Consequently

    E[b(T(tau/n))]/n → 1-exp(-2tau).               (4)

This is an exact expectation scaling theorem: their expected core rank becomes
linear after order-n expected random arc edits. No concentration or typical-trajectory claim is made. It does not say their
whole structure becomes independent or random, and it is not a Chan–Pak
coverage theorem. The mean number of flips by time tau/n is asymptotic to
tau*n/2. No de-Poissonized distribution claim is required.

## 4. Connection and loss ledger

Weighted core -> separator polynomial -> heat response -> expected rank.
The first map retains enough for the specified mean dynamics, the second is
an exact linear representation, and the third evaluates a coefficient. Each
forgets orientation correlations, actual path counts and stabilizers.
The full local rank kernel in the sibling lane retains more one-step data
but still does not select an arbitrary Hamiltonian count. This is a model of
manageable infinite evolution through a closed observable space, not a
substitute for a complete isomorphism grammar.
