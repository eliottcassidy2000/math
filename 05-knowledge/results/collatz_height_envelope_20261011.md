# An exact height envelope and a smaller source-specific ROOT bill

2026-10-11. **PROVED, scoped and self-audited:** C1–C3 below, without
assuming universal convergence. **FINITE-EXACT:** the attached replay.
**OPEN:** an independent lower measurement crossing the new threshold,
universal ROOT entry, and universal Collatz.

## Inheritance and the target

The closest mechanism is A3 of [native alias contexts](context_prime_root_20261011.md):
the source-dependent modulus makes the first alias an actual sibling, reducing
the unweighted r=6 error to 1/(16q^7). The beta-mixture W then had only its
global cap 1/3 applied, giving 1/(48q^7). We improve both the cap and the
sibling multiplier by retaining the *size of the actual source*.

The score is exactly the inherited [adaptive mixture](collatz_adaptive_mixture_flow_20261005.md),
also defined in [source refinement floors](collatz_source_refinement_floor_20261005.md).
The latter's optimal section inequality with constant160/441 is a different,
stronger inherited lower-bound transfer; it is not reproved here. The hostile
is a positive residue class with a deleted target atom. The corrected near
miss would be to mistake a decreasing error allowance for a positive lower
readout. The least-used sidecar is the terminal valuation's parity.

Portfolio: ROOT error control / geometric palette budgets / cubic trace and
elliptic differential. Board: **source height, terminal parity, beta counters,
native siblings, independent lower readout, boundary extension counts**.

For an odd n>1 whose accelerated orbit first reaches1 after valuations
a_1,...,a_tau, set

    A=sum a_i,  K=sum floor((a_i−1)/2),  H=K+tau,
    W(n)=2 K! tau!/(K+tau+1)!.

Set W(1)=1 and W(n)=0 for every unrooted source. In a strict positive
first-hit word, K>=1 and tau>=1. Every intermediate source n_i is at least3.
Let e be the number of even a_i. The terminal valuation is even and at
least4, so

    A=2K+tau+e,                 1<=e<=tau.                    (1)

These assertions are conditional calculations on a supplied finite word;
the resulting upper bounds also hold at unrooted sources because W=0.

## C1. The tail maximum is explicit and attained

For an integer N>=2 define

    h(N)=floor(log_4(3N))+1,      B(N)=2/[h(N)(h(N)+1)].      (2)

Then, unconditionally,

    max {W(n): n positive odd, n>=N}=B(N).                  (3)

In particular this is an exact envelope, not merely an asymptotic estimate.

**Proof.** Multiplying the actual odd-step identities gives

    2^A=n product_i(3+1/n_i)>3^tau n>=3n.

Equation(1) gives A<=2H. Hence 3n<4^H and H>=h(n). Meanwhile

    W(n)=2/[(H+1) binom(H,K)]<=2/[H(H+1)]<=B(n),            (4)

because 1<=K<=H−1 and binom(H,K)>=H. This proves the upper bound for all
rooted n; zero scores satisfy it too. The function B is nonincreasing.

For attainment, put h=h(N) and R_h=(4^h−1)/3. Since 3N is an integer
strictly below4^h, R_h>=N. It is an odd integer with one odd step directly
to1, valuation2h, K=h−1 and tau=1. Therefore W(R_h)=2/[h(h+1)].
This proves(3). The case N=1 instead has maximum W(1)=1. QED.

The envelope is sharp on every threshold interval. The exponent two in
B(N)~2/(log_4 N)^2 cannot be improved in a uniform pointwise upper bound
for this W, because the explicit ROOT ray attains it. This is not a claim
about average weights or typical trajectories. Total mass and pointwise
maxima are distinct observables.

## C2. A source-only quadratic discount along the sibling ray

Write S(n)=4n+1, ell=floor(log_2 n), and j>=0. For every positive odd n,

    W(S^j(n)) <= gamma_j(n) W(n),
    gamma_j(n) = (ell+2)(ell+4) /
                 [(ell+2j+2)(ell+2j+4)].                  (5)

Thus the integrated score has a usable multiplier tending to zero, even
without knowing the source's ROOT word. The earlier inequality with
multiplier1 remains valid but is weaker.

**Proof.** Rootedness is equivalent along siblings. For rooted n>1, S^j
adds2j to the first valuation, so K increases by j and tau is unchanged.
Consequently the exact ratio is

    W(S^j(n))/W(n) = product_{r=0}^{j−1}
                           (K+1+r)/(K+tau+2+r).           (6)

The trajectory identity and n_i>=3 give

    2^A=n product_i(3+1/n_i)<n4^tau<2^(ell+1+2tau).

Since A is integral, A<=ell+2tau. Combining this with e>=1 in(1) gives
2K<=ell+tau−1. Each factor of(6), increasing with K, is at most

    (ell+tau+1+2r)/(ell+3tau+3+2r)
       <= (ell+2+2r)/(ell+6+2r).

For the second inequality differentiate the first ratio in tau: its
derivative has numerator −2ell−4r<=0, and tau>=1. The product telescopes
to(5). At n=1 the strict word is empty, so handle it separately:
W(S^j1)=2/[(j+1)(j+2)]=gamma_j(1)W(1). At an unrooted source both sides
are zero. QED.

No exact *constant* multiplier for W was asserted in the previous note;
the ratio depends on the counters. Equation(5) is a new upper domination,
not an assertion that integrating theta^j preserves an exact multiplier.

## C3. The independently positive ROOT target with the new weights

Let n=2m+1 and j>=1. Retain the native context

    q=(3m+2)(4^j−1)/3,   m+q=(S^j(n)−1)/2,
    nu(a)=(a+1)^−6−(a+2)^−6,
    M_q=sum_{k>=0}nu(m+kq)W(n+2kq).

Let L_q be *any independently justified lower bound* for this M_q. Set

    h_q=h(n+4q),
    A_q=nu(m)+gamma_j(n)nu(m+q),
    E_q=1/[8 h_q(h_q+1) q^7].                             (7)

Then

    W(n)>=max(0,(L_q−E_q)/A_q).                           (8)

In particular **L_q>E_q suffices for ROOT at this specific source**.

**Proof.** C2 bounds the k=1 term by gamma_j(n)nu(m+q)W(n). Every
remaining label n+2kq, k>=2, is at least n+4q, so C1 bounds its score
by2/[h_q(h_q+1)]. The inherited kernel estimate, reproduced here, is

    sum_{k>=2}nu(m+kq)
      <=6q^−7 sum_{k>=2}k^−7
      <=6q^−7(2^−7+integral_2^infinity x^−7 dx)
      =1/(16q^7).

Multiplication gives E_q. Thus M_q<=A_qW(n)+E_q, proving(8). QED.

Relative to1/(48q^7), the threshold improves by the exact factor
h_q(h_q+1)/6. Relative to the prior source coefficient, gamma_j(n)<1
also improves the denominator for j>0. Along j->infinity,

    E_q=O(q^−7/(log q)^2),       gamma_j(n)=O_n(j^−2).

One may use j=1,2,4,8,... for nested moduli, because
q_(2j)=(4^j+1)q_j. A fixed positive gap is therefore not lost merely
through the deterministic error or coefficient allowance. Nothing here
proves that the actual L_q has that gap.

The exact envelope can also be retained alias by alias:

    residual <= sum_{k>=2}nu(m+kq) B(n+2kq),               (9)

which is at most E_q and can be bounded with exact finite sums plus an
integral tail. The uniform envelope is optimal pointwise; summing it
does not assert that all aliases simultaneously attain their individual
maxima. A finite bank of already certified off-source aliases cannot by
itself exceed their valid total capacity and manufacture the missing
source positivity.

## Next discriminating questions

1. Can a native joint measurement supply L_q across the nested moduli,
   with error smaller than the explicit q^−7/(log q)^2 capacity?
2. Can the actual source constrain *several* alias weights jointly more
   strongly than the sum of their sharp individual envelopes? The prior
   fixed-word rigidity excludes a fixed finite common-future explanation
   for every second alias; source-dependent words or a different context
   remain legitimate targets.
3. Can the already proved bounded-distortion bridge transport such a
   lower measurement while retaining its native guard and score multiplier?

These are specific missing estimates, not assumed lemmas. The geometric
analogue in [the boundary note](golden_taxicab_boundary_20261011.md) is to
keep the joint available-colour sets: their individual cardinalities do
not determine whether adjacent hidden vertices can be coloured.

The subsequently received [fixed-universe quantization work](unit_equality_quantized_guards_20261011e.md)
also identifies a boundary of this comparison. Its graph count has a
known positive integer spacing once a chromatic premise is authenticated.
Our height estimate bounds H from below; C2 bounds K relative to tau but
does not bound tau. Thus it gives no source-only bound on the denominator
of a hypothetical positive W(n). An integer-rounding argument still owes
that additional premise. This is consistent with the older
[shell-denominator audit](collatz_shell_denominator_budget_20261005.md),
which already rejects integrality of actual readouts inferred merely from
integer polynomial coefficients.

## Reproduction

[Checker](../../04-computation/experiments/golden_taxicab_boundary_20261011.py)
and [exact output](golden_taxicab_boundary_20261011.json) share a run with
the geometry and pencil notes. The Collatz universe is every odd source
through9999, checked by actual forward first-hit words with cap10000;
every n<1000 has siblings j=1..8 replayed; every integer threshold
2..19999 has its attaining root-ray witness checked. Twenty-four source/
modulus controls compare actual first20 residual aliases with(7). This
finite bank is verification, not a proof of convergence outside the bank.

    python3 -B 04-computation/experiments/golden_taxicab_boundary_20261011.py
    python3 -B -O 04-computation/experiments/golden_taxicab_boundary_20261011.py
