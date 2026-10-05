# Two clocks expose an exact Collatz source-address energy

2026-10-04–05 session. **PROVED:** the source-address partition,
full-Haar orthogonality, odd-Haar nearest-neighbor covariance, and the
weighted difference-square identity below. **FINITE-EXACT:** accompanying
integer/Fraction controls. **OPEN:** the fixed-seed estimate, H1, and Collatz.
These are deductions from finite parity coding; no priority claim is made.

Artifacts: [program](../../04-computation/experiments/source_address_two_clocks_20261004.py),
[exact output](source_address_two_clocks_20261004.out),
[structured controls](source_address_two_clocks_20261004.json).

## Inheritance and the reusable move

The anchor is the source-address reformulation of H1 in
[carry interfaces, section 4](collatz_carry_interfaces_20261004.md).
The niche is the incoming exact random-start second moment; the wildcard
is a discrete energy representation. The board is **binary address / odd-step
clock / halving clock / character / boundary / exceptional seed**.

The incoming batch through `9ba31b1b3b` was read in full: the exact-moment
arguments, two moment programs, integer-seed experiments, HYP-9176 update,
and the explicit correction that heavier q=5 tails do not imply a lower
typical rate. Their second-moment proof identifies distinct path characters.
Their finite ensembles do not supply a bound at the fixed seed 1.
See [Fourier experiments, sections 4k-4l](collatz_fourier_experiments_20261004.md).

The closest proved mechanism is finite inverse-parity coding. The hostile
is transferring an almost-everywhere estimate to a specified integer.
The corrected near miss is forgetting which clock truncates an infinite
sum. The useful sidecar is the pair (odd-step count, halving cost).

The reusable move is: **change the grading of an exact expansion, retain
both truncations, and compute their compatibility before exchanging limits.**
Here the change of grading makes one sum elementary and the other carry H1.

## 1. All compositions at fixed cost partition the odd binary addresses

Fix odd q>=3. For a valuation word w=(a1,...,am), every ai>=1, put

    A=sum ai, Q=2^A, P=q^m,
    B_w=sum_(i=1..m) q^(m-i) 2^(a1+...+a_(i-1)),
    r_w=-B_w*P^(-1) mod Q,  0<r_w<Q.

The r_w are odd. At a fixed A, letting m range from 1 to A gives a
bijection from all compositions of A to all odd residues modulo 2^A.

Proof: use the shortcut map T_q(n)=(qn+1)/2 for odd n and n/2 for even n.
Every A-bit parity string determines one residue modulo 2^A. Inductively,
after choosing its first bit and the unique tail residue t modulo 2^(A-1),
its predecessor is 2t or q^(-1)(2t-1), respectively. Both are unique and
have the required parity. Words w encode exactly the strings

    1 0^(a1-1) 1 0^(a2-1) ... 1 0^(am-1).

Their affine endpoint is (P*n+B_w)/Q, giving the residue r_w above.
There are 2^(A-1) such strings. At this precision the last valuation may
continue beyond the horizon; the extra endpoint-oddness bit modulo 2Q
is deliberately absent. These are prefix addresses, not first-hit receipts.

For q=3 and A=3 the words and addresses are

| Word | Odd-step count | r_w modulo 8 |
|---|---:|---:|
| (3) | 1 | 5 |
| (1,2) | 2 | 3 |
| (2,1) | 2 | 1 |
| (1,1,1) | 3 | 7 |

Thus all finite nonempty words, of all lengths and costs, index each
nontrivial dyadic character exactly once. For q=3 this is the elementary
finite form of the inverse-parity expansion in
[Bernstein and Lagarias, formula 1.6](https://websites.umich.edu/~lagarias/doc/bernstein.pdf).
The proof above specifies the conventions needed here for every odd q>=3.

## 2. The full-Haar sum is universal although its layers depend on q

For R in the 2-adic integers, define chi_(r/2^A)(R)=exp(2*pi*i*r*R/2^A),
using R modulo 2^A. Write

    V_m(R)=sum_(length(w)=m) 2^(-A(w))*chi_(r_w/2^A)(R).

Each sum converges uniformly and absolutely, since its weights sum to one.
At an ordinary integer u, V_m(u) is nu_hat_m(u) of carry interfaces.
The incoming path phase is B_w/2^A; at each fixed m multiplication by
the odd unit -q^(-m) preserves full and odd Haar measure. Its second
moment therefore agrees with ours. Across different m the common source
coordinate R must be retained; those multipliers depend on m.

With normalized Haar measure on Z_2, the characters are orthonormal.
The partition theorem gives

    <V_m,V_n>=0 for m!=n,
    ||V_m||_2^2=sum_w 4^(-A(w))=(sum_(a>=1)4^(-a))^m=3^(-m).    (1)

Their sum converges in L2. It also converges absolutely for almost every R:
sum_m E|V_m| <= sum_m 3^(-m/2)<infinity.
Regrouping its L2 Fourier expansion by halving cost is now legitimate:

    S(R)=sum_(A>=1) 2^(-A) sum_(r odd mod 2^A) chi_(r/2^A)(R)
        =(v2(R)-1)/2, almost everywhere.                         (2)

Indeed for v=v2(R)<infinity the inner sum is 2^(A-1) for A<=v,
-2^(A-1) for A=v+1, and zero afterward. At R=0 the cost partial
sums diverge as A/2; it is a hostile null seed. The cost expansion has
squared L2 tail 2^(-D-1) after cost D, so this identification also holds
in L2. Its squared norm is 1/2, agreeing with sum_m 3^(-m).

Although (2) is q-independent, the allocation of characters to each V_m
depends on q. The universal total has not erased the dynamical question
about its particular layers at a chosen seed.

## 3. Odd seeds give a nearest-neighbor energy

Now use normalized Haar measure on the odd coset 1+2 Z_2, denoted E_o.
The mean of a dyadic character is 1 for the trivial character, -1 for
the character 1/2, and zero for every other character.
Consequently distinct source characters interact only when their
addresses differ by 1/2. Such a pair has the same cost A>=2, and its
residues differ by 2^(A-1).

Their first A-1 shortcut parities agree and their last parities differ.
Therefore their numbers of odd steps differ by exactly one. More explicitly,
if the shorter word ends in am>=2, the longer word is

    (a1,...,a_(m-1),am-1,1).

This proves the exact Gram matrix

    E_o[V_m]= -1/2 if m=1, and 0 otherwise,
    E_o[|V_m|^2]=3^(-m),
    E_o[V_m*conj(V_(m+1))]=-1/(4*3^m),
    E_o[V_m*conj(V_n)]=0 if |m-n|>1.                         (3)

For the adjacent entry, sum the product weights over am>=2:

    -(sum_(a>=1)4^(-a))^(m-1) * sum_(a>=2)4^(-a)
      =-3^(-(m-1))/12=-1/(4*3^m).

This also supplies a hostile to using full-Haar orthogonality on odd seeds:
at q=3 the words (3) and (2,1) have addresses 5/8 and 1/8, and their
characters have odd-Haar inner product -1.

For every finite complex sequence z_1,...,z_M, set z_(M+1)=0.
Expanding (3) and collecting adjacent terms gives the energy identity

    E_o |sum_(m=1..M) z_m V_m|^2
      = |z_1|^2/4 + sum_(m=1..M) |z_m-z_(m+1)|^2/(4*3^m).     (4)

This is a weighted discrete difference energy, plus its left boundary.
Taking z_m=1 through M leaves only the right boundary:

    E_o |sum_(m=1..M) V_m + 1/2|^2 = 1/(4*3^M).              (5)

Equations (3)-(5) are all-depth identities, not fits to an ensemble.
They explain precisely how the cost sum's constant value on odd seeds
emerges from cancellation between adjacent odd-step layers.
Markov's inequality and the first Borel-Cantelli lemma imply that for
every rho>1/sqrt(3), the error in (5) is O_R(rho^M) for almost every
odd R. The implied bound depends on the seed. All ordinary integers form
a countable Haar-null set, so this conclusion certifies none of them
individually.

## 4. The unproved operation is switching clocks at the fixed seed

Retain both clocks by inserting a damping parameter 0<=t<1:

    V_m(t;R)=sum_(length(w)=m) (t/2)^A chi_(r_w/2^A)(R).

The double sum over m and words is absolutely convergent, since
sum_(A>=1) 2^(A-1)*(t/2)^A=t/(2*(1-t)). Hence for odd R,

    sum_(m>=1) V_m(t;R)=-t/2.                              (6)

Letting t tend to 1 proves an Abel-damped value of -1/2 at every odd
seed, including 1. It does **not** by itself identify the undamped
odd-step sum there. Pointwise interchange of those two limits requires
a separate uniform estimate. Equation (2) identifies the limits almost
everywhere, which is a weaker quantifier.

A concrete sufficient new target is to bound the boundary function

    E_M(1)=sum_(m=1..M) V_m(1)+1/2

by C*rho^M for some 1/2<rho<log2(3)-1. Here V_m(1) means t=1 and
seed R=1. Success would bound V_M(1)=E_M(1)-E_(M-1)(1) at the same
exponential rate and, by the O(2^(-M)) phase comparison in carry
interfaces, imply the fixed-frequency H1 bound. This is a sufficient
strengthening: H1 alone does not establish that the limiting boundary
constant is zero. Alternatively, seek the direct V_M(1) bound without
requiring that extra conclusion.

The exact mean-square formula (5) is a benchmark for a pointwise method,
not the missing method itself. A route through binary cylinders must
quantify oscillation within the cylinder containing 1 as its depth grows;
continuity of each individual V_m is not a uniform bound in m. The
fixed-integer Fourier comparison has a factor |u| and is not silently
integrated over Haar-distributed 2-adic frequencies.

Even H1's intended renewal consequence, in
[THM-4519: frequency-one coefficient carries the renewal series](../../01-canon/theorems/THM-4519-frequency-one-coefficient-carries-the-renewal-series.md),
is conditional on H1. Neither that conditional implication nor the
identities here prove universal paid coverage or convergence of all
positive integers. The exceptional-source obligation remains explicit.

## 5. Transfer record and reproduction

| Source | Map | Preserved predicate | Lost information / sidecar |
|---|---|---|---|
| Collatz words | finite inverse-parity address | all nontrivial binary characters, bijectively | odd-step grading must accompany cost |
| Character partition | Haar inner product | exact second moments | the chosen integer seed is averaged away |
| Odd-seed character pairs | flip the last parity bit | adjacent odd-step coupling | source labels needed for deterministic estimates |
| Adjacent covariance | weighted difference squares | exact boundary energy | L2 control does not imply point evaluation |
| Two clock sums | t^A damping | absolute summability for t<1 | uniform endpoint control at t=1 remains OPEN |

Run both:

    python3 -B 04-computation/experiments/source_address_two_clocks_20261004.py
    python3 -O -B 04-computation/experiments/source_address_two_clocks_20261004.py

The exact universe is every composition of every cost A=1,...,12 for
q in {3,5,7,9,11}: 20,475 words. An independent shortcut-map replay
checks every parity string. The program computes every truncated odd-Haar
Gram matrix through cost 12 directly from the paired addresses, then
checks its closed formula, boundary energy, and a nonconstant rational
coefficient probe. Cyclotomic-polynomial arithmetic independently checks
the cost sums through cost 8, at three damping values and thirteen
signed/zero seeds. No floating-point phase cancellation is used.
All 97,842 explicit checks pass identically with and without Python optimization.
The full infinite statements follow from the proofs above.
