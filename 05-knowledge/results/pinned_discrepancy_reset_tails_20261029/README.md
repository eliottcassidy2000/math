# Fractional pin barriers and independent first-reset residual tails

2026-10-29. **PROVED elementary deductions / CITED partial-colouring input /
FINITE-EXACT certificates / OPEN universal Collatz bounds.** Same-agent
proof review, with independent computational checks; no literature-priority
claim and no independent-agent audit. No new canon theorem ID is reserved.

## 0. Inheritance and the research board

Anchor: source-weighted Collatz residual tails. Niche: logarithm-free
completion of pinned cubic-orbit tournaments. Wildcard: the new ramified
Pisano attachment recurrence, tested against spectral growth.

The board is **fractional feasibility / rounding discrepancy / first-reset
fibres / original source weights / residue-conditioned tails / coherent
module growth**. The closest mechanisms are [the previous exact-moment and
signed-residual theorems](../proof_targets_refinement_20261029/README.md),
[orbit covariance P1](../quasirandom_order_lifts_20261029/README.md), and
[Session106's atomic Abel certificate](../session106_collatz_transfer_20261029/spectral/NATIVE_SPECTRAL_BOUNDARY.md).
The hostile is a frequency-aligned pinned bank on one side and a fixed
point in the abstract killed-map test on the other (not a claimed positive
Collatz counterexample). The corrected near miss is
replacing conditioned data by a fresh random law. The least-used sidecars
are a fractional dual separator and the reset endpoint's 3-adic valuation.

The incoming Session108/109 proofs, synthesis and correction were read at
`5ab37f308003`: [operational sigma symmetries](../session108_sigma_symmetry_20261030/README.md),
[Pisano synthesis](../session109_pisano_symmetry_20261030/SYNTHESIS.md), and
[ramified attachment theorem](../session109_pisano_symmetry_20261030/tournaments/results.md).
The corrected Pisano period needs the return-depth sidecar; we do not use
an unconditional no-delay lifting formula. The new coherent attachment
law supplies the exact spectral hostile in §4, rather than a new Collatz
clock. The meta-patterns used are **type the residual** and **canonicalize
mechanics while retaining the predicate**; no new meta-pattern is claimed.

## 1. Logarithm-free completion, with the pin obstruction exposed

### 1.1 The actual Fourier completion problem

Keep the preceding carrier: N>=3 odd; H of odd order r in `(Z/N)^*`;
every nonidentity h has `gcd(h-1,N)=1`. Pair its free nonzero orbits
O_i and -O_i, giving b=(N-1)/(2r) independent signs. The tournament
matrix K has symbol epsilon_i on O_i and -epsilon_i on -O_i.

Let M=(N-1)/2, and define the real M-by-b matrix

    C_(k,i)=2 sum_(x in O_i) sin(2 pi kx/N),  1<=k<=M.

Then `||K||=||C epsilon||_infinity` exactly. Each entry has absolute value
at most L=2r. The inherited covariance gives row square sum N, while
Fourier orthogonality gives `C^t C=Nr I_b`. In particular the orientation
variables are independent coordinates even when their frequency peaks
can cancel. Constants are annihilated; negative frequencies duplicate
absolute values. No Fourier mode is omitted because N is composite.

Prescribe signs on P, let F be the m remaining indices, and set
`d=C_P epsilon_P`. The *fractional pin barrier* is the finite linear program

    t_* = min_{x in [-1,1]^m} ||d+C_F x||_infinity.        (1)

Its exact dual is

    t_* = max_{||y||_1<=1} [<d,y>-||C_F^t y||_1].         (2)

This follows from finite linear-program duality, or exchanging the two
compact convex optimizations. Every y in the displayed ball is a lower
certificate. With p prescribed signs, Parseval also gives

    t_* >= sqrt(N r p / M).                               (3)

Indeed the squared Euclidean norm of the full frequency vector is
`Nr(p+||x||_2^2)`. Free variables can reduce a peak, but cannot erase
the pinned energy. Thus t_*=0 occurs only for the empty pin set.

### 1.2 Rounding theorem and quantitative pin criterion

**Theorem A.** For m>0 there is an integral completion retaining every
prescribed sign and satisfying

    t_* <= min_completion ||K||
        <= t_* + 18 L sqrt(m log(16M/m)),   L=2r.          (4)

Here log is natural. If m=0, the optimum equals t_* exactly. Consequently,
for fixed bounded r, a sequence of pinned carriers has O(sqrt N)
completions **iff** its fractional barriers are O(sqrt N). In particular
the unpinned cubic Sylvester carriers have O(sqrt N) completions; the
logarithm from the previous moment construction is unnecessary for this
existence statement. This is not a bound independent of arbitrary pin
debt, and not a uniform-in-r theorem.

**CITED input.** Lovett–Meka, [*Constructive Discrepancy Minimization by
Walking on the Edges*, Theorem4](https://arxiv.org/html/1203.5747#S2), permits
partial rounding from an arbitrary cube point. For row vectors v_j and
thresholds c_j with `sum exp(-c_j^2/16)<=n/16`, at least half the coordinates
can reach cube faces with row changes at most `c_j ||v_j||_2`.
Their stated approximate version implies the exact existence version by
compactness as its tolerance tends to zero. We use the arbitrary-vector
lemma, not a mistaken substitution of real Fourier entries into a theorem
stated only for incidence sets. The explicit threshold condition is our
input; we do not use the introduction's displayed `log(m/n)` at m=n.

**Proof of (4).** Start at a minimizing fractional point in (1). At a
stage with n unrounded coordinates, restrict each row to those coordinates;
its Euclidean norm is at most L sqrt(n). Choose every threshold as
`4 sqrt(log(16M/n))`. The cited condition holds with equality. Freeze the
at least half integral coordinates and repeat. Previously prescribed
coordinates never move. With n_s<=m/2^s, the total change in any row is at
most

    4L sum_(s>=0) sqrt((m/2^s) log(16M 2^s/m)).

Set a=log(16M/m)>=log16. Concavity of the square root bounds this by
`4L sqrt(ma) sum 2^(-s/2)(1+s/8)`.
The multiplier is `10+(11/2)sqrt2 <18`, proving (4).
This derivation applies to any real matrix with bounded entries and any
fractional starting point. It is the classical partial-colouring move,
made explicit for these pinned Fourier coordinates.

Since `m<=M/r` and `u log(16M/u)` is increasing for 0<u<=M,
the rounding allowance is at most `36 sqrt(r M log(16r))`. For bounded r
this is O(sqrt N). The lower direction of the iff is immediate from
fractional relaxation. No numerical optimization proves the existence
theorem; it follows from the cited lemma and the calculation above.

**Cheap sufficient pin budget.** Taking x=0 in (1) gives t_*<=||d||.
Therefore p pins cost at most2rp, and p=O(sqrt N/r) always permits a
logarithm-free completion for bounded r. Much larger structured pin sets
may pass (1); their actual joint frequency data determines the answer.

**A positive-size free bank does not remove pin debt.** On any frequency,
`sum_i C_(k,i)^2=N` and `|C_(k,i)|<=2r`, so `sum_i |C_(k,i)|>=N/(2r)`.
Pin all but m signs to align with that frequency. Every fractional or
integral completion then has norm at least

    N/(2r) - 4rm.

One factor2rm accounts for discarding the free columns from the aligned
sum, and another for their strongest possible opposing contribution.
Thus m<=N/(16r^2) still forces norm at least N/(4r). This supplies an
all-size hostile with genuinely free variables. Omitting the second
free-column cost would give an invalid lower bound.

**Scope left open.** The rounding theorem is constructive in the classical
randomized partial-colouring framework, but this packet does not implement
that algorithm or claim a new deterministic polynomial-time logarithm-free
construction. The experiment below proposes completions by a numerical LP
and local search, then verifies them independently. Growing r and sharper
constants require additional structure.

## 2. Exact source-preserving first-reset residuals

Write an odd Collatz source as n=2j-1, j>=1, and kill the genuine terminal
j=1. Keep the original law `w(j)=1/[j(j+1)]`; its active mass is1/2.
This is the same source prior as the inherited Abel certificate. A native
valuation-one step is exactly `j -> 3j/2` and occurs iff j is even.

For each active source write `j=2^k ell`, ell odd, k=v_2(j). Its actual
first-reset input and following valuation/endpoint are

    R(j)=3^k ell = y,
    a=1+v_2(3y-1)>=2,
    B(j)=( (3y-1)/2^(a-1) + 1 )/2.                        (5)

There are k valuation-one steps and one reset step. All intermediate
sources are positive and legal. This finite stopping rule exists for
every j independently of the Collatz conjecture. At j=14 (odd source27)
it ends at j=16 (odd source31), so it is not necessarily descent below
the original source.

### 2.1 Exact incoming law and its arithmetic spikes

For 0<z<=1, the discounted residual at the reset input is supported on
odd y>=3 and is exactly

    rho_z(y) = sum_(k=0)^(v_3(y))
                  z^k w(2^k y/3^k).                     (6)

It is zero elsewhere. Every term is a distinct actual original source.
Formula (6) is a finite inverse-fibre sum, not a Haar measure or a law
reinitialized at y. For example

    (3^k)^2 rho_1(3^k) >= (9/4)^k /(1+2^(-k)).            (7)

Thus no constant C gives rho_1(y)<=C/y^2 at all endpoints. A pointwise
quadratic-decay envelope is already false after this one legitimate
prefix operation, even though the original prior has precisely that
decay. Its lost coordinate is v_3(y).

There is nevertheless an independent, all-cutoff tail bound. For real
J>=1 let `P_z(J)=sum_(y>J) rho_z(y)`. Then

    P_z(J) <= 2 / ((1-3z/4)J),
    lim_(J->infinity) J P_z(J) = 1/[2(1-3z/4)].           (8)

**Proof.** Sum over the actual sources 2^k ell, ell odd. For each fixed k,
the tail of `w(2^k ell)` above ell>J/3^k is asymptotic to
`3^k/(2*4^k J)`. Also `w(2^k ell)<=4^(-k)ell^(-2)` and

    t sum_(ell odd, ell>t) ell^(-2) <=2,       t>0.

For t>=1 this follows from the first term plus the integral on a lattice
of spacing2; for t<1 use `sum ell^(-2)<2`. Hence the k contribution to
J P_z(J) is at most `2(3z/4)^k`. Dominated convergence proves the limit,
and summing the bound proves the inequality. At z=1 the exact leading
coefficient is2, rather than the original coefficient1. QED.

### 2.2 Completing the reset restores the coefficient, but changes phase

Let sigma_z be the residual after the *whole* block (5): each source j
contributes `z^(k+1) w(j)` at B(j), unless B(j)=1, which is killed. Put
`S_z(J)=sum_(h>J) sigma_z(h)`. Uniformly for J>=1,

    S_z(J) <= 6z / ((1-3z/4)(4J-1)).                     (9)

For fixed z its exact leading tail coefficient is

    lim_(J->infinity) J S_z(J) = alpha(z)=z/(4-3z).        (10)

The endpoint residue coefficients are respectively

    h=0 mod3: alpha(z)/5,
    h=1 mod3: 4alpha(z)/5,
    h=2 mod3: 0.                                        (11)

**Proof of the bound.** Equation (5) gives `B(j)<=(3R(j)+1)/4`.
Apply (8) at `(4J-1)/3` and multiply by z.

**Proof of the limits.** For fixed k and a>=2, ell lies in the unique
odd class satisfying

    3^(k+1) ell = 1+2^(a-1) mod 2^a.

The endpoint is `(3^(k+1)ell+2^(a-1)-1)/2^a`. Its tail has asymptotic
coefficient

    z^(k+1) * 3^(k+1) / (4^k 4^a).

For a fixed k, the progression-tail estimate
`sum_(ell>L, ell=c mod q) ell^(-2)<=L^(-2)+1/(qL)` and
`L>=(2^a J)/(2*3^(k+1))` give a summable-in-a majorant after multiplying
by J. For the sum over k, the bound used in (8), with the endpoint
inequality from (9), gives the summable majorant `2z(3z/4)^k`.
These justify both limit interchanges, not only a formal geometric sum.
Summing over k>=0 and a>=2 gives (10).

The odd Collatz endpoint is `2h-1=(-1)^a mod3`. Even a therefore gives
h=1, odd a gives h=0, and h=2 is impossible. The sums of4^(-a) for even
a>=2 and odd a>=3 are1/15 and1/60, giving the ratio4:1 in (11). QED.

At z=1 the total tail coefficient is exactly1. Consequently a strict
contraction of that coefficient at every undiscounted reset is impossible
already for this prior. For z<1 the coefficient contracts, but alpha(z)
tends to1 as z tends to1. These facts do not bound survival under
unboundedly many resets. In particular we do not multiply the factor
alpha(z) through later blocks without authenticating their conditional
source law. Its first output is zero on one residue class and has the
arithmetic spikes inherited from inverse fibres; it is not a fresh w.

### 2.3 These are actual Abel trial flows, with an independent residual tail

For 0<z<1, let v_z place `w(j)z^i` at each visited active state before
the whole first-reset block finishes, for each original j>=2. Then

    e_z=w+zQv_z-v_z = sigma_z                              (12)

for the original killed operator Q. This follows by telescoping each
finite native path, including the actual terminal endpoint convention.
The trial flow is summable uniformly in z:

    ||v_z||_1 <= sum_j w(j)(v_2(j)+1) < 7/6.              (13)

Indeed `sum_j w(j)v_2(j)=sum_(h>=1) sum_(ell>=1) w(2^h ell)` is less
than `2 sum_(h>=1)4^(-h)=2/3`; add the active mass1/2.

Thus (9) is an independently proved tail bound for the residual of an
explicit, summable, source-preserving trial flow. It can be used in the
previous signed-residual certificate without simulating unknown complete
orbits. Its scope is one whole reset block, not universal stopping.

The full residual mass at z=1 is determined by the root family below.
It is positive. Therefore (12) alone does not make the Abel bound tend
to zero. The next arithmetic obligation is to control later signed or
positive residuals on their actual conditional fibres, uniformly in both
the number of blocks and z tending to1. Merely retaining the coefficient
in (10) loses the phase information in (11).

## 3. An exact first-block terminal family and a certified mass

The whole block (5) lands at1 exactly at the following source indices:

    j_(k,t)=2^k (2^((2t+1)3^k)+1)/3^(k+1),   k,t>=0,    (14)

excluding j_(0,0)=1, which was already terminal. The original odd integer
is2j_(k,t)-1. In particular the t=0 indices begin1,2,76,13256072.
These are genuine first-reset receipts, not a statement that the repeated
number76 has the same meaning in every earlier construction.

**Proof.** Landing at1 means `3^(k+1)ell=2^(a-1)+1`.
The order of2 modulo3^(k+1) is `2*3^k`, and its unique order-two power
is -1; hence `a-1=(2t+1)3^k`. Conversely (14) is integral, has exactly
v_2(j)=k because its remaining factor is odd, and has final valuation
a=1+(2t+1)3^k. Direct substitution proves the complete native word.
The order formula follows by induction from
`v_3(2^(2*3^k)-1)=k+1`, starting with2^2-1=3. QED.

There are no duplicate sources because the first block has a unique k
and a. This gives the exact first-block terminal mass

    M=sum_(k,t, (k,t)!=(0,0)) w(j_(k,t)).

The checker certifies

    0.258387760297 <= M <= 0.258387760298.                 (15)

Together with the already terminal source of weight1/2 this certifies
total source weight between0.758387760297 and0.758387760298. This is a
closed-form unbounded-family result, not a claim to improve the coverage
of the stronger existing multi-rule atlases.

The truncation is rigorous: increasing t at fixed k at least doubles j;
increasing k at t=0 also at least doubles j, since with a=2^(3^k) its
ratio is `2(a^2-a+1)/3>=2`. Doubling j reduces w by at least a factor3.
Therefore omitted t>T at fixed k costs at most `3w(j_(k,T+1))/2`,
and all k>K cost at most `9w(j_(K+1,0))/4`. The checker uses K=4,T=12,
64 active sources and these independent infinite-tail bounds. It does
not infer completeness from testing a finite range.

## 4. Incoming ramified coherence is a spectral obstruction

Session109 proves that coherent lifts in the ramified metallic ring are
actual lexicographic products: each vertex fibre has p^2 vertices, and
the kernel tournament is regular. If K_A and K_B are regular tournament
matrices of sizes a,b, their product matrix is

    K_(A[B])=K_A tensor J_b + I_a tensor K_B.

The subspace constant on each fibre and its orthogonal complement are
invariant, so the norm is exactly `max(b||K_A||,||K_B||)`.
In the coherent ramified tower, N_k=p^(2k), p>=5. The first term dominates
at each step: `||K_A||>=sqrt(a)` and `||K_B||<=b-1`. Therefore

    ||K_k||=p^(2k-2)||K_1||,
    ||K_k||/N_k=||K_1||/p^2 >0.                          (16)

Thus fully coherent towers are not quasirandom. A compact attachment
recurrence and many free kernel choices do not remove the fixed quotient
mode. This does not exclude low-norm **noncoherent** ramified banks.
Their orbit sizes grow with depth, so the bounded-r conclusion of §1
cannot simply be copied to them. Relative child signs are precisely the
new variables a fractional completion would need to retain.

This explains a real common obstruction: preserving only a parent
summary can freeze the coordinate that needs to cancel. In Collatz that
coordinate is an actual source/residue mass; in the tournament it is a
Fourier mode of a legal sign assignment. We have not identified the two
operators or licensed choosing signs in the deterministic Collatz map.

## 5. Verification and next targets

`check_pins.py` uses SciPy only to propose fractional points, dual vectors
and sign completions. It then recomputes both sides of (2) with rational
witnesses. Sines are enclosed using Machin's formula, an alternating atan
series and a Taylor polynomial with a rational remainder. Fixed-point
interval arithmetic uses192 internal bits; the initial100-bit attempt was
correctly rejected because Horner interval widths exceeded the requested
error. Final coefficients have a proved error at most4r/2^40. Every
reported spectral interval includes the accumulated error. Small carriers
also enumerate every legal sign completion. Diagnostic decimals are not
certificates, and the numerical search is not the cited rounding algorithm.

`check_resets.py` checks the exact native words for20,000 sources,1,002
inverse fibres, and the terminal-family iff on the same source range.
It sums one million original source weights using integer rounding with
explicit total error and the exact omitted mass1/(1000001), at four
cutoffs and two discounts. These finite data test the derived laws; the
proofs above establish their infinite scopes. The root-mass enclosure
uses the separate geometric bounds following (15).

Reproduce both scripts in ordinary and `-O` modes. Dependencies for the
pin script are NumPy and SciPy; the reset script uses only the standard
library. Saved receipts are `pins.json` and `resets.json`.

The logarithm-free **existence** target is closed for bounded orbit size
with the explicit fractional pin criterion. Open extensions are a fast
deterministic constructor, constants useful at moderate N, and growing
orbit sizes. The Collatz next target is a uniform family of conditional
residual-tail estimates through successive resets; the independent
one-block bound is now supplied, with its sharp coefficient obstruction
and the residue information needed for the next step.

The [finite-depth continuation](FINITE_DEPTH.md) sharpens that target:
at every fixed b, `liminf J S_(b,z)(J)>=[z/(4-3z)]^b`. At z=1 no
finite reset depth lowers the critical tail coefficient below1. The
source27 additionally forces C>=13/6 for an all-cutoff, all-depth C/J
envelope. An elementary degrading power bound supplies independent tails
at every finite depth. Uniform tightness of these positive residuals is
equivalent to boundedness of every original orbit; it is itself a major
arithmetic target, not a consequence of a formal renewal equation.

### Reproduction and selected certified comparisons

From the repository root, run each command again with `-O` and compare
its stdout byte-for-byte with the saved receipt:

```bash
OPENBLAS_NUM_THREADS=1 python3 -B 05-knowledge/results/pinned_discrepancy_reset_tails_20261029/check_pins.py
python3 -B 05-knowledge/results/pinned_discrepancy_reset_tails_20261029/check_resets.py
OPENBLAS_NUM_THREADS=1 python3 -B 05-knowledge/results/pinned_discrepancy_reset_tails_20261029/check_words.py
```

All table entries below are outward decimal coarsenings of the rational
certificates in `pins.json`; they are not floating-point LP assertions.

| N | Pinned orbit pairs | Fractional barrier enclosure | A certified completion norm enclosure |
|---|---:|---:|---:|
| 43 | 2 | [5.128571,5.128573] | [6.557438,6.557439] |
| 43 | 6 | [13.771725,13.771727] | [13.771725,13.771727] |
| 133 | 8 | [9.961800,9.961804] | [18.520259,18.520260] |
| 1807 | 30 | [66.761054,66.761074] | [82.651452,82.651453] |
| 1807 | 280 | [574.008645,574.008646] | [574.008645,574.008646] |

The small N=43 cases exhaust every legal completion. The larger cases
give feasible witnesses and lower separators, not an exhaustive optimum
claim. The 280-pin case nearly closes its own interval because the
fractional obstruction is already dominant. The complete review is in
[AUDIT.md](AUDIT.md); `words.json` retains the exact ramified tensor and
source27 receipts.
