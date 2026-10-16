# Collatz break energy must retain the integer clock

2026-10-13, integration of incoming `58b8699d89` / `bd3d32cc31`.
**PROVED:** the complete jump/plateau energy identity and a positive weighted
criterion implying the dyadic energy target. **FINITE-EXACT:** an independent
census through 2^27, cross-checked against every saved parity bit through 2^26.
**REFUTED:** the incoming factor and expectation-to-certificate implication;
also the claim that the finite event word alone determines source energy.
**OPEN:** the required weighted bound for actual Collatz sources.

**Continuation, 2026-10-15:** [the clock-mesh theorem](collatz_clock_mesh_20261015a.md)
now proves a sharp arithmetic refinement allowance and upper certificates
monotone under every physical interval split. [Signed terminal-weight
controllers](collatz_terminal_mesh_20261015a.md) bound the remaining means
on six finite windows without grounding every source. Uniform coarse-mean
and terminal-residual control remains OPEN.

The [incoming break note](collatz_dyadic_energy_breaks_20261013.md) supplies a
useful second coordinate system for the
[source-preserving wavelet energy](collatz_wavelet_energy_20261013d.md).
It writes discrepancy as a cumulative current at sign changes. Retaining
the *integer lengths* between those changes makes this a precise energy
representation. Dropping those lengths changes the counting measure.

## 1. Exact formula on the original source interval

Set `W(n)=0,-1,+1` according as n=0,1,2 mod 3. For any real label sequence t,
summation by parts gives

    H(N)=W(N)t(N)-sum_(n<N) W(n)(t(n+1)-t(n)).               (1)

At every nonzero jump n_i define
`v_i=W(n_i)(t(n_i+1)-t(n_i))/2`, and
`J(N)=sum_(n_i<N) v_i`. Then H(N)=W(N)t(N)-2J(N).
This includes ROOT-supported labels in {-1,0,1}; v_i may be half-integral.
The incoming shortcut `v_i=W(n_i)t(n_i+1)` and the alternating-sign law
are valid on a **binary** labelled interval. They are not general identities
for {-1,0,1}: the formal transition 1->0 already violates
`Delta t=2 t_next [t_next!=t]`. This is a type hostile, not a claimed
nonconvergent Collatz source. All finite sources in the censuses are grounded.

Partition the integer window `[X,2X)` into constant-label plateaus I_i.
Retain the current J_i accumulated from source 1, including jumps before X,
the label s_i on the plateau, and its clipped length L_i. Put
`A_i=sum_(n in I_i)W(n)`, `B_i=sum_(n in I_i)W(n)^2`. Exactly,

    E_X := sum_(X<=n<2X) H(n)^2
         = sum_i [4 L_i J_i^2 - 4 J_i s_i A_i + s_i^2 B_i]. (2)

The periodic coefficients satisfy |A_i|<=1 and B_i<=L_i. Thus the gap
lengths and the incoming current belong to the proof state. Resetting the
current to zero at the start of each octave would change (2).

Equivalently let `P_(n,i)=[n_i<n]` on the same integer window and
`b_n=W(n)t(n)`. The kernel

    K_X(i,j)=max(0,2X-max(X,max(n_i,n_j)+1))= (P^T P)_(i,j)

is positive semidefinite, and

    E_X=4 v^T K_X v-4 v^T P^T b+b^T b.                    (3)

Include every jump before 2X, not just those inside the window. Formula (3)
fits the Gram-certificate transport rule in the wavelet note exactly.
The map from event coordinates to source coordinates is P; the lost
information in an unweighted event summary is K_X and the boundary vector.

## 2. A positive weighted target with an explicit threshold

For |t|<=1 put

    S_X=sum_i L_i J_i^2 = sum_(X<=n<2X) J(n)^2.

The triangle inequality and its reverse in Euclidean norm give

    (max(0,2 sqrt(S_X)-sqrt(X)))^2 <= E_X
                <= (2 sqrt(S_X)+sqrt(X))^2.                 (4)

Consequently, if an independent proof gives `S_X<=alpha X^2` with
`alpha<1/2`, then the desired `E_X<=2X^2` follows whenever

    X >= (sqrt(2)-2 sqrt(alpha))^(-2).                       (5)

This is a deterministic source-preserving sufficient rule. It works also
with zero ROOT labels. If S_X=O(X^2), then
`E_X/X^2=4 S_X/X^2+O(X^(-1/2))`, so retaining this positive statistic costs
only a vanishing normalized boundary correction. Proving (5)'s premise
uniformly at an explicit starting scale remains OPEN. An expected bound
on a different ensemble does not discharge it.

### Same event history, different source energy

Take binary labels initially+1, flipping at the following four positions:

    A: 34,35,59,61;        B: 55,56,59,61.

Both have residues `(1,2,2,1)` and the identical current word
`(+1,+1,-1,-1)`. Hence all unweighted statistics of that finite word agree.
But on the same dyadic window `[32,64)`, direct calculation gives

    E_A=417, E_B=81, E_A-E_B=336.

The shift moves 21 integers from current 0 to current 2, so the leading
change is `4*21*2^2=336`; complete periods make the other changes cancel.
These artificial sequences are **not claimed to obey Collatz**. They
refute the inference from event-word data alone to integer-window energy
and identify the missing coordinate that a Collatz argument must retain.

### A universally valid finite-block compiler, with two positive budgets

For any consecutive integer block A=[a,b), retain its left and right source
ports and the tuple

    L_A=b-a, d_A=J(b)-J(a),
    U_A=sum_(a<=n<b)(J(n)-J(a)),
    V_A=sum_(a<=n<b)(J(n)-J(a))^2.

Its clock cost at incoming current j is the positive quadratic form

    Q_A(j)=L_A j^2+2U_A j+V_A,
    G_A=[[L_A,U_A],[U_A,V_A]] >= 0.

The last jump, from b-1 to b, belongs to d_A even though it is not yet
charged in Q_A. For adjacent A,B, source-preserving concatenation gives

    L_AB=L_A+L_B,                     d_AB=d_A+d_B,
    U_AB=U_A+U_B+L_B d_A,
    V_AB=V_A+V_B+2d_A U_B+L_B d_A^2.              (5a)

Equivalently `G_AB=G_A+R_d^T G_B R_d`, with
`R_d=[[1,d_A],[0,1]]`. This is the integer-clock specialization of the
prefix-moment compiler in the wavelet note, not a new independent
finite-dimensional closure for Collatz. An upper positive-semidefinite
matrix certificate transports through (5a) for **every** incoming j.
Its size is fixed but its entries and source addresses are unbounded;
universal applicability to finite blocks does not prove a uniform budget.

There is a useful positive refinement. Set
`mu_A=U_A/L_A`, `D_A=V_A-U_A^2/L_A>=0`. Weighted variance decomposition gives

    D_AB=D_A+D_B
         +(L_A L_B/(L_A+L_B))*(mu_A-d_A-mu_B)^2.    (5b)

For any recursively split interval, these nonnegative merge terms sum
exactly to D of the whole interval, since point leaves have D=0. Hence

    S_X=X*(J(X)+mu_[X,2X))^2 + sum_(internal nodes) innovation. (5c)

At a balanced split the innovation is
`(|I|/4)*(mean_left J-mean_right J)^2`, the ordinary Haar energy.
Equation (5c) separates two obligations: the actual incoming/mean current
and the accumulation of oscillation across scales. If they are bounded
by `b^2 X^2` and `c X^2` with `b^2+c<1/2`, (5) certifies the target.
Both are deterministic arithmetic obligations, currently OPEN.

This is a concrete use of the Q–Pell and golden controls. The
[Q covariance/Pell map](golden_nested_trees_20261013c.md#6-incoming-q-covariance-work-a-genuine-two-step-pell-mechanism)
suggests retaining the unnormalized matrix before projecting to a
correlation. Here the exact affine shear R_d preserves the incoming port;
there is no asserted Pell invariant. Golden rewriting preserves value
but needs its ordered rewrite history to recover a source; here reordering
blocks changes d_A and the shear. The preserved predicate is the exact
cost for each incoming current, the lost information is the interior
word, and a source-labelled binary tree is its reconstruction sidecar.
The cheapest controls are (5a) at every cut and (5b) on hostile same-word,
different-clock examples; the independent audit below includes them.

## 3. Repair of the incoming stochastic theorem

The old Theorem E omitted the factor 4 from H=-2J+boundary and inferred a
fixed-source amplitude from a model expectation. Its strongest sound
fixed-index version is the following elementary statement.

Let (v_i) be a mean-zero stationary process on a specified probability
space, with covariance gamma(l) and `sum_l |gamma(l)|<infinity`. For a
**deterministic event count m**,

    E[J_m^2] = m gamma(0)+2 sum_(l=1)^(m-1) (m-l)gamma(l),
    E[(-2J_m)^2] = 4m sigma^2+o(m),
    sigma^2=gamma(0)+2 sum_(l>=1)gamma(l).                   (6)

If `gamma(0)=2/3` and `sigma^2=F gamma(0)`, the leading coefficient is
`(8F/3)m`, not `(2F/3)m`. Even the iid model uniform on{-1,0,1} gives
`E[(-2J_m)^2]=8m/3`, providing an exact factor-four control.
An O(1) remainder for (6) follows from the stronger weighted summability
`sum_l l|gamma(l)|<infinity`; absolute summability alone gives o(m).
A bounded endpoint term can add an O(sqrt(m)) cross-term, unless an
additional bound on that cross-moment is supplied.

Additional requirements before applying this to actual source heights are:

- mean zero and a probability space, rather than calling an uncentered
  empirical second moment a covariance;
- a source-clock estimate: B(N) depends on the same gaps and currents,
  so deterministic-m stationarity does not automatically hold conditional
  on m=B(N);
- control of the current inherited from earlier octaves and the covariance
  tail beyond the twelve measured lags;
- a deterministic fixed-source estimate, or a justified almost-sure
  statement about a law actually carried by that source.

Only with an additional deterministic clock `B(n)~d n` and suitable uniform
second-moment control would the heuristic expected energy coefficient be
`4 F d X^2`, about 0.806 at F=.84,d=.24, rather than 0.20. The corresponding
epsilon rms coefficient is about 0.212, rather than 0.105. Neither decimal
is a proved bound on Collatz.

Even if `E[E_(2^k)]<=C 2^(2k)` were established on a suitable probability
space, Markov and Borel–Cantelli give only
`E_(2^k)<=k^2 2^(2k)` eventually almost surely. Fourth moments are not
needed for that implication. This does not give the constant 2 bound at
every dyadic scale beyond the authenticated cutoff 2^34. It also gives no
explicit burn-in or tail budget preventing cancellation of the finite
Mellin numerator. The incoming nonzero-pole conclusion is therefore
withdrawn. The earlier conditional amplitude certificate is unchanged.

A minimal model hostile makes the quantifier issue explicit: let the
energy on every dyadic scale be `3X^2` on one outcome of probability 1/4,
and zero otherwise. Its expectation is only `3X^2/4`, yet the designated
first outcome violates `2X^2` at every scale. These are abstract
nonnegative energy variables, not a Collatz construction. They refute
the expectation-to-uniform-certificate inference on its own.

## 4. Full gap census and the size of the omitted clock

The incoming producer used a histogram restricted to gaps<64, then
normalized it as if it contained all gaps. An independent native 128-bit
descent census gives the following for break positions in `[2^22,2^27)`:

| quantity | exact value |
|---|---:|
| breaks | 31,569,545 |
| consecutive gaps | 31,569,544 |
| gaps at least 64 | 85,315 |
| maximum gap | 1,412 |
| gap residues 0,1,2 mod 3 | 6,083,211; 16,912,343; 8,573,990 |
| total length of complete interior gaps | 130,023,423 |
| length carried by gaps at least 64 | 9,000,105 |
| weighted current sum over complete interior gaps | 1,558,127,929,192,262 |
| contribution from gaps at least 64 | 113,413,856,148,396 |

Thus long gaps are only 0.2702% of events but 6.9219% of the sampled integer
length and 7.2789% of the positive weighted current energy. This is a
finite length-bias effect, not a proof of an asymptotic heavy-tail law.
The corrected full gap probabilities are approximately
`(.1926924,.5357171,.2715906)`.

The old count 31,484,229 was the number of *short gaps*, not breaks. The
lag-product calculation itself used all breaks and its twisted-residue
identity survives. Its fitted Markov helper, however, divided by the
probability that both endpoint residues were nonzero; this is a different
observable from the advertised all-break normalized covariance. The
correct denominator is the constant `sum_c pi_c W(c)^2`. The producer is
repaired and the original output retained under a historical warning;
the independent exact census here supplies the replacement counts.

The source-clock decomposition also identifies the dominant cost in the
largest measured excursion. At X=2^24,

    S_X/X^2 = 0.442190256...,
    X*(mean J)^2/X^2 = 0.414679055...,
    sum innovations/X^2 = 0.027511201....

Thus about 93.78% of the positive energy lies in the octave mean, rather
than its centered fluctuation. The exact H-energy ratio is 1.768761149....
At X=2^26 the two clock components are 0.245573212... and 0.025423223....
These are finite observations, but they sharpen the next test: a good
covariance/variance contraction still needs a bound on the transported
mean and the current entering the interval. Centering it away would
remove most of the difficult octave's energy. The saved census records
all 27 complete octave summaries, with independent Python moment checks
for the first 14.

The unconditional alternative is now specific: use the actual integer
clock and estimate the nonnegative sum S_X, or build an authenticated
upper matrix for K_X coupled to the marked current. Long-gap tails must
be paid in source length, not only in their frequency among events.
This is the same source-retention principle as the Conway denominator
and the unchanged-curve BSD terminal switch, without importing either
problem's quantitative conclusion.

## Reproduction

```sh
clang++ -O3 -std=c++17 04-computation/experiments/collatz_break_clock_20261013d.cpp -o /tmp/break-clock
/tmp/break-clock /tmp/math-amplitude-20261013c/parity26.bits 05-knowledge/results/collatz_break_clock_20261013d.json
python3 04-computation/experiments/collatz_break_clock_audit_20261013d.py --output 05-knowledge/results/collatz_break_clock_audit_20261013d.json
```

The full native census agrees with all 67,108,864 saved parity bits and
checks all twelve raw/twisted lag-product sums exactly. The independent
Python audit checks (2) on all 6,561 ternary words of length 8 on two source
windows (13,122 checks), the fixed-clock hostile, and the iid variance
factor through m=8. It also checks all 729 half-integral jump words of length 6 at every cut (3,645 matrix/variance checks), 7,290 ordered three-block compositions, the corrected incoming Python function bodies through 2^14, and 36 circulant Markov controls against their Fourier eigenvalues. The optional Numba decorator is bypassed only for this independent small audit; no compiled full-producer rerun is claimed. Normal and optimized audit outputs agree.
[Full census](collatz_break_clock_20261013d.json),
[finite controls](collatz_break_clock_audit_20261013d.json).
