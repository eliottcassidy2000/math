# Signed Collatz resonance, exact pair mass, and the qualitative mean

2026-10-13. **PROVED, scoped:** the bounded-antifunction implication,
the finite discrepancy inequalities, the signed model's characteristic
equation, and the forward/backward pair-chain identity. **PROVED from an
inherited input:** the ROOT-supported parity mean is zero using the
almost-sure part of THM-4581, not its sketch-level rate. **FINITE-EXACT:**
the census through `2^34` and pair masses through depth48. **NUMERICAL:**
all fitted frequencies, decays, spectral amplitudes and effective exponents.
Self-audited; no independent peer audit or formalization of this extension.
Universal Collatz and a rigorous mode expansion remain **OPEN**.

## Inheritance and concept board

Closest mechanisms: [THM-4581, affine pair-chain coalescence](../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md)
and [HYP-9230, the class-density resonance model](../hypotheses/HYP-9230-musical-resonance-spectrum-of-collatz-class-densities.md).
The requested targets are [HYP-9258](../hypotheses/HYP-9258-collatz-coalescence-rate-convergent-denominators-and-parity-equidistribution-of-the-total-stopping-time.md).
The hostile is a particular integer on an exceptional infinite branch;
Haar coverage alone does not certify that integer. The corrected near miss
is an observed frequency being identified with a drift approximation, and
a first-descent stopping rule being treated as exactly periodic without
its exceptional-height term. The least-used sidecar is the sign character
of the time shift, together with full branch multiplicities.

Board: deterministic prefix trees; signed scaling; exact absorbing mass;
residue conductor; Pell norm; marked ramification. The last three are
developed in the [companion](tree_resonance_composite_pell_20261013b.md).
The assigned Collatz target is the anchor, composite residue colourings the
niche, and Pell/genus the wildcard. LRC(14) remains open and is not advanced.

## 1. A qualitative result does not need the resonance expansion

Use the **Terras map** `T(n)=n/2` if even, `(3n+1)/2` if odd. Define

    t(n)=(-1)^tau(n) if n eventually reaches1, and t(n)=0 otherwise.

This defines a bounded function on every positive integer without assuming
Collatz. It satisfies `t(Tn)=-t(n)`, including at1 because `1 <-> 2` is an
even cycle. Let

    S(N)=sum_(1<=n<=N) t(n),     A(N)=S(N)/N,
    V(N)=sum_(1<=n<N) |t(n+1)-t(n)|,     delta(N)=V(N)/N.

**P1 (PROVED implication).** Any bounded real function `|t|<=1` satisfying
`t(Tn)=-t(n)` and `V(N)=o(N)` has `A(N)->0`.

Proof with explicit inequalities: set
`E_M=M^-1 sum t(2m)=-A(M)`,
`O_M=M^-1 sum t(2m-1)=-M^-1 sum t(3m-1)`.
Pairing neighbours gives

    |A(2M)+A(M)| <= delta(2M).                         (1)

Grouping three consecutive integers gives

    |M^-1 sum_(m=1)^M t(3m-1)-A(3M)| <= delta(3M).

Consequently

    |A(3M)-A(M)| <= delta(3M)+2delta(2M).              (2)

Fix positive integers `a,b`, with `a` odd. Iteration yields

    |A(2^a M)+A(3^b M)| <= B_(a,b)(M),

    B_(a,b)(M) = sum_(i=0)^(a-1) delta(2^(i+1) M)
                +sum_(j=0)^(b-1)
                   [delta(3^(j+1) M)+2delta(2*3^j M)]. (3)

For any bounded sequence, comparison of its two prefix averages gives

    |A(U)-A(W)| <= 2 |U-W|/max(U,W).

Thus

    |A(3^b M)| <= |2^a-3^b|/max(2^a,3^b)+B_(a,b)(M)/2. (4)

There are arbitrarily close ratios `2^a/3^b` with `a` odd: the sequence
`(2j+1)log(2)/log(3)` is dense modulo1. For any fixed such pair, `B->0`.
Every large N differs by at most `3^b` from a multiple of `3^b`, which
changes its prefix average by `O(3^b/N)`. First let N tend to infinity for
a fixed pair, then improve the pair. Equation(4) proves P1. No rate or
exchange of these two limits was used.

**P2 (PROVED from THM-4581's almost-sure input).** For the ROOT-supported
function above, `V(N)=o(N)` and hence

    #{n<=N: tau(n)<infinity, tau(n) even}
     -#{n<=N: tau(n)<infinity, tau(n) odd} = o(N).      (5)

Indeed let `q_K` be the exact pair-chain nonabsorption probability from
`(0,1)` by depthK. Every absorbed residue modulo `2^K` gives equal-time
affine coalescence for its entire integer cylinder. Antifunction signs
therefore agree on adjacent sources in that cylinder, even if their common
orbit never reaches1. Counting the remaining residues gives

    delta(N) <= 2q_K + 2^(K+1)/N.                     (6)

THM-4581's statement3 supplies `q_K->0`. Fix K, let N tend to infinity,
then let K grow. This proves the required variation limit. This note does
not independently re-prove that input; the implication P1 and the finite
bound(6) remain usable separately if it is re-audited.

If Collatz holds, (5) says the two total-stopping-time parity classes both
have ordinary density1/2. Without it, (5) is an **absolute density
difference**, not a claim that the ROOT basin itself has density1 or even
that its relative even/odd proportions have limits. Its logarithmic signed
average also tends to zero by summation by parts. The difference between
averages at fixed multiplicative endpoints shows the same conclusion on
every fixed-width multiplicative window, without a decay rate.

This separates a newly closed qualitative implication from the still-open
quantitative envelope and the much stronger marked-source entry problem.

## 2. Why neither candidate frequency is the exact model answer

There is an exact finite identity, before any smoothing:

    S(2M)+S(M)+sum_(m=1)^M t(3m-1)=0,

so

    |A(2M)+A(M)/2+A(3M)/2| <= delta(3M)/2.            (7)

Suppressing that error gives the **model**

    h(x)=-h(x-1)/2-h(x+alpha)/2,
    alpha=log_2(3)-1,      x=log_2(n).

Its exponential solutions have

    exp(-lambda)+exp(alpha*lambda)=-2.                (8)

For a half-integer `f=j+1/2`, put `lambda=2 pi i f+w` and
`theta=f log_2(3)-round(f log_2(3))`. Equation(8) becomes

    exp(-w)+exp(alpha*w+2 pi i theta)=2.              (9)

This is exactly HYP-9230's characteristic equation on the **half-integer
frequency lattice**. Its root continued from `w(0)=0` predicts frequency
`f+Im(w)/(2pi)` and decay per octave `-Re(w)`.

| Base f | Predicted frequency | Predicted decay |
|---|---:|---:|
| 2.5 | 2.454923470010307 | 0.212159671655779 |
| 14.5 | 14.468865286662341 | 0.094648348972273 |
| 26.5 | 26.503612094231855 | 0.001241273577568 |
| 79.5 | 79.51046644889884 | 0.010437086290553 |
| 173.5 | 173.48092397028014 | 0.034845518245450 |
| 332.5 | 332.50007587198417 | 0.000000547564541 |

The **linear** expansion `w=2pi i theta/(1-alpha)+O(theta^2)` at base5/2
gives exactly `1/(2-log_2(3))=2.409420839653209`, the half drift.
Thus half drift and5/2 are the first approximation and unshifted lattice
label of the same candidate mode; a superposition of two independent
modes is unnecessary. The nonlinear phase correction is measurable.

The parity rule on near-relations is informative: long-lived model modes
come from `3^q approximately 2^p` with **q odd, p even**, giving base `q/2`.
Examples are `(q,p)=(5,8),(53,84),(665,1054)`. In contrast P1's cancellation
uses near-identity dilations with the exponent of2 **odd**, e.g. `(a,b)=(19,12)`
or `(65,41)`. These are complementary tests of the same character
`chi(2)=-1, chi(3)=+1`, not a coincidence of small numbers.

The error in(7) is only controlled as `o(1)` by P2. That does **not** justify
an expansion into particular decaying modes, let alone its nonzero
coefficients. A rate-sensitive bound for the signed residue discrepancy
in(7) is the remaining analytic target.

## 3. Census through 17,179,869,184 and the aliasing control

The [C++ census](../../04-computation/experiments/collatz_parity_stream_20261013b.cpp)
stores one parity bit per source, about2GiB. Induction on n uses either
the exact halving rule or an even-length Terras jump until it reaches a
smaller positive integer. Passing through the1↔2 cycle preserves parity;
these jumps are not presented as strict first-hit ROOT words.

Two full runs, at jump depths12 and16, agree on **all20,480 bins** at
1,024 bins/octave over `[2^14,2^34)`. The second uses4,096 bins/octave;
aggregating each four gives the first run exactly, including the separated
even-source sums. Both independently check278,563 strict first-hit paths
(every source through2^18,16,384 deterministic pseudorandom samples, all
powers of2). Each bin obeys the even-half identity exactly. The inclusive
prefix sum through2^34 is24,953,018, so its signed mean is about0.00145246.
No convergence beyond this finite universe is inferred.

The [analysis script](../../04-computation/experiments/collatz_parity_modes_20261013b.py)
preselects half-integer bases below512 by `|theta|<=0.06`, fits their
continuous bin integrals, and projects out the other modes before fitting
the low mode's frequency and decay. All inference is **NUMERICAL**; no
independent-noise confidence interval is assigned to deterministic bins.

| Data interval in log2(n) | Fitted low frequency | Fitted decay |
|---|---:|---:|
| 16..26,1024 bins/octave | 2.45548185 | 0.20730294 |
| 20..30,1024 bins/octave | 2.45414787 | 0.21266390 |
| 24..34,1024 bins/octave | 2.45485663 | 0.21247561 |
| 24..34,4096 bins/octave | 2.45484447 | 0.21281597 |

Over24..34 at1024 bins/octave, the residual RMS after fitting the other
model modes is0.007975 for half drift,0.007966 for5/2, and0.006001 for
the signed-root frequency; each candidate is allowed its own fitted decay.
The signed model strongly improves this declared comparison. Training
all fixed modes on20..28 predicts28..34 with RMS0.01277, against raw RMS
about0.2996; this is a holdout in scale, not an independent theorem.

The strongest resolved fine-bin mode is near26.5036, with amplitude about
0.357 at log2(n)=24, much larger than the low mode's0.0319. This makes the
old phrase "the dominant mode is2.455" a statement about the old scan's
bandwidth, not the full data.

The4096-bin run exposes another loss: strong peaks near532.5 and585.5
alias near491.5 and438.5 at1024 bins/octave. Do not identify every coarse
high-frequency fitted coefficient with its labelled root. With only bases
below512 retained, the4096-bin residual is about0.157 over24..34, rather
than0.006; much of the finer signal lies outside that model's bandwidth.
The low mode is stable under this hostile check. The qualitative zero-mean
proof does not assume any modal truncation is complete.

## 4. Exact unmerged mass through depth48

The [pair-chain checker](../../04-computation/experiments/collatz_pair_mass_20261013b.py)
uses the unconditioned start `(k,e)=(0,1)` for `(n,n+1)`. Swapping the two
orbits identifies `(k,e)` with `(-k,-3^(-k)e)`. In integer numerator
coordinates this permits `k>=0`, with `e>=0` at k0. It preserves absorption
and swaps the two fair branch labels; it does not discard their weights.

Let `F_t(s)` count t-bit prefixes ending at state s, including the absorbed
state's **two** self edges. Let `B_r(s)` count r-bit suffixes from s ending
absorbed. Exact finite composition gives

    2^K q_K = 2^K - sum_s F_t(s) B_(K-t)(s).          (10)

The backwards recurrence enumerates every predecessor **with edge
multiplicity**. A split at t31,r17 handles depth48 with1,926,627 forward
states and1,198,189 backward states, instead of enumerating2^48 words.

| K | Unmerged numerator, denominator2^K | q_K | sqrt(K) q_K |
|---:|---:|---:|---:|
| 20 | 626933 | 0.5978899002 | 2.67384 |
| 26 | 38861357 | 0.5790793449 | 2.95274 |
| 32 | 2423829673 | 0.5643418229 | 3.19240 |
| 40 | 602924248450 | 0.5483564095 | 3.46811 |
| 48 | 150619085192045 | 0.5351064843 | 3.70733 |

Controls:7,290 edge/split/mass checks, plus every one of2^16 residues
compared against independent large integer orbit pairs, retaining odd-step
counts to distinguish an affine identity from an isolated coincidence.
Normal and optimized Python outputs agree. The depth20 value reproduces
THM-4581; depth26 reproduces the incoming census.

Local fitted exponents `-Delta log(q)/Delta log(K)` are0.1218 on20..26,
0.1242 on26..32,0.1288 on32..40,0.1342 on40..48. These finite transients
do **not** determine the asymptotic exponent; the existing large-time
Monte Carlo and its heuristic constant concern a much later range.

The first-hit numerators are `h_K=2 a_(K-1)-a_K`, where `a_K=2^K q_K`.
They are the weights of a deterministic prefix-free accepting tree:
`sum_(k<=K) h_k/2^k=1-q_K`. Its quotient states are a compression of
proof histories, not a substitute for their multiplicities or the source
address. This is precisely why a measure-one result cannot certify every
fixed source.

## Reproduction and next proof targets

Saved exact data: [1024-bin census](collatz_parity_bins_20261013b.json),
[4096-bin census](collatz_parity_bins4096_20261013b.json),
[pair masses](collatz_pair_mass_20261013b.json).
Numerical analyses: [1024](collatz_parity_analysis_20261013b.json),
[4096](collatz_parity_analysis4096_20261013b.json).

```sh
python3 04-computation/experiments/collatz_parity_modes_20261013b.py --edges /tmp/parity-edges.txt
clang++ -O3 -std=c++17 04-computation/experiments/collatz_parity_stream_20261013b.cpp -o /tmp/parity-stream
/tmp/parity-stream 34 /tmp/parity-edges.txt /tmp/parity.json 12
OPENBLAS_NUM_THREADS=1 python3 04-computation/experiments/collatz_parity_modes_20261013b.py --input /tmp/parity.json --output /tmp/parity-analysis.json
python3 04-computation/experiments/collatz_pair_mass_20261013b.py --output /tmp/pair48.json
```

For the second full run use `--bits 4096`, and
jump depth16. Endpoints are explicitly stored integers;90- and120-digit
Decimal computations agree, but a separate interval proof of each floor
is not claimed or needed to define the census universe.

Next useful obligations: bound the **signed** residue error in(7), rather
than bounding it by all surviving pair mass; prove renewal estimates with
the actual `(k,e)` return law; and attach an independently decreasing
height to the marked integer source. None of the present averaged bounds
supplies that last coordinate.
