# Collatz amplitudes as boundary currents: all-prefix discrepancy and a conditional pole

**2026-10-13, renewal_shapes_oct11 continuation.**
**Status:** PROVED elementary identities and implications below; FINITE-EXACT
integer census through `2^34`; interval-VERIFIED finite coefficient bounds;
CONDITIONAL noncancellation of one Dirichlet pole; NUMERICAL amplitude estimates.
The required infinite discrepancy estimate, a full modal expansion, independent
positive ROOT coverage at every source, and universal Collatz remain **OPEN**.
Self-audited; no independent reviewer or Lean verification is claimed.

## Inheritance and the live board

The anchor is the signed mod-3 error, the niche is golden spectra with their
source ensemble specified, and the wildcard is Conway's nested recurrence.
The closest mechanism is P4 of
[the signed-resonance note](collatz_signed_resonance_pair_mass_20261013b.md):
the exact Dirichlet defect identity. The hostile is its infinitely many
characteristic zeros in `Re(s)>1`, which necessarily cancel. The corrected
near miss is treating a spectral denominator as an amplitude. The least-used
sidecar is the finite terminal forcing, already retained for clock residues
in P3, and the positive adjoint boundary method of
[THM-3588](../../01-canon/theorems/THM-3588-amm-r512-truncated-adjoint-pascal-repair-horizons.md).

The incoming four commits through `2e6d7c5940` were read before deriving.
[THM-4629](../../01-canon/theorems/THM-4629-twisted-transfer-law-step-count-characters-of-collatz-orbits-hear-the-near-returns-by-congruence-class.md)
assigns near-return frequencies to characters and computes the spectra of
finite 3-adic compressions. It does not prove the infinite regularity needed
for the integer sequence. The present work supplies a finite projection
formula and a rigorously quantified sufficient tail condition.

The five live concepts are: marked terminal phase; signed residue current;
Mellin boundary flux; an ordered tree with a decreasing index; and two real
embeddings with the source counting measure retained. The companion
[golden and nested-tree note](golden_nested_trees_20261013c.md) carries the last
two. The repository cards used are “Keep observer, recurrence class, and
finite head,” “Retain local incidence through counting,” and “Type every
analogy and every implication.” No new general meta-pattern is promoted.

## 1. Terminal forcing for every character

Use the Terras map `T(n)=n/2` or `(3n+1)/2`. On the basin of 1 let `tau`
be the strict first-hit length and `omega` its odd-step count. Put

    chi(n) = zeta^(alpha tau(n)+beta omega(n)),
    chi(n) = 0 outside that basin,
    a=zeta^alpha, b=zeta^(alpha+beta), kappa=1-ab.

Then, exactly for every positive integer,

    chi(2m)=a chi(m),
    chi(2m-1)=b chi(3m-1)+kappa 1_(m=1).             (1)

For `n>1` this is first-step decomposition. At `n=1`, the right-hand branch
without the correction is `b chi(2)=ab`, whereas `chi(1)=1`.
Thus the terminal defect vanishes iff `2 alpha+beta=0 mod m`.
It vanishes for `(-1)^tau`, but is **2** for both `(-1)^omega` and
`(-1)^(tau+omega)`. This is the two-letter return word of the terminal cycle,
not a negligible change to an amplitude calculation. Away from the root it
does not change THM-4629's homogeneous characteristic polynomial.

The [exact group-algebra audit](../../04-computation/experiments/character_boundary_audit_20261013c.py)
checks strict paths through `2^16`, all characters modulo 2 through 8 on its
declared local range, and finite rational Dirichlet identities for every
character modulo 2 through 5: **6,655,589 checks**. Its minimal hostile is
`chi(1)=1` versus the uncorrected prediction `-1` for odd-step parity.

For `Re(s)>1`, define

    F(s) = sum chi(n)n^(-s),
    E(s) = sum (1_(n=2 mod3)-1/3)chi(n)n^(-s),
    R(s) = sum_(n=2 mod3) chi(n)[(n-1/2)^(-s)-n^(-s)],
    D(s) = 1-a 2^(-s)-b (3/2)^s/3.

**Lemma B1 (PROVED).**

    D(s)F(s)=b(3/2)^s[E(s)+R(s)]+kappa.              (2)

Split the even and odd terms of `F` and use (1). The `R` series is
holomorphic in `Re(s)>0`, since its terms are bounded by a compact-dependent
constant times `n^(-Re(s)-1)`. This extends P4, including its previously
zero forcing, to all step-count characters.

## 2. A finite shell, rather than unstable forward extrapolation

Write `F_N(s)=sum_(n<=N) chi(n)n^(-s)`, and truncate `E,R` at the same N.
Let `C_N=E_N+R_N`, `L=floor((2N-1)/3)`, `J=floor(L/2)`, and `N>=2`.

**Lemma B2 (PROVED, with all cutoffs retained).**

    b(3/2)^s C_N+kappa-D(s)F_N
      = -(F_N-F_L)+a 2^(-s)(F_N-F_J).                (3)

Indeed the odd sources whose image is at most N are precisely the odd
sources at most L. Their sum is `F_L-a 2^(-s)F_J`, including the forcing
`kappa` at source 1. Rearrangement gives (3). The off-root rational checks
are important: an initial draft had the wrong sign on `D(s)F_N`; this was
caught before publication, whereas testing only at `D(s)=0` would hide it.

At any zero `s0` of D the right side is the entire finite numerator of the
residue. It only uses sources in `(J,N]`, approximately `(N/3,N]`.
For `t(n)=(-1)^tau` on ROOT and zero elsewhere, `a=b=-1`, so

    C_N(s0)=(2/3)^s0 [(F_N-F_L)+2^(-s0)(F_N-F_J)],
    r_N=-(3/2)^s0 C_N(s0)/D'(s0).                   (4)

This formula is an exact projection at a characteristic zero; it is not a
least-squares sinusoid fit. It carries both the cutoff and the source phase.
If `Re(s0)>1`, each shell tends to zero using just `|chi|<=1`. Thus the
whole numerator of (2) necessarily vanishes there. The old growing-root
hostile is now visible directly as a boundary estimate. A characteristic
root in the slow half-plane still needs a tail estimate and noncancellation.

The AMM adjoint analogy is precise at this level: apply a dual functional
that annihilates the interior recurrence, leaving a boundary value and an
explicit remaining budget. In THM-3588 the budget is positive and comes
from admissible correction boxes; here it is a signed arithmetic tail and
is not yet paid. The map is between proof methods, not between their source
populations or their numerical constants.

## 3. All-prefix discrepancy, rather than selected dyadic samples

Set

    H_r(N)=sum_(n<=N) (3 1_(n=r mod3)-1)t(n), H=H_2.

The [packed census](../../04-computation/experiments/collatz_discrepancy_20261013c.cpp)
inspects every prefix through `2^34`, retaining all three H_r, their exact
extrema and locations, squared sums, and exact maxima of `H_r(n)^2/n` on
each dyadic interval. Integer comparisons use 128-bit arithmetic.

**FINITE-EXACT:**

    |H(N)| < (3/2) sqrt(N), 1024 <= N <= 2^34.       (5)

The greatest ratio in this range is about `1.446665`, attained at
`N=28,337,297` with `H=7701`. The last octave is appreciably smaller.

| Prefix interval `(2^(k-1),2^k]` | max \|H\| | location of that max | max \|H(n)\|/sqrt(n) |
|---:|---:|---:|---:|
| k=20 | 784 | 1,042,632 | 0.767804 |
| k=25 | 7701 | 28,337,297 | 1.446665 |
| k=30 | 24467 | 1,036,800,457 | 0.800696 |
| k=32 | 51751 | 3,968,079,202 | 0.821540 |
| k=34 | 65907 | 10,298,900,696 | 0.649436 |

The last column maximizes the ratio separately from the absolute maximum.
The exact positive energy also satisfies
`integral_X^(2X) H(floor(x))^2 dx < 2 X^2` for every dyadic
`2^10<=X<=2^33`. The largest normalized energy is about1.768761 near
`X=2^24`, and the last is about0.095384. Endpoint conversion from the
saved inclusive sums is exact.
At `2^34`, `H=26869` and `S=sum t=24953018`, agreeing with the old census.
Two full runs at jump depths 12 and 16 agree on all 35 dyadic rows for all
three residues and on the packed prefix through `2^26`. Each run also checks
278,563 strict, single-step first-hit paths. The fixed even-length jumps can
pass through the terminal two-cycle without changing parity; they are not
claimed to be strict first-hit certificates themselves. Overflow and
descent caps abort rather than certify an unsupported source.

The [data](collatz_discrepancy_20261013c.json) are a complete finite-prefix
test, not a proof that any square-root or square-root-logarithm bound holds
forever. A constant multiple of sqrt(N) is also stronger than needed and
could fail at rare scales. The conditional result below permits logarithmic
growth in that multiplier.

There is a useful boundary-current identity. For
`w_n=3 1_(n=2 mod3)-1`, its prefix primitive W is `0,-1,1` according as
`n=0,1,2 mod3`. Discrete integration by parts gives

    H(N)=W(N)t(N)-sum_(n<N) W(n)[t(n+1)-t(n)].       (6)

The unsigned variation bound forgets the orientations and residue positions
of these cuts. The all-prefix experiment measures what survives their signed
cancellation. A proof must retain these signs or supply another reason they
cancel; multiplying the surviving pair mass by a square-root heuristic is
not a derivation.

## 4. An interval-certified conditional nonzero pole

Let s0 be the simple zero of the parity D in the disk of radius `10^-30`
about

    0.9982092207652568136445639331675041376460
      +240.2478301406912443253007903875287137674 i.

It corresponds to frequency `26.5036120942` per octave. The
[Arb verifier](../../04-computation/experiments/collatz_pole_certificate_20261013c.py)
certifies that this disk contains exactly one zero using Rouche's theorem,
with `|D'|>0.1449` throughout the enclosing rectangle. Its center residual
is less than `2.8e-46`; the second derivative bound pays the disk remainder.

**Theorem B3 (CONDITIONAL; the implication and finite interval calculation
are verified).** Suppose, for every real `x>=2^34`,

    |H(x)| <= (1/4) sqrt(x) log_2(x).                (7)

Then F continues meromorphically to `Re(s)>1/2`, and has a **nonzero simple
pole** at s0. If r is its residue, the interval computation proves

    0.1193 < 2 |r| 2^(29(Re(s0)-1)) < 0.5908.       (8)

Under the stronger hypothesis `|H(x)|<=1.5 sqrt(x)` for every `x>=2^34`,
the corresponding interval sharpens to

    0.3165 < 2 |r| 2^(29(Re(s0)-1)) < 0.3936.       (9)

There is a second, positive-budget sufficient condition: for every dyadic
`X>=2^34`, assume

    integral_X^(2X) H(floor(x))^2 dx <= 2 X^2.       (9a)

The same interval verifier then bounds the normalized residue between
0.3215 and 0.3886. This hypothesis allows
isolated larger pointwise errors; it asks for a bound on each whole block.
It remains OPEN, just as (7) does.

These are bounds on a Dirichlet residue in the normalization used by the
spectral model. Turning that residue into a full expansion of the integer
sequence still needs inversion and growth estimates. Neither (7) nor its
stronger alternative is proved. There is no unconditional pole claim here.

**Proof and the precise debt.** Partial summation gives, for `sigma>1/2`,

    |E(s)-E_N(s)| <= |H(N)| N^(-sigma)/3
       + |s|/3 integral_N^infinity |H(x)| x^(-sigma-1) dx. (10)

Under (7), writing `eta=sigma-1/2`, the integral term is at most

    |s|/(12 log 2) N^(-eta)[log(N)/eta+1/eta^2].      (11)

Under (9a), Cauchy--Schwarz and a geometric sum instead bound the
integral in (10) by

    N^(-eta) sqrt(2(1-2^(-2sigma-1))/(2sigma+1))
                  / (1-2^(-eta)).                 (11a)

This proves analytic continuation via the Mellin integral of the step
function H. It does not by itself assert convergence of the original
Dirichlet series throughout that half-plane. The integral agrees with
that series where it originally converges and therefore supplies the
needed continuation and the same finite-cutoff tail formula.

An unconditional bound on the other remainder is

    |R(s)-R_N(s)| <= |s|/(2 sigma)
       (1-1/(2N))^(-sigma-1) N^(-sigma).            (12)

The bounds hold uniformly on compact subsets of the stated half-plane,
which gives holomorphic E there. Apply (2), (4), and (10)--(12).

For the finite shell, the census stores exact signed zeroth and first
moments on intervals of length `2^17`. For a block `[l,l+B)`, center
`c=l+(B-1)/2`, use

    sum t(n)n^(-s) = c^(-s)[sum t(n)-(s/c)sum(n-c)t(n)] + error.

The error is at most

    |s(s+1)|/2 * B(B^2-1)/12 * l^(-sigma-2).

Partial starting blocks are bounded separately, and the endpoint `n=2^34`
is included with its known positive sign. All powers, root enclosures,
quadrature errors and tail arithmetic use 160-bit Arb balls. The finite
numerator's quadrature disk error is less than `8.29e-6`. The two-complex-
coordinate rectangle deliberately overestimates that disk.

The [1.1 MB exact moment file](collatz_boundary_moments_20261013c.npz) has
SHA256 `54ca9bef99c8774bda6149b61ad45147a0256d4e8abecf5aaac4e04ce3d3d323`.
Its total agrees with the full census; 256 blocks are independently checked
from raw bits. The [certificate output](collatz_pole_certificate_20261013c.json)
records every final interval. The sign of the finite residue cannot be
cancelled by the allowed future tail, proving (8). No probabilistic model is
used in this conditional implication.

**Consequence without a modal expansion.** Under (7), `S(X)` cannot be
`O(X^beta)` for any `beta<Re(s0)`: such a bound would make F holomorphic at
s0 by partial summation. Thus the signed mean cannot have a power-law bound
`O(X^(-delta))` with `delta>0.0017907792347432`. This is a conditional lower
obstruction on the rate of cancellation. Since the ROOT count B(X) satisfies
`|S(X)|<=B(X)`, it also gives `limsup log(B(X))/log(X)>=Re(s0)` under the same
hypothesis. B1--B3 use neither THM-4581 nor its sketch-level quantitative
rate. This is an independent conditional population bound; it does not
certify every source, a positive natural density, or a uniform lower bound
at every height.

## 5. What the amplitude calculation predicts numerically

The direct finite sum `E_N+R_N` and the independently evaluated shell (4)
agree to between `1e-14` and `3e-13` in ordinary floating arithmetic at the
five tested roots. This is a numerical control distinct from the Arb proof.
For the Mercator root, no fitted coefficient is used:

| data through N | `2|r_N|2^(29(sigma-1))` | phase of r_N |
|---:|---:|---:|
| 1024 | 0.279900 | 1.572341 |
| 2^14 | 0.321482 | 1.366818 |
| 2^18 | 0.357228 | 1.471835 |
| 2^22 | 0.349774 | 1.460455 |
| 2^26 | 0.354674 | 1.454037 |

This agrees with the independently fitted roughly 0.3545 line near `2^29`.
It also qualifies the earlier phrase “the amplitude is fixed by 2^14”:
a visible line is already present then, but simply truncating the forcing
at that point does not yet give its final amplitude to high accuracy.
At `2^26`, the same formula gives 0.011045, 0.064722, 0.103267 and 0.020287
for the bases 2.5, 14.5, 79.5 and 332.5, respectively. In particular, a
nearer resonance need not have a larger residue. The last weak line is not
certified nonzero by the tail budgets used here.

## 6. Reroutes for the actual open estimate

1. **A weaker positive target:** dyadic energy instead of every prefix.
   If `sum_(X<n<=2X) H(n)^2 = O(X^2 (log X)^A)` for all dyadic X, then
   Cauchy--Schwarz on each block makes E holomorphic in `Re(s)>1/2`.
   Indeed the block contribution to `integral |H(x)|x^(-sigma-1)dx`
   is `O(X^(1/2-sigma)(log X)^(A/2))`. The saved exact squared sums
   supply the consequence object for this target. To certify a residue,
   explicit constants and the future energy budget are still required.
2. **Keep a signed boundary controller:** use (6), retaining W, the phase,
   the actual pair endpoint, and its source address. Absorption removes a
   cut exactly. An absolute-value estimate before cancellation returns to
   the much larger old pair-mass bound. An all-depth bound for this signed
   operator remains open.
3. **A source-aware compression:** store interval endpoints, its terminal
   character defect, signed prefix extrema, and its first two Mellin
   moments. Concatenation adds moments with an explicit translation and
   updates prefix extrema; forgetting either the endpoint or the root phase
   destroys (3). Conway's copying tree is a control in which the required
   rank is already the index. A native decreasing rank for every Collatz
   source has not been constructed.

The common gain is a concrete boundary budget. The remaining barrier has
not been relabelled as a proof: (7), or an effective substitute, is exactly
the new outstanding arithmetic obligation.

## 7. Integration of the incoming fast-mode and covariance work

The final fetch brought `bba55f5d6f` and `be97fc0dad`. The former's
[fast-mode discussion](meta_fibonacci_one_sided_walks_20261013.md)
samples 64 prefixes per octave and fits a roughly -0.452 exponent for its
normalized mod-3 discrepancy over the stated finite range. That is
consistent with a square-root heuristic; it does not prove a bound at
every prefix or future scale. Our complete census through 2^34 retains
the intervening extrema and exact squared sums. These are complementary
observations, with different strength and no change to the infinite
hypothesis in section 4.

The incoming Q covariance model prompted a new exact calculation:
the correlation update `(r-1)/(3-r)` has a two-generation Pell matrix
and an explicit variance law; see
[golden/nested section 6](golden_nested_trees_20261013c.md#6-incoming-q-covariance-work-a-genuine-two-step-pell-mechanism).
This suggests a precise helper question for the positive energy target:
can the signed-cut refinement preserve a finite covariance state while
charging all omitted, source-conditioned correlations to a summable
remainder? The Q calculation works because the model specifies exactly
which correlations survive. Its stationarity and sampling assumptions
have not been proved for Q, much less for Collatz. A numerical covariance
closure would not certify the future energy bound used here.

## Reproduction

All commands run from the repository worktree; scratch files stay outside it.

```sh
clang++ -O3 -std=c++17 04-computation/experiments/collatz_discrepancy_20261013c.cpp -o /tmp/discrepancy
/tmp/discrepancy 34 12 /tmp/d12.json /tmp/parity26
/tmp/discrepancy 34 16 /tmp/d16.json /tmp/parity26-control /tmp/moments34.bin
python3 04-computation/experiments/collatz_amplitude_20261013c.py --bits /tmp/parity26.bits --output /tmp/amplitudes.json
python3 04-computation/experiments/collatz_pole_certificate_20261013c.py --raw-moments /tmp/moments34.bin --moments /tmp/moments34.npz --bits /tmp/parity26.bits --census /tmp/d12.json --output /tmp/pole.json
python3 04-computation/experiments/collatz_pole_certificate_20261013c.py --moments 05-knowledge/results/collatz_boundary_moments_20261013c.npz --census 05-knowledge/results/collatz_discrepancy_20261013c.json --output /tmp/pole-replay.json
python3 04-computation/experiments/character_boundary_audit_20261013c.py --output /tmp/character.json
```

The first calculation needs about 2 GiB of packed source bits. Requirements:
C++17; Python numpy, scipy, mpmath, python-flint for the analyses. Censuses
abort on unsupported overflow/descent, and their source universe is finite.

Audit manifest: [research_audit_20261013c.json](research_audit_20261013c.json)
records exact replay comparisons, normal/optimized checks, moment controls,
source/data hashes, and the finite energy extrema.
