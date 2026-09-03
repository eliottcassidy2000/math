# LRC(14): adaptive clocks, Brownian current, and the next sharp wall problems

Status: **SYNTHESIS OF PROVED / FINITE-EXACT RESULTS; LRC(14) OPEN.**

## Inheritance pass and live concept board

The closest proved mechanisms were THM-2066's labelled two-lift owner word,
THM-4330's anchored `2+12` minority reduction, THM-4335's exact renewal
owners, and THM-638's signed pair law.  The canonical hostiles were the
degree-twelve `420|h` minority, the divisor-complete bodies that defeat scalar
entry gates, and the tight half-turn current row

```text
(1,3,5,7,9,11,13,15,17,19,21,45).
```

The corrected near misses were fixed-bank universality, address-free pair
reuse, independent `1/7` deletion of an anchor strip, and the hope that cubic
occupation moments force a free sheet.  The least-used sidecars were the two
physical lift labels, the centered wall residue, the absolute-current
histogram, and divisor transport between clocks.

The session kept five live concepts:

1. labelled owner clocks and their inverse system;
2. prefix-envelope interval renewal and third-tooth survival;
3. the even-anchor residual wall;
4. signed half-turn current and its stochastic filtration;
5. finite occupation moments versus realizable runner geometry.

Every successful move came from retaining one coordinate that a tempting
scalar quotient had discarded.

## What became theorem

### 1. Bounded entry and adaptive completeness

THM-4347 proves, relative to the existing seam theorems and an exact dual
implementation, that every eleven-element body `H` with `max(H)<=40` makes
`2H union {a,b}` safe for every pair of distinct positive odd tails.  The new
exact layers `31,...,40` scan `2,257,174,140` bodies and leave no owner-word
survivor.

The same fixed clock bank does not extend uniformly.  The divisor-complete body

```text
H=(1,5,11,13,17,19,23,37,41,70,72)
```

escapes every clock `15,...,43`; clock `47` certifies it.  This is a hostile to
the proof bank, not an unsafe LRC row.

THM-4349 repairs the quantifier rather than enlarging the fixed bank.  For a
fixed eleven-body `H`, universal safety against all distinct odd tail pairs is
equivalent to emptiness of one owner relation.  One deliberately coarse
complete choice is

```text
N_H=28(42h+1)^2 lcm(1,...,12h),       h=max(H).
```

If every clock is nonempty, inverse-limit compatibility cannot hide a
profinite-only obstruction: the LRC(12) safe arc stabilizes both tail residues
to actual odd integers below `12h`.  Thus the owner-clock method is complete
on this seam.  The remaining problem is to force emptiness, not to justify the
certificate language.

### 2. Third teeth became wall arithmetic

THM-4348 turns farthest-reach renewal into an exact prefix-record recurrence
and gives a centered-residue successor test with strict endpoints and fixed
ties.  The important unexpected reduction is its wall shadow.  If a selected
speed is `u`, a residual speed `s` kills exactly

```text
gcd(u,s) [C_q(S)+C_q(-S)]
```

of the `2u` signed `u`-walls, with `q=u/gcd(u,s)` and
`S=s/gcd(u,s)`.  Hence it kills at most

```text
2 gcd(u,s) ceil(u/(7 gcd(u,s))).
```

One residual cannot cover off total resonance.  Two residuals can reach the
capacity boundary only in one classified quotient-three pattern.  In the
primitive `420|h` branch this closes banks of at most two residuals whenever
`7u` does not divide `h`.  The exact next wall is therefore

```text
7u | h,
```

not an undifferentiated request for “more third-tooth competition.”

### 3. The stochastic-process analogy became literal

THM-4346 proves that the half-turn current covariance of two odd speeds is

```text
[d_14(U+V)-d_14(U-V)]/(7UV),
```

after gcd reduction.  Its nonzero residue numerator has rank three and is,
after the order `1,5,3`, twice the covariance `min(i,j)` of a three-step
Brownian path.  For an arbitrary odd family the Fourier coefficients are a
Dirichlet convolution, and the energy splits into orthogonal reverse-martingale
differences indexed by the `7`-adic valuations of the speeds.

This makes the earlier Green-current language exact rather than metaphorical,
but it also identifies its limit.  A twelve-speed `9`-power family has total
current energy below the second-moment threshold.  Two abstract exchangeable
sheet laws agree on every labelled occupation tensor through degree three yet
have different free-sheet mass; their first separation is quartic.  They are
not runner-realizable, so the conclusion is precise: cubic moment data alone
does not imply freedom without an arithmetic or geometric sidecar.

THM-4345 supplies that sidecar at an even anchor.  Restriction of any
half-periodic observable along an odd dilation to the anchor-danger strip is
exactly a Euclidean quotient plus one residual-wall integral.  The induced
fourteen-state radius/phase transducer is associative.  On the canonical
hostile, even the exact positive part of the quadratic current certificate is
negative, while the sharp nonlinear envelope `integral g(|C|)` is positive.
The useful stochastic object is therefore the anchor-conditioned current
histogram, not merely its variance.

## Connection contracts

| Source | Target | Exact map | Preserved | Lost / needed sidecar |
|---|---|---|---|---|
| owner relations over clocks | actual odd tail pairs | inverse limit, then safe-arc stabilization | strict danger, parity, both lifts | practical clock size |
| greedy tooth cover | signed wall arithmetic | selected tooth to centered successor residue | endpoints, ties, reachability | multiplicity when three or more residuals overlap |
| signed current | Brownian / martingale energy | mod-14 kernel and `7`-adic Fourier projections | full quadratic energy | anchor location and nonlinear depth |
| anchor restriction | residual-wall automaton | `D=7q+r` radius/phase update | exact restricted integrals | joint wall/current state if absolute values are taken |
| degree-three sheet moments | free-sheet event | factorial-moment projection | all cubic labelled tensors | quartic/arithmetic realizability |

## Next sharp problems

1. **Resonant wall `7u|h`.** Change the distinguished selected speed or anchor
   and derive the transition law between the two wall quotients.  The target
   is a two-anchor lemma showing that a row cannot be totally resonant in both
   charts unless it enters a previously closed AP/divisor family.

2. **Three-residual wall overlap.** Replace the union-capacity bound by an
   exact intersection or additive-energy formula for three centered residue
   shadows.  The quotient-three equality pattern is the hostile control.  A
   useful theorem must measure overlap, not just sum individual capacities.

3. **Compress the complete adaptive clock.** THM-4349 proves existence but its
   `lcm(1,...,12h)` clock is intentionally enormous.  Seek a polynomial-size
   clock selected from the finite wall denominators, or prove a small-clock
   dichotomy after divisor completeness.  The max-72 fixed-bank escape must be
   included as a negative control.

4. **Force the anchor-conditioned current histogram.** Express
   `integral g(|C|)` through tail probabilities `P(|C|>=j)` and propagate them
   through the residual-wall transducer and `7`-adic martingale shells.  A
   variance-only argument is already refuted by the exact canonical hostile.

5. **Locate the quartic realizability boundary.** Determine which quartic
   occupation tensors are forced by the arithmetic current kernel and by
   wall geometry.  The abstract cubic-blind pair is a hostile; the desired
   separator must vanish on its nonrealizable directions while remaining
   computable on physical rows.

The strongest immediate route is the conjunction of problems 1 and 2: it
attacks the exact remaining wall rather than searching for a global scalar.
Problems 3 and 4 are orthogonal complete-certificate routes and should be
kept alive as cross-checks.  None of these results closes the anchored seam or
LRC(14).
