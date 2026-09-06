# Independent norm-five projection profiles and compact global bases

**Status: INDEPENDENT ANALYTIC RECONSTRUCTION + FINITE-EXACT;
INDEPENDENT AUDIT PASS.** See the [profile referee](overnight5_20260906_lrc_norm5_audit.md)
and [generic reduction referee](overnight5_20260906_lrc_nonadditive_reduction_audit.md). The resulting sharp norm-five conclusions are already proved in
incoming
[THM-4441 / signed (1,2,2) sharp ray closure](../../01-canon/theorems/THM-4441-lrc14-signed-122-sharp-ray-closure.md),
read directly at `origin/main` after incoming `058a8ded9`:

```text
min_i E_i(w)<=46/665, equality iff w=(2,19,20),
mu(F_w)<=51/770,      equality iff w=(1,11,20).          (1)
```

This is independent concurrent corroboration, not a new closure or theorem
namespace. The independent route records all three continuum projections,
their selectors, a saturated-plane proof of the physical integral, and
smaller complete bases for the two **global** constants: 174 rows through
height 50 for the network maximum and 131 through height 45 for the physical
maximum. It does not replace THM-4441's larger bases for its finer per-cone
constants.

The final section separately checks an exact c=535 nonadditive count cutoff
from the already audited 747-pattern box. Incoming proved
[THM-4437 / all-parity reduction to three low circuits](../../01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md)
now supplies a stronger generic all-height result, so that cutoff is an
independent audit route, not a missing canonical premise.

## 1. Full carrier scope and the four sorted sectors

Assume `w=(a,b,c)` consists of distinct positive sorted ternary-unit integers,
is primitive, and admits a signed primitive relation v with coefficient
magnitudes `(1,2,2)`. No oddness assumption is made. In raw `y=3x` coordinates,
the physical error cube has radius `r=3/14`. The exact owner residue rule
for a full-support relation allows only `delta=0 mod3`, while

```text
delta=v dot n is integral, |delta|<r||v||_1=15/14.
```

Hence delta is zero. The primitive line identity `v cross C=delta w`
therefore confines the **complete** carrier set to `C=k v`. All coefficients
are units modulo three, so the owner-live word is exactly `3` not dividing k.
The strict live cutoff must retain all three roofs:

```text
Lambda(w)={k v: k in Z, 0<|k|<B, 3 does not divide k},
B=min_i r(w_j+w_k)/|v_i|.                              (2)
```

The parity-free physical address and branch normalization were independently
derived in the [earlier additive audit](overnight3_20260906_lrc_additive_audit.md)
and checked in the [all-parity native comparison](overnight4_20260906_lrc_parityfree_native.md).
No odd-only ceiling is used here.

After setting `t=a/c` and normalizing c to one, positivity and sorting give
the following exhaustive sectors. Let `kappa=B/(rc)`.

| Sector | Speed identity | v | b/c | t range | kappa |
|---|---|---|---|---|---|
| I | `c=2(a+b)` | `(2,2,-1)` | `1/2-t` | `(0,1/4)` | `1/2` |
| II | `2c=a+2b` | `(1,2,-2)` | `1-t/2` | `(0,2/3)` | `(2+t)/4` |
| III | `2c=2a+b` | `(2,1,-2)` | `2-2t` | `(1/2,2/3)` | `(2-t)/2` |
| IV | `2b=2a+c` | `(2,-2,1)` | `t+1/2` | `(0,1/2)` | `(1+t)/2` |

To see completeness, isolate a signed coefficient on one side of its
positive speed relation. A coefficient of magnitude one on c yields I or
IV. A coefficient of magnitude two on c yields II or III. The other cases
contradict positivity or sorting. The possible II/IV overlap forces a
multiple of `(2,5,6)` and violates the ternary-unit filter; the other
intersections are improper. The complete eligible finite bases therefore
have a unique sector per row.

In THM-4441's notation these are F1, F3, F4, F2. Its F4 parameter is b/c,
which equals `2-2t` in our sector III. Thus its crossing at `b/c=7/8` is
exactly our `t=9/16`; the two formula systems agree.

## 2. Exact projection integrals with the complete cutoff

Write `z=w/c`. For coordinate i, let `P_i=z_j z_k`, `S_i=z_j+z_k`, and
`m_i=|v_i|`. The actual raw projection profile, restricted to the full live
line, is

```text
h_i(u)=(r/c) min(2,(S_i-m_i u/(rc))/P_i),  0<=u<rc kappa,
h_i(u)=0,                                              u>=rc kappa.
```

Its value at zero is `q=2r/c=3/(7c)`. Some profiles have a downward jump at
the common cutoff because another roof, rather than this one, vanishes
there. Removing that cutoff would count nonexistent carriers.

Define

```text
J_i(t)=integral_0^kappa min(2,(S_i-m_i s)/P_i) ds.
```

Then `integral_0^infinity h_i=r^2 J_i`. Since the owner word retains two
thirds of positive and negative integers, the continuum projection is

```text
A_i(t)=(4/3)r^2 J_i(t)=(3/49)J_i(t).                  (3)
```

Every roof starts above the common cap. If
`g_i=m_i kappa-S_i+2P_i`, exact plateau/trapezoid integration gives

```text
J_i=2 kappa-(max(g_i,0))^2/(2m_i P_i).                 (4)
```

The resulting **three** integrals are:

| Sector | J_a | J_b | J_c |
|---|---|---|---|
| I | `(2t+7)/8` | `1-t/4` | `2t^2-t+1` |
| II | `(9t+14)/16` | `1-t/16` | `(t^2-t+2)/2` |
| III | `(t+7)/8` | `(16-9t)/8` | `2t^2-3t+2` |
| IV | `(16t+7)/(8(2t+1))` | `1` | `1+t-t(4t-1)_+^2/[4(2t+1)]` |

Here `x_+=max(x,0)`. The positive part in IV is essential: J_c remains
`1+t` until t=1/4. The other cap-switch numerators are positive on their
whole parameter intervals. The standalone verifier checks the formulas
symbolically and checks every required rational inequality over its full
open interval using exact real-set arithmetic.

## 3. The exact selector and its uniform continuum ceiling

In I, J_b is never the minimum. The other difference factors as

```text
J_c-J_a=(2t-1)(8t-1)/8.
```

Thus J_a wins until t=1/8, and J_c thereafter. The former increases, while
the latter decreases on this range, so the maximum of their minimum is
`29/32`, attained at the crossover.

In II, the exact selector again changes from J_a to J_c at t=1/8:

```text
J_c-J_a=(t-2)(8t-1)/16,
J_b-J_c=t(7-8t)/16>0.
```

Before the crossover J_b also lies above J_a. Afterward J_c is convex, so its
maximum on the remaining closed parameter interval occurs at an endpoint;
its value at t=2/3 is `8/9<121/128`. Therefore the sector maximum is
`121/128`, at t=1/8.

In III, J_b is greater than one throughout. The relevant difference is

```text
J_c-J_a=(t-1)(16t-9)/8.
```

The increasing J_a and decreasing J_c cross at t=9/16, again with value
`121/128`. In IV, J_a is below one, J_b equals one, and J_c is above one.
Its selected integral is J_a, which increases to the boundary supremum
`15/16` as t approaches 1/2.

Consequently

```text
min_i A_i(t) <= (3/49)(121/128)=363/6272.              (5)
```

The two real-ratio maximizers correspond to primitive triples `(2,15,16)`
and `(9,14,16)`, which have a speed divisible by three. They still give a
valid closed-domain bound. No equality statement about the discrete network
is inferred from a continuum maximizer.

The selector in (5) chooses **one projection for the given speed triple**.
It does not exchange a minimum over coordinates with the sum over carriers.
That distinction is necessary at `(10,11,16)`:

```text
E=(17/176,9/140,3/55),
mu(F_w)=331/6160,  min E-mu(F_w)=1/1232.
```

At k=1 the second roof wins, and at k=2 the third wins. There is no fixed
physical roof on this entire norm-five line.

## 4. Independent physical integral by a saturated error plane

The constant physical integral can be obtained without listing its piecewise
envelope. Choose an integer Bezout vector q with `q dot w=1`, and set
`n_1=v cross q`. Then `w cross n_1=v`. The linear map

```text
(y,u) -> e=u n_1-y w
```

parametrizes the plane `v dot e=0` with area Jacobian `||v||`. Its fibre
length in y is exactly the real physical carrier length `L_w(u v)`. Thus

```text
integral_R L_w(u v)du
 =area({e in [-r,r]^3:v dot e=0})/||v||.               (6)
```

There is no longitudinal lattice index: w and v are primitive and the
displayed integer lift realizes the step v exactly. This is the same
physical raw-y normalization used in (2), before applying the owner word.

Eliminate the coefficient-one error coordinate. In the other two coordinates
the section is a square `[-r,r]^2` cut by `|e_j+-e_k|<=r/2`. The coordinate
area factor cancels the `||v||` in (6). Two omitted corner triangles each
have legs `3r/2`, giving

```text
integral_R L_w(u v)du
 =4r^2-2*(1/2)*(3r/2)^2
 =(7/4)r^2=9/112.                                    (7)
```

The physical continuum is therefore `(2/3)(9/112)=3/56`. Root supplied this
saturated-plane route independently; the fifth-round referee also integrates
the full lower envelope independently in all four sectors. It agrees with
THM-4441's existing piecewise physical identity. The companion independently
checks the complete envelope at rational controls, including both sides of
the genuine sector-II owner switch.

## 5. Strict sampling error and compact complete bases

For T>0, let `R_<(T)=#{k:1<=k<T, 3` does not divide `k}`. Writing
`M=ceil(T)-1` gives

```text
R_<(T)=M-floor(M/3)<2T/3+2/3.
```

Every h_i and the full physical profile is nonnegative and nonincreasing,
including its downward cutoff jump. Layer-cake integration gives

```text
2 sum_(k>=1,3 not dividing k) h_i(k)
 <(4/3) integral_0^infinity h_i(u)du +(4/3)h_i(0).
```

Hence the exact tail estimates are

```text
min_i E_i <363/6272+4/(7c),
mu(F_w)  <3/56+4/(7c).                               (8)
```

The requested strict `11/140` bound already follows from the first line at
c>=28; its real cutoff is `17920/649`. The complete norm-five H27 base has
48 rows, all below the target. The stronger global sharp bounds (1) need
only these cutoffs:

```text
4/[7(46/665-363/6272)]=340480/6731<51,
4/[7(51/770-3/56)]=1760/39<46.                        (9)
```

An independent speed-universe generator constructs every primitive sorted
ternary-unit norm-five triple through height 50: exactly 174 rows, with
sector I/II/III/IV counts `25/67/32/50`. Of those, 131 have c<=45. The exact
complete multiplier sum reconstructs all three E values and the physical
mass, then compares every row to the frozen H63 table previously audited by
root's independent native six-sheet engine. Thus the literal data are an
explicit inherited finite premise, not a newly claimed large replay.

Both finite equality loci are singletons:

```text
(2,19,20): E=(48/665,11/140,46/665), mu=173/2660;
(1,11,20): E=(51/770,3/35,3/35),    mu=51/770.
```

Together with the strict tails in (9), these establish exactly the global
sharp conclusions of THM-4441 by compact bases. This does not claim the
stronger per-sector maxima in that theorem from an insufficient H50 base.

## 6. Separate general nonadditive cutoff certificate

The physical/network count gate for `11/140` is `N<=11c/60`. Recompute the
already audited parity-free coefficient universe `0<=p_1<=p_2<=p_3<=18`,
primitive, at least two nonzero coordinates, at most one zero residue modulo
three, excluding `(0,1,1)`, `(1,1,1)`, `(1,1,2)`, and `(1,2,2)`.
There are exactly 747 patterns. The companion reuses the audited fractional-
knapsack width algorithm but recomputes every exact record and checks the
complete frozen semantic digest.

For pattern p let alpha_p be its normalized slope and B_p its exact
open-interval intercept. An empty defect list is discharged as N=0;
strict counting is only used for a nonempty list. The finite compiler gives

```text
max_p B_p/(11/60-alpha_p)=35280/199<178,
attained at p=(7,13,18), alpha_p=17/147, B_p=12.        (10)
```

Thus the entire small-coefficient generic sector is count-safe at c>=178.
For `M>=19`, the inherited zonotope estimate is

```text
alpha=6/49+4/(7*19)=142/931,
Delta=11/60-alpha=1721/55860.
```

Together with `B_v<=2S/7+4/3`, the sufficient cutoff is

```text
c>=A S+B,
A=15960/1721, B=74480/1721.                           (11)
```

The parity-free short relation `S<4sqrt(c/3)` pays (11) for every c>=535:

* For S<=53, `A S+B<=920360/1721<535`.
* For S>=54, `g(S)=3S^2/16-A S-B` is positive at 54, where
  `g(54)=18547/6884`, and its derivative is positive thereafter.

This independently validates the c=535 analytic reduction. It does not
repeat root's separate full H535 native computation, and that larger finite
result is not a premise of the norm-five theorem (1). After the incoming
promotion, THM-4437 already gives every generic projection `<=6/77<11/140`
at all heights. Combining that current canon with THM-4441 and the separate
norm-four sharp closure is the shorter proof route to a global nonadditive
`11/140` theorem. The new H535 transcript is distinct independent evidence,
not a reason to hide the stronger incoming result.

The companion also checks two exact subsets requested for the parent's
combined classification of the old `6/77` boundary. Its complete additive
H61 subset has 146 rows; only `(1,4,5)` fails to exceed the old ceiling,
with both targets `1/28`. Its complete norm-four H34 subset has 88 rows;
only `(2,11,20)` exceeds the old ceiling, and only `(1,5,11)` equals it.
The exact tail controls are
`3/49+4/(7*35)=19/245<6/77` and
`9/98-6/(7*62)=237/3038>6/77`. These are finite-input and arithmetic checks
for a combined corollary using current canon, not a new independent
classification theorem claimed by this note.

## 7. Verification and exact scope

The companion uses SymPy for exact rational-function identities and real
inequalities, not floating point extremum searches. It checks all coefficient
records against semantic digest

```text
cf808062354debbefc1d8ead8ad0d10e9da5427cb42b8f083b6af24d0059c87c.
```

The required H63 data file has SHA256
`c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa`.
No repository mathematics is imported. The inherited width algorithm is
included directly and explicitly identified as reused, not counted as a
new independent implementation.

```text
python -B 04-computation/overnight5_20260906_lrc_norm5_profiles.py --head-table 05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv
python -B -O 04-computation/overnight5_20260906_lrc_norm5_profiles.py --head-table 05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv
```

Normal and optimized outputs match byte-for-byte on **161,225** explicit
gates. The frozen LF-byte hashes are

```text
source 38165ecdd4925b4015cf716ae0f741297428b3e06211a7a8403cbb95d4361f52
output cb5e46db7a5a6ea0bf84fbc9414be4c6e87eff3a92845169653fcea5fa6a959f
```

The finite bases,
literal-input lineage, norm-five switched-roof hostile, and empty-defect
branch are all explicit. No theorem ID, shared navigation entry, or Git
mutation is made in this lane. No arbitrary body floor, actual entry,
synchronization, or LRC(14) conclusion follows from these local results.

**Filing note.** Root filed this completed audit after checkpoint `07b2d91b2`.
Local paths were made repository-relative and output line endings normalized;
normal and optimized verification was rerun with matching output. The proof
and finite universes are unchanged.
