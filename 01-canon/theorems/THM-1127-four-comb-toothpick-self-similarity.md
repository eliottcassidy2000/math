---
id: THM-1127
title: Exact toothpick self-similarity on the (0,4,5,9) four-comb ray
status: PROVED / COMPUTER-ASSISTED — exact all-K tail on one fixed core and translated four-comb ray; not a uniform r=5 closure
source: codex-2026-07-18-S73 (r5-toothpick-function)
depends_on:
  - THM-1101
related:
  - THM-1097, MISTAKE-164
verification:
  - 04-computation/lrc14_r5_toothpick_self_similarity_exact_codex_S73.py
  - 05-knowledge/results/lrc14_r5_toothpick_self_similarity_exact_codex_S73.out
---

# THM-1127 — exact toothpick self-similarity on one four-comb ray

Fix

```text
P=(1,2,4,5,7,9,11,12),
(k1,k2,k3,k4)=(K,K+4,K+5,K+9).
```

Write `h=1/14`, let

```text
G_P = {t in R/Z : ||p t|| >= h for every p in P},
E_K = G_P minus (D_K union D_(K+4) union D_(K+5) union D_(K+9)),
```

and denote the component count, measure, and longest component length of
`E_K` by `N(K)`, `mu(K)`, and `L(K)`.  As in THM-1101, set

```text
Tmass(K) = N(K)/(6 mu(K)),
Tcomp(K) = 1/(3 L(K)),
T(K)     = min(Tmass(K),Tcomp(K)),
R(K)     = T(K)/(K+9).
```

The clustered range begins at `K>13 max(P)=156`.

## The theorem

### 1. Fixed torus carrier

For `A={0,4,5,9}` define the fixed rational polygon

```text
H = {(t,x) in (R/Z)^2 : ||x+a t|| >= 1/14 for every a in A}.
```

Then, exactly,

```text
E_K = {t in G_P : (t, Kt mod 1) in H}.                 (1)
```

Thus the four tooth combs do not form four unrelated moving objects.  They
are one fixed polygon sampled by the slope-`K` winding line.  This is the
precise toothpick self-similarity law.

### 2. Exact longest-component law

For every integer `K>=771`, put

```text
m_K = floor((97K-107)/126).
```

Then

```text
                 14m_K - 2K + 47
L(K) = -----------------------------------------.       (2)
       14 (K+4)(K+5)
```

Equivalently, with

```text
d_K = (97K-107) mod 126,       0 <= d_K < 126,
```

the period-126 phase form is

```text
                 79K + 316 - d_K
L(K) = -----------------------------------------.       (3)
       126 (K+4)(K+5)
```

The threshold is structurally sharp for this formula: its dominant strip is
not fully contained at `K=770`, where the containment margin is `-1` and

```text
L(770)=17/21672,
formula_candidate(770)=6781/8397900.
```

### 3. The sharp final-comb target holds at every legal scale

For every integer `K>=157`,

```text
7(K+9)L(K)>1.                                           (4)
```

The exact finite bank `157<=K<=770` has minimum

```text
min 7(K+9)L(K) = 695/346,       attained at K=173.
```

Formula (3) proves the entire tail.  Thus this ray satisfies the sharp
`L>1/(7k4)` target directly relevant to a fifth comb, even at scales where
the older `R<1` target fails.  This distinction is essential: the `R`
crossing is a failure of a coarser sufficient mechanism, not evidence that a
fifth comb covers this ray.

### 4. The last `R` crossing on this ray

Exact endpoint subtraction over `157<=K<=770`, combined with (2) for the
tail, proves

```text
R(K) < 1  for every K>=364.                             (5)
```

The boundary is genuine:

```text
R(363)=242/237 > 1,
R(364)=633696/759055 < 1.
```

There are exactly `27` clustered values `K>=157` with `R(K)>=1`; the largest
is `363`.  The maximum in the clustered scan is

```text
R(173)=2422/2085.
```

This closes the older `1/(3 k4)` component target on this ray only.  It does
not quantify over a different core or a different four-killer shape.

### 5. Exact drift and phase laws

Let `Q=194040`.  For every `K>9`, and hence throughout the clustered range,

```text
N(K+Q) = N(K) + 121304.                                (6)
```

Moreover,

```text
mu(K)              -> 59151097/627525360,
N(K)/K             -> 15163/24255,
K L(K)             -> 79/126,
Tmass(K)/K         -> 65382856/59151097,
Tcomp(K)/K         -> 42/79,
R(K)               -> 42/79.                           (7)
```

The two threshold slopes are separated:

```text
42/79 < 65382856/59151097.
```

Consequently the component branch is the eventual branch.  On each residue
class modulo `Q`, `mu(K)` is an exact rational function obtained by summing
affine-wall intersection endpoints.  Hence `Tmass` is quasirational with
phase period `Q`; (3) makes `Tcomp` rational with phase period `126`; and the
actual `R` is their lower envelope.  This is the functional form behind the
nonmonotone finite crossings.

## Proof

### Step 1: from four combs to one polygon

Put `x=Kt mod 1`.  For each `a in A`,

```text
||(K+a)t|| = ||x+a t||.
```

Taking the four inequalities simultaneously gives (1).  Notice also the
additive parallelogram

```text
(K+9)+(K)=(K+4)+(K+5).
```

The offset set is `{0,4}+{0,5}`.  The apparent four-comb motion is therefore
a rank-two torus arrangement, not four independent erosion events.

### Step 2: exact rational polygon atlas

The endpoints of a vertical danger arc are

```text
x=-a t +/- 1/14  (mod 1).
```

Two such walls cross only when

```text
(b-a)t = integer, integer+1/7, or integer-1/7.
```

There are `53` distinct wall times.  Their denominator lcm is `1260`.
Adding the endpoints of the `20` components of `G_P` gives `91` wall times,
`28` cells lying in `G_P`, and total denominator lcm

```text
Q=194040.
```

Exact Fraction integration over those `28` cells gives

```text
integral_(G_P) |H_t| dt       = 59151097/627525360,
integral_(G_P) #components(H_t) dt = 15163/24255.       (8)
```

The largest vertical safe gap is

```text
M=79/126,
```

attained exactly at `t=29/126` and `t=97/126`.  These are the outer
endpoints of the reflected core intervals

```text
I_-=[29/126,13/56],       I_+=[43/56,97/126].
```

Outside `I_- union I_+`, every vertical safe gap is at most

```text
M_2=67/154.                                             (9)
```

All of these assertions are finite rational comparisons in the referee.

### Step 3: solve the dominant strip

On `I_+`, the four danger arcs leave one vertical safe strip.  In a compatible
lift of the `x`-circle its walls are

```text
4-5t+1/14 < x < 4-4t-1/14,
```

of vertical width `t-1/7`.  Intersecting it with the lift `x=Kt-m` gives

```text
J_m = ((14m+57)/(14(K+5)), (14m+55)/(14(K+4))).        (10)
```

The right endpoint of `J_m` is at most `97/126` precisely when

```text
126m <= 97K-107.
```

Therefore `m=m_K` is the rightmost full component.  Its left endpoint lies in
`I_+` precisely when

```text
B(K):=56m_K-43K+13 >= 0.                               (11)
```

The thirteen initial values `B(771),...,B(783)` are

```text
12,25,38,51,8,21,34,47,4,17,30,43,0.
```

Also `B(K+13)-B(K)` is `1`, except when the floor remainder wraps, when it is
`57`.  Hence (11) holds for every `K>=771`.  By contrast, `B(770)=-1`.

The length of (10) at `m=m_K` is exactly the right side of (2).  Earlier full
components are shorter because their numerator drops by `14` each time.  To
check the only possible component clipped at `97/126`, write

```text
d=97K-107-126m_K.
```

It exists only for `d>=47` and has, over the common denominator
`126(K+4)(K+5)`, numerator `(d-47)(K+4)`.  The full component exceeds it by

```text
(126-d)K + 504 - 5d > 0       (K>121).                 (12)
```

Reflection handles `I_-` and gives the same maximum.

It remains to exclude the other core cells.  Their vertical gaps obey (9),
and every wall has slope of absolute value at most `9`.  Following a winding
line through such a labelled strip therefore gives

```text
component length <= (67/154)/(K-9).                    (13)
```

Indeed the line gains `K(v-u)` in the lifted `x` coordinate while either
strip wall can drift by at most `9(v-u)`.  From (3), using `d<=125`, the
dominant component is larger than

```text
79K/(126(K+5)^2).
```

The comparison with (13) is the increasing quadratic inequality

```text
79*154*K*(K-9) > 67*126*(K+5)^2,
```

which already holds at `K=771`.  Thus no other cell wins, proving (2) and
(3).

The same lower bound first gives the sharp target

```text
7(K+9)L(K) > 79K(K+9)/(18(K+5)^2) > 1
```

for every `K>=771`.  The exact `157<=K<=770` bank has minimum `695/346`,
proving (4).  It also gives the stronger old target

```text
L(K)>1/(3(K+9))   for K>=771.
```

Hence `R(K)<=Tcomp(K)/(K+9)<1` on the tail.  The exact finite bank through
`770` then proves (5) and the stated last crossing.

### Step 4: count recurrence and limits

On each of the `28` polygon cells, every safe strip is bounded by affine walls
of integer slope in `{0,-4,-5,-9}`.  Because `K>9`, the winding line crosses
each wall monotonically.  Replacing `K` by `K+Q` adds `Qt` to its phase.  At
every cell endpoint `u`, `Qu` is an integer, so the endpoint phase and all
gluing data are unchanged.  A strip over `[u,v]` gains exactly `Q(v-u)`
crossings.  Summing over strips and using (8) gives

```text
Q * 15163/24255 = 121304,
```

which proves (6), including the boundary-gluing terms.

The same affine-wall intersections are Riemann sums for the vertical widths.
Their limit is the first integral in (8).  This proves the first two limits
in (7).  Formula (3) proves `KL ->79/126`; the remaining threshold limits
follow by substitution.  Since the component limit is strictly smaller than
the mass limit, it is eventually the minimum and `R->42/79`.

## Verification

The standalone referee uses `fractions.Fraction` for every construction and
decision.  It checks the `53/91/28` rational wall atlas, all `1351` exact
endpoint rows `150<=K<=1500`, the `157..770` legal finite bank, the complete
period-126 tail guard, and four direct high-slope checks of (6).  Ordinary and
`PYTHONOPTIMIZE=1` executions are byte-identical to the frozen output.

```text
04-computation/lrc14_r5_toothpick_self_similarity_exact_codex_S73.py
SHA-256 d4024238821019506fd3b77a70efb244ebbb8e0de85f6c3cbe820614b7173b9e

05-knowledge/results/lrc14_r5_toothpick_self_similarity_exact_codex_S73.out
SHA-256 5f19804c0b821e668b6708f1adba422707f170d2c79ff0a835733ce0224b168d
```

## Kakeya/tournament reading

The useful Kakeya-like object is the family of slope lines

```text
Gamma_K={(t,Kt mod 1)}
```

through one fixed labelled torus polygon.  The hard transient is a wall-phase
problem: a line can graze a core endpoint or change which affine strip owns
the longest section.  That explains both the isolated late crossing at
`K=363` and why translating one tuple does not justify monotonicity for all
four-killer shapes.

For Tournament Analysis, runners, combs, core sections, section boundaries,
wall events, residues, cover arcs, Fourier modes, and proof obligations were
all challenged as vertices.  Exact wall-coordinate order gives a transitive
`91`-vertex tournament with no directed cycles, `91` singleton SCCs, score
multiset `{0,...,90}`, and one Hamiltonian path after ties are coalesced.  That
plain order is not faithful: it destroys wall slopes, affine owners, strip
adjacency, lengths, and the core mask.  The faithful carrier is the labelled
polygon plus the slope line.  It preserves the exact LRC safe-set predicate;
the naked tournament preserves only wall order.

## Scope and next use

This theorem proves a complete all-scale result for one offset shape.  It is
not a reduction from arbitrary

```text
k1<k2<k3<k4
```

to `(K,K+4,K+5,K+9)`, and it does not repair uniform `r=5` in THM-1101.
Its reusable contribution is the fixed-polygon method: for any bounded offset
shape, one can build a rational wall atlas, read exact quasiperiods, isolate
maximal vertical gaps, and replace unproved tail monotonicity by a finite
phase proof.  Uniform `r=5` now asks whether arbitrary offset shapes admit a
finite compactification or a separate estimate once the offsets grow with
the base scale.
