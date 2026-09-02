---
id: THM-4330
title: "LRC(14) affine two-adic root types and anchored clock/pool entry sieves"
status: >
  PROVED RELATIVE TO CITED LRC THROUGH THIRTEEN TOTAL RUNNERS,
  THM-366/771/2060/4150/4156/4191/4326/1042/4331/4332 + VERIFIED-EXACT.  The graph of pairwise
  differences of minimum two-adic valuation is the complete bipartite graph
  between the two primitive parity classes.  Its degree at a prescribed
  runner is exactly that runner's number of odd primitive relatives.  The
  singleton root type is safe for every prescribed runner, so the 12+2 split
  is the first unresolved dyadic type.  For a collision-free degree-twelve
  minority anchor with row {2h} union W, the inherited small-denominator
  sieve proves safety unless 420 divides h.  In that residual, the clocks
  1/2 plus-or-minus 1/(28h) still give an exact residue-band certificate and
  prove safety when max(W)<=12h.  Exact doubled-denominator unit banks and an
  adaptive THM-2060 divisor-capacity test give further sufficient criteria.
  A unit equality-wall safe-band lemma closes dC union T when 2<=d<=7,
  |T|=d, every tail is a d-unit, and max(T)>=(14-d)max(C); its d=5 consequence
  is new here, while d=2,3,7 recover prior cones.  Two explicit safe rows defeat the
  half-turn grid, all unit banks, and every represented capacity scale
  simultaneously; one is then closed by THM-771's owner-endpoint defect.
  MC12 separately defeats every represented capacity scale and rational
  denominator through forty.  These are anchored certificates, not a closure
  of the root type.  In the degree-two lane, midpoint
  collisions reduce to settled lower dimension and every collision-free row
  is literally an eleven-even/two-odd seam.  Re-referencing changes the
  distinguished runner and is not an entry proof.  A positive rational
  projective refactorization into at most two labels outside the fixed pool
  closes the seam; hence a counterexample in this branch has at least three
  outsiders after every such refactorization.  The actual primitive tail
  ratio gives the sharper sufficient test mu(G_H)>=mu(C_(p,q)). The
  displayed `H={1,...,11}`, tails `(1,9)` row fails all displayed projective,
  mass, component, and literal-pool gates but is safe with clearance `1/12`, so
  the gate union is not complete. LRC(14) and arbitrary-row entry remain open.
source: root + parity_entry + entry_corpus_audit + minority_clock_atlas + minority_residue_cover + minority_owner_quotient + owner_overlap_scout / LRC14 continuation session, 2026-09-01
depends_on:
  - LRCUpTo13
  - "T. Sungkawichai and T. Trakulthongchai, Eleven, twelve, and thirteen lonely runners, arXiv:2604.23906v1 (preprint)"
  - THM-366-lrc-small-denominator-divisibility-sieve
  - THM-771-sheet-endpoint-defect-and-reduced-winding-pierce
  - THM-2060-crt-tail-coset-saturation
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
  - THM-1042-component-length-obstruction-to-additive-certificates
  - THM-4331-lrc14-safe-component-endpoint-denominator-odd-wall-escape
  - THM-4332-lrc14-fixed-pool-single-constraint-implication-rigidity
related:
  - THM-639-hamiltonian-path-frame-for-runner-families
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2105-small-clock-affine-carrier-gate
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-4062-lrc-divisor-star-affine-intercept-obstruction
  - THM-4070-lrc14-d2-complete-small-denominator-sieve-and-affine-ray-firewall
  - THM-2888-eight-body-first-apex-global-pair-cap-atlas
  - THM-566-lrc14-covering-sets-have-no-uniform-bounded-denominator-witness
  - THM-4203-fixed-pool-seventeen-body-depth-eight-haar-completion
  - THM-4329-lrc14-complete-thirty-label-fixed-outsider-and-thirty-two-label-pascal-chart
  - THM-4333-lrc14-rank-three-surplus-and-cofinal-third-tail-completion
script: 04-computation/lrc14_entry_parity_affine_classification_probe.py
output: 05-knowledge/results/lrc14_entry_parity_affine_classification_probe.out
script_sha256: 2f026c4a09e3d44bf699f9cb92e4bcfd4658bb1ba82d6c6c4bff77bf3353e5f5
output_sha256: c29058584edf4c15bbb750e56a5193650cc92bac532865c5eab67ac23f2314da
safe_band_script: 04-computation/lrc14_unit_equality_wall_safe_band_probe_20260901.py
safe_band_output: 05-knowledge/results/lrc14_unit_equality_wall_safe_band_probe_20260901.out
safe_band_script_sha256: 122f95ffddcaece455d59f57c4bfb1bed895020e513574cd48901b4d18278c7e
safe_band_output_sha256: 697563661fd263560c8a52b613acfe556d795e08398475ba29f6a9c759b76dc0
joint_hostile_script: 04-computation/lrc14_minority_anchor_joint_clock_capacity_hostile_probe.py
joint_hostile_output: 05-knowledge/results/lrc14_minority_anchor_joint_clock_capacity_hostile_probe.out
joint_hostile_script_sha256: 6aca59304b31c623d1a44c631eaa8a345723f5fcbc5ef09ecd953bafc9251e9f
joint_hostile_output_sha256: b7ed216c305b6935d6b27b6eca84e91d753ae9603c52dc06cc5b353e99f3c0c1
endpoint_pierce_script: 04-computation/lrc14_r7_minority_endpoint_pierce_hostile_probe_20260901.py
endpoint_pierce_output: 05-knowledge/results/lrc14_r7_minority_endpoint_pierce_hostile_probe_20260901.out
endpoint_pierce_script_sha256: 423076c1c1f29cf5403ae732be29686661cfcebaae37cfd66c3122b6edb2f8ec
endpoint_pierce_output_sha256: 6c4d9ded0e7cf02689681085f9a52babaf15f1bc621cd661a6cd090f9b8e9c9d
audit: >
  PASS.  The graph theorem is proved symbolically below.  The exact probe
  checks all seven root-size profiles, all 3,060 fourteen-subsets of
  {-2,...,15}, all 858 configurations in a structured 12+2 universe, the
  AP and one-reference hostiles, the exact minority-anchor denominator,
  half-turn, unit-bank, adaptive-capacity, and stopping-hostile controls,
  literal and projective pool-entry controls, and 90 affine normalization
  images.  All 1,155,876 assertions pass.  Normal, optimized, and hash-seeded
  invocations produce the frozen transcript.  Three focused exact companions
  independently audit the unit safe band and both joint method hostiles; their
  normal, optimized, and hash-seeded transcripts also match byte-for-byte.
---

# THM-4330 -- affine two-adic root types and anchored clock/pool entry

**PROVED RELATIVE TO CITED LRC THROUGH THIRTEEN TOTAL RUNNERS,
THM-366/771/2060/4150/4156/4191/4326/1042/4331/4332 + VERIFIED-EXACT.  LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance pass

Let

```text
V={v_0,...,v_13} subset Z
```

be fourteen distinct labelled velocities.  Put

```text
g=gcd_(i<j)|v_i-v_j|,                 s=nu_2(g),          (1)
```

and form a labelled graph `Gamma(V)` on the runners by

```text
ij in E(Gamma(V))  iff  nu_2(v_i-v_j)=s.                 (2)
```

For a prescribed distinguished runner `r`, call

```text
u_i=(v_i-v_r)/g,                      i!=r,               (3)
```

its signed primitive relatives.  The values `|u_i|` form the primitive,
scale-equivalent positive-speed row after repetitions are identified.  If
`A={|u_i|:i!=r}`, the original relative row is `gA` and

```text
G_(gA)=m_g^(-1)(G_A),                    m_g(x)=gx.       (3a)
```

Thus the original and normalized rows have equivalent nonemptiness; the
prescribed anchor is unchanged.

> **Affine two-adic root theorem.**  The graph `Gamma(V)` is exactly the
> complete bipartite graph between the two parity classes of the normalized
> velocities.  Thus
>
> ```text
> Gamma(V)=K_(m,14-m)                                    (4)
> ```
>
> for some `1<=m<=13`, and
>
> ```text
> #{i!=r:u_i odd}=deg_Gamma(r).                          (5)
> ```
>
> Writing `k=min(m,14-m)`, the seven unordered root types and their edge
> counts are
>
> ```text
> k                 1   2   3   4   5   6   7
> |E(Gamma)|       13  24  33  40  45  48  49.           (6)
> ```

The root type `k=1` is safe at the LRC(14) threshold for **every prescribed
runner**.  Consequently a counterexample configuration must have `k>=2`.
This packages an old one-tail mechanism into the first intrinsic affine
classification; it is not a new proof of the one-tail lemma.

At `k=2`, the twelve majority-class runners have graph degree two and the two
minority-class runners have degree twelve.  For a prescribed degree-two
runner exactly one of the following holds:

1. two signed relatives have the same absolute value, equivalently the
   prescribed runner is the midpoint of two other velocities; settled
   lower-dimensional LRC closes the row;
2. the thirteen absolute relatives are distinct and have the literal form

   ```text
   2H union {a,b},       |H|=11,       a!=b positive odd. (7)
   ```

Globally, `k=2`, `|E(Gamma)|=24`, and the existence of a degree-two runner are
equivalent.  At least one degree-two runner is collision-free.  That last
existence statement is only diagnostic: changing to it changes the runner
whose loneliness is being proved.

For a collision-free prescribed degree-twelve runner, the primitive row has
the complementary form

```text
{2h} union W,      h>=1,      |W|=12,      W positive odd and distinct. (7a)
```

The THM-366 small-denominator sieve now has a sharp specialization:

> **Minority-anchor denominator sieve.**  The anchored row `(7a)` is safe
> whenever
>
> ```text
> 420 does not divide h.                                (7b)
> ```
>
> Consequently every counterexample in this anchored branch has `420|h`.
> It must also satisfy, for each `m in {9,11,13}`,
>
> ```text
> m|h or m|w for some w in W.                          (7c)
> ```

Condition `(7b)` is exact only for the inherited universal even-denominator
bank: when `420|h`, each of `4,6,8,10,12,14` divides `2h`, so its reciprocal
clock puts the even speed at zero.  This is a certificate failure, not an
unsafe row.

There is a complementary fixed-clock certificate even inside that residual.
Put `r_w` for the least positive residue of `w mod 28h`.  The two clocks

```text
t_+=1/2+1/(28h),              t_-=1/2-1/(28h)           (7d)
```

have identical distance vectors and witness the anchored row `(7a)` **if
and only if**

```text
r_w<=12h or r_w>=16h              for every w in W.      (7e)
```

In particular, `max W<=12h` is sufficient.  Since the members of `W` are
odd, that hypothesis gives `max W<=12h-1`; every odd relative then has
strict clearance at least `1/14+1/(28h)`, while `2h` is tight at `1/14`.
Twelve distinct odds force `h>=2`.  The cutoff is sharp only as the initial
monotone band for these clocks: `w=12h+1` is their first failure, but their
safe residues resume at `16h+1`.  This is an elementary fixed-clock
certificate for one prescribed minority anchor, not a closure of root type
`k=2`.

The closest proved mechanisms are THM-366's small-denominator sieve,
THM-2060's one-tail sheet dodge, and THM-2061's exact two-tail dyadic seam.
The canonical hostile is
`V={0,1,...,13}`, whose graph is `K_(7,7)`.  The corrected near miss is the
forbidden inference "some reference has degree two, therefore the original
runner enters (7)"; THM-2888 already isolates this change-of-observer error.
The least-used relevant sidecars are the prescribed runner label and, at a
minority anchor, the normalized same-class half-gap `h`.

## 2. Proof of the root theorem and the singleton closure

Choose any base velocity `v_*` and put `z_i=(v_i-v_*)/g`.  The gcd of all
`z_i-z_j` is one, so both parities occur.  Moreover

```text
nu_2(v_i-v_j)=s
 iff (v_i-v_j)/g is odd
 iff z_i and z_j have opposite parity.                  (8)
```

This proves `(4)`.  With `v_*=v_r`, equation `(8)` says that the neighbors of
`r` are precisely its odd primitive relatives, proving `(5)`.  Formula `(6)`
is `k(14-k)`.

Now suppose `k=1`.  If the prescribed runner is the unique minority vertex,
all thirteen primitive relatives are odd, and the clock `x=1/2` gives every
relative distance `1/2`.

If the runner lies in the thirteen-vertex majority class, its primitive row
has the form

```text
2C union {w},                    w odd,                 (9)
```

with at most twelve distinct body speeds.  Cited LRC through thirteen total
runners supplies a phase `y` at which every member of `C` has distance at
least `1/13>1/14`.  The two lifts

```text
x_0=y/2,                    x_1=(y+1)/2                 (10)
```

preserve every doubled-body distance.  Since `w` is odd, its two phases
differ by `1/2`; the two open danger arcs of radius `1/14` cannot contain
both.  One lift therefore handles the tail.  This is THM-2060's `q=2`
one-tail mechanism and proves the complete `k=1` closure.

For any prescribed runner, and distinct `i,j` with `i,j!=r`, repeated
absolute relatives occur exactly when

```text
|v_i-v_r|=|v_j-v_r|  iff  v_i+v_j=2v_r.                (11)
```

Such a row has at most twelve distinct positive constraints, so cited LRC
through thirteen total runners gives the stronger threshold `1/13`.  If a
degree-two runner has no collision, `(5)` leaves two distinct odd absolute
relatives and eleven distinct even ones, proving `(7)`.

If instead the prescribed runner has degree twelve and no collision, `(5)`
leaves twelve distinct positive odd relatives and one positive even relative.
Writing the latter as `2h` proves `(7a)`.

Suppose first that `420=lcm(2,3,4,5,6,7)` does not divide `h`.  Choose
`d in {2,3,4,5,6,7}` with `d` not dividing `h`, and put `m=2d`.  No odd
member of `W` is divisible by the even integer `m`, while `m` does not divide
`2h`.  THM-366 therefore makes `t=1/m` a lonely witness at threshold `1/14`.
Conversely, all six even denominators `m in {4,6,8,10,12,14}` divide `2h`
exactly when every `d in {2,...,7}` divides `h`, equivalently `420|h`.
This proves `(7b)` and the exact boundary of this particular clock bank.

THM-366 also says that every counterexample row must contain a speed divisible
by each `2<=m<=14`.  Once `420|h`, the speed `2h` handles every such `m`
except possibly `9,11,13`.  For each of those three odd denominators, a
divisible speed is either `2h`, equivalently `m|h`, or an odd member of `W`.
This proves `(7c)`.

Finally assume `Gamma=K_(12,2)`.  Take the least and greatest velocities in
the twelve-vertex parity class.  Neither can be the midpoint of two vertices
from that same class.  A collision can therefore use only the two minority
vertices, which have one midpoint; at most one of the two extremes equals
it.  Hence at least one extreme is collision-free.  The example

```text
V={-1,1,0,2,4,...,22}                                  (12)
```

has exactly one collision-free degree-two reference, namely `22`, so the
lower bound one is sharp.  Equation `(12)` does not license a change from a
different prescribed runner.

It remains to verify the half-turn clocks.  At either clock in `(7d)`,

```text
||2h t_+||=||2h t_-||=1/14.                            (12a)
```

For odd `w`, reduction modulo `28h` gives

```text
w t_+ = 1/2 + r_w/(28h)                    (mod 1),
w t_- = 1/2 - r_w/(28h)                    (mod 1),
||w t_+||=||w t_-||=|r_w-14h|/(28h).                  (12b)
```

Thus the odd constraint is safe exactly when
`|r_w-14h|>=2h`, which is `(7e)`.  Also `t_-=1-t_+`, explaining the identical
distance vectors.  The initial-band corollary and its strict margin follow
immediately.  No coprimality assumption is used.

## 3. Minority-anchor residual atlas and stopping hostiles

Retain the primitive collision-free row

```text
S_h={2h} union W,                    |W|=12 odd.          (MC1)
```

The complete safe set of the even speed alone is

```text
E_h=union_(j=0)^(2h-1)
    [(14j+1)/(28h),(14j+13)/(28h)]       in R/Z.         (MC2)
```

Indeed, this is the full preimage under multiplication by `2h` of the closed
arc `[1/14,13/14]`.  The following two sieves sample `(MC2)` in complementary
ways.

The faithful continuous quotient is already implicit in THM-4062 and in the
two-tail fold of THM-2061.  Put `x=2t`, choose a real representative of `x`,
and label its two lifts by

```text
t_epsilon=(x+epsilon)/2,                   epsilon in Z/2. (MC2a)
```

Both lifts preserve the anchor exactly:

```text
||2h t_epsilon||=||h x||.
```

For odd `w`, if `||wx||<1/7`, let `N_w(x)` be its unique nearest integer.
Then `w` is strictly dangerous at exactly the lift

```text
epsilon=N_w(x) mod 2;                                  (MC2b)
```

if `||wx||>=1/7`, it is dangerous at neither lift.  Equality at `1/7`
becomes equality at `1/14` and is safe.  Consequently, with

```text
G_h={x:||hx||>=1/14},
Sigma_2(W)={x:both lift labels occur among the active owners}, (MC2c)
```

the exact obstruction is

```text
M(S_h)<1/14                 iff                 G_h subset Sigma_2(W). (MC2d)
```

Changing the representative `x` to `x+1` swaps both lift labels and every
owner parity, so coverage in `(MC2d)` is gauge invariant although an
individual owner bit is not.  The clock and capacity tests below quotient
this labelled location data; `(MC2d)` records what their stopping hostiles
say must be restored.

### 3.1 Doubled-denominator unit bank

Fix `8<=p<=14` with `p` not dividing `h`, and put

```text
U_p=(Z/(2p)Z)^x/{+1,-1},
R_p(W)={[w]:w in W and gcd(w,2p)=1}.                    (MC3)
```

Then the bank

```text
t_a=a/(2p),                    a in (Z/(2p)Z)^x         (MC4)
```

contains a witness for `(MC1)` **if and only if** `R_p(W)!=U_p`.

To prove this, first

```text
||2h t_a||=||ah/p||>=1/p>=1/14,
```

because `a` is a unit and `p` does not divide `h`.  For odd `w`, a bad phase
has an odd circular numerator `r` satisfying `7r<p`.  In the range
`8<=p<=14`, this forces `r=1`.  Thus `w` is bad at `t_a` exactly when

```text
aw=+1 or -1                         (mod 2p).            (MC5)
```

A nonunit `w` deletes no clock, while a unit `w` deletes exactly the sign
class `[w^(-1)]`.  This proves the iff criterion.  On the `420|h` wall, the
only potentially useful `p` are `8,9,11,13`, and only those not dividing the
actual `h`; `10,12,14` already divide `h`.

### 3.2 Adaptive divisor capacity

Let `a>=2` divide at least one speed of `S_h`.  Put

```text
D_a={v in S_h:a|v}=aC_a,
T_a=S_h\D_a,
q_a(v)=a/gcd(a,v)                    for v in T_a.       (MC6)
```

Primitivity makes `D_a` proper, so `C_a` has at most twelve distinct positive
speeds.  Cited lower-dimensional LRC supplies a `1/14`-safe core phase.
THM-2060's `a` lifts preserve `aC_a`, and tail `v` occupies at most the
fraction `ceil(q_a(v)/7)/q_a(v)` of them.  Consequently

> **Adaptive divisor criterion.**  The row `(MC1)` is safe whenever
>
> ```text
> sum_(v in T_a) ceil(q_a(v)/7)/q_a(v) < 1.             (MC7)
> ```

Hence a counterexample on the `420|h` wall must reverse `(MC7)` for every
integer `a>=2` represented as a divisor of at least one row speed.  Some
immediate necessary restrictions are

```text
r_3>=3,             r_5>=5,             r_7>=7,         (MC8)
```

where `r_p=#{w in W:p does not divide w}`.  For the remaining THM-366 gates,
write

```text
A=#{w:3 does not divide w},
B=#{w:3|w and 9 does not divide w},
d_p=#{w:p|w}.
```

Then a counterexample must satisfy

```text
9 does not divide h:  d_9>=1 and 2A+3B>=6;
9 divides h:                     2A+3B>=9;
11 does not divide h: 1<=d_11<=7;    11 divides h: d_11<=6;
13 does not divide h: 1<=d_13<=6;    13 divides h: d_13<=5.  (MC9)
```

These are necessary capacity and carrier constraints, not a classification.

There is a sharp repair on the equality walls in `(MC8)`.  It retains a
tail's **location** on the common sheet deck rather than only its capacity.

> **Unit equality-wall safe-band lemma.**  Let `2<=d<=7`, and let
>
> ```text
> S=dC union T,       |C|=13-d,       |T|=d,            (MC9a)
> ```
>
> be a set of thirteen distinct positive speeds with `gcd(d,w)=1` for every
> `w in T`.  Put `R=max C` and `mu=M(C)`.  Then `S` is safe whenever some
> `w in T` satisfies
>
> ```text
> w >= dR/[14(mu-1/14)].                                (MC9b)
> ```
>
> In particular, cited lower-dimensional LRC gives `mu>=1/(14-d)`, so the
> scale-free condition
>
> ```text
> max T >= (14-d)R                                     (MC9c)
> ```
>
> is sufficient.

To prove this, choose a maximizer `theta_0` of the core lower envelope.  For
the selected tail `w`, define the closed no-danger band

```text
P_w={theta:{w theta} in [d/14,1-d/14]}.                (MC9d)
```

Every complementary component of `P_w` has length `d/(7w)`, so some
`theta in P_w` lies within `d/(14w)` of `theta_0`.  The core is `R`-Lipschitz;
thus `(MC9b)` gives

```text
min_(c in C)||c theta|| >= mu-dR/(14w) >= 1/14.        (MC9e)
```

On the `d` lifts `t_j=(theta+j)/d`, every speed in `dC` keeps this clearance.
Because `w` is a unit modulo `d`, its phases are the translated `d`-grid;
membership in `(MC9d)` says none is strictly dangerous.  Each of the other
`d-1` unit tails occupies at most one lift, since the open danger arc has
length `1/7<=1/d`.  At least one lift remains safe.  Closed endpoints and
the weak inequality in `(MC9b)` are load-bearing.

At `d=2`, condition `(MC9c)` is exactly THM-4052 `(12)`'s coarser
`max T>=12R` cone.  The same proof likewise recovers already-known cones at
`d=3` and `d=7`; among the minority-anchor equality branches
`p in {3,5,7}`, its genuinely new specialization is `p=5`.

For the minority anchor on the `420|h` wall, `(MC9c)` says that a
counterexample on an equality branch `r_p=p`, `p in {3,5,7}`, must obey

```text
max T_3 < 11 max C_3,
max T_5 <  9 max C_5,
max T_7 <  7 max C_7,                                  (MC9f)
```

where `pC_p` is the `p`-divisible core and `T_p` its complement.  The
`p=3` cone is already subsumed by THM-4052, and `p=7` is THM-771 `(11)`.
The unified proof exposes the previously unused `p=5` equality-wall cone.

### 3.3 Joint method hostiles and the owner-endpoint repair

The preceding three marginal certificates can fail simultaneously.  Put

```text
h=420,
D={11+1680k:0<=k<=6},
C={525,945,1365,1575},
W_P=D union C union {P},              P in {1287,9009}. (MC10)
```

Both rows `{840} union W_P` are primitive and THM-366-complete: the anchor
handles every denominator through fourteen except `9,11,13`, while
`1287=9*11*13` and `9009=7*9*11*13` each handle those three.

Neither row has a witness in the complete integer half-turn grid

```text
t_j=1/2+j/11760,                         j in Z.         (MC11)
```

Indeed, `14|j` puts the anchor at zero.  If `7` does not divide `j`, the
seven products of `D` are

```text
11j+1680jk,                       0<=k<=6,
```

and traverse the seven `1680`-blocks modulo `11760`; the within-block
residue is nonzero because its vanishing would force `1680|j`.  Exactly one
product is therefore strictly in the bad block `(5040,6720)`.  If `j=7s`
with `s` odd, write

```text
C=105{5,9,13,15}.
```

After division by `735`, danger is `sv=+7 or -7 (mod 16)`.  The four values
of `v` represent the four odd-unit sign classes modulo sixteen, so exactly
one member of `C` is strictly bad.  These assigned blockers avoid the
endpoints; other, harmless endpoint incidences do occur.

The common set `D union C` also represents every unit sign class needed for
`p=8,9,11,13`, while the anchor kills `p=10,12,14`.  Thus both rows defeat
every bank `(MC4)`.  Exact divisor enumeration gives the sharper joint
stopping data

```text
P=1287: 71 represented scales, min(MC7)=8/7
        exactly at a=7,21,35,105;
P=9009: 78 represented scales, min(MC7)=1
        exactly at a=7,21.                              (MC11a)
```

Neither is an unsafe row.  For `P=1287`, the clock `t=6/17` has minimum
clearance `2/17>1/14`.  The `P=9009` row shows why owner location is the
right repair.  At the seven-sheet core phase `theta=1/22`,

```text
C_7={v/7:7|v}={120,75,135,195,225,1287}
```

has least residue numerator `3 mod 22`, hence is strictly core-safe.  Tail
`11` lies on a THM-771 endpoint and owns no strictly bad sheet; the other
six tails own sheets

```text
0,2,2,4,6,6.
```

Thus THM-771's exact defect is `F=Q+Omega=1+2=3`, with free sheets
`1,3,5`.  The middle lift `t=67/154` has minimum clearance
`3/22>1/14`.  Capacity forgot both the endpoint slack and the two owner
collisions.  The two focused exact companions freeze both rows, including
normal, optimized, and hash-seeded replay.

### 3.4 A deeper finite-clock and capacity hostile

Put

```text
h=180180,
W={5,7,11,13,17,23,29,31,49,75,111,135}.               (MC12)
```

Here `2h=360360=lcm(1,...,15)`, so every THM-366 denominator is already paid.
An exact reduced-rational scan proves that every clock of denominator at most
`40` which is safe for `2h` is strictly killed by `W`.  After the symmetry
`t -> 1-t`, the complete anchor-safe universe has `112` clocks, at

```text
q/count = 16/4,17/7,19/8,23/10,25/10,27/9,29/12,
          31/13,32/8,34/7,37/16,38/8.                  (MC13)
```

The adaptive capacity `(MC7)` also fails for every one of the `202` integers
`a>=2` dividing at least one displayed speed; its exact minimum is `10/7` at
`a=7`.  Yet `(MC12)` is safe: `t=1/41` has minimum clearance `5/41>1/14`.
Thus denominator `41` is its exact first reduced-rational witness.  This
freezes the stopping reason for both clock-polishing and divisor-capacity
routes; it is not an LRC counterexample.

## 4. Exact projective fixed-pool entry for degree-two anchors

Retain the fixed pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (13)
```

Let `(7)` be an already anchored, collision-free row.  Suppose there is a
positive reduced rational

```text
lambda=p/q,       gcd(p,q)=1,                             (14)
```

such that

```text
H'=lambda H subset Z_(>0),       |H'\P|<=2.             (15)
```

Then `(7)` is safe.  Equivalently, a counterexample in this literal dyadic
branch must satisfy

```text
|lambda H\P|>=3                                           (16)
```

for every lawful `lambda` in `(14)` for which `lambda H` is integral.  In
particular, `|H\P|>=3` when `lambda=1`.  More generally, an already displayed
row `2cH_0 union {a,b}` cannot be a counterexample when `|H_0\P|<=2`.

### Proof

Integrality in `(15)` and `gcd(p,q)=1` imply `q|h` for every `h in H`; write
`H=qK`, so `H'=pK`.  Multiplication by any positive integer preserves Haar
measure on the circle.  Hence

```text
mu(G_H)=mu(G_K)=mu(G_(H')).                             (17)
```

If `H'` has zero outsiders, THM-4156 gives `mu(G_(H'))>4/63`; with one
outsider THM-4191 gives `mu(G_(H'))>=4/63`; with two outsiders THM-4326 gives
the same bound.  Equation `(17)` transfers that bound to the original body
without changing either tail, and THM-4150 closes `(7)`.

Conversely, every rational rescaling available to this Haar-invariance move
comes from an equality `pH=qH'` of positive integer multiples and therefore
has `H'=(p/q)H`.  Thus `(14)--(15)` exactly characterize **entry through this
projective fixed-pool mass template**, not safety itself.

The test is finite.  If `(15)` holds, at least nine image labels lie in `P`.
For one of them, `lambda=u/h` with `u in P` and `h in H`.  It is therefore
enough to inspect at most `30*11` candidate rational scales, reduce them,
check integrality, and count pool hits.  No tail coordinate is rescaled, so
the numerator parity is irrelevant.

## 5. Parity entry is strictly weaker than pool entry

Consider

```text
V={0,2,4,...,22} union {1,3}.                            (18)
```

The anchors `0` and `22` are collision-free degree-two references.  At either
one the half-body is

```text
H_0={1,2,...,11}.                                       (19)
```

Thus `(18)` passes the affine parity gate.  It does not pass the projective
pool gate.  Since `1 in H_0`, every integral reduced scale in `(14)` has
`q=1`; the exact finite candidate scan finds at most six labels in `P`, with
six attained only at scale `10`.  The required number is nine.  This is a
certificate failure, not an unsafe row.

The rational freedom is nevertheless real.  Put

```text
H_1={1,2,4,5,8,10,15,20,21,30,40}.                     (20)
```

There is no literal fixed-pool entry at the displayed anchor, but the even
numerator scale `lambda=2` gives

```text
2H_1={2,4,8,10,16,20,30,40,42,60,80},                  (21)
```

with nine pool labels and two outsiders.  This is the boundary control that
prevents replacing `(14)` by an integer body divisor test.

## 6. Pair-adaptive Haar transfer

For an eleven-body `H`, put

```text
G_H={y in R/Z:min_(h in H)||hy||>=1/14}.                (22)
```

Write the distinct odd tails as

```text
a=pt,       b=qt,       0<p<q,       gcd(p,q)=1,        (23)
```

so `p,q,t` are odd.  Let `C_(p,q)` be THM-4150's open two-sheet cross-comb.
Then

> **Pair-adaptive criterion.**  If `G_H` is nonempty and
>
> ```text
> mu(G_H)>=mu(C_(p,q)),                                 (24)
> ```
>
> the row `2H union {a,b}` is safe.  Hence every counterexample in this
> branch satisfies the strict, ratio-sensitive obstruction
>
> ```text
> mu(G_H)<mu(C_(p,q))<=4/63.                            (25)
> ```

For an eleven-body, the nonemptiness in `(24)` is automatic: cited
lower-dimensional LRC supplies a point with clearance at least `1/12`, and
continuity supplies a nonempty open subset at the weaker threshold `1/14`.

Indeed, failure gives exactly the THM-4150 containment

```text
G_H subset m_t^(-1)(C_(p,q)),             m_t(y)=ty.    (26)
```

The right side is a proper open set of measure `mu(C_(p,q))`, while `G_H` is
nonempty and compact.  Strict measure excess contradicts inclusion.  At
equality, the open difference would have measure zero and hence be empty, so
the two sets would form a nonempty proper clopen subset of the circle, again
impossible.  This proves `(24)--(25)`.

The exact threshold is executable from THM-4150:

```text
mu(C_(p,q))
 =2/49+2[B_2({1/2+(q-p)/14})-B_2({1/2+(q+p)/14})]/(pq), (27)
```

where `B_2(u)=u^2-u+1/6` on `[0,1)`.  Its maximum is `4/63`, uniquely at
`(p,q)=(1,9)`; it is zero at `(1,3)` and `(1,5)`.  Thus `(24)` is strictly
stronger than the uniform THM-4150 test on every submaximal ratio.

### 6.1 A combined-gate survivor

The present sufficient gates are not complete, even in a highly structured
degree-two row. Take

```text
V_0={0,2,4,...,22} union {1,9}                         (28)
```

with prescribed anchor `0`. Its body is `H_0={1,...,11}` and its odd-tail
ratio is `(1,9)`. Section 4's exact scan shows that every lawful projective
image of `H_0` has at most six pool labels, below the required nine. Also

```text
mu(G_(H_0))=10931/194040<4/63=mu(C_(1,9)),             (29)
```

so `(24)` fails. THM-1042 gives maximum positive component width `1/77`.
Every endpoint denominator is divisible by `14`, hence for `b=9`

```text
bW<=9/77<3/14<=2/7-1/min(Q_L,Q_R),                    (30)
bW<=9/77<1/7<=2/7-1/Q_L-1/Q_R.
```

so no THM-4331 component passes either. THM-4332 separately prevents an
unscaled pool subset from pointwise implying the outside labels needed by
this body.

Nevertheless the row is strongly safe: at `x=5/24`, direct reduction gives

```text
min_(v in V_0\{0})||vx||=1/12>1/14.                   (31)
```

Thus failure of the combined projective, mass, component, and literal-pool
certificates is not danger. A bounded residual-clock sidecar is genuinely
missing from this gate atlas.

## 7. Connection contract and scope firewall

```text
source:       labelled affine velocity configuration with prescribed runner
target:       minimum-v2 cut graph, then either one-even/twelve-odd or
              eleven-even/two-odd anchored rows
map:          divide affine content; take the prescribed anchor star;
              identify equal absolute relatives only after recording them
preserved:    anchor-relative LRC predicate, parity, common dilation class
destroyed:    magnitudes and signs after the parity quotient; the graph alone
              also forgets midpoint collisions and fixed-pool membership
sidecar:      anchor label, signed relatives, midpoint flag; either h, W,
              residue/unit classes, sheet owners and divisor capacities,
              or H and p/q
decisive test: anchor degree; midpoint equation; h mod 420; exact clock/owner
               deck and capacity; finite u/h pool scan; actual-ratio Haar test
```

Scope firewalls:

1. Re-referencing does not preserve the distinguished runner.  The existence
   of one collision-free degree-two anchor proves nothing for another anchor.
2. The `k=1` result is a complete root-type closure, but the theorem does not
   close `k=2,...,7`.  In type `k=2`, `(7b)--(7e)` certifies a cofinite
   divisibility branch and one explicit residue subcone for a degree-twelve
   anchor; it does not close the `420|h` residual.  Those anchors are not in
   the two-tail lane.
3. The unit-bank and adaptive-capacity tests `(MC3)--(MC9)` are sufficient,
   not necessary.  The two safe rows `(MC10)--(MC11a)` defeat all three
   marginal banks simultaneously; MC12 separately defeats every reduced
   denominator through forty.  They are not evidence of an LRC
   counterexample.  THM-771 closes one after restoring endpoint-owner data:
   endpoint slack closes the row, while two owner collisions increase the
   free-sheet count to three.
4. Conditions `(14)--(15)` characterize entry through the stated rational
   fixed-pool mass template.  Their failure is not danger, and pool labels
   are not a projective invariant.
5. The scalar condition `(24)` is sufficient, not necessary.  When it fails,
   THM-2061's component locations and dyadic owner words remain load-bearing.
6. No arbitrary-row entry, owner/arrival transfer, or proof of LRC(14) follows.

THM-4332 supplies the exact literal boundary at scale `lambda=1`:
`G_P subset G_h` holds only for `h in P`. This does not contradict the
whole-body positive-rational Haar refactorization in Section 4. THM-4333's
rank-three surplus is likewise a fixed residual-pair mass theorem, not an
affine entry map.

The successful research move is to pay the affine rank first and then audit
the physical action: the quotient graph exposes the exact parity type, while
the anchor and projective-scale sidecars prevent a change-of-observer or
confusion between body-mass equivalence and physical tail rescaling.

## 8. Exact probe and replay

The probe is corroborative; the complete-bipartite proof above is
symbolic.  It freezes the following hostile and positive controls:

- all seven root sizes and edge counts in `(6)`;
- all `binom(18,14)=3,060` bounded configurations in `{-2,...,15}`;
- all `13*binom(12,2)=858` configurations in the structured split universe;
- the `K_(7,7)` AP hostile, the sharp one-reference example `(12)`, and the
  parity-without-pool example `(18)`;
- the even-denominator sieve through four complete `420`-periods, its exact
  `420|h` wall, and the residual small denominators `9,11,13`;
- three exact periods of the minority-anchor residue formula for every
  `1<=h<=40`, its first nonvacuous twelve-odd positive row, the first failed
  residue, and resumed safety;
- the doubled-denominator iff criterion for every `8<=p<=14` and an
  additional complete `h=420` grid/unit-bank hostile with essential-blocker
  controls;
- focused exact replays of both joint rows `(MC10)`, all their represented
  capacities, and the `P=9009` THM-771 owner defect;
- all reduced rational clocks through denominator `40` and all `202`
  represented adaptive divisor scales for `(MC12)`, plus its exact first
  witness `1/41`;
- literal fixed-pool and even-numerator projective positive entries;
- translation and nonzero common-dilation invariance on 90 images.

Replay:

```bash
python3 04-computation/lrc14_entry_parity_affine_classification_probe.py \
  | diff -u \
  05-knowledge/results/lrc14_entry_parity_affine_classification_probe.out -

python3 -O 04-computation/lrc14_entry_parity_affine_classification_probe.py \
  | diff -u \
  05-knowledge/results/lrc14_entry_parity_affine_classification_probe.out -

PYTHONHASHSEED=1729 \
python3 04-computation/lrc14_entry_parity_affine_classification_probe.py \
  | diff -u \
  05-knowledge/results/lrc14_entry_parity_affine_classification_probe.out -
```

All three commands pass with the assertion count printed in the frozen
transcript.  The focused companions replay by the same pattern:

```bash
python3 04-computation/lrc14_unit_equality_wall_safe_band_probe_20260901.py \
  | diff -u 05-knowledge/results/lrc14_unit_equality_wall_safe_band_probe_20260901.out -
python3 04-computation/lrc14_minority_anchor_joint_clock_capacity_hostile_probe.py \
  | diff -u 05-knowledge/results/lrc14_minority_anchor_joint_clock_capacity_hostile_probe.out -
python3 04-computation/lrc14_r7_minority_endpoint_pierce_hostile_probe_20260901.py \
  | diff -u 05-knowledge/results/lrc14_r7_minority_endpoint_pierce_hostile_probe_20260901.out -
```

Optimized and hash-seeded versions also byte-match.  All hashes are recorded
in the frontmatter.
