# Three aligned / four drift compression — scratch proof note

**Status:** finite-exact scratch result under audit; not canonical and not an
LRC(14) closure.

## 1. Scope and inherited object

Fix a literal six-speed body

```text
F subset {1,...,14}, |F|=6,
L=14 lcm(F),
J_F = body-safe 1/L-cell word.
```

Assume three tail combs are aligned to the body ruler and four are genuine
drifts.  A denominator-one drift is another aligned comb, so this note treats
only

```text
2 <= d_1 <= d_2 <= d_3 <= d_4,
D=lcm(d_1,d_2,d_3,d_4) | L.
```

This exclusion uses the now-proved mixed four-aligned / three-drift closure
in THM-2928; it is not hidden inside the computation.

Put `S_D=J_F mod D` and `S=|S_D|`.  THM-2928's body projection and the safe
floors used by THM-2941 give

```text
mu(aligned safe carrier) >= u_3=55/91,
four-drift union capacity <= 1-u_4=625/1183,
S/D <= (1-u_4)/u_3 = 125/143.                 (1)
```

The exact inherited universe has `26,970` body/divisor rows on `217`
resolving denominators.  Divisor-lattice inversion gives `694,921,995`
denominator four-multisets and `21,357,714,101` row-occurrences.  The point
of this note is to compress that universe without iterating those
occurrences.

## 2. Exact divisor-lattice generating function

For `d|D`, let

```text
C_D(d)=(D/d) ceil(d/7).
```

The arity-four multiset series is

```text
G_(D,E)(y,x,z)
 = product_(d|E,d>1)
     (1-y x^(C_D(d)) z_2^[d=2] z_3^[d=3] z_4^[d=4])^(-1)
       mod y^5.
```

Its lcm-exact part is

```text
H_D = sum_(E|D) mu(D/E) G_(D,E).                       (2)
```

Divisors with an identical feature are grouped before multiplying.  If a
feature occurs on `g` divisor symbols, selecting it `m` times has coefficient
`binom(g+m-1,m)`.  Thus the coefficient of `y^4` in `(2)` is stored only by

```text
(number of d=2 labels,
 number of d=3 labels,
 number of d=4 labels,
 total ambient capacity of all d>4 labels).            (3)
```

State `(3)` is sufficient for every common-phase screen below, and summing
its multiplicities recovers the full `694,921,995` shapes exactly.  It is a
compression of the computation, not a quotient assertion about realizable
drift packets.

## 3. Common-\(u\) activity and the equality-strengthened obstruction

Let `K_A` be the normalized safe carrier of the three aligned combs.  It is
compact and

```text
mu(K_A) >= 55/91.                                      (4)
```

Write a denominator-`d` drift numerator as `c`, with `gcd(c,d)=1`.  For
`d<7`, its fixed-`u` address mask is nonempty exactly when

```text
A_i(u): ||c_i u|| < d_i/14.
```

When active it occupies one lifted residue class, and

```text
mu(A_i)=d_i/7.                                         (5)
```

Only `d_i<=4` has activity mass below `(4)`.  Denominators five and six can
be synchronized on a set already larger than `55/91`, so arbitrary-coupling
one-marginal transport alone cannot exclude them.

More generally, if `u_k` is the safe floor for `k` aligned labels, the
denominators whose single-label activity marginal is at most `u_k`—and
hence are detectable when that label is essential by this one-marginal
screen—are exactly

```text
d <= 7u_k.
```

For `k=2,3,4,5,6` this gives, respectively,

```text
{2,3,4,5}, {2,3,4}, {2,3}, {2}, empty.               (5a)
```

Here equality is included: compact-safe/proper-open separation also kills
the endpoint `d/7=u_k`.  (For the displayed `k=2,...,6` floors no integer
lies at an endpoint, so the lists themselves are unchanged.)  This explains
the small-denominator alphabets found independently in the
five-, four-, and three-drift sectors.  Larger denominators are **not**
proved viable; they are merely invisible to this one-marginal test.  Indeed,
when every marginal is strictly above `u_k`, the arbitrary-coupling
relaxation permits a joint law synchronizing all statuses on mass at least
their smallest marginal, so no nonempty upward event can have maximum mass
at most `u_k`.

For each `d<=4`, one may grant either its ambient load `D/d`, or the stronger
row-specific bound

```text
lambda_d(S_D)=max_(a mod d) |S_D intersect (a mod d)|. (6)
```

Grant every `d>4` mask its full ambient capacity and call their total `B`.
For certified small-mask loads `C_i`, define the upward threshold family

```text
A_R={E subset J: sum_(i in E) C_i >= R},
R=S-B,       J={i:d_i<=4}.                             (7)
```

Any cover forces the actual common-phase status into `A_R` at every
`u in K_A`.

Let `H=min(A_R)`.  The general upward-status theorem gives, for arbitrary
dependence with the exact marginals `(5)`,

```text
mu(status in A_R)
 <= tau
 := min(1, min {
        sum_i (d_i/7) w_i:
        w_i>=0, sum_(i in H0)w_i>=1 for every H0 in H
      }).                                               (8)
```

The usual measure-only contradiction kills when `tau<55/91`.  Here equality
also kills.  The event in `(7)` is open, while `K_A` is compact.  If
`tau<=55/91<1`, the event is proper.  A compact subset of a proper open
subset of the circle has strictly smaller measure, so

```text
55/91 <= mu(K_A) < mu(A_R) <= tau <= 55/91,
```

a contradiction.  Therefore the exact necessary condition is

```text
tau > 55/91.                                           (9)
```

Capacity equality is nevertheless included in `(7)`.  For ambient pattern
`(d_1,d_2)=(2,3)` and residual `R/D=1/3`, either singleton can pay, so
`tau=5/7` and the endpoint survives.  Immediately above `1/3`, only the
denominator-two label can pay and `tau=2/7`, so it dies.  Replacing `>=` in
`(7)` by `>` would falsely kill the endpoint.

Across all ambient patterns from `{2,3,4}` of length at most four, no cover
value equals `55/91`.  The closest killing value is `4/7`, with gap `3/91`;
the closest surviving value is `13/21`, with gap `4/273`.  Local
support-aware loads `(6)` can in principle meet the equality seam, which is
why the topological upgrade in `(9)` is retained.

## 4. Discrete \(q\)-fibre status is a distinct screen

The common-`u` event above must not be conflated with the divisor-fibre
pushforward.  For `d,q|D`, put

```text
ell=ceil(d/7), g=gcd(d,q), ell=A g+r, 0<=r<g,
H=D/lcm(d,q).
```

THM-2928's exact quotient/remainder law says that an enlarged
denominator-`d` needle has load

```text
H(A+X_i(b))
```

in the `q`-fibre over `b`, where `X_i` is true on exactly

```text
R_i=q r/g
```

fibres.  For four masks, sum the constant baselines.  At any target load
threshold `T`, the fibres capable of paying `T` form an upward event in the
nonconstant `X_i`.  The unnormalized fractional-cover theorem gives

```text
maximum real-relaxed heavy fibres
 =min(q, min_w sum_i R_i w_i).                         (10)
```

If the exact support histogram has more `T`-heavy fibres than `(10)`, the
denominator profile is impossible.  Unlike `(9)`, this finite count has no
compact/open strictness upgrade: equality in `(10)` is retained.

The bounded profile scout aggregates a denominator multiset by

```text
(ambient capacity,
 sum of q-fibre baselines,
 multiset of (high increment H, exact high count R_i)).
```

This is sufficient to reconstruct every threshold in `(10)`.  Different
`q` screens are currently tested independently; taking their rowwise minimum
is an honest upper bound on the all-`q` intersection, not an exact
intersection count.

Two complete bounded controls are already terminal after the common-phase
compression:

```text
D=56:  54 no-small shapes, 15 rows;
       scalar leaves 14 shapes / 110 occurrences;
       q=8 (local modulus 7) leaves 0.

D=84: 246 no-small shapes, 65 rows;
       scalar leaves 20 shapes / 20 occurrences;
       q=4 (local modulus 21) leaves 0.
```

Thus every `D<=100` occurrence surviving the support-aware common-phase
screen is killed: that screen leaves only `110+20=130` no-small occurrences,
and the two discrete quotient tests remove all of them.

## 5. Full aggregate census and replay

The full ordinary and optimized runs are recorded at:

```text
.scratch/lrc14_k3_four_drift_divisor_status_gf.ordinary.out
.scratch/lrc14_k3_four_drift_divisor_status_gf.opt.out
```

The full ordinary and optimized outputs are byte-identical.

The ordinary full census gives:

```text
stage             occurrences      shapes      rows   bodies  divisors
raw                21,357,714,101  694,921,995  26,970  3,003   217
scalar             13,796,122,287  694,631,606  26,954  3,003   213
ambient status     13,285,494,474  694,279,329  19,798  3,003   192
support status     13,280,722,299  694,254,050  18,599  3,003   186
```

The semantic digests of the row/pattern aggregate consequence objects are:

```text
raw             e5bf013929a772c94c62605fa42900360308bc660edc6b0cf7915396ed93f5a8
scalar          a13036637365abd6e27a56ef3ede32b0933e2cc65666472a1490a6b9b69d888a
ambient status  9f642b2c4332078af4b545c081d5c07093f71a414db76949c3e51cffb148f549
support status  01267d3e6e359eb632727b63dfa463a067f3534d1c957364c014fb9dd33c2982
```

Thus scalar capacity removes `7,561,591,814` occurrences, fractional-cover
activity removes another `510,627,813`, and the exact small-class sidecar
removes another `4,772,175`.  The full compression removes
`8,076,991,802`, or `37.8177%`, of the raw occurrence ledger.  This is a
material aggregate reduction, not an emptiness result.

Reproduction:

```bash
python3 .scratch/lrc14_k3_four_drift_divisor_status_gf.py \
  > .scratch/lrc14_k3_four_drift_divisor_status_gf.ordinary.out
python3 -O .scratch/lrc14_k3_four_drift_divisor_status_gf.py \
  > .scratch/lrc14_k3_four_drift_divisor_status_gf.opt.out
cmp .scratch/lrc14_k3_four_drift_divisor_status_gf.ordinary.out \
    .scratch/lrc14_k3_four_drift_divisor_status_gf.opt.out
```

The current-source replays took `655.77` and `651.30` wall seconds under
concurrent machine load.  Their source and common output SHA-256 are

```text
source  77c4271698350a67ed26aa32c27b07d5e31f24d56d143a3c416ee80dbe54776c
output  f60cf5039485bf6721e0b84a4ddb5399d44e78e29cfc7725cbdb1bb1a87e3440
```

The bounded `q`-profile scout and its byte-identical ordinary/optimized
outputs have hashes

```text
source  7a2ff7c9bee0974cc40adb14326b7e0fdc775151b4f1f14c6c1af36d6463765c
D=56   457a6d58bb1d364cc093ffe8a7a39663da2ef70931b2923396e8ed81b4779124
D=84   44ce06552b2bf9875b411d3a8c3944a82cf0de2ade0859b2da9dd2b7c09a229f
```

The independent checks retained in both modes are:

1. exact body/divisor universe and support cutoff `(1)`;
2. weighted generating-function equality against brute force for every
   `D<=60`;
3. unweighted Möbius totals for every resolving denominator;
4. exact fractional-cover vertex enumeration over `Q`;
5. union/intersection, equality, post-equality, and closest-margin controls;
6. monotonicity of raw, scalar, ambient-status, and support-status ledgers;
7. a literal four-needle positive control through every `q|84`;
8. no-small profile totals independently reconstructed by Möbius inversion.

## 6. Loss ledger and next bottleneck

The common-phase compression preserves denominator lcm, small-denominator
multiplicity, exact one-marginals, the threshold predicate, and (in its
strong form) the largest `S_D` residue load for `d=2,3,4`.  It forgets
numerators and all joint phase geometry, deliberately granting arbitrary
coupling.

The `q`-profile compression restores exact quotient baselines, increments,
and high counts at one chosen quotient, but still forgets numerator/AP
locations and correlations between different `q`.  A survivor is therefore
only a survivor of upper relaxations.

The structural stopping signal is now exact.  Of the
`13,280,722,299` support-status survivors,

```text
12,852,450,428 = 96.7752%
```

have no denominator in `{2,3,4}`.  Another `428,271,871` contain a small
denominator only because the other three masks already receive enough
capacity in this relaxation; common-phase activity cannot charge an
inessential label.  Progress beyond that wall must use the discrete
divisor-fibre tree `(10)`, simultaneous multi-`q` profiles, local AP
locations, or the physical ordered-drift suffix bank from THM-2941.
Repeating a common-phase marginal calculation cannot see further.

The aggregate denominator GF does **not** itself impose THM-2941's current
ordered physical cap `z_1<=378`: `z_1` is the smallest actual ordered drift
speed, not a denominator coordinate.  That cap is an orthogonal physical
sidecar to be rejoined only after the quotient relaxation has been pushed
as far as possible.
