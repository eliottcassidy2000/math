# Two aligned / five drift compression for LRC(14)

**Status:** `FINITE-EXACT`; the hostile census `D<=100` is empty.  A general
expected-spike/Mobius screen followed by an exact one-spike refinement reduces
the full `k=2` quotient ledger from `951,545,890,235` to
`317,793,841,816` occurrences.  Retaining the deterministic floor and every
optional spike bit then reduces it further to **`200,389,247,292`**
occurrences.  This does **not** close the global `k=2` sector or LRC(14).

## Inheritance and concept board

- **Closest proved mechanism:** THM-2928's support-transfer ladder, enlarged
  unit-needle fibre law `(37o)--(37p)`, and upward-status/fractional-cover
  theorem.
- **Hostile example:** the raw `k=2,p=5` universe contains
  `50,874,159,718` denominator five-multisets and
  `951,545,890,235` body-row occurrences.  Literal tuple enumeration is the
  wrong object.
- **Corrected near miss:** a phase-free scalar or largest-fibre capacity is
  only a necessary relaxation.  It loses the locations of exceptional fibres,
  the common arithmetic status table, and the local unit-AP shape.
- **Least-used sidecar:** the strict topology of the common aligned phase.
  The aligned safe carrier is compact while every proper activity threshold
  event is open.  This makes equality at the safe-floor mass impossible.
- **Live concepts:** support transfer; common-`u` activity alphabet;
  fractional blocker; divisor-Mobius GF; discrete `q`-profile; located local
  needle.

## 1. Support reduction and the exact activity alphabet

Let `F` be one of the `3003` literal six-body supports, let `D|L_F`, and let
`S_D` be its projected safe-cell support.  With two aligned tails and five
drifts, THM-2928 gives

```text
mu(R_A) >= u_2 = 66/91,
mu(union of five drift dangers) <= 1-u_5,
u_5 = 478/1365.
```

Consequently every cover requires

```text
|S_D|/D <= (1-u_5)/u_2 = 887/990.                  (1)
```

Exactly `27,163` body/divisor rows on `219` resolving denominators pass `(1)`.

For a drift denominator `d<7`, the enlarged mask is either empty or one
residue class modulo `d`, with common-phase activity event

```text
||cu|| < d/14,
```

of Haar mass `d/7`.  Since

```text
5/7 = 65/91 < 66/91 < 6/7,
```

the exact alphabet ruled out by the one-marginal compact/open test is

```text
d in {2,3,4,5}.                                    (2)
```

This is the next rung of the general one-marginal stopping rule

```text
d/7 <= u_k.
```

For such `d`, no single activity set can contain the whole compact carrier
`R_A`.  Denominators `d>=6` are merely **not excluded by this test**; they are
initially granted their full ambient enlarged capacity

```text
C_D(d)=(D/d)ceil(d/7).                              (3)
```

## 2. Exact multiplicity-pattern blocker

Let

```text
m=(m_2,m_3,m_4,m_5),       sum m_d <=5
```

be the multiplicities of the small denominators.  Label the corresponding
small masks by `i`.  Give label `i` load `a_i` and exact one-marginal

```text
r_i=d_i/7.
```

Two sound choices of load are used:

```text
ambient: a_i=D/d_i;
support-aware: a_i=max_(b mod d_i)|S_D intersection {x=b mod d_i}|.
```

If the `d>=6` masks have total ambient capacity `C_large`, coverage for every
`u in R_A` forces the upward threshold event

```text
A_T={x in {0,1}^s: sum_i a_i x_i >= T},
T=|S_D|-C_large.                                    (4)
```

Let `H_T` be the clutter of inclusion-minimal true states of `A_T`.  The
upward-status theorem gives the exact real relaxation

```text
max Pr(A_T)
 =min(1, min {sum_i r_i w_i:
              w_i>=0, sum_(i in H)w_i>=1 for all H in H_T}).  (5)
```

The inequality needed here is **strict**:

```text
max Pr(A_T) > 66/91.                                (6)
```

Indeed, the actual `R_A` is compact and is contained in the open activity
event `A_T`.  If that event is proper, compact containment gives positive
measure room; equality in `(6)` cannot contain a carrier of mass at least
`66/91`.  If the event is the whole circle, its mass is already one.

Define `beta_D(m;a)` as the largest attainable active load for which the
left side of `(6)` exceeds `66/91`.  Then every cover satisfies

```text
C_large + beta_D(m;a) >= |S_D|.                    (7)
```

For ambient weights, `beta_D/D` depends only on `m`; call it `beta(m)`.
There are `126` patterns.  Exactly seven have `beta=0`:

```text
(0,0,0,0);
(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1);
(2,0,0,0),(1,1,0,0).                               (8)
```

The pure-type thresholds are:

| denominator | copies `0,1,2,3,4,5` | `beta` |
|---|---|---|
| `2` | `0,1,2,3,4,5` | `0,0,0,1/2,1/2,1/2` |
| `3` | `0,1,2,3,4,5` | `0,0,1/3,1/3,2/3,2/3` |
| `4` | `0,1,2,3,4,5` | `0,0,1/4,1/2,3/4,3/4` |
| `5` | `0,1,2,3,4,5` | `0,0,1/5,2/5,3/5,4/5` |

The exact script generates and scale-checks all `126` mixed thresholds.  Their
value histogram is:

```text
0:7, 1/5:6, 1/4:5, 1/3:5, 2/5:6, 9/20:5,
1/2:11, 8/15:5, 7/12:5, 3/5:6, 13/20:5,
2/3:3, 7/10:10, 11/15:5, 3/4:10, 47/60:5,
4/5:4, 5/6:10, 17/20:3, 13/15:2, 9/10:4,
11/12:1, 14/15:1, 19/20:2.                         (9)
```

This exact table is the first scalable replacement for iterating the five
small drift phases.

## 3. Divisor-poset generating function

For `E|D`, group every divisor `d>1` of `E` by the feature

```text
(1_(d=2),1_(d=3),1_(d=4),1_(d=5),
 C_D(d)1_(d>=6)).
```

Equivalently, truncate at degree five the multiset series

```text
G_(D,E)(y,z_2,z_3,z_4,z_5,x)
 = product_(d|E,d>1)
   (1-y z_2^(1_(d=2)) ... z_5^(1_(d=5))
        x^(C_D(d)1_(d>=6)))^(-1).                   (10)
```

If `N_D(m,c)` counts five-multisets with lcm exactly `D`, small pattern `m`,
and large capacity `c`, divisor-lattice inversion gives

```text
N_D(m,c)
 =sum_(E|D) mu(D/E)
   [y^5 z^m x^c]G_(D,E).                            (11)
```

The coefficient within each feature class of alphabet size `g` is
`binom(g+j-1,j)`, so `(11)` propagates only proof-relevant states.  Summing
over `(m,c)` recovers independently

```text
sum_(E|D) mu(D/E) binom(tau(E)+3,5),
```

and the implementation brute-compares every weighted coefficient against
literal multisets for all `D<=60`.

## 4. Exact hostile census through `D<=100`

Only seven resolving denominators occur:

```text
D=14,28,42,56,70,84,98,
row counts=1,6,14,15,23,65,23.                     (12)
```

The exact ledgers are:

| relaxation | occurrences | shapes | rows | bodies | divisors |
|---|---:|---:|---:|---:|---:|
| raw lcm-exact | 174,448 | 3,680 | 147 | 78 | 7 |
| scalar capacity | 117,415 | 3,369 | 147 | 78 | 7 |
| ambient status | 9,972 | 1,618 | 110 | 73 | 5 |
| support-aware status | 6,089 | 1,097 | 53 | 38 | 5 |
| support status + all-`q` largest-fibre cap | 5,059 | 874 | 53 | 38 | 5 |
| support status + full all-`q` two-level profile | **500** | **120** | **21** | **15** | **3** |

Thus the common-`u` blocker removes `94.8%` of scalar survivors in this
bounded hostile universe.  Every shape containing `d=5` is eliminated.

The phase-free largest-fibre cap uses, for each `q|D`,

```text
ell=ceil(d/7), g=gcd(d,q), H=D/lcm(d,q),
cap_q(d)=H ceil(ell/g).                              (13)
```

It removes another `1,030` occurrences but no additional row.  The exact
two-level profile retains

```text
ell=A g+r,
baseline=H A,
increment=H,
high-fibre count=(q/g)r.                            (14)
```

For every target load threshold, `(5)` with marginals
`((q/g)r)/q` bounds the number of simultaneously heavy fibres.  This cuts

```text
6,089 -> 500,
1,097 -> 120,
53 -> 21.                                           (15)
```

The final divisor occurrence counts are

```text
D=28:3, D=56:401, D=84:96.                          (16)
```

The `D=42,70` ledgers are empty.  Among the final `500`, only `60` have no
small denominator; the full profile, rather than the scalar fibre cap, is
what removes the earlier no-small bulk.

## 5. Equality wall, integer audit, and missing coordinate

There are two distinct equality rules:

1. **Common phase `u`:** equality with `66/91` is impossible because a compact
   carrier must sit inside a proper open status event.  This justifies the
   strict `>` in `(6)`.
2. **Discrete `q` fibres:** equality with the fractional-cover value is
   allowed by the real counting relaxation.  It may still fail integer
   realizability, but topology supplies no strict gap.

After thresholds already covered by the baseline are omitted, every one of
the `500` survivors has zero slack at a nontrivial `q`-profile threshold.  An
exact integer feasibility model was therefore run before invoking AP
location.  Its variables are the counts

```text
n_x in Z_(>=0), x in {0,1}^s,
```

with total `q`, the exact high-fibre marginals, and every nested target-heavy
count imposed simultaneously.  All `500` real-profile survivors are integer
feasible at every `q|D`: flooring was not the remaining obstruction in this
bounded ledger.

The minimal equality survivor is

```text
F=(1,2,3,4,6,12), L_F=168, D=28, |S_D|=22,
(d_1,...,d_5)=(4,4,4,4,e), e in {7,14,28}.          (17)
```

At `q=4`, the baseline is one and the four denominator-four traces each have
increment seven on exactly one fibre.  The target has two fibres of load at
least six and all four of load at least five.  Equality forces the integer
law

```text
n_0001=n_0010=n_0100=n_1000=1,                     (18)
```

all other status counts zero: the four exceptional classes partition
`Z/4Z`.  This law is locally realizable, so the audit stops honestly at
located lifted unit-AP compatibility.

At this intermediate stage the `500`-occurrence ledger forgets:

- one **common integer status table** across all target thresholds and all
  divisors `q`;
- the **locations** of each exceptional fibre, not only their marginals;
- compatibility of those locations under quotient maps in the divisor poset;
- numerator/phase coupling and the terminal **unit arithmetic-progression
  shape** inside a selected fibre.

The next two sections supply exactly the missing common-`u` coupling.

## 6. Composing common phase with every `q` profile

For a fixed one of the `500` occurrences, let `E` be the subset of
denominators in `{2,3,4,5}` that are active at the common phase `u`.  Remove
the whole `q`-profile of every inactive small mask and grant all large masks
their complete enlarged profiles.  Call `E` feasible when this same subset
passes every target threshold for every `q|D`.  Feasibility is upward.

Thus the feasible subsets form one upward event `A` on the common phase, not
one separately selected event for each quotient.  Coverage forces

```text
R_A subset {u: active-small-set(u) in A}.            (19)
```

Apply the exact fractional-cover theorem with small-mask marginals `d_i/7`
and the compact/open strict inequality.  The exact result is

```text
500 occurrences / 120 shapes / 21 rows
 -> 60 occurrences / 4 shapes / 15 rows.            (20)
```

The kills split as

```text
D=28: 3,  D=56: 341,  D=84: 96.
```

Every occurrence containing a denominator at most five dies.  The four
surviving shapes are

```text
(7,8,8,8,8),
(8,8,8,8,14),
(8,8,8,8,28),
(8,8,8,8,56),                                      (21)
```

each on the same fifteen `D=56` rows.  The semantic digest of `(20)` is
`be9d41336a4be0a62ef980305b6260c8e68184fa025198c5f0251013251fdc03`.

## 7. The `D=56` exact-mean/eight-spike terminal

Put `q=8`.  In every shape `(21)`, the fifth mask contributes at most one
point uniformly to every `q`-fibre.  Each denominator-eight mask hits
`Y_i(u)` fibres.  Exact averaging gives

```text
integral_T Y_i(u)du=8/7.                             (22)
```

All fifteen target histograms have eight fibres of load greater than one:

```text
support 48: (6^8),                                  ten rows;
support 44: (5^4,6^4),                              four rows;
support 40: (4^2,5^4,6^2),                          one row.  (23)
```

Consequently coverage forces

```text
sum_(i=1)^4 Y_i(u)>=8
```

throughout `R_A`.  This is a proper open event.  Markov and `(22)` give

```text
mu{sum_i Y_i>=8} <= (4*(8/7))/8 =4/7<66/91.         (24)
```

This contradiction closes all `60`.

There is an independent tooth-level explanation.  A denominator-eight mask
has one guaranteed seven-point `q=8` spike and a second spike exactly when

```text
{cu} in (3/7,4/7),                                  (25)
```

an open event of mass `1/7`.  The fifth mask gives capacity one in every
fibre, and all eight target fibres need a spike.  Hence the four base spikes
and all four extra spikes are required; `(25)` again contradicts the aligned
safe floor.  Exact controls include the positive state `c=1,u=1/2`, whose
two spike classes are `{0,7}`.

Therefore

```text
the complete k=2 quotient ledger with D<=100 is empty. (26)
```

## 8. Full expected-spike Mobius compression

In fact every one of the `27,163` support-transfer rows has `7|D`.  Put
`q=D/7`.  A denominator `d|D` has one of two roles:

- if `d` does not divide `q`, its enlarged trace contributes one uniformly
  to every `q`-fibre;
- if `d|q`, its literal mask hits `Y_d(u)` `q`-fibres, with the exact mean

```text
integral_T Y_d(u)du
 =(q/d) integral_T #{active classes mod d}du
 =(q/d)(d/7)=q/7.                                  (27)
```

Suppose `c` of the five masks are uniform and `m=5-c` are spike-type.  Let

```text
N_c=#{b mod q: lambda_q(b)>c}.
```

Coverage forces `sum_i Y_i>=N_c` on `R_A`.  If `N_c>0`, Markov plus
compact/open strictness gives the exact necessary tariff

```text
N_c < m q/(7u_2)=13m q/66.                          (28)
```

This needs only the number `c`, so the full denominator universe again has a
degree-five exact-lcm generating function.  For `E|D`, let

```text
U_D(E)=#{d|E:d>1,d not dividing q},
V_D(E)=tau(E)-1-U_D(E).
```

Equivalently, if `D=7^a q_0` with `7` not dividing `q_0`, then

```text
U_D(E)=tau(E/7^a) if v_7(E)=a, and 0 otherwise.
```

Writing `Mult(n,j)=binom(n+j-1,j)`, with `Mult(n,0)=1`, the exact number of
lcm-`D` five-multisets with `c` uniform symbols is

```text
A_(5,c)(D)
 =sum_(E|D)mu(D/E)Mult(U_D(E),c)Mult(V_D(E),5-c).   (29)
```

Summing `(29)` over `c` recovers `(37lM)` exactly; `A_(5,0)(D)=0`, since a
full septimal exponent must occur.

Applying `(28)` rowwise gives

| ledger | shapes | occurrences | rows | bodies | divisors |
|---|---:|---:|---:|---:|---:|
| raw | 50,874,159,718 | 951,545,890,235 | 27,163 | 3,003 | 219 |
| expected-spike survivors | 36,962,285,549 | 320,011,786,356 | 4,592 | 2,999 | 149 |

Thus `(28)` kills `631,534,103,879` occurrences and `22,571` rows without
enumerating a denominator tuple.  The surviving occurrences by `c` are

```text
c=1:  13,975,689,268
c=2: 143,128,657,438
c=3: 119,436,903,764
c=4:  41,980,227,729
c=5:   1,490,308,157.                              (30)
```

The smallest survivor is

```text
D=168, q=24, F=(1,2,3,4,6,12), |S_D|=88,
c=4, m=1, N_4=4.                                   (31)
```

Its `1,605` denominator shapes already reveal the next coordinate.  If the
lone spike denominator is `d=2,3,4,6,8,12,24`, the exact masses of
`{Y_d>=4}` are respectively

```text
2/7, 3/7, 4/7, 6/7, 1/7, 5/7, 3/7.                (32)
```

Only `d=6` survives the aligned floor; it accounts for `204` of the `1,605`
shapes.  The next scalable quotient should therefore retain the two-level
hit law of each spike type, rather than only its expectation.  Located AP
phases and compatibility across quotient scales remain later coordinates.

## 9. Global exact `c=4` one-spike law

The suggested coordinate scales exactly.  In the `c=4` stratum there is a
unique denominator `d|q`; the other four masks contribute one point to every
`q`-fibre.  Write

```text
d=7a+r, 0<=r<7.
```

Because the `d` residue centres are equally spaced, the number of active
classes at common phase `u` is exactly `a+X_d(u)`, where `X_d` is Boolean and

```text
mu{X_d=1}=r/7.
```

Thus the number of hit `q`-fibres is the two-level variable

```text
Y_d(u)=(q/d)(a+X_d(u)).                            (33)
```

If `N_4=#{b:lambda_q(b)>4}`, coverage requires `Y_d>=N_4` throughout the
compact aligned-safe carrier.  There are only three cases:

```text
N_4 <= (q/d)a:       event mass 1;
(q/d)a < N_4 <= (q/d)(a+1): event mass r/7;
N_4 > (q/d)(a+1):    event mass 0.                 (34)
```

The middle event is proper open, so its mass must be **strictly** greater
than `u_2=66/91`.  Since `r` is integral, only `r=6` passes.  This is an exact
one-mask law, not a Markov estimate.

Retaining `d` does not require tuple enumeration.  If `U_D(E)` is the number
of divisors of `E` that do not divide `q`, then the number of exact-lcm-`D`
five-multisets with four uniform symbols and lone spike `d` is

```text
A_(4,d)(D)
 =sum_(E|D,d|E) mu(D/E) Mult(U_D(E),4),  d|q.      (35)
```

Summing `(35)` over `d|q,d>1` recovers the complete `c=4` coefficient of
`(29)`.  A literal multiset comparison holds for every septimal `D<=100`, and
an independent divisor-downset recurrence reproduces every weighted
coefficient for all `219` resolving moduli.

The exact global result is:

| ledger | shapes | occurrences | rows | bodies | divisors |
|---|---:|---:|---:|---:|---:|
| incoming `c=4` | 3,246,262,178 | 41,980,227,729 | 3,609 | 2,940 | 134 |
| exact one-spike survivors | 3,208,397,602 | 39,762,283,189 | 3,373 | 2,871 | 133 |

The refinement kills `2,217,944,540` occurrences and removes `D=588`
entirely from this stratum.  Of the survivors, `39,336,054,399` are already
covered at the lower level of `(33)`; the remaining `426,228,790` use the
proper `r=6` extra-tooth event.  Leaving the other four `c`-strata unchanged,
the full relaxation becomes

```text
36,924,420,973 shapes / 317,793,841,816 occurrences. (36)
```

The minimal killed group is `(31)` with lone `d=2`: its lower/upper hit
counts are `0,12`, and the needed event has mass `2/7`; its coefficient is
`184`.  The minimal surviving group has the same `D`, body, support and
`N_4=4`, but lone `d=6`: its hit counts are `0,4`, its proper event has mass
`6/7`, and its coefficient is `204`.  This freezes the boundary mechanism:
the next useful refinement must couple two or more spike events, or retain
their locations, because a lone `r=6` event genuinely has enough measure.

Frozen audit hashes:

```text
full_q7_expected_spike_gf.py:
  3e26985308834c4c0dcece6f17a214fa7622c981a9a0e9f085930512347f39b7
ordinary/-O output:
  7ce4080f899f6e319b24977fc3b4f13a30290452f2b8113bb8ddb1f60e136d30
full_q7_one_spike_refinement_gf.py:
  8f655eb8d398cf309f4982392eb761c7ce8d45409ca41cd03bcfae16390bcd99
ordinary/-O output:
  05e123c4c01beebfb161de8a87d59c2507ffaaa3cfa2512a308e14a22bb33a38
survivor semantic digest:
  14c9e777f98ca3d44ec3747eeb08153f4edffd23e0d74ff8f4abf4f5195440fa
```

## 10. All-arity floor/exception blocker

### Why every support-hard modulus is septimal

The fact `7|D` has a short structural proof.  Write

```text
L=14 lcm(F),  M=L/D,
```

and suppose `7` does not divide `D`.  In the fibre `j=r+Dt`, a body speed
`f` has period `P=L/f`.  Since `P` is divisible by seven while `D` is not,
the orbit size

```text
T=P/gcd(P,D)
```

is divisible by seven.  The half-open danger block has length `P/7`, and
`gcd(P,D)` divides `P/7`; hence exactly `T/7` points of each orbit lie in
that block.  The orbit repeats `M/T` times in the fibre, so speed `f` removes
exactly `M/7` of its cells.  Six body speeds remove at most `6M/7`, leaving a
safe cell in every residue class.  Thus

```text
7 not dividing D  ==>  S_D=Z/DZ.                    (37q7a)
```

The exact census checks `(37q7a)` on all `96,235/96,235` nonseptimal
body/divisor rows.  Since the support-transfer cutoff is strictly below one,
every support-hard row is therefore septimal for a reason, not by accident.

### Exact floor plus optional-bit law

Fix `q=D/7`.  Suppose `c` denominator symbols do not divide `q`; they each
pay one point in every `q`-fibre.  For each of the remaining `m=5-c` symbols
write

```text
d_i=7a_i+r_i,  w_i=q/d_i,  0<=r_i<7.
```

At common phase `u`, the exact number of fibres hit by that mask is

```text
Y_i(u)=a_i w_i+w_i X_i(u),
Pr(X_i=1)=r_i/7.                                    (37q7b)
```

Bits with `r_i=0` are absent.  Put

```text
B=sum_i a_i w_i,
N_c=#{b mod q:lambda_q(b)>c}.
```

Coverage forces the upward event

```text
sum_i w_i X_i >= N_c-B                              (37q7c)
```

throughout the compact aligned-safe carrier.  Let `beta_u2(w,r)` be the
largest subset-sum threshold whose exact fractional-cover maximum is
strictly greater than `u_2=66/91`.  The complete necessary condition is

```text
N_c <= B+beta_u2(w,r).                              (37q7d)
```

This simultaneously extends `(28)` and the one-spike law `(33)--(34)`.  It
is stronger than the mean screen: the exact audit finds no survivor of
`(37q7d)` rejected by `(28)`.

### Exact-lcm completion GF

There is still no need to enumerate a full denominator five-tuple.  For a
positive-residue spike multiset `P`, let `s=|P|`, `ell=lcm(P)`, and let `z`
be the number of `r=0` spike symbols, so `c+z+s=5`.  In a divisor downset
`E|D`, write

```text
U_D(E)=#{d|E:d>1,d not dividing q},
Z_D(E)=#{d|gcd(E,q):d>1, 7|d}.
```

The exact number of uniform/zero-residue completions of `P` is

```text
K_(D,c,z)(ell)
 =sum_(ell|E|D) mu(D/E)
    Mult(U_D(E),c) Mult(Z_D(E),z).                  (37q7e)
```

Every such completion has allowance

```text
A(P,z)=zq/7
       +sum_(d in P)(q/d)floor(d/7)
       +beta_u2((q/d)_(d in P),(d mod 7)/7).        (37q7f)
```

Grouping `(37q7e)` by `(c,A)` gives an exact suffix query for every row.
Summing over `A` recovers the full c-refined Mobius coefficient for every
one of the `219` resolving moduli.

The global ledgers are:

| stage | shapes | occurrences | rows | bodies | divisors |
|---|---:|---:|---:|---:|---:|
| raw support-transfer | 50,874,159,718 | 951,545,890,235 | 27,163 | 3,003 | 219 |
| expected-spike `(28)` | 36,962,285,549 | 320,011,786,356 | 4,592 | 2,999 | 149 |
| exact floor/exception `(37q7d)` | **26,908,162,790** | **200,389,247,292** | **4,414** | **2,977** | **149** |

Thus the exact blocker removes another `10,054,122,759` shapes and
`119,622,539,064` occurrences, respectively `27.20%` and `37.38%` of the
mean survivors.  Relative to the raw quotient ledger, `78.94%` of all
occurrences are now gone.

The arity split makes the leverage visible:

| `c` | exact shapes | exact occurrences | rows | minimum surviving `D` |
|---:|---:|---:|---:|---:|
| 1 | 3,457,885 | 3,524,756 | 45 | 13,860 |
| 2 | 12,908,211,669 | 46,320,209,782 | 779 | 840 |
| 3 | 10,562,966,029 | 112,812,921,408 | 2,425 | 336 |
| 4 | 3,208,397,602 | 39,762,283,189 | 3,373 | 168 |
| 5 | 225,129,605 | 1,490,308,157 | 2,409 | 392 |

In particular, `c=1` falls from `13,975,689,268` mean-surviving
occurrences to only `3,524,756`.  The minimum global survivor is still the
row `(31)` with lone `d=6`, allowance `4`, and coefficient `204`; the all-c
engine independently reproduces the entire one-spike `c=4` census.

### Audit and remaining coordinate

Two exact engine builds (`-O2` and `-O3`) are byte-identical.  Literal
exact-lcm enumeration through `D<=300` checks `52,925` five-multisets and
`1,880` independent suffix queries.  Ordinary and `python3 -O` referee runs
are byte-identical.  Frozen hashes are:

```text
full_q7_all_spike_floor_exception_gf.py:
  4f42009b972a552b112e9523c4c7bbedbe41430d91c20ce1cd005e426221fc13
full_q7_all_spike_floor_exception_engine.cpp:
  664d0df36d104d959279605c8ea8539d61ab595b155e5157fa7d0433f1b7944c
ordinary/-O output:
  f711376bdfa0064f70d76e42505cf1eb89dadc4c66bcacd90d985b6641c2cd75
engine query output:
  49cbf8fb160a78125d426cf7348c3c35e4859f110f493bb4b03acfc2f0c92125
survivor semantic digest:
  2eec9a97f02a7b8f8e36e50f747d53186ff5b84234a14fc3ac64818b54033675
```

The refinement still optimizes an arbitrary joint law of the optional bits.
It forgets their actual common-phase interval locations, their locations
inside the `q` fibres, compatibility across quotient scales, and numerator/
lifted-AP coupling.  These are now the honest missing coordinates.

For one bit, the extra event itself is a translated wide comb: if `d=7a+r`
and the reduced numerator is `h`, then

```text
X_d=1 iff {hu} lies in the open interval centered at -a/2
         with radius r/14.                            (37q7g)
```

The odd-`a`, `r=6` branch is genuinely hostile: `A=(h,14h)` gives a literal
containment already at `(h,a)=(1,1)`, so a measure-only kill cannot work.
The even-`a` center-zero branch remains `OPEN`.  A tempting component proof
was rejected: its subclaim that no `D_b` can contain one component of the
half-shifted narrow complement is false (`H_3`'s middle component lies inside
the central tooth of `D_2`).  Exact scouts found no full center-zero
containment in the tested window, but no theorem is claimed.

## Reproduction

```bash
python3 .scratch/lrc14_k2_five_drift_pattern_compression_20260729/scout.py --max-D 100
python3 .scratch/lrc14_k2_five_drift_pattern_compression_20260729/qfiber_hostile.py --max-D 100
python3 .scratch/lrc14_k2_five_drift_pattern_compression_20260729/qprofile_hostile.py --max-D 100
python3 .scratch/lrc14_k2_five_drift_pattern_compression_20260729/composed_common_u_qprofile.py --max-D 100
python3 .scratch/lrc14_k2_five_drift_pattern_compression_20260729/d56_eight_spike_terminal.py
python3 .scratch/lrc14_k2_five_drift_pattern_compression_20260729/full_q7_expected_spike_gf.py
python3 .scratch/lrc14_k2_five_drift_pattern_compression_20260729/full_q7_one_spike_refinement_gf.py
python3 .scratch/lrc14_k2_five_drift_pattern_compression_20260729/full_q7_all_spike_floor_exception_gf.py
```

Each program uses `RuntimeError` guards, not Python `assert`; ordinary and
optimized runs are required to be byte-identical before handoff.
