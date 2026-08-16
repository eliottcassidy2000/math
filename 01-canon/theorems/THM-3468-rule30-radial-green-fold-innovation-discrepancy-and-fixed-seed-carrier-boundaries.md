---
id: THM-3468
title: "Rule 30 radial Green fold, innovation discrepancy, and fixed-seed carrier boundaries"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A centered Green
  restart exposes a canonical valuation-parity backbone and folds every later
  collision into radial source/transport grades.  Every finite correction
  band is eventually 2-power-periodic, while the backbone plus any finite band
  is nonperiodic; hypothetical global cancellation must escape all finite
  bands.  The innovation phase atlas extends to a Haar--Bernoulli coordinate
  homeomorphism, while an exact zero-intersection current controls spatial
  bias and a Baire hostile blocks average-to-origin transfer.  Fixed-seed
  query, Hankel, and binary-address ranks give finite exact restricted-carrier
  lower bounds.  None of the three Rule 30 prizes is claimed.
source: root-rule30-crossfrontier-20260815
audit: >
  Independent centered-restart, radial-fold, valuation-clock, two-slack,
  period-tariff, innovation-atlas, transition-current, discrepancy, Baire,
  query/Hankel/digit-carrier, model-scope, hash, and documentation audits;
  ordinary and optimized executions reproduce the stored transcript exactly.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
related:
  - THM-2227-sharp-parity-three-checkpoint-bellman-profile-exclusion
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3459-rule30-ternary-intersection-factorial-truth-lift-and-keller-boundaries
  - THM-3466-factorial-face-stokes-and-keller-boundary-current
external:
  - Stephen Wolfram, "Announcing the Rule 30 Prizes", https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/ (2019; CITED for the problem statements only)
  - Wolfram Rule 30 Prizes, https://rule30prize.org/ (CURRENT OFFICIAL LISTING checked 2026-08-15; active listing and submission status only)
script: 04-computation/rule30_radial_innovation_carriers_thm3468.py
output: 05-knowledge/results/rule30_radial_innovation_carriers_thm3468.out
script_sha256: dd0ac8c4e4a298e665b162e12e5e3333d714e04261ce00e552e1c21dc2d6ea17
output_sha256: 7ac34b56cfae942d00f6873c2dae4a65a0004775b04705970522aaf8a812430f
hash_basis: raw bytes
---

# THM-3468 -- Rule 30 radial Green fold, innovation discrepancy, and fixed-seed carrier boundaries

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem sharpens the three open Rule 30 fronts without claiming a
solution.  Its common mechanism is a filtration invoice: a nonperiodic or
balanced carrier is already present, but the distinguished center is obtained
only after an unbounded graded cancellation or a marked evaluation.

## 1. Inheritance and conventions

The closest proved inputs are THM-3456's left-permutive trace torsor,
THM-3458's right-edge odometer and moving observer, THM-3459's typed
intersection/lift boundaries, and THM-3463's Mealy, Green-current, period-lift,
innovation-atlas, and promise-complexity theorem.  The canonical hostiles are
Rule 60 for free trace structure, the width-six equal-state/different-center
collision, and THM-3463's full twelve-bit dyadic-kernel word universe.  The
corrected near misses are:

1. a fixed edge offset is periodic, but the center offset moves;
2. Haar balance in phase does not control one marked phase;
3. a universal-input or seed-neighborhood lower bound does not lower-bound the
   one fixed seed; and
4. a nonperiodic summand proves nothing until cancellation after the radial
   grading is controlled.

Let `a_t(j)` be Rule 30 from one seed,

```text
a_0(0)=1,       a_0(j)=0 for j!=0,
f(l,c,r)=l+c+r+cr                 over F_2.             (1)
```

Use centered Laurent rows and their nonlinear bond sources

```text
A_t(y)=sum_(j in Z) a_t(j)y^j,
Q(y)=y^(-1)+1+y,
e_t(j)=a_t(j)a_t(j+1),
E_t(y)=sum_j e_t(j)y^j.                                (2)
```

Then the local rule is exactly

```text
A_(t+1)=Q A_t+E_t.                                     (3)
```

For `n,q,d>=0`, define the ternary coefficient and centered Rule-150 Green
shell

```text
K(n,q)=[x^q](1+x+x^2)^n,
H_d(n)=[y^d]Q(y)^n=K(n,n+d).                           (4)
```

and put `H_d(n)=0` when `d>n`.  All additions below are in `F_2` unless a
real density or signed bias is explicitly displayed.

## 2. A canonical nonperiodic backbone after the time-two restart

The second Rule 30 row is

```text
A_2=y^(-2)+y^(-1)+y^2.                                 (5)
```

Restarting (3) at time two and taking the center coefficient gives, for every
`t>=2`,

```text
c_t:=a_t(0)
 =H_1(t-2)
  +sum_(s=2)^(t-1) sum_(j in Z)
       e_s(j) H_|j|(t-1-s).                            (6)
```

Indeed the two distance-two terms in (5) have equal Green weights and cancel;
the lone distance-one term survives.  Formula (6) is a genuine distinguished-
seed identity, not a random-input or linearized approximation.

For a radial fold, put

```text
F_s(0)=e_s(0),
F_s(d)=e_s(d)+e_s(-d)             (d>=1).              (7)
```

Palindromy of `Q^n` turns (6) into

```text
c_t=H_1(t-2)
    +sum_(s=2)^(t-1) sum_(d=0)^min(s,t-1-s)
       F_s(d)H_d(t-1-s).                               (8)
```

This folds only the value of the oriented bond indicator: `e_s(d)` is the
bond `(d,d+1)`, while `e_s(-d)` is its half-integer reflected bond
`(-d,1-d)`.  Rule 30 itself is not reflection-symmetric, and the orientation
and injection time remain load-bearing sidecars.

### 2.1 The first radial shell is valuation parity

THM-3463's Green digit recursion specializes to

```text
H_(2e)(2m)       =H_e(m),
H_(2e+1)(2m)     =0,
H_(2e)(2m+1)     =H_e(m),
H_(2e+1)(2m+1)   =H_e(m)+H_(e+1)(m).                  (9)
```

Since `H_0=1`,

```text
H_1(2m)=0,
H_1(2m+1)=1+H_1(m),                                   (10)
```

and hence

```text
boxed: H_1(n)=nu_2(n+1) mod 2.                         (11)
```

Thus the canonical Rule 30 backbone in (6) is the valuation-parity, or
period-doubling, word.  It is not eventually periodic.  If `p>0` were an
eventual period, put `r=nu_2(p)` and choose `M>r`, past the transient, with
`M` and `r` of opposite parity.  For `n=2^M-1`,

```text
H_1(n)=M mod 2,
H_1(n+p)=nu_2(2^M+p) mod 2=r mod 2,                   (12)
```

a contradiction.  Its natural density of ones is

```text
sum_(q>=0) 2^(-(2q+1)-1)=1/3.                         (13)
```

Its ordinary generating function over `F_2` obeys the Mahler equation

```text
G(z)=z/(1-z^2)+zG(z^2).                                (13a)
```

Thus (11) is simultaneously an `O(log n)` evaluator, a two-automatic
sequence, and a nonrational obstruction.

More generally,

```text
H_(2^q)(n)=H_1(floor(n/2^q)),                          (14)
```

so every power-of-two radial shell contains a block-inflated nonperiodic
valuation clock.  Equation (9), or the two-carry digit compiler for (4),
evaluates every fixed shell from the binary digits of `(n,d)`; the nonlinear
source array `F_s(d)`, not the linear transport, is the unresolved object.

### 2.2 Source grading and the raw-uniqueness boundary

For every `d>=1`, the first possible folded source satisfies

```text
F_d(d)=1.                                              (15)
```

At time `d`, the right extreme bond `(d,d+1)` meets the exterior zero, while
the leftmost two cells `(-d,-d+1)` are both one.  Thus every radial grade is
actually activated; it is not a formal zero summand.

A literal unique-event argument is nevertheless blocked.  At target `t=3`,
the unfurled THM-3463 current has distinct nonzero summands `(s,k)=(1,2)` and
`(1,3)`.  For every `t>=4`, the summands `(1,2)` and `(3,4)` both have the
all-target weight `H_0=1`.  Thus raw uniqueness occurs only at `t=2` among
nonlinear targets.  The time-two restart produces a unique **linear
backbone**, not a unique total collision event; graded noncancellation is the
viable target.

## 3. Every finite source/transport band is periodic

For a term in (8), write

```text
d=s-u,
n=t-1-s=d+v,                                           (16)
```

where `u,v>=0` are respectively source-edge depth and Green slack.  Equivalently,

```text
t=2s-u+v+1,
H_d(d+v)=[x^v](1+x+x^2)^(d+v).                        (17)
```

Fix `(u,v)`.  Write the right and left edge offsets as

```text
b_r(s)=a_s(s-r),
ell_r(s)=a_s(-s+r),       b_(-1)=ell_r=0 for r<0.      (18)
```

For `s>u`, the folded source is exactly

```text
F_s(s-u)=b_u(s)b_(u-1)(s)+ell_u(s)ell_(u+1)(s).       (18a)
```

At the single head value `s=u`, the radial distance is zero and
`F_s(0)=e_s(0)` instead.  The right offsets in (18a) are purely periodic by
THM-3458.  For the left offsets, the local rule gives

```text
ell_r(s+1)
 =ell_(r-2)(s)+ell_(r-1)(s)
  +(1+ell_(r-1)(s))ell_r(s).                           (19)
```

After the finite startup, the boundary bits satisfy `ell_0=ell_1=1`; induct
on `r`.
Once the two lower drivers have common 2-power period `P`,
the monodromy on the one bit `ell_r` is either constant, if some multiplier
in (19) vanishes, or a translation of `F_2`.  It therefore has eventual
period `P` or `2P`.  Every fixed left offset, and hence the source factor in
(18a), is eventually 2-power-periodic.  More precisely, `ell_0` and `ell_1`
have eventual period one.  If the lower drivers for `ell_r` have common
period dividing `2^(r-1)`, the same monodromy argument gives period dividing
`2^r`.  Hence `ell_r` has eventual period dividing `2^r`.  The right factors
`b_u,b_(u-1)` have common period dividing `P_(u+1)<=2^u`, so the complete
source factor in (18a) has eventual period dividing `2^(u+1)`.

Choose `2^M>v`.  Frobenius factors of
`(1+x+x^2)^(d+v)` belonging to binary places at least `M` have degree at least
`2^M` and cannot affect `[x^v]`.  The Green factor in (17) therefore depends
only on `d+v mod 2^M`.  It too is 2-power-periodic.

For a target `t`, put `s=(t-1+u-v)/2` and define

```text
J_(u,v)(t)=F_s(s-u)K(s-u+v,v)                         (19a)
```

when `s` is integral, `2<=s<t`, and `s>=u`; put it equal to zero otherwise.
Then the later-collision current in (8) is the coefficientwise finite sum

```text
R_t=sum_(u,v>=0) J_(u,v)(t).                          (19b)
```

Consequently every fixed `(u,v)` channel is eventually 2-power-periodic in
the target time.  More explicitly, its eventual period divides

```text
2^[1+max(u+1,ceil(log_2(v+1)))].                      (19c)
```

The extra factor two records the required parity class in
`t=2s-u+v+1`.

Let `B` be any finite set of `(u,v)` channels and let `R_B(t)` be their sum.
Then

```text
C_B(t)=H_1(t-2)+R_B(t)                                 (20)
```

is not eventually periodic: otherwise subtracting the eventually periodic
`R_B` would make `H_1` periodic.  Thus **every finite-band approximation to
the distinguished center is already nonperiodic**.

There is also a sharp conditional escape tariff.  Suppose the true center
were eventually periodic.  For any finite `B`, take a common eventual period
`P` of `c_t` and `R_B(t)`, and choose `M>nu_2(P)` with the opposite parity and
large enough that both `t=2^M+1` and `t+P` are past their transients.
At

```text
t=2^M+1,                                               (21)
```

the two sequences `c` and `R_B` repeat after `P`, while (12) shows that the
backbone flips.  Hence the complementary sum over channels outside `B` must
flip.  Under a periodicity hypothesis, cancellation therefore escapes every
finite rectangle or finite set of source/transport grades: at least one of
`u` or `v` must become unbounded.  No finite radial band can carry it.

### 3.1 Three forced boundary grades, not three independent spines

The first three fixed source-edge depths have closed forms.  The two leftmost
row bits are `11`, the third is zero from time two onward, while on the right
edge the packed bit recurrence gives

```text
b_0(s)=1,
b_1(s)=b_2(s)=s mod 2.                                 (21a)
```

On the zero-slack Green boundary, `K(d,0)=1` is the unique
all-one-direction ternary path.  Hence

```text
F_s(s)=1             (s>=1),       target t=2s+1,
F_s(s-1)=s mod 2     (s>=1),       target t=2s,
F_s(s-2)=s mod 2     (s>=2),       target t=2s-1.       (21b)
```

These are three exact forced source grades, but not three independent scalar
spines: the `u=0` and next-row `u=2` contributions already cancel on one odd
target class.  Likewise, for `n>=2`, the middle Green coefficient `H_0(n)=1`
is odd after many ternary ancestry words rather than coming from one word.
This is the minimal boundary between the genuine three-branch Green tree and
a Berggren ancestry interpretation.

Let `R_t` denote the later double sum in (8), as in (19b), and let
`rho_(H_1 R)` denote the density of `H_1(t-2)R_t`.  If this joint density and
`rho_R` both exist, center balance would require the exact correlation invoice

```text
1/2=1/3+rho_R-2rho_(H_1 R),
rho_R-2rho_(H_1 R)=1/6.                               (22)
```

The missing `1/6` is a signed correlation, not merely a collision density.

## 4. The innovation atlas extends to a Haar--Bernoulli coordinate chart

Retain THM-3458's edge columns and periods,

```text
b_k(t)=a_t(t-k),
P_(k+1)=2^(epsilon_k)P_k.                              (23)
```

Each fixed `b_k` is purely periodic.  Let `btilde_k:Z_2 -> F_2` be its
continuous periodic extension and define

```text
T_h(k)=btilde_k(h+k),             h in Z_2.            (24)
```

For integer `h` with `h+k>=0`, this is the physical cell `a_(h+k)(h)`.
When `h+k<0`, (24) is only the odometer-periodic extension; it asserts no
negative-time Rule 30 cell.

Let

```text
I_K={1<=k<=K:epsilon_k=1},
m_K=|I_K|,
P_(K+1)=2^m_K,                                         (25)
```

and enumerate the infinite innovation set as
`kappa_1<kappa_2<...`.  It is infinite because the periods are unbounded, and

```text
P_(kappa_r)=2^(r-1),       P_(kappa_r+1)=2^r.          (26)
```

THM-3463 proves that the compatible finite maps

```text
Gamma_K:Z/P_(K+1)Z -> F_2^(I_K),
h |-> (T_h(k))_(k in I_K)                              (27)
```

are bijections.  Reduction of phase and projection of innovation coordinates
commute as `K` grows.  Passing to inverse limits gives a homeomorphism

```text
Gamma:Z_2 -> F_2^N,
h |-> (T_h(kappa_r))_(r>=1).                           (28)
```

Every cylinder of `m` innovation bits has one preimage residue modulo `2^m`,
so `Gamma` sends Haar measure to Bernoulli product measure.  Consequently:

1. under Haar phase, `T_h(kappa_(r+1))` is a fresh fair bit independent of
   every earlier depth readout;
2. each readout before `kappa_(r+1)` is a Boolean function of the first `r`
   innovation coordinates;
3. Haar-almost every phase is normal, hence balanced, along innovation depths;
4. under uniform phase in `Z/P_(K+1)Z`, the base-two Shannon entropy of the
   full depth-`K` phase word is exactly `m_K`; and
5. the prize center is the one marked address

```text
Gamma(0)=(c_(kappa_r))_(r>=1).                         (29)
```

This is a topological and measure-space coordinate theorem, not a group or
dynamical Bernoulli conjugacy.  The finite source groups are cyclic whereas
`F_2^m` has exponent two for `m>=2`; addition by one becomes a nonlinear
odometer on the cube and still has zero entropy.  At `K=4`, with innovations
`1,3,4`, the encoded phase permutation is

```text
7,4,1,6,3,0,5,2.                                      (30)
```

Thus even an innovation readout which is linear in `Gamma` coordinates need
not be linear in ordinary binary phase bits.

## 5. The exact spatial transition current and discrepancy bound

Write `T_k(h)=T_h(k)`, with `T_k=0` for `k<0`.  The packed recurrence gives,
for `k>=2`,

```text
T_k(h+1)+T_k(h)=Q_k(h),
Q_k(h)=T_(k-1)(h+1) or T_(k-2)(h+2).                  (31)
```

For integer `h` with `h+k>=0`, both parents in (31) are physical edge bits at
the common time `h+k`.  For general `h in Z_2`, they are the corresponding
common-phase periodic extensions.  The identity holds everywhere by
continuity.  The current has period dividing `P_k`.  Put

```text
D_k(t)=b_k(t+P_k)+b_k(t).
```

Its lower-bit driver is `P_k`-periodic, so `D_k` is constant in `t`; at zero
it is `b_k(P_k)+b_k(0)=epsilon_k`.  Telescoping one lower-width cycle therefore
gives

```text
xor_(0<=h<P_k) Q_k(h)=epsilon_k,
T_k(h+P_k)=T_k(h)+epsilon_k.                           (32)
```

If `epsilon_k=1`, (32) is anti-periodicity and `T_k` is exactly balanced over
its `2P_k=P_(k+1)` phase cycle.

Assume now `epsilon_k=0` and put `p=P_k=P_(k+1)`.  Let

```text
B_k=sum_(h=0)^(p-1)(-1)^(T_k(h)).                      (33)
```

If `Q_k=0`, then `T_k` is constant, `B_k=(-1)^c_k p`, and (35) below is an
equality.  Otherwise the event count is `2r>0`.  If its cyclic gap lengths are
`ell_1,...,ell_(2r)`, then integrating (31) from the marked value
`T_k(0)=c_k` gives

```text
B_k=(-1)^c_k sum_(h<p)(-1)^(xor_(u<h)Q_k(u)),
|B_k|=|ell_1-ell_2+...-ell_(2r)|.                      (34)
```

The event locations determine the bias magnitude; the marked initial owner
`c_k` determines its sign.  Since each bit value has at least `r` nonempty
runs,

```text
|B_k|<=p-|supp Q_k|.                                  (35)
```

Define the aligned simultaneous-zero set

```text
Z_k={h:T_(k-1)(h+1)=T_(k-2)(h+2)=0},
rho_k=(1/p)|{h:T_k(h)=1}|.                             (36)
```

Because `Q_k` is the OR of the two parents, (35) is the exact useful bound

```text
boxed: |rho_k-1/2| <= |Z_k|/(2p).                     (37)
```

For aligned parent spins `U,V`, the missing statistic is

```text
|Z_k|/p=(1+E U+E V+E(UV))/4.                          (38)
```

Here `E` is uniform averaging over `h in Z/pZ`.
Thus separate parent densities do not close the bias: their aligned
zero-intersection correlation is the required sidecar.

If either parent depth is innovative, `|Z_k|<=p/2` and
`rho_k in [1/4,3/4]`.  More sharply, if

```text
epsilon_(k-2)=epsilon_(k-1)=1,                        (39)
```

put `q=P_(k-2)`, so `P_k=4q`.  At each base phase `t<q`, write the aligned
parent pair in the first block as `(A,B)` and its first coordinate in the
second block as `A'`.  The four blocks give

```text
(A,B), (A',B+1), (A+1,B), (A'+1,B+1),
```

which contain every element of `F_2^2` exactly once.  Thus the aligned parent
coordinates are uniform.  THM-3463's no-`111` law forces `epsilon_k=0`,
`|Z_k|=p/4`, and

```text
boxed: 3/8 <= rho_k <= 5/8.                           (40)
```

Every period-lift motif `110` therefore has a quantitatively near-balanced
spatial child.  The following are all such instances in the companion's
audited range `k<=30`:

```text
k       p       |supp Q_k|       wt(T_k)
5       8            6              3
8       32          24             16
17      256        192            128
26      1024       768            516.                (41)
```

The lower endpoint in (40) is attained at `k=5`.  This is a phase-density
theorem at one fixed depth; it gives no distribution statement for the marked
values `c_k=T_k(0)` as depth varies.

## 6. A dense cylinder hostile from the light-cone boundary

Fix a depth cutoff `K`, a phase cylinder `h=r mod P_(K+1)`, and a length `L`.
Choose `q>max(K,L)` with

```text
q=-r mod P_(K+1).                                     (42)
```

Then the phase `h=-q` lies in that cylinder, so its readouts through depth `K`
match those of `r`.  At depths `q<=k<=2q`, however, the time `k-q` is
nonnegative and the readout is physical:

```text
T_(-q)(k)=a_(k-q)(-q).                                (43)
```

For `q<=k<2q`, the site is outside the light cone and (43) is zero.  At
`k=2q` it is the extreme-left front and equals one.  Every phase cylinder
therefore contains a phase word with a later exact block

```text
0^q 1.                                                (44)
```

For the center cylinder take `r=0` and `q` divisible by `P_(K+1)`.  Agreement
at depths `k<=K<q` uses the odometer extension and would correspond to
negative times at site `-q`.  The block at depths `k=q+t`, `0<=t<=q`, is
physical and occurs at times `t=0,...,q`.

For each `q`, define the clopen set

```text
V_q=intersection_(k=q)^(2q-1){h:T_h(k)=0}
    intersection {h:T_h(2q)=1}.                       (44a)
```

It contains a phase cylinder around `-q` modulo `P_(2q+1)`.  The construction
above shows that `union_(q>=L)V_q` meets every phase cylinder, so it is open
and dense.  Baire's theorem therefore gives a residual set of extended
2-adic phases with arbitrarily long zero runs.  A general 2-adic phase is not
a physical lattice site.  Long zero runs are compatible with density `1/2`;
this is a decisive hostile to bounded-run or uniform-cylinder transfer, not a
refutation of balance.

## 7. Fixed-seed finite carrier lower bounds

The next statements concern the actual single-seed lookup `n |-> c_n`.  They
are exact on the displayed finite universes and are not extrapolated.

### 7.1 Binary-index decision queries

For `m>=1`, let

```text
f_m:{0,1}^m -> F_2,
f_m(x_0,...,x_(m-1))=c_[x]_2,                         (45)
```

where `[x]_2` is the padded `m`-bit index.  Exact Boolean decision-tree
recursion gives

```text
D_query(f_m)=m                   for 2<=m<=18.          (46)
```

For standard positive `m`-bit indices, put

```text
g_m(x)=c_(2^(m-1)+[x]_2),       x in {0,1}^(m-1).     (47)
```

Then

```text
D_query(g_m)=m-1                for 5<=m<=16.          (48)
```

These are named-seed bit-query lower bounds, but only `Omega(log n)` in the
numeric index.  They make no linear-in-`n` time claim.

### 7.2 Unary homogeneous-state obstruction

Let

```text
Hankel_d=(c_(i+j))_(0<=i,j<d)                         (49)
```

over `F_2`.  Exact elimination gives

```text
d       16  32  64  128  256  512  1024  2048  4096  8192
rank    15  31  64  128  255  511  1023  2047  4096  8192. (50)
```

In particular the `8192` matrix is nonsingular.  Any representation

```text
c_n=lambda(T^n v)                                     (51)
```

over `F_2` which matches through `c_16382` has state dimension at least
`8192`.  Bounded all-prefix Hankel rank would be a rational OGF and hence an
eventually periodic binary sequence; proving unbounded rank would settle the
first prize.  Formula (50) is a finite lower bound only.

### 7.3 Binary-address tensor obstruction

Flatten the padded eighteen-bit lookup after `k` low bits:

```text
M_k(r,q)=c_(r+2^k q),
0<=r<2^k,       0<=q<2^(18-k).                        (52)
```

The complete `F_2` rank profile for `k=1,...,17` is

```text
2,4,8,16,32,64,128,256,510,256,128,64,32,16,8,4,2.   (53)
```

Every cut has maximal possible rank except the `9/9` cut, whose deficiency is
two in characteristic two.  The central `512 x 512` matrix has nonzero
determinant modulo three, so its rank over `Q` is `512`.  Since the other cuts
are already maximal modulo two, the rational tensor-train/Schmidt-rank profile
is maximal at every digit cut.

Restrict the high nine-bit word to `256<=q<512`, so the index has standard
eighteen-bit form.  The resulting `256 x 512` matrix has `F_2` rank `256`, all
`256` high-prefix rows are distinct, and all `512` low-prefix columns are
distinct.  Therefore any exact DFAO for the full center sequence on
standard binary representations needs at least

```text
256 states when read MSD first,
512 states when read LSD first.                        (54)
```

Equations (53)--(54) also lower-bound finite digit-matrix and tensor-train
carriers, and OBDDs in the corresponding MSD-first or LSD-first contiguous
block orders.  They do not lower-bound a general algorithm.  Equality of two binary
halves has maximal communication rank yet an `O(log n)`-bit comparison, and
factorial sequences have large homogeneous Hankel carriers but simple
time-varying recurrences.  Those are the canonical hostiles to a general-time
inference.

### 7.4 Three distinct carrier invoices

The tests above separate three operations which must not be conflated:

1. ordinary Hankel rank tests a unary homogeneous iterate `T^n`;
2. the dyadic flattenings test bounded matrix products along the binary address
   of `n`, as in an automatic or fixed Berggren-address compiler; and
3. THM-3463's Mealy action retains the real Rule 30 observer, but the required
   section address is `0^(n-1)`, of length `n`, rather than the `O(log n)`
   binary address.

The Green kernel already has an `O(log n)` digit compiler.  A lawful remaining
state-growth target suggested by these invoices is a compressed, time-varying
cocycle for the nonlinear source/current array.

## 8. Cross-frontier connection ledger

| source -> target | exact preserved object | destroyed coordinate | needed sidecar / stopping boundary |
|---|---|---|---|
| ternary Green paths -> `H_d(n)` | endpoint-fibre parity | branch order and path multiplicity | source time, radial grade and bond orientation |
| `H_1` -> valuation clock | `nu_2(n+1) mod 2`, nonperiodicity, density `1/3` | all later Rule 30 sources | full folded source current |
| radial current -> scalar center | exact sum after setting every grade weight to one | which shell cancels which | `(u,v)` grading; every finite band is insufficient |
| phase -> innovation cube | topology, cylinders, Haar law and innovation joint distribution | no information; coordinatewise XOR is not the native cyclic law | `tau=Gamma o(+1)o Gamma^(-1)` for dynamics |
| Haar atlas -> prize center | spatial phase averages | Dirac evaluation at `h=0` | marked address plus every noninnovation readout |
| negative phase `-q` -> physical ray segment | depths `k>=q`, with time `t=k-q` | depths `k<q` are periodic-extension values only | flag `h+k>=0` |
| fixed-depth phase density -> temporal center density | Haar average at one depth | marked evaluation and cross-depth dependence | `delta_0`, readout family and depth discrepancy |
| `T_k -> Q_k` | transition positions | additive constant / owner sign | `c_k=T_k(0)` |
| `Q_k -> |supp Q_k|` | event count | cyclic spacing and alternating gaps | ordered event set or gap word |
| parent marginals -> zero current | separate densities | simultaneous-zero intersection | aligned Walsh correlation `E(UV)` |
| center prefix -> Hankel/digit ranks | exact ranks and restricted-carrier lower bounds | arbitrary time-varying algorithms | charged model, advice and preprocessing policy |

The valuation parity in (11) is a lawful base-two analogue of the hidden
valuation-parity time coordinate in THM-2227; it is not an LRC phase transfer.
The radial variable is the same kind of load-bearing grading whose removal
breaks direct factorial-moment transfers.  The three Green branches form a
commutative endpoint compiler, a quotient of the ternary word tree rather than
a Berggren ancestry tree: folding keeps endpoint parity and loses words.  No
intrinsic pairwise orientation appears, so a tournament would be cosmetic.
THM-3459's characteristic and representative ledger remains the Jacobian
boundary; nothing here supplies a planar Keller map.

## 9. Exact verification and open frontier

Reproduce the exact audit with

```bash
python3 04-computation/rule30_radial_innovation_carriers_thm3468.py
python3 -O 04-computation/rule30_radial_innovation_carriers_thm3468.py
```

The companion uses two independent center recurrences and a frozen prefix
hash.  It checks the centered restart and Green recurrences, period lifts,
transition currents, the complete motif-`110` table through `k=30`, physical
cylinder blocks, exact decision trees, Hankel ranks, and every eighteen-bit
flattening in the stated universes.  Universal claims are proved above; finite
ranges audit the implementation or are explicitly labelled finite exact.

The three prize-facing probes are now precise:

1. prove that the unbounded radial tail cannot cancel the valuation-parity
   backbone into an eventually periodic word;
2. control the aligned zero-intersection drift and then the marked address,
   rather than fitting another unmarked phase density; and
3. compress the nonlinear collision current in a fixed uniform binary-index
   cost model, or prove that its required state grows.

The official prize listing was last checked on `2026-08-15`; on that dated
evidence the repository treats all three questions as open.  This theorem
claims no prize solution and no literature novelty.
