---
id: THM-3516
title: "Rule 30 marked van der Put carry and power-section bridge"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The marked
  center at every positive depth is an explicitly calibrated high binary
  digit of an odd all-n normalized displacement.  That displacement is both
  a telescoped sum of equal-scale units and a van der Put path sum over the
  set bits of the marked time.  Exact prefix fibers discard precisely the
  shells above the current innovation count, while a physical depth-six
  example proves that the surviving target digits do not combine by XOR:
  the lower binary carry is load-bearing.  At powers of two the path
  collapses to one shell unit and one deep orbit-signalizer section.  No Rule
  30 prize consequence or fixed-seed complexity lower bound is claimed.
source: root/rule30-normalized-displacement/marked-carry/2026-08-16
audit: >
  PASS (2026-08-16).  The auditor independently rederived the all-n
  normalization, period owner, active-shell and carry formulas, projective
  gauge, signalizer cancellation, and exact unary future-output quotient;
  regenerated conventional Rule 30 centers through depth 512; and confirmed
  ordinary, optimized, and stored transcript equality with no assertion gates.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3493-rule30-dyadic-wrap-atlas
  - THM-3507-rule30-normalized-dyadic-displacement-sibling-trace-and-assouad-spectrum
  - THM-3511-rule30-orbit-signalizer-gap-renormalization-and-shallow-portrait-hostile
  - THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary
related:
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-3489-rule30-packed-restart-and-pointed-pascal-face
script: 04-computation/rule30_marked_van_der_put_carry_thm3516.py
output: 05-knowledge/results/rule30_marked_van_der_put_carry_thm3516.out
script_sha256: 64f394aa3556c19394f598e0f2c2a850d42417b2d33b502a4d56a7c7d27f03c3
output_sha256: 140c473ee75f0b0a94bc90e4c5205f22799e33e4f00f1f57cc8553a0483edf02
hash_basis: raw bytes
---

# THM-3516 -- Rule 30 marked van der Put carry and power-section bridge

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3512 turns THM-3507's normalized displacements into literal van der Put
coefficients, but deliberately stops before the moving marked center.  The
missing operation is now explicit: follow the calibrated binary phase path,
add its oriented shell increments by ordinary integer addition, and retain
the carry into one moving high digit.  At a power-of-two time the path has
one edge, so THM-3511 turns the bit into one deep section activity.

The inheritance pass is:

1. closest proved mechanism: THM-3512's shell coefficients and exact prefix
   fibers;
2. canonical hostile: at depth six the XOR of the two visible target digits
   is wrong because their lower residues carry into the target;
3. corrected near miss: projective ratios recover a marked value only after
   one signed amplitude gauge is calibrated; and
4. least-used sidecars: the low phase representative, ordinary carry, and
   zero-ray section at depth `2^m-v_m`.

The live concepts are the all-n displacement, van der Put path, phase fiber,
carry parity, projective gauge, and orbit-signalizer zero ray.  There is no
intrinsic pairwise relation and hence no tournament.

## 1. Notation and the moving marked bit

Retain

```text
R_0=1,
R_(t+1)=Phi(R_t),             Phi(x)=x xor ((2x) or (4x)),
iota(t)=R_t for t>=0,
c_n=bit_n(R_n),
J(t)=(iota(t)-1)/2.                                  (1)
```

Let

```text
P_(k+1)=2^(epsilon_k)P_k,     E_k=log_2 P_k,
kappa_1<kappa_2<...={k:epsilon_k=1},
v_m=kappa_(m+1),              q_m=2^m,
d_m=v_(m+1)-v_m,                                    (2)

U_m(t)=(iota(t+q_m)-iota(t))/2^(v_m).                (3)
```

Every `U_m(t)` is an odd 2-adic unit.  At nonnegative integer phases it is
an odd integer.  Since `R_n=1+2J(n)`, the physical center is exactly

```text
boxed: c_n=bit_(n-1)(J(n)),             n>=1.        (4)
```

No fixed-width quotient is being used in (4): the requested output digit
moves with `n`.

## 2. The all-n normalized displacement

For `n>=1`, put

```text
alpha=nu_2(n),
W_n=(R_n-1)/2^(v_alpha).                             (5)
```

THM-3512's radial distance law, applied to phases `0` and `n`, gives

```text
nu_2(R_n-1)=v_alpha.                                 (6)
```

Consequently `W_n` is odd.  Combining (4)--(6) gives the first marked
bridge, including both index edge cases:

```text
boxed:
c_n=0                              if v_alpha>n,
c_n=bit_(n-v_alpha)(W_n)           if v_alpha<=n.    (7)
```

If `v_alpha=n`, the requested digit is digit zero and is therefore one.  If
`v_alpha>n`, the factor `2^(v_alpha-1)` in `J(n)` starts beyond digit `n-1`,
so the center is zero.

There are two exact decompositions of `W_n`.  First, split `[0,n)` into
blocks of length `q_alpha=2^alpha`.  Equation (3) telescopes to

```text
boxed:
W_n=sum_(j=0)^(n/2^alpha-1) U_alpha(j2^alpha).       (8)
```

The number of summands is odd, which also proves directly that `W_n` is odd.

Second, for every set bit `m` of `n`, put

```text
r_m(n)=n mod 2^m.                                    (9)
```

Recursively deleting the highest set bit in the van der Put expansion gives

```text
J(n)=sum_(m:bit_m(n)=1) C_(2^m+r_m(n))(J)
    =sum_(m:bit_m(n)=1) 2^(v_m-1)U_m(r_m(n)).        (10)
```

The least set bit is `alpha`, and `v_m` is strictly increasing, so division
by `2^(v_alpha-1)` is integral and yields

```text
boxed:
W_n=sum_(m:bit_m(n)=1)
        2^(v_m-v_alpha)U_m(r_m(n)).                  (11)
```

Thus a center value uses at most `popcount(n)` dyadic block owners.  This is
an exact representation count, not a running-time bound: the needed unit
digits occur at precision proportional to `n`.

The first edge cases are literal.  At `n=1`, `alpha=0`, `v_0=1`, and
`W_1=3`, so (7) reads `c_1=bit_0(3)=1`.  At `n=2`, `alpha=1` and `v_1=3>2`,
so (7) reads `c_2=0` even though `W_2=3` is odd.

## 3. Exact calibrated phase fiber and active shells

Put

```text
e=E_(n+1)=log_2 P_(n+1),
p=P_(n+1)=2^e,
rho_n=n mod p,               0<=rho_n<p.             (12)
```

Both `n` and `rho_n` lie below `2^n`.  THM-3512's exact prefix-fiber law
therefore gives

```text
boxed:
J(n)=J(rho_n) mod 2^n,
c_n=bit_(n-1)(J(rho_n)).                             (13)
```

The owner `rho_n`, rather than an unpointed phase orbit or Haar average, is
the calibrated phase coordinate.  Its indexing is exact:

```text
boxed: m<e  iff  v_m<=n.                             (14)
```

Indeed, `e` counts the innovations `kappa_j` at depths at most `n`, while
`v_m=kappa_(m+1)`.  Since `rho_n<2^e`, applying (10) to `rho_n` gives

```text
boxed:
J(n)=sum_(0<=m<e, bit_m(rho_n)=1)
       2^(v_m-1)U_m(rho_n mod 2^m)       mod 2^n.    (15)
```

Equivalently, start with the full sum (10): every term with `v_m>n` is zero
modulo `2^n`, and the remaining set bits are exactly those of `rho_n`.  The
exact prefix fiber and van der Put typing discard the same information.

The support bound supplies a check on (13).  The row `R_r` is supported
through bit `2r`, so `J(r)` is supported through bit `2r-1`.  Hence
`2rho_n<n` implies `c_n=0`; if the formal equality `2rho_n=n` occurs, the
persistent extreme-left bit gives `c_n=1`.  THM-3493's dyadic period floor
rules out that equality on the actual period tower.  These statements
recover, rather than strengthen, its wrap-prefix atlas.

For example, write `n=2^m+r`, `0<=r<2^m`.  Then

```text
J(n)=J(r)+2^(v_m-1)U_m(r).                           (16)
```

The base `J(r)` has no digit at `n-1`.  If `v_m>n`, (16) gives `c_n=0`; if
`v_m=n`, it gives `c_n=1`.  This is precisely the zeros-then-one part of a
nonempty wrapped dyadic block.

## 4. The marked bit is an ordinary carry, not an XOR

Assume `v_alpha<=n` and put

```text
L=n-v_alpha,
A_n={m:bit_m(n)=1 and v_m<=n}.                        (17)
```

For `m in A_n`, retain exactly the necessary unit precision and define the
least nonnegative residue

```text
T_(m,n)=2^(v_m-v_alpha)
         [U_m(r_m(n)) mod 2^(n-v_m+1)]
         mod 2^(L+1).                                (18)
```

Then (11) gives

```text
W_n=sum_(m in A_n)T_(m,n) mod 2^(L+1).               (19)
```

Split each residue at the target digit,

```text
b_(m,n)=bit_L(T_(m,n)),
lambda_(m,n)=T_(m,n) mod 2^L.                        (20)
```

Ordinary binary addition gives the exact carry law

```text
boxed:
c_n=
 xor_(m in A_n)b_(m,n)
 xor [floor(sum_(m in A_n)lambda_(m,n)/2^L) mod 2]. (21)
```

For `L=0`, every `lambda` is zero and (21) is the odd leading digit.  For
`L>0`, the second term is load-bearing in general.

### 4.1 The minimal physical carry hostile

At `n=6=4+2`,

```text
R_6=6409,              J(6)=3204,
alpha=1, v_alpha=3,    W_6=801.                      (22)
```

The two shell units and both normalizations are

```text
U_1(0)=3,              U_2(2)=399,

J(6)=12+3192,          (12,3192) mod 64=(12,56),
W_6=3+2*399=3+798,     (3,798) mod 16=(3,14).        (23)
```

The target in `W_6` is digit `L=3`.  The two target digits XOR to one, but
their lower residues are `3` and `6`; these sum to `9` and carry once across
`2^3`.  Hence

```text
boxed: c_6=0 xor 1 xor 1=0.                          (24)
```

The companion checks depths one through five first; depth six is the
smallest positive physical example where the carry correction in (21) is
one.  Coefficient parities, individual target digits, a mod-two Haar detail,
or a Cartier shadow of those digits cannot replace the ordinary carry.

### 4.2 Sharp adaptive precision in the ambient unit model

For `v_m<n`, knowing an active unit only modulo `2^(n-v_m)` is uniformly
insufficient.  Replacing that odd unit by

```text
U_m(r_m)+2^(n-v_m)                                   (25)
```

preserves its oddness and every lower retained digit, but changes (18) by
exactly `2^L` and toggles `c_n`.  Thus precision `n-v_m+1` is sharp when the
shell units are treated as independent odd inputs.  Actual Rule 30 units
obey cross-scale and cross-phase relations, so (25) is not a fixed-seed
complexity lower bound.

## 5. Projective gauge and the power-of-two bridge

Recall THM-3512's projective sibling quantities

```text
G_m(t)=-U_m(t+q_m)/U_m(t),
Z_m(t)=(1-G_m(t))/2^(d_m)=U_(m+1)(t)/U_m(t).         (26)
```

At the marked origin, `U_0(0)=3`, so

```text
boxed:
U_m(0)=3 product_(j=0)^(m-1) Z_j(0).                 (27)
```

More generally,

```text
U_m(t)=U_0(t) product_(j=0)^(m-1)Z_j(t),
U_0(t)=3(-1)^t product_(s=0)^(t-1)G_0(s).            (28)
```

Thus the complete projective field plus the signed gauge `U_0(0)=3`
recovers every integer-phase shell owner.  The gauge is necessary.  The
abstract isometry `J -> -J` preserves all valuations, exact prefix fibers,
and every `G_m,Z_m`, but changes the gauge from `3` to `-3`; already
`J(3)=55` has marked bit one whereas `-55 mod 8=1` has marked bit zero.  This
is not a second Rule 30 orbit; it identifies what the quotient forgets.

Let now `q=2^m`.  Its binary expansion has one set bit, so

```text
W_q=U_m(0).                                          (29)
```

Equations (7), (27), and (29) give

```text
boxed:
c_(2^m)=0,                                           if v_m>2^m,

c_(2^m)=bit_(2^m-v_m)(U_m(0))
         =bit_(2^m-v_m)(3 product_(j<m)Z_j(0)),      if v_m<=2^m.
                                                            (30)
```

The equality case requests digit zero of an odd unit and gives one.

## 6. Orbit-signalizer section and its cancellation

Let `s_m` be THM-3511's exact marked orbit signalizer, let
`w_m in {A,B,C}^(2^m)` be its right-action word, let `tau_0` take a section
at input child zero, and let `activity(w)` be the root permutation bit.  Put

```text
Zray(s)=[(0^infinity)^s]_2.                           (31)
```

THM-3511 gives `U_m(0)=Zray(s_m)`.  Digit `ell` of the image of the zero ray
is the root activity after taking the input section at `0^ell`.  Therefore

```text
boxed:
c_(2^m)=0,                                           if v_m>2^m,

c_(2^m)=activity(tau_0^(2^m-v_m)(w_m)),              if v_m<=2^m.
                                                            (32)
```

At equality the section depth is zero and `activity(w_m)=1`.  When
`v_m>2^m`, the exponent would be negative and is deliberately not defined.

There is a sharper cancellation.  In THM-3511's convention, `A` fixes the
first seed bit and `A|_1=B`, while the marked cylinder is
`u_m=1 0^(v_m-1)`.  Hence, as an exact word identity,

```text
boxed: w_m=tau_0^(v_m-1)(B^(2^m)).                   (33)
```

If `v_m<=q=2^m`, substituting (33) into (32) cancels the innovation depth:

```text
tau_0^(q-v_m)(w_m)=tau_0^(q-1)(B^q).                 (34)
```

If `v_m>q`, the first `q` tail bits of `B^q(0^infinity)` are zero, so the
same final activity is zero.  Thus both branches unify as

```text
boxed:
c_q=activity(tau_0^(q-1)(B^q))       for every q=2^m. (35)
```

In fact (35) holds for every positive integer `q`: after the first seed bit,
the section of `A^q` is `B^q`, and its tail output digit `q-1` is physical
packed digit `q`.

The exact finite control through `m=12` is

```text
m:       0   1   2   3   4   5   6    7    8    9    10   11    12
v_m:     1   3   4   6   7   9  15   16   24   25    27   29    34
section: 0   -   0   2   9  23  49  112  232  487   997 2019  4062
c_2^m:   1   0   1   1   1   0   1    1    1    0     1    0     1
                                                               (36)
```

The dash is the `v_1>2` branch, not a missing computation.

### 6.1 Exact unary future-output quotient lower bound

This same word orbit supports a precisely restricted complexity statement.
For any integer `q>=1`, define

```text
T=tau_0,
W_(q,k)=T^k(B^q),
a_k=activity(W_(q,k)).                                (37)
```

The calibrated tail identity is

```text
boxed: a_k=bit_(k+1)(R_q).                            (38)
```

The isolated row `R_q` has its unique extreme top bit at `2q` and vanishes
above it.  Thus

```text
a_(2q-1)=1,               a_k=0 for k>=2q.           (39)
```

Moreover,

```text
boxed: W_(q,2q)=A^q,       T(A^q)=A^q.               (40)
```

For (40), consider the section of each of the `q` ordered `B` factors after
the input prefix `0^(2q)`.  Before the `i`th factor, the transformed zero ray
is the tail of the finite row generated by `B^i`, whose support ends before
depth `2q`; the two last input bits seen by that factor are `00`, and the
corresponding `B` section is `A`.  Hence every factor has become `A`.

Declare future-output equivalence only on this unary Moore system:

```text
u equiv v iff activity(T^r u)=activity(T^r v) for every r>=0. (41)
```

Then

```text
boxed:
W_(q,0),W_(q,1),...,W_(q,2q)
have exactly 2q+1 distinct future-output classes.     (42)
```

If `i<j<=2q` were equivalent, the tail `(a_k)_(k>=i)` would be periodic with
period `j-i`.  Its eventual zero in (39) would then force every term from
`i` onward to be zero, contradicting `a_(2q-1)=1`.  Equation (40) makes the
last class fixed, so (42) is exact.

Consequently any deterministic unary Moore quotient which intertwines `T`
and preserves all future activities on this one orbit needs at least `2q+1`
states.  In particular the marked target `k=q-1` occurs before any raw-state
repeat.  This excludes raw cycle detection and a fixed finite
future-output syntactic quotient.  It is not a lower bound for a
target-only, nonuniform, arithmetic, circuit, word-RAM, or general Rule 30
algorithm.

### 6.2 A sharp Prize-1 subtarget

If the full center were eventually periodic with period `T=2^a u`, `u` odd,
then the power subsequence in (32) would be eventually periodic in `m`.
For `m>=a`, the residues `2^m mod T` have period `ord_u(2)` when `u>1` and
are constant when `u=1`.  Consequently,

```text
the right side of (32) is not eventually periodic
    ==> the Rule 30 center is not eventually periodic.         (43)
```

This is a one-way reduction to an exact deep-section diagonal.  Nothing here
proves its nonperiodicity.  The powers of two have density zero, so it gives
no center-balance conclusion.

## 7. Interpretation and preservation ledger

Equation (15) is a nonstationary Bratteli path integral in literal algebraic
terms: `rho_n` is the low phase address, a set bit selects an edge, and the
van der Put coefficient is its oriented cocycle.  Innovation heights say
where each edge becomes visible in the output tree.  Incidence or entropy
retains those heights; the marked center additionally asks for the carry of
the oriented path integral at moving output level `n-1`.

| source -> target | preserved | destroyed | required sidecar / hostile |
|---|---|---|---|
| `n -> rho_n` by prefix fiber | state prefix through the marked digit | high phase bits | calibrated least phase representative |
| van der Put path -> `W_n` | exact all-n displacement | nothing | ordinary integer addition |
| shell target digits -> `c_n` | direct digit contributions | lower carry | carry parity; physical hostile `n=6` |
| shell unit -> fixed residue | digits below cutoff | moving target digit | precision `n-v_m+1` |
| sibling pair -> `G_m,Z_m` | gaps and amplitude ratios | signed gauge | `U_0(0)=3`; hostile `J` versus `-J` |
| `U_m(0) -> s_m` | marked power-of-two owner | nothing | zero ray and right-action convention |
| shallow portrait -> center | bounded initial zero-ray digits | digit at depth `2^m-v_m` | deeper section / overflow |
| unary word -> future-output class | every later raw section activity | other computations of one target | model restricted by (41) |

The carry is ordinary binary carry, not XOR holonomy.  The zero ray is an
ordered input path, not an unpointed portrait.  The powers-of-two reduction
is a sparse exact target, not a balance statistic.

## 8. Finite verification and no-prize boundary

The companion has the frozen universes

```text
packed recurrence:        exact untruncated R_t, 0<=t<=4096;
innovation scales:        0<=m<=12;
all-n path/fiber audit:    1<=n<=32;
seed-period controls:      exact orbit and first return modulo 2^(n+1)
                           for every 1<=n<=32;
block telescoping:         342 equal-scale shell terms;
van der Put path:          80 marked shell terms;
projective origin:         12 exact rational identities;
signalizer bridge:         independent right-action words, 0<=m<=12;
unary future output:       every 1<=q<=256 and all 2q+1 classes. (44)
```

The packed path and signalizer recurrence are distinct finite
representations.  The companion checks `n=1,2`, freezes the depth-six carry
data, verifies the period owner and both decompositions at every declared
all-n depth, compares (32) and (35) with the direct center, and checks
(38)--(42) throughout the declared unary universe.

Run

```bash
python3 04-computation/rule30_marked_van_der_put_carry_thm3516.py
python3 -O 04-computation/rule30_marked_van_der_put_carry_thm3516.py
```

and compare both byte-for-byte with
`05-knowledge/results/rule30_marked_van_der_put_carry_thm3516.out`.  The
companion contains no Python assertion gates.

This theorem proves no center nonperiodicity, no limiting center density,
no balance, no bounded or unbounded innovation gaps, and no general
computation lower bound.  The `popcount(n)` shell count is not an algorithm
until the adaptively precise phase-owned units are computed.  The ambient
precision hostile (25) is not a physical Rule 30 perturbation.  The unary
lower bound is only for the explicitly declared future-output quotient.  The
finite table is not extrapolated.  No Rule 30 prize or literature novelty
claim is made.
