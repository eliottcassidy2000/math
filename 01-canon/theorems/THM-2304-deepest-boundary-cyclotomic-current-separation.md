---
id: THM-2304
title: "Deepest-boundary cyclotomic current separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Let a rational step
  function have jumps only on the guard, five unit, and three blocker comb
  boundary banks of a strict first-depth-one row. At the exact THM-2276
  shallow-pair frequency K, every nondeep endpoint phase lies in the
  cyclotomic field of thirteen-conductor 13^(b-1), while deepest c_3
  endpoint phases lie in an extension of degree R=13^(c-b). The endpoint
  derivative current therefore has one lower-field coordinate and R-1
  independently separated deepest-current coordinates. If the Fourier
  coefficient at K vanishes, all R-1 nontrivial currents vanish separately
  and the trivial deepest current balances the lower current. On the 120
  live interior profiles R-1 is at least 168. Recursive shell peeling gives
  the same R-1 constraints and leaves only a virtual c_3/13^(c-b)
  boundary current at the middle depth. A full deepest comb is the sharp
  zero-current control, so this is a structural reduction rather than a
  profile closure; LRC(14) remains open.
source: codex-2026-07-25-deepest-boundary-current
depends_on:
  - THM-2276-shallow-owner-residue-aligned-crossing
related:
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2303-terminal-component-phase-current-and-defect-rank
script: 04-computation/lrc14_deepest_boundary_cyclotomic_current_thm2304.py
output: 05-knowledge/results/lrc14_deepest_boundary_cyclotomic_current_thm2304.out
script_sha256: 3410dca88824fb8c558cbf3425e78737f5a976f0a1d14c6ce6490f5084d04099
output_sha256: eaf292433beb1a20b4e604ae787c578cea541cafd40cbd0d9b1977ae50c11afd
hash_basis: working-tree bytes (LF)
---

# THM-2304 -- full scalar boundaries split the component phase current

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2303 identifies the missing LRC datum as a complex current on terminal
handoff components. Its hostile carrier is deliberately local: the component
endpoints may be cut wherever the desired cancellation requires. A full
scalar exclusive-owner or current-service set has much less freedom. Every
jump belongs to one of the nine actual comb boundary banks.

At THM-2276's shallow-pair frequency, the deepest blocker is the only bank
whose endpoint phases reach the top thirteen-cyclotomic conductor. Linear
independence over the lower conductor field separates those phases before
they can cancel. The faithful global sidecar is therefore not one more
unsigned root energy. It is a vector of signed endpoint currents:

```text
full component current at K
  =
lower/trivial current
  + 13^(c-b)-1 independently separated deepest currents.       (1)
```

The theorem does not prove that one of these currents is nonzero. It proves
that any exact cancellation must kill every one of them separately.

## 1. Comb boundaries and endpoint current

For a positive integer `a`, put

```text
D_a={x in R/Z:||a x||<1/14}.                         (2)
```

Its boundary bank is

```text
partial D_a
 ={(14r+epsilon)/(14a) mod 1:
      0<=r<a, epsilon in {-1,+1}}.                  (3)
```

The guard-safe set

```text
C_H={x:||Hx||>1/7}
```

has endpoints `(7r+epsilon)/(7H)`, hence also has denominator dividing
`14H`. Boolean combinations of the guard, five unit danger combs, and three
blocker danger combs are finite step sets with no other jump locations.

Use a right-continuous representative of a periodic step function `f` and
write

```text
Delta f(x)=f(x+)-f(x-).                              (4)
```

Distributional differentiation gives, for every nonzero integer `K`,

```text
2*pi*i*K f_hat(K)
 =sum_(x in Jump(f)) Delta f(x)exp(-2*pi*i*Kx).     (5)
```

Thus the desired Fourier coefficient is exactly an oriented boundary
current. Formula (5), rather than the component magnitudes, is the endpoint
version of THM-2303's component-current sum.

At coincident scalar boundaries, (4) is the actual total jump. If one wants
bank provenance, the one-sided product rule assigns the coincident terms;
their sum remains (4). A nontrivial deepest conductor class cannot coincide
with a lower bank, so no convention affects the separated coordinates below.

## 2. Full-gap cyclotomic separation

The algebraic mechanism works in the following slightly wider scope. Let

```text
p=13,

nu_p(K)=1,

a_*=p^c u_*                         (p does not divide u_*),

2<=b<c,

nu_p(a)<=b<c
```

for every other boundary-bank coefficient `a`. Let `N` be any positive
integer prime to `p` which contains every prime-to-`p` endpoint denominator.
For the scalar application one may take `N` to be the prime-to-thirteen part
of

```text
lcm(14H,14q_1,...,14q_5,14u_1,14u_2,14u_3).        (6)
```

Define

```text
L=Q(zeta_(N p^(b-1))),

M=Q(zeta_(N p^(c-1))),

xi=zeta_(p^(c-1)),

R=p^(c-b).                                          (7)
```

Because `N` is prime to `p`,

```text
[M:L]
 =phi(Np^(c-1))/phi(Np^(b-1))
 =p^(c-b)
 =R.                                                (8)
```

Moreover

```text
xi^R=zeta_(p^(b-1)) in L.
```

The minimal polynomial of `xi` over `L` consequently has degree `R`, so

```text
1,xi,...,xi^(R-1)                                  (9)
```

is an `L`-basis of `M`.

Every nondeep endpoint `x` has phase

```text
exp(-2*pi*i*Kx) in L.                               (10)
```

Indeed, multiplication by `K` removes one factor of thirteen from its
denominator, leaving thirteen-conductor at most `p^(b-1)`.

For a deepest endpoint

```text
x=(14r+epsilon)/(14p^c u_*),                        (11)
```

the phase has conductor dividing `Np^(c-1)`. Splitting its prime-to-thirteen
and thirteen-primary parts, then reducing the exponent modulo `R`, gives a
unique expression

```text
exp(-2*pi*i*Kx)=alpha_x xi^(s_x),

alpha_x in L,             0<=s_x<R.                 (12)
```

For `1<=s<R`,

```text
nu_p(s)<=c-b-1,
```

so `alpha xi^s` has thirteen-conductor at least `p^b`. It cannot be a
lower-bank phase, whose conductor is at most `p^(b-1)`.

Group the actual deepest-boundary jumps by their class in (12):

```text
C_s
 =sum_(x deepest jump, s_x=s)
      Delta f(x)alpha_x in L,           1<=s<R.     (13)
```

Let `B in L` be the sum of every nondeep jump and every deepest jump with
`s_x=0`. Equation (5) becomes the unique basis expansion

```text
2*pi*i*K f_hat(K)
 =B+sum_(s=1)^(R-1) C_s xi^s.                       (14)
```

The basis (9) proves the exact equivalence

```text
f_hat(K)=0

iff

B=0 and C_1=C_2=...=C_(R-1)=0.                     (15)
```

These are `L`-valued signed currents, not raw endpoint counts. Equation
(15) is much stronger than saying that the total deepest contribution
vanishes: every conductor/address class must vanish separately. It is also
stronger than a relative trace, which sees only the trivial coordinate.

The hypotheses `nu_13(K)=1` and “all jumps lie in the named scalar banks”
are load-bearing. An arbitrary cut endpoint can introduce another conductor,
and a multiplier divisible by thirteen changes the base field and the number
of separated coordinates.

## 3. Application to the 120 interior LRC rows

On a strict first-depth-one profile, use

```text
c_1=13u_1,

c_2=13^b u_2,

c_3=13^c u_3,

3<=b<=c-2,                    5<=c<=19.             (16)
```

THM-2276 supplies a shallow pair with

```text
K=plus or minus m c_1,

1<=m<=757,                    13 does not divide m. (17)
```

Thus `nu_13(K)=1`. Every endpoint from the guard, the five unit masks,
`c_1`, or `c_2` has phase in the field `L` of (7); only `c_3` can reach
`M minus L`.

Let `f` be the indicator of any full Boolean scalar stratum formed from
those nine masks. In particular, one may take the strict `c_1`-private set,
a prescribed current-service set, or an exact blocker-word stratum. Then
(15) applies with

```text
R=13^(c-b)>=13^2=169.                               (18)
```

Hence a zero exact pair coefficient forces at least

```text
R-1>=168                                             (19)
```

separate nontrivial deepest-boundary currents to vanish. Across the complete
strict profile bank,

```text
2<=c-b<=16,

168<=R-1<=13^16-1=665416609183179840.               (20)
```

This is the first constraint in the current phase program that uses the
global scalar boundary atlas and is unavailable to arbitrary local
components. It does not require choosing a tournament. The intrinsic object
is a signed current vector indexed by relative conductor/address class.

The same separation remains valid for step amplitudes in `L`. In particular,
the thirteenth-root coefficients of a rooted word lie in `L` because
`b-1>=2`. This observation does not by itself turn the fractionally gauged
current of THM-2299 into the exact pair atom; the Fourier frequency and target
restriction still have to match.

## 4. Recursive shell peeling

The full-gap basis has an equivalent recursive description. Let

```text
f=L_0 1_(D_a^c),                 a=13^d u,          (21)
```

where `L_0` is a lower-bank step mask. At a noncoincident boundary
`(14r+epsilon)/(14a)`, the `D_a^c` jump has orientation `epsilon`, multiplied
by the appropriate one-sided value of `L_0`.

For collisions, fix the right-sided product convention

```text
Delta(L_0 g)(x)
 =L_0(x+)Delta g(x)+g(x-)Delta L_0(x).              (21a)
```

Define the bank-provenance current explicitly by

```text
J_a
 =sum_(x=(14r+epsilon)/(14a))
    L_0(x+) epsilon exp(-2*pi*i*Kx),                 (21b)
```

and let `U_a` be the sub-sum with
`13` not dividing `14r+epsilon`. The second term in (21a) belongs to the
lower-bank current. Thus coincident boundaries are neither duplicated nor
dropped, and the definition agrees with the ordinary jump contribution
away from collisions.

If

```text
13 divides 14r+epsilon,
```

then `r=13r'-epsilon` and

```text
14r+epsilon=13(14r'-epsilon).                       (22)
```

The same point is therefore the opposite-side boundary

```text
(14r'-epsilon)/(14(a/13))
```

of the virtual comb `D_(a/13)`. Its orientation is reversed, while the
right-sided multiplier `L_0(x+)` in (21b) is unchanged. Therefore

```text
J_a=U_a-J_(a/13).                                   (23)
```

Peeling (23) from `a=c_3` down to valuation `b` says:

```text
J_(c_3)
 =U_(c_3)-U_(c_3/13)+...+
   (-1)^(c-b-1)U_(c_3/13^(c-b-1))
   +(-1)^(c-b)J_(13^b u_3).                         (24)
```

At the first peel, the twelve nonzero last digits give

```text
12*13^(c-b-1)
```

`L`-coordinates. The successive trivial-digit branches contribute

```text
12*13^(c-b-2),...,12.
```

Their total is

```text
12(1+13+...+13^(c-b-1))
 =13^(c-b)-1
 =R-1.                                              (25)
```

Thus recursive current descent and the one-shot basis (9) are the same
decomposition in opposite orders. If all nontrivial currents vanish, the
only deepest contribution left is the virtual boundary current

```text
(-1)^(c-b)J_(13^b u_3)                              (26)
```

at the same thirteen-depth as the actual middle blocker `c_2`. This turns an
unstructured component cancellation into an exact final interaction:

```text
actual lower-bank current
  versus
virtual-descended deepest current at middle depth.                (27)
```

Closing (27), or proving one colored `U` current nonzero from the scalar
cover, is the residual mathematical obligation.

## 5. Sharp controls and stopping boundary

### 5.1 A single deepest tooth fires a separated current

Take one component interval of `D_(c_3)`. Its two endpoint exponents differ
by a thirteen-unit multiple of `2`. They occupy distinct classes modulo
thirteen and therefore distinct classes modulo `R`; at least one class is
nontrivial. The corresponding `C_s` cannot cancel with the other endpoint
or with any lower bank. Formula (15) proves its coefficient at `K` is
nonzero.

### 5.2 The full deepest comb is the exact hostile control

For the whole comb,

```text
(1_(D_a))_hat(n)=0 unless a divides n.               (28)
```

Here `nu_13(c_3)=c>1=nu_13(K)`, so

```text
(1_(D_(c_3)))_hat(K)=0.                             (29)
```

The full deepest bank has every boundary but all its separated currents
cancel by complete comb uniformity. Thus “a deepest endpoint is active” is
not enough. A future positive theorem must prove a colored defect from that
uniform current, not merely boundary incidence.

The multiplier-four carrier of THM-2299 is also untouched. Its two
components use freely cut local endpoints rather than the full nine-bank
Boolean boundary, and its pair cancellation lies in the lower/trivial
coordinate. It remains the correct proof that root energy and discrete
labels alone do not determine phase.

## 6. Connection and loss ledger

```text
source:
  a full Boolean scalar stratum on a strict profile, THM-2276's exact
  shallow-pair frequency, and THM-2303's component-current requirement;

map:
  replace components by their oriented endpoint derivative current,
  factor endpoint phases into prime-to-thirteen and thirteen-primary
  parts, and expand over the full relative cyclotomic basis;

preserved:
  exact frequency, actual scalar boundary banks, jump signs and
  amplitudes, deepest address/conductor class, and total Fourier current;

destroyed by total current alone:
  which of the R-1 deepest classes supplied or cancelled a contribution;

restoring sidecar:
  (C_1,...,C_(R-1)) together with the lower/trivial balance B;

cheapest next tests:
  prove a nonuniform top-digit occupancy current on one full private set,
  or close the final actual-c_2/virtual-c_3 interaction (27).           (30)
```

This connection is exact:

```text
THM-2303 component U(1) phase current
  -> endpoint derivative current
  -> THM-2304 cyclotomic address vector.
```

What it loses is component pairing away from the named frequency. What it
gains is the full scalar-cover endpoint restriction absent from every local
phase no-go.

## 7. Exact verification

The companion uses exact integer arithmetic. It checks:

- all `120` strict pairs `(b,c)`;
- the conductor gap and relative degree `13^(c-b)`;
- the complete recursive identity
  `12(1+13+...+13^(c-b-1))=13^(c-b)-1`;
- the minimum `168` and maximum `13^16-1` current counts;
- separation of the two endpoints of every single-tooth residue control;
- the strict conductor gap between every nontrivial deepest class and every
  lower-bank phase;
- the full-comb and arbitrary-cut stopping controls; and
- that rooted thirteenth-root amplitudes remain inside the base field.

Reproduce with

```bash
python3 04-computation/lrc14_deepest_boundary_cyclotomic_current_thm2304.py
python3 -O 04-computation/lrc14_deepest_boundary_cyclotomic_current_thm2304.py
```

Both runs are byte-identical to the stored output, with every load-bearing
check active under optimized Python. The exact field-degree and endpoint
current arguments are the proof above rather than delegated numerical
conclusions. No scalar profile is excluded. QED.
