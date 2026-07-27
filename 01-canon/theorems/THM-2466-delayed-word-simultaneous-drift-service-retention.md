---
id: THM-2466
title: "Delayed-word simultaneous drift and service retention"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For any
  fixed positive rational Boolean
  word, one sufficiently delayed base-thirteen copy simultaneously
  retains a nonzero Hilbert-valued endpoint drift and positive
  nonnegative common-root service. After disintegrating both densities
  onto the same oriented root base, the BV errors are at most
  C_Q V_b/13^(k-1) and C_Q V_s/13^(k-1), where
  C_Q=min(mu(Q)(1-mu(Q)),Var(Q)/12). Applied after
  THM-2459, the same delayed semantic word retains endpoint drift at
  least mu(Q)^2 D_0/252004 and service at least
  mu(Q) M_0/32768 in a union of at most four enhanced atoms. This is
  conditional on a supplied owner-supported packet, common oriented
  physical root section, and fixed pre-delay densities. Without owner
  support the result retains only a terminal-stratum filter, not a
  semantic owner-to-word coupling. It does not identify the canonical
  owner, preserve a prescribed first-expiration clock, or exclude a
  scalar row.
source: codex-2026-07-26-delayed-word-joint-retention
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
  - THM-2459-four-atom-drift-and-root-service-coarsening
  - THM-2460-idempotent-semantic-word-copy-and-word-block-cosupport-boundary
related:
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2442-delayed-word-septimal-source-completion
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2464-linked-blocker-depth-threshold-and-clock-cell-universality
script: 04-computation/lrc14_delayed_word_joint_retention_thm2466.py
output: 05-knowledge/results/lrc14_delayed_word_joint_retention_thm2466.out
script_sha256: 0940c48febdf5490cc2814234801a68b84fbc011d543b6f4a1b5e287c689e1b0
output_sha256: 3e1ad096ac6a1d6f8a4ef93103eaaaa69ce193355fff8a34614935b1a6461f65
hash_basis: working-tree bytes (LF)
---

# THM-2466 -- delay decorrelates one word from both surviving densities

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2460 isolates the exact remaining word obstruction at a fixed clock:
selected endpoint drift may live in one semantic word block while all
common-root service lives in another.  That hostile is sharp when the three
word blocks are frozen.  It is not stable under an arbitrarily late common
base-thirteen delay.

The missing coordinate is the order of quantifiers.  First fix the endpoint
drift density and the common-root service density.  Then dilate one fixed
positive terminal word.  The expanding map makes that same word
simultaneously independent of both fixed densities:

```text
fixed root-base drift b, fixed root-base service s, then k -> infinity

  Q(13^(k-1)y)b(y)   -> mu(Q)b(y) in integrated drift,

  Q(13^(k-1)y)s(y)   -> mu(Q)s(y) in integrated service.       (1)
```

The one-power offset is physical rather than cosmetic.  On the oriented
root section `x=(y+u)/13`, first disintegrate the endpoint table over its
thirteen inverse branches.  The endpoint word at clock `k` then becomes the
root-constant base word at clock `k-1`, so drift and service are functions
on the **same** base coordinate.

This theorem proves the resulting conditional bridge with an exact BV rate.
It does not supply the common root section or identify it with the canonical
THM-2305 owner.  Those typing statements remain the structural frontier.

## 1. A two-density dilation lemma

Let `T=R/Z` carry normalized Haar measure.  Let `H` be a finite-dimensional
real or complex Hilbert space.  Fix on one common base coordinate `y`:

```text
Q:T->{0,1},                 a rational circle step word,

mu=integral_T Q>0,

b:T->H,                    a fixed Bochner-integrable root-base drift density,

s:T->R_(>=0),              a fixed integrable service density.          (2)
```

Put

```text
q=integral_T b,            D=||q||_H^2,

M=integral_T s.                                                   (3)
```

Assume

```text
D>0,                       M>0.                                  (4)
```

For every integer `k>=1`, put `N_k=13^(k-1)` and define

```text
q_k=integral_T Q(N_k y)b(y)dy,

D_k=||q_k||_H^2,

M_k=integral_T Q(N_k y)s(y)dy.                           (5)
```

Then

```text
q_k -> mu q,

M_k -> mu M.                                                  (6)
```

In particular, there is `k_0` such that for every `k>=k_0`, the **same**
delayed copy of `Q` satisfies

```text
D_k>0,                       M_k>0.                         (7)
```

The statement remains true on either prescribed parity class of `k`, or on
any other cofinal subsequence.  This matters because

```text
13^k mod 7=(-1)^k mod 7.                                    (8)
```

Thus a lawful THM-2442 septimal parity may be frozen before the delay is
chosen.  Equation (7) by itself does not preserve a previously selected
septimal Fourier colour; that requires the phasewise THM-2442 construction.

### Proof of qualitative convergence

The functions `Q(N_k y)` converge weak-star in `L^infinity(T)` to the
constant `mu`.  One elementary proof is enough here.  For a trigonometric
monomial `e(my)`, the Fourier support of `Q(N_k y)` is contained in
`N_k Z`.  Hence

```text
integral_T Q(N_k y)e(my)dy=0
```

for every nonzero fixed `m` once `N_k>|m|`; for `m=0` the integral is
`mu`.  Approximate each scalar coordinate of `b` and the scalar function
`s` in `L^1` by trigonometric polynomials.  Since `0<=Q<=1`, the
approximation errors are uniform in `k`.  This proves (6).  The two limits
in (6) are respectively nonzero and positive by (4), so one common tail of
clocks proves (7).  Restricting that tail to a cofinal subsequence proves
the parity assertion.  QED.

## 2. Exact BV threshold

Suppose now that `b` and `s` have periodic bounded variation.  Use Hilbert
total variation for `b`:

```text
V_b=Var_H(b),              V_s=Var(s),

V_Q=Var(Q),

C_Q=min(mu(1-mu),V_Q/12).                                   (9)
```

Then the uniform estimates are

```text
||q_k-mu q||_H
 <=C_Q V_b/N_k,                                             (10)

|M_k-mu M|
 <=C_Q V_s/N_k.                                             (11)
```

Consequently, if

```text
N_k
 >=2 C_Q V_b/(mu sqrt(D)),

N_k
 >=2 C_Q V_s/(mu M),                                      (12)
```

then

```text
D_k>=mu^2 D/4,

M_k>=mu M/2.                                               (13)
```

Using only the Boolean-discrepancy half of `C_Q`, the sufficient thresholds
in (12) simplify to

```text
N_k>=2(1-mu)V_b/sqrt(D),

N_k>=2(1-mu)V_s/M.                                         (12a)
```

If only nonvanishing is needed, the factor `2` in (12) or (12a) may be
replaced by `1`, with strict inequalities.

### Proof of the two BV constants

Put `h=Q-mu` and choose the periodic primitive

```text
G(t)=integral_0^t h(z)dz,                 G(0)=G(1)=0.       (14)
```

For every `t`, the amount of a set of measure `mu` which can lie in
`[0,t]` is between `max(0,t-(1-mu))` and `min(t,mu)`.  Therefore

```text
||G||_infinity<=mu(1-mu).                              (15)
```

The function `G(Ny)/N` is a periodic primitive of `h(Ny)`.  Periodic
Riemann--Stieltjes integration by parts, also valid for Hilbert-valued BV
functions by duality, gives

```text
norm(integral_T h(Ny)f(y)dy)
 <=mu(1-mu) Var(f)/N.                                  (16)
```

For a scalar or Hilbert-valued circle-BV function `f`, integration by
parts gives

```text
||f_hat(n)||<=Var(f)/(2 pi |n|),             n!=0.          (17)
```

The Fourier covariance of a dilation is

```text
integral_T Q(Ny)f(y)dy-mu integral_T f
 =sum_(n!=0) Q_hat(n) f_hat(-nN).                          (18)
```

The right side is absolutely summable by (17).  Therefore

```text
norm of (18)
 <=V_Q Var(f)/(4 pi^2 N) sum_(n!=0)1/n^2
 =V_Q Var(f)/(12N).                                       (19)
```

Take the better of (16) and (19), then set `(f,N)=(b,N_k)` and
`(s,N_k)` to get (10)--(11).
The reverse triangle inequality and (12) give

```text
||q_k||>=mu sqrt(D)/2,

M_k>=mu M/2,
```

which is (13).  QED.

For the LRC step functions the hypotheses are automatic after the physical
speeds and the Boolean atom union are fixed.  If `a,f:T->[0,13]` are the
two root-count functions and `s=af`, then the useful side estimate

```text
V_s<=13 Var(a)+13 Var(f).                                 (20)
```

follows from the product rule for BV functions.

## 3. Endpoint disintegration and the common root word

Use the common oriented physical section

```text
phi(y,u)=(y+u)/13,                  u in F_13.             (21)
```

At endpoint clock `k>=1`, put

```text
W_k(x)=Q(13^k x).
```

Periodicity gives the exact identity

```text
W_k(phi(y,u))
 =Q(13^(k-1)y+13^(k-1)u)
 =Q(13^(k-1)y)
 =:w_k(y).                                                (22)
```

Thus `W_k` is constant on every physical `C_13` root fibre.  It can be
inserted on both semantic root packets, where idempotence makes the common
service filter

```text
w_k(y)^2=w_k(y).                                         (23)
```

Let `B(x)` be the finite-dimensional pointwise endpoint-table integrand
after the whole Boolean atom union has been formed and before smoothing or
Fourier expansion.  Define its root-base disintegration

```text
b(y)=(1/13)sum_(u in F_13) B(phi(y,u)).                    (24)
```

Then transfer along the thirteen inverse branches and (22) give

```text
integral_T W_k(x)B(x)dx
 =integral_T w_k(y)b(y)dy.                                (25)
```

This is the drift density used in Sections 1--2.  It and the service
density must use this same oriented root base.  A separately averaged or
unoriented density is not an input to the theorem.  Equation (25) also
explains why the endpoint clock `k` has mixing scale `N_k=13^(k-1)`.

The word remains target-neutral under the THM-2350/2365 target co-shift:
the transported shift is multiplied by `13^k/13`, an integer.  A source
shift is different.  One must freeze a lawful parity and use the shifted
word family of THM-2442; holding an unshifted word fixed under a source
translation is not licensed by this theorem.

## 4. Application to one THM-2459 atom union

Assume the complete common-root hypotheses of THM-2457 and THM-2459.
Also supply one THM-2305 owner `j` and assume that both the fixed endpoint
packet and the semantic `A/F` packets to which the atom bank is applied are
supported on its positive exclusive-owner event `E_j` before terminal
transport.  Thus there are `128` complete Boolean atoms, an aggregate
endpoint drift

```text
D_0>0,                                                       (26)
```

and a positive total directed root-service mass

```text
M_0>0.                                                       (27)
```

THM-2459 supplies a Boolean union `I` of at most four atoms for which

```text
D_I>=D_0/63001,

M_I>=M_0/16384.                                             (28)
```

Freeze that union before choosing a new clock.  Form the whole Boolean
endpoint packet first, apply the same THM-2365 noncirculant/H-drift
projection to its pointwise table integrand, and disintegrate it as in (24).
This gives a full-circle finite-dimensional Hilbert-valued rational step
function `b_I` on the same base as service, satisfying

```text
integral b_I=q_I,                  ||q_I||^2=D_I.           (29)
```

On the common root base put

```text
a_I(y)=sum_u A_I(y,u),

f_I(y)=sum_u F_I(y,u),

s_I(y)=a_I(y)f_I(y),               integral s_I=M_I.        (30)
```

THM-2457 allows its rational finite-interval base `Y` to be a proper subset
of the circle.  Extend `a_I,f_I,s_I` by zero off `Y`; then `b_I,s_I` are
honestly full-circle functions in the lemma.  They are fixed rational step
functions once `I` and the physical root section are fixed.  The
owner-support hypothesis on both endpoint and semantic packets is not used
by the BV estimate; it is used to interpret the later terminal stratum as
the outgoing word of this source packet.

Let `Q=1_(Q_(j,sigma))` be **any** positive rational THM-2305 terminal
word for a supplied owner `j`.  Apply Sections 1--3 to `b_I,s_I,Q`.
For every sufficiently large clock `k`, the enhanced union

```text
P_(I,k)=P_I W_k                                             (31)
```

still consists of at most four **enhanced** atoms, all carrying the same
semantic owner-to-word label `(j,sigma)`, and obeys

```text
D_(I,k)>=mu(Q)^2 D_0/252004,

M_(I,k)>=mu(Q) M_0/32768.                                  (32)
```

Here

```text
252004=4*63001,

32768=2*16384.                                              (33)
```

This is stronger than selecting a largest word block at one frozen clock.
Every fixed positive word eventually receives its asymptotic share of both
fixed densities.

If the owner-support hypothesis is removed, all analytic statements and
bounds remain true, but (31) is only a common terminal-stratum filter.  It
must not be called the THM-2305 coupling `(j,sigma)`; multiplying by
`Q_(j,sigma)` does not manufacture the source event `E_j`.

The common word filter is multiplied into the whole Boolean packets before
any Poisson--Abel smoothing or Fourier expansion.  Idempotence (23) then
preserves the THM-2457 disjoint-root hypotheses.
Applying its root-current bounds inside (31) gives

```text
sum_(r!=0)|J_(I,k)(r)|^2
 >=mu(Q)^2 M_0^2/368005682823168,

max_(r!=0)|J_(I,k)(r)|
 >=mu(Q) M_0/66453504.                                     (34)
```

The denominators are

```text
368005682823168=32768^2*342732,

66453504=32768*2028.                                       (35)
```

Qualitatively, every nonzero root colour supplied by THM-2457 survives.
Independently, because `b_I` used the identical THM-2365 drift projection,
positive drift in (32) feeds that theorem and permits a fresh
eligible target/deep colour and then a fresh exact relation address with a
`91`-unit deepest harmonic.  The two selections now occur inside one
semantic-word-filtered Boolean packet.  They are not asserted to be the
same Fourier colour or the same relation coordinate.  A lawful retained
septimal source colour still requires THM-2442's phasewise shifted word
family, as stated in Sections 1 and 3.

## 5. What this repairs in THM-2460

THM-2460's split-base hostile has the quantifiers

```text
fix one clock and its three word cells;
then place drift in one cell and service in another.               (36)
```

This theorem has the opposite, physically available order:

```text
fix the drift and service densities;
then choose one common word clock sufficiently late.               (37)
```

Equation (6) forbids a fixed positive density from avoiding the same
positive word along every sufficiently late clock.  The hostile therefore
does not refute delayed simultaneous retention.  It does prove that the
delay and the fixed-density hypothesis are load-bearing.

The result also sharpens the clock-cell viewpoint of THM-2464.  A strict
open clock cell gives pointwise realizability after a mesh threshold.
THM-2466 instead gives positive integrated drift and service even when no
single preselected phase cell is retained.  The two mechanisms solve
different problems:

```text
THM-2464: preserve finitely many prescribed pointwise truth values;

THM-2466: preserve two fixed nonzero integrated observables.        (38)
```

If the terminal word and a blocker requirement share one physical factor
or one clock, THM-2464's joint-cell condition still applies.  BV mixing
cannot make an empty pointwise cylinder nonempty.

## 6. Sharp boundaries

### 6.1 No prescribed-clock conclusion

Fix any finite dilation `N` and any nonconstant Boolean word `Q`.  Take

```text
b=s=1-Q(Ny).                                                   (39)
```

Both unfiltered integrals are positive, while

```text
integral Q(Ny)b(y)dy=integral Q(Ny)s(y)dy=0.                   (40)
```

After `b,s` are frozen, later dilations recover (6).  Hence no theorem
from positivity alone can replace "sufficiently delayed" by the prescribed
first-expiration clock.

### 6.2 The densities must be fixed before the delay

If one is allowed to choose

```text
s_k(y)=1-Q(13^(k-1)y),                                      (41)
```

then every `s_k` has the same positive mean but `M_k=0` for every `k`.
The same construction works for drift.  Thus (6) is not uniform over
clock-dependent densities, even though each individual density is a
rational step function.

### 6.3 Positive marginals and a common root gauge are essential

If `D=0`, `M=0`, or `mu(Q)=0`, the corresponding conclusion is false.
If the endpoint packet and semantic descendants are placed on unrelated
root sections, (22)--(25) are unavailable: an endpoint word need not be
the common filter appearing in the service integral.  This theorem does
not construct or identify that section.  All statements ignore the finite
null set of step-function seams, and all variations are periodic circle
variations.

### 6.4 A later word is not a first word

For a genuine prescribed THM-2305 arm, the first expiration is

```text
k=lambda_j+1.                                               (42)
```

The threshold in (12) may be much larger and depends on the variations
and positive margins of the selected packet.  Its semantic label is
honest, but its time is a delayed restart.  Renaming it the prescribed
first word would erase the temporal coordinate isolated by THM-2461.

## 7. Exact gain and remaining frontier

The theorem closes the abstract analytic part of the one-word coupling:

```text
supplied common root section
+ supplied positive owner-supported packet E_j
+ fixed nonzero endpoint drift density
+ fixed positive root-service density
+ any fixed positive semantic word

  -> one sufficiently delayed copy of that word retains all three. (43)
```

What remains is no longer an asymptotic word-selection problem.  The
missing LRC data are structural:

1. identify a canonical owner and common oriented root section for the
   actual selected endpoint packet and its semantic descendants;
2. prove that the resulting service density is positive, or force an
   empty semantic joint clock cell;
3. reconnect the retained root current and endpoint drift to one lawful
   relation-current sidecar without confusing their Fourier colours;
4. if first-expiration minimality is required, supply it separately -- BV
   delay cannot do so.

No scalar valuation profile is excluded.  The ledger remains `165`, and
LRC(14) remains open.

## 8. Exact companion

Run

```text
python3 04-computation/lrc14_delayed_word_joint_retention_thm2466.py
python3 -O 04-computation/lrc14_delayed_word_joint_retention_thm2466.py
```

The dependency-free `Fraction` companion verifies:

- exact root constancy and the one-power clock offset;
- representative scalar drift/service dilation integrals and both the
  Boolean-discrepancy and `1/(12N)` BV inequalities;
- a frozen-density first-clock hostile followed by delayed recovery;
- the clock-dependent-density hostile;
- the THM-2459/2457 quantitative denominators; and
- target neutrality and the two cofinal septimal parity classes.

Every truth-bearing executable check uses `require`; there are no `assert`
statements or floating-point truth tests.  Normal and optimized runs must
reproduce the stored transcript byte for byte after LF normalization.

QED.
