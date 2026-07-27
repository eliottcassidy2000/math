---
id: THM-2461
title: "Temporal blocker-word cocycle and diagonal-polarized repair boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. At the
  untwisted all-safe section, the complete local atom bank contains
  the source owner and only one pure target word; the other pure word
  and fork require the fixed deepest-safe bit to change. The correct
  semantic object is therefore a prescribed-clock atom-to-word
  coupling, not a state label. An exact strict-profile local hostile
  has three points in the same complete source atom with the same
  two-digit owner address but terminal words fork, first-pure, and
  second-pure. The THM-2305 predecessor root and THM-2365 deep-probe
  root are distinct typed coordinates. Same-time idempotence cannot
  be transported through an expanding clock except for null/full
  atoms. Common co-shifts keep a danger/safe pair identically zero,
  whereas THM-2379 repair is necessarily polarized and off-diagonal.
  Only the two graft units have nonzero target covectors; after
  q_*=k_b is fixed, only k_a among the five first-failure roles is
  target-active. No semantic/root intertwiner, scalar-row exclusion,
  or LRC(14) consequence is proved.
source: codex-2026-07-26-temporal-word-repair-boundary
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
  - THM-2350-owner-pivot-dual-dipole-normal-form
  - THM-2379-anchored-guard-unit-deletion-factor-repair-current
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
related:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2456-two-root-replica-uniform-offset-boundary
script: 04-computation/lrc14_temporal_word_polarized_repair_boundary_thm2461.py
output: 05-knowledge/results/lrc14_temporal_word_polarized_repair_boundary_thm2461.out
script_sha256: 6ba42d0c613378c0f55254dbaaa4af01a2b9d617814958d302935be1baa683ba
output_sha256: 16564f61e439311131670b0425242b63e893e95e3dcebd1809a8542caf7ba479
hash_basis: working-tree bytes (LF)
---

# THM-2461 -- blocker words are temporal and repair is polarized

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2452 copies a complete Boolean atom from one endpoint copy to the
other at the same physical time. THM-2305's word, by contrast, is read
after a prescribed expanding clock, and THM-2379's repair compares a
failed factor with a separately shifted safe complement.

These are three different operations:

```text
same-time state:       P_tau(x)P_tau(x)=P_tau(x);

temporal handoff:      E_j(x)Q_(j,sigma)(T^(k_j)x);

polarized repair:      d_f(x)g_f(x-s/13).                    (1)
```

The theorem proves that neither of the last two objects can be recovered
from the first by relabelling. It replaces the attempted state label by a
directed temporal coupling, separates two root coordinates that both happen
to be indexed by `F_13`, and gives the exact target-covector gate on repair.

## 1. The untwisted local truth table

Use THM-2305 notation. Let `j` be the source blocker, let `a` be the
other nondeep blocker, and write `c=c_3`. Put

```text
g_v=1-d_v.
```

At the untwisted section where the guard and five ordinary roles are all
safe and the deepest blocker is fixed safe, the four local blocker atoms
are

```text
L_(epsilon_j,epsilon_a)
 =1_(A_0)
   d_j^(epsilon_j)g_j^(1-epsilon_j)
   d_a^(epsilon_a)g_a^(1-epsilon_a)
   g_c.                                                   (2)
```

The scalar cover and the definitions of the exclusive owners give

```text
L_00=0,

L_10=1_(E_j),

L_01=1_(E_a)=1_(Q_(j,{a})),

L_11=1_(A_0 intersection D_j intersection D_a
            intersection D_c^c).                         (3)
```

Every complete atom with a guard/unit failure is disjoint from `A_0`.
The other two THM-2305 word strata are

```text
Q_(j,{c})   subset D_c,

Q_(j,{a,c}) subset D_c.                                (4)
```

They cannot occur in the bank (2), which fixes `g_c`. Thus the same-time
bank contains:

- the source state `E_j`;
- the pure target state `Q_(j,{a})`;
- no local copy of `Q_(j,{c})` or the fork `Q_(j,{a,c})`.

Equation (3) is deliberately untwisted. Freezing it while the target
action moves one of its factors would repeat MISTAKE-266.

## 2. The correct object is a prescribed-clock coupling

Let `{P_tau}` be any complete Boolean partition of the source circle.
For

```text
sigma in {{a},{c},{a,c}},
```

define

```text
W_(k_j)(tau,sigma)
 =mu(P_tau intersection E_j
          intersection T^(-k_j)Q_(j,sigma)).               (5)
```

This is a nonnegative directed atom-to-word matrix. THM-2305's word
partition gives the exact mass identity

```text
sum_(tau,sigma) W_(k_j)(tau,sigma)=rho_j.             (6)
```

For the complete local bank in Section 1, `E_j` lies in the unique source
atom

```text
tau_src=(all guard/unit safe, j dangerous, a safe).
```

Therefore

```text
W_(k_j)(tau,sigma)=0              for tau!=tau_src,

W_(k_j)(tau_src,sigma)=rho_(j,sigma).                (7)
```

THM-2305 proves

```text
max_sigma W_(k_j)(tau_src,sigma)>=rho_j/3,            (8)
```

but it does not make the row one-sparse. A terminal word is an outgoing
edge label of the source atom, not a function of that atom.

THM-2306 already retains the corresponding exact current:

```text
f_sigma(y)
 =1_(Q_(j,sigma))(y) P^(k_j)1_(E_j)(y),

integral f_sigma=rho_(j,sigma).                       (9)
```

Equation (5) is the state-resolved version of (9). It is the minimal
semantic replacement for an invalid map `tau -> sigma`.

## 3. Two roots with different types

Write

```text
lambda=lambda_j,

G_j=P^lambda 1_(E_j).
```

The THM-2305/2306 predecessor-root packet is

```text
v_(r_owner)(y)
 =G_j((y+r_owner)/13),

P G_j(y)
 =(1/13)sum_(r_owner in F_13)v_(r_owner)(y).          (10)
```

Here `r_owner` labels the last inverse branch of the source ancestry over
the terminal base point `y`. The word current is

```text
f_sigma(y)
 =(1/13)1_(Q_(j,sigma))(y)
       sum_(r_owner)v_(r_owner)(y).                  (11)
```

The THM-2365/2452 root has a different definition:

```text
Delta_(r_deep)(x)
 =d(c_3x-r_deep/13).                                 (12)
```

It labels a translated deepest-probe phase at the integration point `x`.
Equations (10) and (12) have the same index set only because both use the
prime thirteen. Their base variables, maps, and observables differ:

```text
r_owner:
  inverse ancestry sheet over y;

r_deep:
  deepest danger translate at x.                    (13)
```

No theorem identifies them. A lawful bridge must retain both on one
physical fibre product or prove a base-dependent affine intertwiner.
Matching their printed labels would repeat the anchor error in
MISTAKE-260.

## 4. Same-time idempotence cannot be copied through time

Let

```text
T(x)=13x mod 1
```

on the circle, and let `P` be a Boolean measurable function. Suppose for
some `k>=1` that

```text
P(x)P(T^k x)=P(x)                                  (14)
```

almost everywhere. Then

```text
P<=P composed_with T^k.
```

Multiplication by `13^k` preserves Haar measure, so the two sides have
the same integral. They are therefore equal almost everywhere:

```text
P=P composed_with T^k.                              (15)
```

The map `T^k` is ergodic. For example, if `q=13^k`, Fourier invariance
gives

```text
P_hat(n)=0                    when q does not divide n,

P_hat(qn)=P_hat(n).
```

Iteration, or square summability along every nonzero `q`-ray, forces every
nonzero coefficient to vanish. Hence `P` is constant almost everywhere,
and

```text
mu(P) in {0,1}.                                     (16)
```

Thus the THM-2452 identity

```text
P_tau(x)^2=P_tau(x)
```

cannot be replaced by

```text
P_tau(x)P_tau(T^k x)=P_tau(x)
```

for a nontrivial atom. This does not forbid a positive temporal overlap;
it proves that the diagonal endpoint copy does not supply that overlap
for free.

## 5. Exact three-word hostile with one source atom and address

The failure is already visible on one exact strict local profile. Take

```text
H=1,

(q_1,...,q_5)=(4,2,3,6,10),

(c_1,c_2,c_3)=(13,13^3,2*13^5)
              =(13,2197,742586).                    (17)
```

The nine speeds are positive, distinct, primitive, `H` is odd, every
`q_i` is a thirteen-unit, and

```text
(nu_13(c_1),nu_13(c_2),nu_13(c_3))=(1,3,5).         (18)
```

Select `j=1`; its prescribed clock is

```text
k_j=2,

T^(k_j)x=169x mod 1.                                (19)
```

With common denominator `2222`, put

```text
x_1=862/2222,
x_2=863/2222,
x_3=865/2222.                                      (20)
```

At all three points the guard and five ordinary roles are safe, `c_1`
is dangerous, and `c_2,c_3` are safe. If `q_*=q_1`, their seven-bit
complete source mask is

```text
(0,0,0,0,0,1,0).                                   (21)
```

Their first two base-thirteen digits are also identical:

```text
(floor(13x_i),
 floor(13{13x_i}))=(5,0).                           (22)
```

But their terminal residues are

```text
169*(862,863,865)
 = (1248,1417,1755) mod 2222,                       (23)
```

and the exact THM-2305 terminal words are

```text
x_1: {c_2,c_3},

x_2: {c_2},

x_3: {c_3}.                                        (24)
```

The minimum normalized distance from any source or terminal inequality
to its strict boundary is respectively

```text
177/15554,

111/7777,

83/7777.                                           (25)
```

Thus small rational intervals around the three points preserve the
patterns: this is not an endpoint artifact.

The owner address in (22) still does not determine the deep-probe roots.
The roots `r` for which

```text
||c_3x_i-r/13||<1/14
```

are respectively

```text
{11,12}, {1,2}, {6,7}.                              (26)
```

This directly exhibits the type distinction in Section 3.

The hostile is intentionally scoped as a **local-profile** obstruction,
not a global scalar-cover row. Indeed,

```text
x_0=159/1111
```

is guard/ordinary safe and safe for all three blockers. On denominator
`2222` its nine centered distance numerators are

```text
(318,950,636,954,314,958,310,938,702).              (27)
```

So (17) does not satisfy the global cover hypothesis. Equations
(20)--(26) refute a local state-to-word rule; they do not exclude a scalar
profile or change the `165`-row ledger. THM-2457 independently gives a
live-row pure/fork hostile with fixed deep labels.

## 6. Common shifts are diagonal; repair is polarized

Use THM-2379 notation:

```text
d_L(y)=1_(||y||<L/14),

u_L=1-d_L,                  L in {1,2}.              (28)
```

A common shift of a failed factor and its safe complement is always
diagonal, pointwise:

```text
d_L(y-s/13)u_L(y-s/13)=0
```

for every `s`, almost everywhere. This is precisely the algebra of a
matched same-time endpoint copy.

THM-2379 instead fixes the failure at shift zero inside a nonnegative
weight `w` and translates only the displayed complement:

```text
support(w) subset {d_1(cx-alpha)=0}
                    intersection
                  {d_L(vx-beta)=1},

K^+_(r,s)
 =integral w(x)
    d_1(cx-alpha-r/13)
    u_L(vx-beta-s/13) dx.                           (29)
```

The failure anchor gives

```text
K^+_(r,0)=0                       for every r.       (30)
```

Yet the exact translated counts are

```text
sum_r d_1(y-r/13)=2-d_1(13y),

sum_s u_L(y-s/13)=13-2L+d_L(13y).                  (31)
```

If

```text
rho=integral w>0,
```

then (29)--(31) give

```text
sum_(r,s)K^+_(r,s)
 >=(13-2L)rho>0.                                   (32)
```

All this mass is off the plane `s=0`. On a shift cell containing the
unshifted failure, the number of available nonzero safe translates is at
least

```text
11 for L=1,

9 for L=2,                                         (33)
```

and both constants are attained.

Thus a repair is not another diagonal Boolean copy. It is a polarized
edge from a fixed danger anchor to an independently shifted safe
complement. THM-2379 converts that edge into its positive mixed Fourier
coefficient; THM-2452's common co-shift does not create the edge.

## 7. Only paired graft units are target-active

THM-2350 gives the two lawful target covectors in the unique omitted-unit
gauge `ell_(u_0)=0`:

```text
eta=e_a-e_(k_a),

ell=e_c-e_(k_b).                                    (34)
```

For a role `f`, let

```text
lambda_f=(eta_f,ell_f).
```

Modulo the common source-translation/scalar gauge, the exact table is

```text
lambda_a    =( 1, 0),
lambda_c    =( 0, 1),
lambda_(k_a)=(-1, 0),
lambda_(k_b)=( 0,-1),

lambda_f=(0,0)
 for the omitted unit, three ungrafted units, and source j.       (35)
```

Therefore among the six guard/unit roles, exactly `k_a,k_b` are
target-active.

Here `k_a,k_b` are distinct and neither is the omitted unit, as required
by the owner pivot. THM-2445 chooses the deepest graft

```text
q_*=k_b
```

and makes the five other guard/unit roles its first-failure bank. By
(35), only

```text
k_a
```

among those five failures has a nonzero target covector. The other four
are target-neutral. A generic THM-2445 first-failure label therefore
cannot be declared a lawful target repair.

Even for `k_a` or `k_b`, THM-2379 translates the displayed unit complement
while holding its paired blocker inside the weight. A lawful target action
in (34) co-translates that blocker in the opposite direction. The
THM-2379 coefficient is consequently an auxiliary factor-repair current,
not automatically the THM-2334/2365 target current.

There is nevertheless an exact algebraic pullback when a physical packet
does realize the repair shift through a nonzero target covector. Let

```text
V=F_13^2,

lambda_f:V->F_13
```

be nonzero and define

```text
Ktilde(r,delta)=K^+(r,lambda_f(delta)).              (36)
```

With normalized transforms and the same sign convention as THM-2379,
finite character orthogonality gives

```text
Ktilde_hat(a,b lambda_f)=K^+_hat(a,b),

Ktilde_hat(a,xi)=0
        when xi notin span(lambda_f).                (37)
```

Indeed, every fibre of `lambda_f` has thirteen points; the kernel-character
sum is thirteen on the target line and zero off it, exactly cancelling the
extra normalization. If `lambda_f=0`, equation (30) instead makes
`Ktilde` identically zero.

Equations (36)--(37) are a conditional Fourier lift, not a claim that the
canonical packet has the form (36). They make the obstruction sharp:
target-neutral failures cannot carry repair charge, while target-active
failures still need the paired blocker to transform covariantly.

## 8. Exact stopping boundary

The proved chain is

```text
same-time complete atom
  -/-> prescribed terminal word;

same printed F_13 label
  -/-> common owner/deep root;

common danger/safe co-shift
  -/-> polarized repair edge;

generic first-failure role
  -/-> lawful target direction.                    (38)
```

The strongest faithful survivor is the joint packet

```text
temporal coupling W_(k_j)(tau,sigma)
  + owner predecessor root r_owner
  + deep-probe root r_deep
  + target covector lambda_f.                       (39)
```

A closure theorem must either:

1. physicalize (37) on one common base and prove positive co-support;
2. produce a base-dependent affine intertwiner between the two roots;
3. coarsen to a semantic THM-2457 edge while retaining positive drift;
   or
4. realize the polarized repair together with its paired blocker dipole.

Local truth bits, root-cardinality coincidence, or same-time
idempotence do not supply any of these sidecars.

## 9. Exact companion

Run

```text
python3 04-computation/lrc14_temporal_word_polarized_repair_boundary_thm2461.py
python3 -O 04-computation/lrc14_temporal_word_polarized_repair_boundary_thm2461.py
```

The dependency-free companion uses integers and `Fraction` arithmetic
with explicit `require` checks. It:

- verifies distinctness, primitivity, parity, and the strict `(1,3,5)`
  valuation profile;
- checks every source and terminal distance in the three-word hostile,
  the common mask/address, exact margins, and deep-probe root sets;
- checks the explicit uncovered point and therefore the local-only scope;
- verifies the translated danger/safe counts on every rational shift
  cell and the sharp off-diagonal counts `11,9`;
- exhausts the thirteen-cycle analogue of the temporal idempotence
  no-go;
- checks a sharp three-word temporal coupling row; and
- verifies the complete target-covector table, the one-of-five
  first-failure gate, and all `168` nonzero covector pullbacks.

Both runs reproduce

```text
05-knowledge/results/lrc14_temporal_word_polarized_repair_boundary_thm2461.out
```

byte-for-byte.

An independent semantic audit rederived the local truth table, temporal
coupling, strict rational hostile, local-cover caveat, root typing, and
ergodic no-go. A separate repair audit rederived the diagonal/polarized
distinction and target-covector gate. No scalar row is excluded, the
ledger remains `165`, and LRC(14) remains open.

QED.
