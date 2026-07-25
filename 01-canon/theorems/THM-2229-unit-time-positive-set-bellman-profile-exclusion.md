---
id: THM-2229
title: "Unit-time positive-set Bellman sieve and 218-profile exclusion"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. In the scalar five-unit/
  three-blocker branch, retain the positive unit residual itself as one
  clause and every available even divided-blocker union as another. Resolving
  time one factor of thirteen at a time gives the pointwise root caps 10/13
  for the positive residual and 2/13 for each unit danger atom. An exact
  arbitrary-coupling Bellman relaxation closes 436 profiles: respectively
  106, 107, 119, and 104 at first depths two through five. Its overlap with
  the ten shallow exact exclusions is exactly (3,3,4), and it contains all
  217 THM-2226 exclusions. Relative to the pre-THM-2227 ledger it independently
  removes 218 profiles and reduces 458 to 240; because THM-2227 already closes
  the same six high-first profiles, its additive current decrement is 212,
  from 452 to 240. Exact primal and dual vertex enumerations agree, and
  ordinary and optimized outputs are byte-identical. The remaining profiles
  are 238 low-first rows and (4,6,8),(5,7,9). LRC(14) remains open.
source: klein-2026-07-24-unit-time-positive-set-bellman
depends_on:
  - THM-2219-scalar-depth-four-sparse-tail-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2224-transfer-owner-word-temporal-union-bound
  - THM-2226-three-checkpoint-bellman-sieve-and-eight-profile-residue
  - THM-2227-sharp-parity-three-checkpoint-bellman-profile-exclusion
related:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
script: 04-computation/lrc14_unit_time_positive_set_bellman_thm2229.py
output: 05-knowledge/results/lrc14_unit_time_positive_set_bellman_thm2229.out
script_sha256: c57eb84a7e8e3374482e4006be78c1fe25b1b32e7546f3d15a82e4a9b1df730c
output_sha256: d163befa29b7424f69fd4f6c03165d18737bc3fdaca988cf0c8491f60e50f4e1
hash_basis: working-tree bytes (LF)
---

# THM-2229 -- the unit-time positive-set Bellman sieve

The multiplication-by-`169` clause sieve of THM-2226 sees blocker
checkpoints but forgets the set whose mass supplies the contradiction.
Retaining that positive set, and resolving each step as two separate
multiplications by thirteen, closes almost half of the remaining scalar
ledger.

## 1. The same positive set at every even checkpoint

On `T=R/Z`, put

```text
D_a={x:||ax||<1/14},
C_H={x:||Hx||>1/7}.                                  (1)
```

In the scalar `5+3` branch write the three blocker coefficients as

```text
c_j=13^(lambda_j)u_j,        13 does not divide u_j,
1<=lambda_1<=lambda_2<lambda_3<=19.                  (2)
```

For the five unit terminal coefficients `q_i`, let

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),

A_+={R>0}
   =C_H\union_(i=1)^5D_(q_i).                        (3)
```

THM-2198 gives

```text
measure(A_+)>=delta_5=961/6930.                      (4)
```

For `0<=k<=lambda_1`, set

```text
U_k=union_(j=1)^3 D_(c_j/13^k).                     (5)
```

The transfer-parity theorem THM-2222 gives the almost-everywhere inclusions

```text
A_+ subset U_k                     for every even k<=lambda_1. (6)
```

The point is that all the sets in (6) live at the same point `x`. They are
not unrelated images of `A_+`.

## 2. Resolve one factor of thirteen at a time

Let

```text
S(x)=13x mod 1,
X_(j,t)(x)=1_(D_(u_j))(S^t x).                       (7)
```

Exact divisibility gives

```text
1_(D_(c_j/13^k))(x)
 =X_(j,lambda_j-k)(x).                              (8)
```

Thus for every available even checkpoint,

```text
U_k={x: OR_(j=1)^3 X_(j,lambda_j-k)(x)=1}.          (9)
```

This is the literal unit-time schedule. A coefficient does not merely have
a parity label after a `169`-step; its unit danger atom occurs at the exact
integer time `lambda_j-k`.

There are two pointwise root-fibre caps. For a thirteen-unit `u`, THM-2198's
root count is

```text
sum_(Sx=y)1_(D_u)(x)=2-1_(D_u)(y)<=2.               (10)
```

Also `A_+ subset C_H`, and the guard-root law gives

```text
sum_(Sx=y)1_(A_+)(x)
 <=sum_(Sx=y)1_(C_H)(x)
 =10-1_(C_H)(y)
 <=10.                                              (11)
```

These statements hold away from the finite endpoint set. Consequently an
`A_+` atom occupies at most `10/13` of every root fibre, while each unit
danger atom occupies at most `2/13`. Equation (11) is pointwise in the
terminal point; it is not an averaged measure estimate.

## 3. The positive-set clause Bellman relaxation

Fix a profile with `2<=lambda_1<=5`, and let

```text
E={0,2,...,2 floor(lambda_1/2)}.                     (12)
```

The clause set is

```text
Omega={star} union E.                               (13)
```

Clause `star` asks for membership in `A_+` at time zero. Clause `k in E`
asks that at least one of the three atoms at times

```text
lambda_1-k, lambda_2-k, lambda_3-k                  (14)
```

occur. By (6) and (9),

```text
A_+
 =A_+ intersection intersection_(k in E)U_k         (15)
```

almost everywhere.

For a function `v:2^Omega -> R`, eliminate one time slice as follows. If
the active clauses have marginal caps `q_r`, put

```text
(mathcal E_t v)(B)
 =max sum_(Z subset Omega)p_Z v(B union Z),          (16)
```

where

```text
p_Z>=0,
sum_Z p_Z=1,
sum_(Z:r in Z)p_Z<=q_r.                             (17)
```

Inactive clauses never occur in `Z`. At time zero, clause `star` has cap
`10/13`. Every active blocker atom contributes `2/13` to its checkpoint
clause; if several atoms serve the same clause at one time, the union bound
adds their caps. There are at most three, so every such cap is at most
`6/13`.

Start from

```text
v_(-1)(B)=1_(B=Omega)                               (18)
```

and apply (16) successively over all unit times. The backward-root induction
of THM-2224 applies verbatim. For a fixed terminal point, the actual law of
the set of clauses hit by its thirteen roots satisfies (17), by
(10)--(11). Conditional on one chosen root, the earlier chain is bounded by
the previous Bellman value. Averaging and then enlarging to every joint law
allowed by (17) proves

```text
measure(
  A_+ intersection intersection_(k in E)U_k
 )
 <=B(lambda_1,lambda_2,lambda_3),                   (19)
```

where `B` is the exact rational output of the finite recursion.

No independence is assumed. Coincident atoms, shared roots, arbitrary
cross-core correlation, and every owner tie are included in the admissible
joint laws. Combining (4), (15), and (19), every profile satisfying

```text
B(lambda_1,lambda_2,lambda_3)<961/6930              (20)
```

is empty.

## 4. Exact profile classification

The exact recursion gives

```text
first depth 2:  106/153 profiles close,
first depth 3:  107/136 profiles close,
first depth 4:  119/120 profiles close,
first depth 5:  104/105 profiles close.             (21)
```

At first depth two, the closed profiles are exactly

```text
lambda_2=3, lambda_3>=6,

or

lambda_2>=5, lambda_3 != lambda_2+2.                (22)
```

At first depth three, they are exactly

```text
lambda_2=3 and (lambda_3=4 or lambda_3>=6),

or lambda_2=4 and lambda_3>=7,

or lambda_2>=6 and lambda_3 != lambda_2+2.          (23)
```

At first depths four and five only

```text
(4,6,8),        (5,7,9)                             (24)
```

fail the strict Bellman test. The worst passing values at first depths
two through five are, respectively,

```text
625252/4826809,
49000/371293,
2896/28561,
74441360/815730721,                                 (25)
```

all strictly below `961/6930`. The least failing values are

```text
4480/28561,
824660/4826809,
60792/371293,
59340/371293.                                       (26)
```

Thus the strict boundary is not a rounding artifact.

The raw set in (21) has `436` profiles. Its exact overlap with the ten
already closed profiles of deepest depth at most four is

```text
{(3,3,4)}.                                          (27)
```

It contains every one of THM-2226's `217` exclusions and is disjoint from
THM-2224's `455` first-depth-at-least-six exclusions. Among THM-2226's
eight residual profiles, it additionally closes

```text
(4,6,9), (4,7,8), (4,7,9),
(5,7,10), (5,8,9), (5,8,10).                       (28)
```

Therefore, relative to the `458`-profile ledger before THM-2227, this
theorem independently closes

```text
436-1-217=218                                       (29)
```

new profiles and leaves `240`. THM-2227 is an independent hidden-state
proof of the six profiles in (28), so after that theorem's promotion the
additive decrement here is `212`, from `452` to the same exact residue

```text
240 = 238 low-first profiles + {(4,6,8),(5,7,9)}.   (30)
```

The companion freezes the complete profile sets and their SHA-256 digests,
so (29)--(30) are checked as set identities, not obtained by subtracting
possibly overlapping headline counts.

## 5. Why the private-owner image pump does not itself give checkpoints

Two tempting refinements explain why the successful state in Section 3 is
`A_+`, rather than an unlabelled owner-deleted image.

First, let `a` be a thirteen-unit. The centered danger function

```text
g_a=1_(D_a)-1/7                                    (31)
```

satisfies

```text
L g_a=-g_a,                                         (32)
```

where `L` is the unnormalized thirteen-root transfer operator. Since
`LR=-R`, the centered owner-deleted residual

```text
F_a=R-g_a=R-1_(D_a)+1/7                             (33)
```

also lies on the `-1` eigenline. Its top level is exactly

```text
{F_a=8/7}=A_+\D_a.                                  (34)
```

Indeed `R<=1`, equality means membership in `A_+`, and every point outside
(34) has `F_a<=1/7`.

But the indicator-relevant threshold is

```text
h_a=F_a-1/7=R-1_(D_a),
{h_a>0}=A_+\D_a.                                    (35)
```

It no longer lies on the eigenline:

```text
L^m h_a
 =(-1)^m(h_a+1/7)-13^m/7.                           (36)
```

For every positive even `m`, the right side is already negative at the top
value `h_a=1`. Thus centering identifies the correct owner-deleted level
set, but thresholding introduces a growing constant and does not reproduce
THM-2222's same-set parity inclusion.

Second, THM-2198's image-pump inclusions for the post-owner remainder are
genuine image statements, but pulling them back does not multiply clauses.
For a remaining blocker `b` divisible by `13^u`,

```text
(S^u)^(-1)D_(b/13^u)=D_b.                           (37)
```

Hence every pulled-back inclusion over the shallow-owner gap is the same
original two-blocker cover. It is not a sequence of distinct checkpoints.
The quantitative remainder floor `2593/69300` remains useful, but advancing
it needs the private label/carry, a post-deletion expansion law, or another
sidecar. Equations (36)--(37) are the first precise obstruction to treating
the image pump as a free Bellman recursion in this carrier.

## 6. Connection and loss ledger

The carrier is

```text
source:
  one scalar 5+3 cover and its positive unit residual A_+;

map:
  write each divided blocker as a unit danger observed at its exact
  multiplication-by-thirteen time, and make A_+ plus every available
  even checkpoint into labelled clauses;

preserved:
  the same positive set, every even parity inclusion, exact integer time
  offsets, clause labels, and pointwise root-fibre marginal caps;

destroyed:
  parent-bit dependence inside the exact 1-versus-2 and 9-versus-10 root
  counts, cross-core phase, owner identity within a clause, winding, and
  all higher joint correlations;

needed sidecar:
  a labelled root bit/carry, the integral correlation phase of THM-2218,
  or a genuine post-owner floor/expansion law;

cheapest hostile controls:
  (4,6,8) and (5,7,9), whose exact relaxed bounds in (26) remain above
  the target.                                             (38)
```

The useful refinement over THM-2224/2226 is not a sharper static cap. It is
the change of scale from `169` to `13` together with the otherwise discarded
positive-set clause. The toothpick schedule is thereby resolved at its
natural digit scale.

## 7. Exact audit and scope

Run

```bash
python3 04-computation/lrc14_unit_time_positive_set_bellman_thm2229.py
python3 -O 04-computation/lrc14_unit_time_positive_set_bellman_thm2229.py
```

For every marginal LP the companion independently enumerates every primal
basic feasible solution and every dual vertex, requiring exact equality.
All arithmetic is rational. Ordinary and optimized transcripts are
byte-identical to the stored output.

The script checks:

```text
all 1,140 legal valuation profiles in the ledger;
all 514 profiles with first depth two through five;
the four per-depth closure counts in (21);
the exact classifications (22)--(24);
the worst-passing and least-failing fractions (25)--(26);
the ten-profile, THM-2224, THM-2226, and THM-2227 intersections;
the independent 218-profile and additive 212-profile decrements;
the complete 240-profile residue and two SHA-256 set digests. (39)
```

This theorem does not close any first-depth-one profile, does not close the
two rows in (24), and does not prove LRC(14). QED.
