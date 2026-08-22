---
id: THM-3692
title: "Ordinary normalized two-by-three shift-root peeling obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the normalized
  planar chart P=x plus exactly
  two distinct nonlinear monomials and Q=y plus exactly three distinct
  nonlinear monomials, suppose both P supports have positive x exponent and
  all three Q supports have positive y exponent.  Then Jac(P,Q) is not
  constant, with no hypothesis on the six cross determinants.  A
  degree-unbounded shift-root peel forces one of six scalar towers, each with
  an active terminal singleton.  Axis branches and JC(2) remain open; no
  counterexample is claimed.
source: jc-sparse-direct-search / 2026-08-22
audit: >
  PASS -- root independently checked bottom uniqueness under the additive
  order, self-edge activity at 2s, every Case I--III interval-forcing step,
  all six terminal multiplicity profiles, and the leading-edge parallel gate.
depends_on: []
related:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3587-consecutive-keller-fibre-factor-toric-and-coefficient-span-gates
  - THM-3689-fully-transverse-two-by-two-sparse-support-gate
  - THM-3690-complete-normalized-two-by-two-sparse-support-closure
script: 04-computation/jc2_ordinary_two_by_three_shift_root_peeling_thm3692.py
output: 05-knowledge/results/jc2_ordinary_two_by_three_shift_root_peeling_thm3692.out
script_sha256: 5b5e40a35a11fbd2fb0fc8b9076c46e4ecbd72f1770def28ee3682d739d759d8
output_sha256: 27ed345b91bcb0790229bf70756be4d5eeeb88c4f7da67e5bcee57b73649040e
semantic_sha256: 48c059b4da5666623bf5f82c7b99b7c92f5156ba05efab16f38788b50f14222e
hash_basis: raw LF bytes for files; forced scalar-tower debt profiles and two bounded support censuses for semantic hash
---

# THM-3692 -- ordinary normalized two-by-three shift-root peeling obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is the first
step beyond the complete two-by-two closure of THM-3690.  It replaces a large
activity-partition atlas by an acyclic translation grammar for the Jacobian
debts.

Work over a characteristic-zero field `K`.  For an exponent vector `u`, write
`X^u=x^(u_1)y^(u_2)`, and put `e_1=(1,0)`, `e_2=(0,1)`.

## 1. Statement

Let `p_1,p_2,q_1,q_2,q_3` be nonnegative integer exponent vectors such that

```text
p_1 != p_2,                         q_1,q_2,q_3 pairwise distinct,
|p_i|,|q_j| >= 2,
(p_i)_1 > 0 for i=1,2,              (q_j)_2 > 0 for j=1,2,3.       (1)
```

For arbitrary nonzero coefficients `a_i,b_j in K*`, set

```text
P=x+a_1 X^p_1+a_2 X^p_2,
Q=y+b_1 X^q_1+b_2 X^q_2+b_3 X^q_3.                       (2)
```

Then

```text
Jac(P,Q) is not constant.                                   (3)
```

No cross-transversality is assumed: any subset of
`det(p_i,q_j)` may vanish.  The hypotheses in the last line of `(1)` say only
that all five divergence labels are active.  The cases where a `p_i` is pure
`y`, or a `q_j` is pure `x`, are not proved by this theorem.

## 2. The shift-root debt grammar

Translate the supports by the differentiated linear terms:

```text
s_i=p_i-e_1,                         r_j=q_j-e_2.           (4)
```

By `(1)`, the two shifts `s_i` and three roots `r_j` are distinct within their
own families and belong to `N^2 minus {0}`.  Direct expansion gives

```text
Jac(P,Q)-1
 = sum_i a_i (p_i)_1 X^s_i
 + sum_j b_j (q_j)_2 X^r_j
 + sum_ij a_i b_j det(p_i,q_j) X^(s_i+r_j).               (5)
```

Thus the eleven potential debts are represented by

```text
two seed labels at S={s_1,s_2},
three root labels at R={r_1,r_2,r_3},
six edge labels at R+S.                                    (6)
```

Every seed and root coefficient in `(5)` is nonzero.  An edge may be inactive,
but only when its displayed determinant vanishes.  Since every label in `(5)`
has positive total degree, a constant Jacobian would equal `1`; every active
debt location would then need at least two labels.

Fix graded lexicographic order on `N^2`.  It is a translation-invariant strict
total order with `0<u` for every nonzero `u`.  Consequently

```text
r+s > r and r+s > s                                       (7)
```

for every root and shift.  This makes `(6)` an acyclic two-shift recurrence.

## 3. Bottom forcing

Let `s=min S` and `r=min R`.  If `s<r`, the seed at `s` is the unique least
active debt: the other seed and all roots are larger, while every edge is
larger than its shift and root.  If `r<s`, the root at `r` is similarly unique.
Therefore no-singleton cancellation forces

```text
min S=min R=s.                                             (8)
```

There is exactly one seed and one root at `s`.  The edge obtained from this
pair lies at `2s`.  It is unconditionally active, because its original
supports are `p=s+e_1`, `q=s+e_2`, and

```text
det(s+e_1,s+e_2)=s_1+s_2+1 != 0.                          (9)
```

Write `S={s,t}` and `R={s,u,v}`, with `u<v`.  Any other edge is the sum of a
root at least `s` and a shift at least `s`; translation invariance shows that
it can equal `2s` only when both summands equal `s`.  Hence the active self-edge
at `2s` needs a seed or root mate:

```text
t=2s or u=2s (or both).                                   (10)
```

Whenever the peel below forces a shift `ks` and a root `ell s`, the
corresponding edge is active.  Indeed

```text
det(ks+e_1, ell s+e_2)=k s_1+ell s_2+1 != 0               (11)
```

for positive integers `k,ell`.  This identity is the key guard against hiding
a terminal singleton on a parallel-exponent boundary.

## 4. Exhaustive three-case peel

All comparisons below use the fixed additive order.  A location strictly
between two displayed multiples cannot be rescued by a later edge, because an
edge is larger than each of its summands.

### Case I: `t=2s` and `u!=2s`

If `u<2s`, its root is singleton.  Thus `u>2s`.  The edge `s+2s=3s` is active
by `(11)`, so comparison with the next root forces `u=3s`.  The active edge
`u+s=4s` then forces `v=4s`.  The remaining locations have profile

```text
1s:2, 2s:2, 3s:2, 4s:2, 5s:2, 6s:1.                    (12)
```

The terminal `6s` edge is active by `(11)`, a contradiction.

### Case II: `u=2s` and `t!=2s`

A seed `t<2s` would be singleton, so `t>2s`.  The third root also satisfies
`v>2s`.  There is an active edge at `3s` from `2s+s`.  If `t` and `v` shared a
location strictly below `3s`, they could cancel each other, but the `3s` edge
would then remain singleton; if only one lay below `3s`, it would already be
singleton.  Therefore neither lies below `3s`, and at least one equals `3s`.

There are exactly three subcases:

```text
t=3s, v>3s:       the next edge forces v=4s;  7s is singleton;
v=3s, t>3s:       the next edge forces t=4s;  5s is singleton;
t=v=3s:                                      5s is singleton.            (13)
```

Every terminal in `(13)` is one of the positive factors `(11)`.

### Case III: `t=u=2s`

The two edges at `3s` already meet.  If `v<3s`, its root is singleton.  If
`v=3s`, the edge at `5s` is terminal.  For `v>3s`, the active edge at `4s`
forces `v=4s`; then the edge at `5s` is terminal.  The two profiles are

```text
v=3s:  1s:2, 2s:3, 3s:3, 4s:2, 5s:1;
v=4s:  1s:2, 2s:3, 3s:2, 4s:2, 5s:1, 6s:1.              (14)
```

This contradicts no-singleton cancellation in every branch of `(10)`, proving
`(3)`.

## 5. Leading-edge scout and controls

There is a second useful structural reduction.  Under graded lex order, the
sum of the leading nonlinear `P` support and leading nonlinear `Q` support is
the unique largest cross bucket and is larger than every divergence bucket.
Therefore any normalized two-by-three support survivor, including an axis
branch, must make that determinant zero: the two leading supports are
parallel.

The companion contains two bounded checks:

```text
direct, no structural filter, total degree <=7:
  33 monomials, 2,880,768 support-pair/triple pairs, 0 survivors;

leading-parallel filtered, total degree <=10:
  63 monomials, 77,555,583 full pairs,
  3,739,003 necessary-filter tests, 0 survivors.          (15)
```

The first census constructs all active debt buckets directly, with no use of
the shift-root proof.  The second uses only the unique-leading-edge necessary
condition to search farther.  Neither bounded computation is used to infer the
degree-unbounded theorem.

Support and axis drops are genuine.  The tame map

```text
P=x+alpha y^2,                  Q=y+delta P^2
```

has one nonlinear term in `P`, three in `Q`, and Jacobian `1`.  Likewise
`P=x+alpha y^2+gamma y^3, Q=y` has two pure-`y` nonlinear supports and
Jacobian `1`.  These do not inhabit `(1)`, but show why its boundaries should
not be silently discarded.

An exact collision is also a separate obligation.  The two-by-three map

```text
P=x+a x^2-x^3,
Q=y+b x^2+d x^4+e x^6
```

collides at `(1,0)` and `(-1,0)`, yet

```text
Jac(P,Q)=1+2ax-3x^2,
```

so the Keller debt survives.

## 6. Consequence and reproduction

Run

```bash
python3 -B 04-computation/jc2_ordinary_two_by_three_shift_root_peeling_thm3692.py
python3 -O -B 04-computation/jc2_ordinary_two_by_three_shift_root_peeling_thm3692.py
```

The exact companion verifies `(11)`, all six forced tower profiles in
`(12)--(14)`, both censuses in `(15)`, the tame and collision controls, and the
absence of inactive Python `assert` statements.

Combining this theorem with THM-3690 gives a sharp next sparse frontier: an
exact normalized two-sided candidate with only two nonlinear terms in `P`
needs at least three in `Q`, and at the two-by-three level it must enter an
axis branch

```text
(p_i)_1=0 for some i, or (q_j)_2=0 for some j.            (16)
```

Those axis cells remain open here.  No Keller counterexample or nonproperness
certificate is produced.  **QED.**
