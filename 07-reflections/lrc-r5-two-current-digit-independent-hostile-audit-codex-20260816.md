# Independent hostile audit of two current inverse digits

**Status: FINITE-EXACT INDEPENDENT HOSTILE AUDIT ACCEPTS THE SCOPED
TWO-STATIC-CURRENT-DIGIT PROJECTIVE CYLINDER AND REPAIRS THE RECURRENCE
CLAIM; LRC(14) remains OPEN.**

The audit did not import the submitted two-digit implementation.  Before
reading it, it fixed the current-sheet radix convention

```text
a = r_0 + 13 r_1 + 169c,             0 <= r_0,r_1 < 13,
```

and independently rebuilt the THM-2471 source packet from the accepted
one-digit audit's canonical dependencies.  It then folded by `13^3`, pulled
the `r_1` window first, the `r_0` window second, and the collision-root
window last.  The resulting finite object agrees with every submitted
source, response, rank, support, and spectral digest.

## Radix, endpoints, and reflection

A generic half-open interval routine, written independently of the
candidate, proves literal profile-list equality for all `169` digit pairs:

```text
window_(r_0)(window_(r_1)(H))
  = window_(r_0+13r_1; modulus 169)(H).
```

This is equality with the same `[left,right)` convention, including every
endpoint; it is not merely an almost-everywhere mass comparison.  The source
reconstruction has `120234` primitive pieces, `547` folded pieces, `715`
two-digit pieces, `2743` root pieces, and `34` joint boundaries.  Every
active `(state,u)` cell has exactly `132` of the `169` current addresses.

The independently derived involution is

```text
(y,u,r_0,r_1,state)
  -> (1-y,12-u,12-r_0,12-r_1,state XOR 2).
```

It has no borrow because every digit of `13^5-1` is twelve.  The audit checks
the involution on every joint interval and on the full address-mass table;
the latter has digest
`4869b7db8c3a713cdb7dabc91949d0c37d4220560b640db3d88e853e3ce2f49a`.

## Exact marginal chain and endpoint audit

The first source marginal holds pointwise on every source interval:

```text
sum_(r_1) U_(u,r_0,r_1) = U_(u,r_0).
```

After multiplying the literal endpoint factors, integrating, and inverting
characters, the full chain is exact:

```text
T_2(state,r_0,r_1,t)
  --sum r_1--> independently audited T_1(state,r_0,t)
  --sum r_0--> audited Boolean square T_0(state,t).
```

The first arrow reproduces the accepted one-digit tensor digest
`2a195fac7fbd60a4bbd2597bf34baf87302ead16c7c0e8c67c0667b0d320dfba`;
the second reproduces the Boolean-square digest
`5f4b9609faaa5f148d112a7cde5cfba0ab2c1385b4c53ea9c4bcfc6e93d106fc`.
All thirteen per-character gamma digests also match.  Three literal guard
controls pass, the same-root sector is pointwise zero, and four independently
coded direct `2197`-term inversions equal the compressed result.  Thus the
sparse endpoint aggregation is an exact Fubini regrouping, not a change of
normalization, sheet, sign, or phase.

## Rank obstruction and its exact boundary

In axis order `(state,r_0,r_1,relation,combined-address)`, the actual and
flat-hostile ranks are

```text
actual:        (4,13,13,6,4)
flat:          (4, 4, 1,6,4)
pure actual:   (3,12,12,6,4)
pure flat:     (0, 0, 0,0,0).
```

For each fixed `r_0`, the actual raw and contrast ranks of the
`r_1 x (state,relation)` flattening are

```text
(4,3,3,4,3,3,4,3,3,4,3,3,4),
```

while the same-parent flat hostile has raw ranks one and contrast ranks
zero.  All `5184=3*12^3` pure four-coordinate Fourier modes fire.  The
combined-address flattening nevertheless has rank only four, so the two
digit axes do not create a `169`-dimensional response sector.

These conditional ranks are amplitude ranks: rows are `r_1` and columns are
the combined `(state,relation)` output.  By contrast, commit `2d52215a3`
separates a rank-four channel-to-output amplitude map from a rank-six
relation-row carrier.  The two-current tensor has already marginalized the
pointed root/tail coordinate, so it cannot test equality with that carrier or
define an operator on it.  The repeated value four is compatible with an
amplitude bottleneck, but without a coordinate map it proves neither a shared
amplitude quotient nor a recurrence.

The rank test has one precise consequence.  It rules out every scalar,
state/relation-independent factorization

```text
T_2(state,r_0,r_1,t)=K(r_0,r_1)T_1(state,r_0,t),
```

including arbitrary signed and noncirculant `K`, because every fixed-`r_0`
flattening of such a product has rank at most one.  It does **not** rule out
a hidden-state recurrence, a matrix-valued carrier transition, a
state/relation-dependent or nonlinear rule, or any temporal recurrence.
The former headline “not an autonomous recurrence” was therefore too broad
and is repaired to a scalar-kernel obstruction.

## Connection contract and boundary

| field | audited answer |
|---|---|
| source | THM-2471 current sheet above the audited owner Boolean square |
| target | `V_4 x F_13(r_0) x F_13(r_1) x F_13(relation)` |
| map | radix split of the same current inverse index, exact nested half-open windows, literal endpoint multiplication, integration, character inversion |
| preserved | common base, current sheet, both current digits, collision root through endpoint selection, state, guards, phase, source weights, reflection |
| destroyed | three older current digits, source time, root difference, pointed tail, deep label, grouped exact address, chronology |
| positive gate | exact two-step marginals plus fixed-`r_0` conditional ranks `3/4` against scalar bound `1` |
| hostile | uniform `r_1` lift with the same parent, conditional rank `1`, contrast rank `0`, and no four-way interaction |
| cheapest successor | retain the audited six-pointed/root-difference carrier and test a matrix transition; a third digit alone is still static |
| prohibited inference | complete address, `U_clock` action, arrival ancestry, physical current, row exclusion, or LRC(14) |

The accepted phrase is therefore **two-static-current-digit projective
cylinder with a scalar state/relation-blind factor obstruction**.  It is not
a current chronology theorem and does not close or reduce LRC(14).

The concurrently integrated source/current-digit probe therefore sharpens
the proposed successor: a matrix carrier is the right level to test.  It
crosses `b_source` with `r_0`, however, not `r_1` with `r_0`; it neither
supplies nor refutes the missing two-current-digit carrier transition.

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_independent_audit_20260816.py
```

Normal and optimized runs have identical complete transcripts.  Independent
script/output/semantic LF-normalized SHA-256:

```text
126be106a34a990f22e66b26d79ea3568ee4f394419936ed47f1bfbe0656788f
e12b1eecf6c7b4cb2a25398a6fdfcf7d3261693e810a3104b9ca1811d54ca15e
f2af726e9d5abd1487e841623ce2f62ca647c86e5a1a68e41eda4d9dda6c81ac
```
