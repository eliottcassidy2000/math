# Two current inverse digits form a projective tree with a scalar-kernel obstruction

**Status: FINITE-EXACT POSITIVE TWO-DIGIT CURRENT-ANCESTRY OBJECT ON ONE
`r=5` OWNER BASE; INDEPENDENTLY HOSTILE-AUDITED WITH SCOPE REPAIR;
LRC(14) remains OPEN.**
Retaining the two lowest base-thirteen digits of the THM-2471 current inverse
branch gives an exact `V_4 x F_13^2 x F_13` response tensor.  Its
`r_1` marginal is the pinned one-digit `r_owner` candidate row by row,
and the next marginal is the audited Boolean square.  Nevertheless, the
second digit is not generated from the first by any scalar
state/relation-independent one-step factor: all thirteen conditional ranks
are three or four, whereas
every such scalar lift has rank at most one.  This does not obstruct
hidden-state, matrix-valued, state-dependent, nonlinear, or temporal
recurrences.

This is a static two-cylinder ancestry statistic.  It is not a complete
THM-2334 address, a `U_clock` transition, an arrival/source intertwiner, a
physical current, a row exclusion, or a proof of LRC(14).

## Inheritance and the exact question

The closest current-leg object is the candidate at
`18d5f423a`.  It retains

```text
r_0 = r_owner = a mod 13
```

on the THM-2471 sheet

```text
X_(u,a)=(w_u+a)/13^5,       w_u=(y+u)/13.
```

Its branch marginal recovers the audited owner-visible Boolean square, its
actual three-way ranks are `(4,4,6)`, and its pure interaction ranks are
`(3,4,6)`.  The canonical hostile is still MISTAKE-417: full character
support can be a separable lift, so exact marginal recovery and a rank gate
must precede any spectral interpretation.

The least-used lawful sidecar is the immediately older digit on this same
current leg.  Write

```text
a = r_0 + 13 r_1 + 13^2 c,
0 <= r_0,r_1 < 13,
0 <= c < 13^3.                                      (1)
```

Unlike the source-time high digit of the left `P_169 e` packet, both
`r_0` and `r_1` belong to the same `a` in `X_(u,a)`.  Equation
(1) supplies an explicit support map; no equality of unrelated
`F_13` labels is being assumed.

## Lawful fold/window factorization

Let

```text
H(z) = sum_(c=0)^(13^3-1) f((z+c)/13^3).
```

This is the numerator profile of the `13^3` Perron fold.  Window
composition is contravariant, so the higher retained digit is taken first:

```text
H_(r_1)(x)       = H((x+r_1)/13),
H_(r_0,r_1)(w)  = H_(r_1)((w+r_0)/13)
                 = H((w+r_0+13r_1)/13^2),
U_(u,r_0,r_1)(y)= H_(r_0,r_1)((y+u)/13).             (2)
```

Substituting (1) into (2) gives exactly

```text
U_(u,r_0,r_1)(y)
 = sum_c f((w_u+r_0+13r_1+13^2c)/13^5).             (3)
```

All interval lists use the convention `[left,right)`.  The extraction
window for a digit `r` is
`[rT/13,(r+1)T/13)`; the right endpoint belongs to the next window.
The script checks profile-list equality, not merely almost-everywhere
agreement:

```text
window_(r_0)(window_(r_1)(H))
  = window_(r_0+13r_1, modulus 169)(H)               (4)
```

for all `169` pairs.  Thus no endpoint can be double-counted or omitted.
At the outer event sweep, masks are toggled at the left endpoint and every
source boundary is inserted before integration, preserving the same
half-open convention.

## The exact marginal chain

Put `c'=r_1+13c`.  Summing (3) over `r_1` gives

```text
sum_(r_1) U_(u,r_0,r_1) = U_(u,r_0),                 (5)
```

the existing one-digit profile.  Summing `r_0` next gives

```text
sum_(r_0,r_1) U_(u,r_0,r_1) = U_u.                  (6)
```

The program checks (5) pointwise on every joint source interval and then
checks both response marginals after endpoint multiplication, integration,
and character inversion:

```text
T_2(state,r_0,r_1,t)
  --sum r_1--> T_1(state,r_0,t)
  --sum r_0--> T_0(state,t).                         (7)
```

The first arrow reproduces the pinned one-digit gamma digest
`344fc6f4...d9e72` and tensor digest `2a195fac...fba`.  The second
reproduces the audited Boolean-square tensor digest
`5f4b9609...06fc`.  Three literal guard controls pass before and after
the first marginal, and the same-root sector is pointwise zero at every
two-digit fibre.

Equation (7) is a genuine projective cylinder law.  It is not yet a temporal
recurrence: no second clock or action carrying one cylinder level to another
has been identified.

## Geometry and reflection

The apparent `13^2` enlargement remains geometrically small.  The exact
source decomposition has

```text
13^3-fold profile pieces:          547
first-window pieces by r_1:
  (7,55,53,54,54,55,3,55,54,54,53,55,7)
all mod-169 window pieces:         715
all root-window pieces:          2743
joint source boundaries:           34.
```

On every owner-visible source cell, each active root has exactly `132`
of the `169` addresses.  All thirteen `r_0` values survive after the
`r_1` marginal, so the older digit is not a cosmetic duplication of the
one-digit branch.

The exact sheet reflection is

```text
1-X_(u,a)(y) = X_(12-u,13^5-1-a)(1-y).
```

Because `13^5-1` has all base-thirteen digits equal to twelve, no borrow
occurs.  Hence

```text
(y,u,r_0,r_1,state)
  -> (1-y,12-u,12-r_0,12-r_1,state XOR 2).           (8)
```

The program checks (8) on every joint interval and checks the full
`169`-entry address-mass reflection.  This is a cylinder involution, not a
tournament on the address labels.

## Sparse exact implementation

A direct loop over every weighted endpoint segment, all thirteen guard
shifts, and all `169` addresses would request roughly

```text
186244 * 13 * 169 = 409178068
```

address updates.  The exact integrand is constant on only 33 intervals
between the 34 joint source boundaries.  The implementation therefore:

1. sweeps the literal endpoint events without changing their masks or guard
   semantics;
2. accumulates the scalar endpoint contribution by
   `(source cell, selected U-root mask)`;
3. expands the `169` address weights only after that accumulation; and
4. sums the beta and alpha character axes before the final thirteen-point
   relation transform.

The complete run performs `1,248,644` scalar coefficient updates and only
`11,086` post-aggregation keys across all guard shifts.  Four nonzero
entries are also recomputed by a separately coded direct `2,197`-term
character inversion, agreeing with the compressed route.  This compression
uses exact constancy of the source profiles; it does not coarsen a coordinate
or replace the endpoint factors.

## Conditional-rank hostile: projective is not a scalar Markov lift

Let `T_1` be the one-digit parent.  Any scalar state/relation-independent
one-step digit factor, even a completely arbitrary non-circulant one, has
the form

```text
T_K(state,r_0,r_1,t)
 = K(r_0,r_1) T_1(state,r_0,t).                      (9)
```

For fixed `r_0`, flatten (9) with rows indexed by `r_1` and columns by
`(state,t)`.  It is an outer product, so

```text
rank T_K[r_0] <= 1.                                  (10)
```

The concrete flat hostile `K=1/13` has the same `r_1` marginal as the
actual tensor.  Its thirteen raw conditional ranks are all one and its
digit-contrast ranks are all zero.

The exact tensor instead gives

```text
raw conditional ranks =
contrast ranks =
(4,3,3,4,3,3,4,3,3,4,3,3,4).                       (11)
```

Thus the second digit carries state/relation-dependent response information
at every first digit.  Equations (7) and (11) together are the main result:
the two levels are projectively consistent but cannot be explained by a
scalar state/relation-independent first-order factor.  No conclusion about
general autonomous recurrences follows from this rank test.

The period-three appearance in (11) is recorded but not promoted to a
recurrence.  A recurrence claim would require the same typed operator at a
third cylinder level or a lawful clock transport.

## Rank and spectral record

In the order
`(state,r_0,r_1,relation,combined-address)`, the flattening ranks are

```text
actual:             (4,13,13,6,4)
flat hostile:       (4, 4, 1,6,4)
pure four-way:      (3,12,12,6,4)
flat pure four-way: (0, 0, 0,0,0).                  (12)
```

The combined address flattening has rank only four.  This is an important
limitation: two full-rank thirteen-valued coordinate axes do not create a
169-dimensional response sector.

The actual four-dimensional Fourier tensor has all
`4*13^3=8788` coefficients nonzero.  After removing every lower-order
marginal, exactly the

```text
3*12*12*12 = 5184
```

pure four-coordinate modes remain, all nonzero.  The flat hostile has no
`r_1`-bearing mode and zero pure interaction.  Fourier fullness is
secondary here; the marginal chain and conditional-rank obstruction are the
decisive gates.

At relation `(1,0,6)`, the actual state rows have support
`(132,133,132,133)`, state rank four, and the same conditional pattern
(11).  The flat hostile has `169` entries in every state row but conditional
rank one.  Dense support is therefore not being mistaken for ancestry
information.

## Typed boundary and next tests

The coordinate hierarchy now contains several distinct thirteen-valued
objects:

```text
source-time b=floor(n/13) in P_169 e,
current r_0=a mod 13,
current r_1=floor(a/13) mod 13,
collision root u,
source-root difference u-q,
deep label r_deep,
THM-2334 grouped exact address.
```

Only `r_0,r_1` have been joined here by an explicit radix map on the same
current sheet.  Nothing in this computation identifies them with any other
line.

The concurrently arrived current-branch/root-difference crossing adds a
useful constraint on the recurrence question.  With only `r_0` retained,
its statewise branch row spaces equal the six independently audited pointed
source-tail row spaces.  Thus a true finite-state ancestry recurrence should
be tested on that six-carrier bundle.  The present probe marginalizes the
root pair and pointed tail, so its conditional ranks detect failure of a
scalar autonomous digit kernel but do not yet prove or refute a six-by-six
pointed-carrier transition.

The cheapest next decisive tests are:

- an independent reconstruction of (2)--(7) that does not import the
  one-digit candidate except for final digest comparison;
- the joint pointed/root-difference tensor retaining both `r_0` and
  `r_1`, testing whether the second digit acts on the audited six-carrier
  row bundle;
- a third current digit, testing whether the conditional-rank pattern is
  stable under another projective marginal; and
- a lawful clock action carrying one typed two-cylinder stalk to another.

The third digit alone would still be static ancestry.  Chronology requires
the last item.

## Connection contract

| field | exact answer |
|---|---|
| source | THM-2471 current sheet `X_(u,a)` above the audited owner-node Boolean square |
| target | `V_4 x F_13(r_0) x F_13(r_1) x F_13(relation)` |
| map | write `a=r_0+13r_1+169c`; fold by `13^3`; take `r_1`, `r_0`, then root windows; multiply literal endpoint factors; integrate and invert |
| preserved | common base, current sheet, collision root until endpoint selection, Boolean state, literal guards, source weights, endpoint factors, both current digits, reflection |
| destroyed | three older current digits, source-time branch, root difference, absolute source tail, deep label, grouped exact address, chronology |
| mandatory chain | `sum r_1 ->` pinned one-digit current tensor; `sum r_0 ->` audited Boolean square |
| positive gate | every fixed-`r_0` conditional rank is `3` or `4`, versus the scalar state/relation-blind factor bound `1` |
| concrete hostile | uniform `r_1` lift with identical one-digit parent, rank-one conditionals, and zero four-way interaction |
| boundary | finite candidate on one owner base; no current, row theorem, exclusion, or LRC(14) conclusion |

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_probe_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_two_digit_current_ancestry_probe_20260816.py
```

The pinned semantic digest is
`61743457afc0cff984c87affa7f2e67bf3a21e08a401ea69b679319f2f51e826`.
