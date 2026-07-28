---
id: THM-2635
title: "Half-tooth opposite-graph unit section and reversed-digit closure"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  Split every translated deep-probe tooth in THM-2616's complete common-x
  carrier into its two literal incident halves before normalization.  The
  full refined carrier still has global content 26.  On THM-2629's imposed
  graph r=-h-1 and THM-2600's canonical rail selector, either fixed half is
  positive on the same nine future labels in all 84 base cells.  The right
  half epsilon=0 has a uniform fixed-half unit at h=9; the left half
  epsilon=1 has uniform fixed-half units at h=3,8,10.  Thus one literal
  physical edge preserves a nonzero seven-clock class uniformly.  In a
  pointed THM-2620 determinant sector, r=-h-1 is exactly minus the outgoing
  endpoint coordinate and the future carrier is an eleven-edge path in the
  parabolic C13 cycle.  THM-2630's carry formula has only two reversed-digit
  closures on that path; only (epsilon,h,kappa,j,r)=(1,3,0,4,9) also belongs
  to the canonical uniform unit set.  This is an imposed positive/unit
  section, not a successor checksum.  It does not align the THM-2625
  allocation gauge or clock, transport its endpoint current, identify two
  semantic roots, exclude a row, or prove LRC(14).
source: endpoint-pair-audit-2026-07-28-half-tooth-unit-section
depends_on:
  - THM-2600-constant-six-middle-rail-common-x-atlas-and-uniform-bockstein-section
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2620-endpoint-pair-parabolic-transvection-and-translation-gauge-boundary
  - THM-2629-fixed-deep-affine-graph-spectrum-and-puncture-cancellation-boundary
  - THM-2630-old-wall-affine-clutching-and-successor-sector-no-go
related:
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2634-endpoint-pair-two-carry-cospan-and-single-carry-no-go
script: 04-computation/lrc14_half_tooth_unit_section_thm2635.py
output: 05-knowledge/results/lrc14_half_tooth_unit_section_thm2635.out
script_sha256: 9835f86aae0ce4401cc05d0e7a7ea8a90a0d00a7497be34bde0ef29ab1160f6b
output_sha256: 5f7bea76659c6462c2def139027bfa28e8ed6a79f6718803fd90f2e2f6a73ca4
hash_basis: LF-normalized bytes
---

# THM-2635 -- one literal later edge retains the opposite-graph unit

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2629 finds that the coefficient graph

```text
r=-h-1 mod 13                                               (1)
```

is uniquely optimal for the unsplit seven-clock unit score.  THM-2630 then
proves that (1) is not a physical successor checksum: even after retaining
the future-owner clock and the incident half-edge, every affine graph leaves
large future-digit fibres.

There is nevertheless a sharper positive fact between those two results.
If (1) is imposed as a section, one *fixed literal half* of the later tooth
already carries a nonzero globally primitive Bockstein in every base cell.
This removes the possible obstruction that all of THM-2629's noncancellation
comes from adding the two physical incident edges.  The resulting edge has an
exact parabolic dictionary and one reversed adjacent-digit closure.  Neither
fact supplies the still-missing endpoint-current/semantic-root gluing.

All residues below lie in `F_13`.

## 1. Split the complete carrier before normalization

Use THM-2616's full nonnegative common-`x` tensor

```text
A_(e,q;ell,r,h),

e: one of 162 middle rails,
q,h in F_13,
ell in F_7,
r in F_13^*.                                               (2)
```

Here `q` is the present target section, `h` is the delayed physical future
digit, and `r` labels the actual translated deep probe

```text
Delta_r(x)=d_1(2*13^5 x-r/13).                             (3)
```

In the exact integer interval coordinates of THM-2616/2630, split the tooth
in (3) into

```text
epsilon=0: [14r,    14r+13]       (right half),
epsilon=1: [14r-13, 14r   ]       (left half).             (4)
```

Let `A^epsilon` be the numerator after applying (4), before any reduction
modulo 13.  Endpoint conventions on null interval boundaries do not change
the measure, and the exact half-open implementation gives

```text
A^0+A^1=A                                                    (5)
```

on every labelled digit entry.

The global primitive content is recomputed over the **entire** refined
carrier, not on the graph (1) and not separately on a rail or half:

```text
g_half=gcd{A^epsilon_(e,q;ell,r,h): all labels}=26.          (6)
```

Thus it is exactly the same as THM-2616's full-carrier content, but this is a
new theorem rather than an assumed divisibility.  The companion checks

```text
1,201,816 positive half entries,
2,299,752 digitwise identities in (5).                      (7)
```

Consequently

```text
a^epsilon=A^epsilon/26                                     (8)
```

is an honest globally primitive integer refinement.  In particular, the
mod-13 unit test below does not illegally divide by `26=0 mod 13` after a
half-carrier changed the content.

## 2. Canonical rail and fixed-half unit classes

For a base cell `c=(s,ell_4)`, use THM-2600's canonical rail selector

```text
theta_c=1 iff s in {6,11} or (s,ell_4)=(8,2),
theta_c=0 otherwise.                                      (9)
```

This chooses one of the already proved positive rails in every one of the
`12*7=84` base cells.  For `h in {1,...,11}`, put

```text
r_h=-h-1 !=0                                               (10)
```

and define the fixed-half seven-clock vector

```text
Y^epsilon_ell(c,h)
 =a^epsilon_(e_c,h;ell,r_h,h) r_h^(-1) mod 13.             (11)
```

As in THM-2616/2629, subtract the last coefficient and regard (11) as an
element of

```text
F_13[z]/(Phi_7(z)).                                        (12)
```

Call `(c,h,epsilon)` a fixed-half unit when this class is invertible.  No
row, graph, or half is reprimitivized.

Let `P_c^epsilon` be the `h` for which at least one coefficient in (11) is
positive before reduction, and let `U_c^epsilon` be its unit set.  Exact
computation gives the common positivity law

```text
intersection_c P_c^0 = intersection_c P_c^1
 =H_9={1,3,4,5,6,7,8,9,10}.                              (13)
```

More strongly,

```text
intersection_c U_c^0={9},
intersection_c U_c^1={3,8,10}.                            (14)
```

The complete cell-size histograms are

| half | positivity sizes | unit sizes |
|---|---|---|
| `epsilon=0` | `10^21,11^63` | `5^1,9^18,10^15,11^50` |
| `epsilon=1` | `9^2,10^82` | `5^1,8^2,9^21,10^60` |

Therefore the single literal left incident edge `epsilon=1`, the single
canonical rail selector (9), and each one of `h=3,8,10` give a positive
common-`x` carrier and a nonzero seven-clock unit in all `84` base cells.
The realizing rail does not vary with `h` or with the choice of unit witness.

If either of the one or two already proved rails over a cell may be used,
the common unit sets enlarge to

```text
epsilon=0: {3,5,6,7,9},
epsilon=1: {1,3,4,6,8,10}.                               (15)
```

Equation (15) is recorded as a control; the canonical theorem uses (9).

## 3. Recombination is a hostile control, not a monotonicity law

Recombine the two vectors in (11) entrywise before testing unitness.  The
canonical-rail common set and size histogram are

```text
{3,4,5,8,9,10},                   7^1,10^22,11^61.        (16)
```

This recovers the unsplit fixed-`r` coefficient class on the same rail.
Notice that (16) is neither the union nor the intersection of the two sets
in (14).  Cyclotomic unitness is nonlinear under addition: a half may create
or lose a unit by changing cancellation.  The meaningful new statement is
the direct globally normalized half computation, not inheritance from the
unsplit graph.

The complete half-graph coefficient bank has digest

```text
6827925d5fd66cc96d92db80e61b7b3f0ed275cc9a45708015d697f71dba9eb9. (17)
```

## 4. Exact endpoint dictionary: an eleven-edge parabolic path

Fix one nondegenerate THM-2620 sector

```text
q_vec!=0,                     Delta!=0.                    (18)
```

Choose `R_0` with `det(q_vec,R_0)=Delta` and put

```text
R_h=R_0+h q_vec,
L_h=R_0+(h+1)q_vec.                                    (19)
```

Then

```text
L_h-R_h=q_vec,              det(L_h,R_h)=Delta,           (20)
```

and the thirteen pairs `(L_h,R_h)` are exactly the transvection fibre.  The
projective parabolic sends `[R_h]` to `[L_h]`.

Because `(q_vec,R_0)` is a basis, there is a unique covector `lambda` with

```text
lambda(q_vec)=1,                 lambda(R_0)=0.            (21)
```

It reads

```text
lambda(R_h)=h,                   lambda(L_h)=h+1.          (22)
```

Hence THM-2629's opposite graph has the exact pointed-sector form

```text
r_h=-h-1=-lambda(L_h).                                   (23)
```

The full affine chart is the oriented cycle

```text
0->1->...->12->0.                                        (24)
```

The available future carrier `h=1,...,11` retains precisely

```text
1->2, 2->3, ..., 11->12,                                 (25)
```

so it is an eleven-edge directed path.  Its two missing adjacent cycle edges
are

```text
12->0  (would have r=0),          0->1 (future h=0 absent). (26)
```

This explains both the functional form `r=-h-1` and THM-2629's puncture
alignment.  It is an exact abstract endpoint dictionary, not yet a physical
endpoint allocation.

The choice of `R_0`, and therefore the zero of `lambda`, is gauge data.
Replacing `R_0` by `R_0+c q_vec` preserves `(q_vec,Delta)` but translates
the affine coordinate `h`.  Nothing in THM-2625 identifies its left-deep
allocated endpoint origin with the physical digit-zero convention in (21).

## 5. One reversed adjacent-digit closure

THM-2630 supplies the literal later-scale carry law.  If

```text
j=floor(13^6 x) mod 13,
u={13^6 x},
kappa=floor(2u),
h=floor(13u),                                              (27)
```

then on half `epsilon`

```text
r=2j+kappa+epsilon,
j=7(r-kappa-epsilon).                                     (28)
```

Thus `j` and `h` are adjacent base-thirteen digits on the same physical
ancestry, while the abstract outgoing coordinate in (23) is

```text
-r=h+1.                                                   (29)
```

On `h=1,...,11`, including both possible `kappa` values in the split cell
`h=6`, equations (28)--(29) have exactly two solutions with

```text
j=-r=h+1:

(epsilon,h,kappa,j,r)=(1,3,0,4,9),
(epsilon,h,kappa,j,r)=(1,7,1,8,5).                        (30)
```

The mechanism is a fixed-point equation, not an unexplained finite hit.
The pointed-sector condition is `r=-j`, while the physical tooth law is
`r=2j+kappa+epsilon`; hence

```text
3j=-(kappa+epsilon).                                      (30a)
```

The three possible sums `kappa+epsilon=0,1,2` give respectively

```text
j=0,4,8.                                                  (30b)
```

The first is the punctured `r=0`, `h=12` edge.  The middle value has two
formal bit decompositions, but `j=4` forces `h=j-1=3`, whose physical half
has `kappa=0`; therefore only `(epsilon,kappa)=(1,0)` is consistent.  The
last value gives the second solution `(epsilon,kappa,h)=(1,1,7)`.  This
explains both closures and why only the left half occurs.

Only the first lies in the canonical uniform fixed-half unit set (14):

```text
(epsilon,h,kappa,j,r)=(1,3,0,4,9).                        (31)
```

At (31), the positive/unit parabolic label edge is

```text
h -> -r : 3 -> 4,                                        (32)
```

while physical digit time is

```text
j -> h : 4 -> 3.                                         (33)
```

It is therefore a **reversed** adjacent-digit closure.  The second closure
`7->8` is positive but not a uniform canonical unit.  Conversely the other
uniform units `h=8,10` do not close to the adjacent predecessor digit.  This
three-way hostile separates positivity, unit noncancellation, and physical
chronology.

## 6. Exact stopping boundary

The strongest conclusion is now:

```text
one canonical physical half-edge
 + one imposed opposite-affine section
 + one seven-clock unit in every base cell
 + one reversed adjacent-digit closure.                  (34)
```

It is not a positive successor theorem.  THM-2630 proves that retaining the
edge, future-owner clock, and any affine graph does not make `h` a function
of the physical probe data.  Equation (1) is externally selected from a
dense positive relation; it is not inferred from `x` as a checksum.

Nor does (34) transport THM-2625's endpoint current.  The two constructions
have different typed events:

```text
THM-2625: owner c1, word {a}, delay k=2 (R=13^2),
          complex left-deep allocated current J(L,R);

here:     THM-2616's present/delayed {a}->{b} common-x tensor,
          later probe at speed 2*13^5, delayed digit at R=13^6. (35)
```

THM-2625 proves nonzero current on every abstract edge in every
nondegenerate sector, including both orientations after changing sector
signs.  It does **not** identify its allocation covector with `lambda`, put
its current on the half-carrier (4), or show that the two separately nonzero
cyclotomic classes have nonzero product on a common atom.

Finally, `h` is a proved delayed physical root digit, but `j` in (27) is at
present a predecessor/carry digit, not a proved semantic root of an adjacent
word event.  The class (11) is a seven-clock aggregate, not one clock sheet.
The first missing map is therefore a common-ancestry cospan which

1. realizes `j` as a semantic predecessor root on an adjacent clock;
2. fixes the pointed-sector origin and reversal gauge in (21)--(33);
3. transports one THM-2625 allocated endpoint coefficient to the same
   half-tooth atom while preserving noncancellation; and
4. keeps both endpoint carries rather than attaching one carry only after
   Fourier aggregation.

The last item is the unproved target reserved by THM-2634; that stub is not a
dependency.  No scalar row is excluded and no LRC(14) conclusion follows.

## 7. Exact companion

Run

```bash
PYTHONPATH=04-computation \
  python3 04-computation/lrc14_half_tooth_unit_section_thm2635.py

PYTHONPATH=04-computation \
  python3 -O 04-computation/lrc14_half_tooth_unit_section_thm2635.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_half_tooth_unit_section_thm2635.out.
```

The companion independently rebuilds the full THM-2616 interval carrier,
splits every translated tooth, recomputes the global content over all labels,
checks every digitwise additive partition, tests every canonical and
either-rail fixed-half unit, recombines the halves, checks the complete bank
digest, and exhausts the adjacent-digit closure equations.  Every logical
decision is an exact integer or finite-field optimized-mode guard.

QED (candidate; independent hostile audit pending).
