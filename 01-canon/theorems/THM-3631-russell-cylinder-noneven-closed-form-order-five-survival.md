---
id: THM-3631
title: "Russell-cylinder non-even closed-form survival through order five"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For the fixed non-even hostile polynomial
  Q_h of THM-3627, a closed polynomial target two-form pulls back to the common
  normalized unit through total source degree five on all three collision
  branches.  The exact combined system has rank 92 in 168 unknowns and a
  76-dimensional affine solution space.  THM-3627 obstructs even arbitrary
  target two-forms at degree six, so six is the first failure degree for this
  fixed Q_h both before and after imposing closedness.  Local formal Darboux
  gives target-pair germs, not global polynomial functions or a JC(2)
  counterexample.
source: kps-s189 / THM-3627 de Rham leak continuation, 2026-08-21
audit: >
  PASS -- an independent hostile reconstruction using the THM-3627 Fraction
  sparse-series backend matched every pullback and de Rham matrix, recovered
  the complete rank table and a 49-term rational witness with
  (P(0),K(0),R(0))=(0,4,0), checked the exterior-derivative sign and the
  common nonzero constant, and verified byte-identical normal, optimized, and
  stored transcripts after the declared LF normalization.  The audit also
  repaired two wording/provenance imprecisions without changing the claim.
depends_on:
  - THM-3624-russell-cylinder-noneven-fold-weighted-cokernel-boundary
  - THM-3627-russell-cylinder-noneven-hostile-degree-six-closure
related:
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
script: 04-computation/jc2_russell_cylinder_noneven_closed_form_order_five_thm3631.py
output: 05-knowledge/results/jc2_russell_cylinder_noneven_closed_form_order_five_thm3631.out
script_sha256: 9eb341a3c6c451865ecc1891faee4db7aa96e210133375aec0cd7277d56440be
output_sha256: a37e8b1a7d057096f1f53a55e1dfd7ce132817b574bef9687b8619156f0680bc
semantic_sha256: f8ce8a774ffca01430a00519e6672d30c7458e51829037e1daf731163c138796
hash_basis: LF-normalized bytes for files; canonical JSON for semantic ledger
---

# THM-3631 -- non-even closed-form survival through order five

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
THM-3627 found the first arbitrary-form
obstruction for one explicit non-even fold.  This refinement shows that the
positive side immediately below that obstruction survives the necessary de
Rham condition and therefore represents genuine local target-pair germs.

## 1. Fixed hostile and target two-forms

Use the THM-3627 polynomial

```text
Q_h(x)=1/5408(
  44069x^11+112059x^10-154749x^9-406377x^8
 +188107x^7+513081x^6-82835x^5-230931x^4
 +5408x-4056).                                           (1)
```

At the three source points `x=-1,0,1`, the fold

```text
(x,t) |-> (x,Q_h(x)+t^2,t)                              (2)
```

has the common smooth Russell-cylinder target point used in THM-3624.  Write
its regular local coordinates as

```text
(c,epsilon,w)=(c,e+3,w).                                (3)
```

For `N>=0`, let `V_N` be the space of triples of target polynomials of total
degree at most `N`, and write the corresponding two-form as

```text
Omega=P dc wedge d epsilon
     +K dc wedge dw
     +R d epsilon wedge dw.                             (4)
```

Thus

```text
dim V_N=3 binom(N+3,3).                                 (5)
```

Pull `(4)` back along the three branch germs of `(2)`, retain the coefficient
of `d xi wedge dt`, and truncate at total source degree `N`.  Requiring the
same normalized constant `12` on all branches defines the THM-3627 linear
system

```text
P_N Omega=tau_N.                                        (6)
```

## 2. Closedness is one exact leak equation

In the regular coordinates `(3)`, direct differentiation gives

```text
d Omega=(P_w-K_epsilon+R_c)
          dc wedge d epsilon wedge dw.                  (7)
```

For coefficients of degree at most `N`, equation `d Omega=0` is therefore the
complete finite system

```text
P_w-K_epsilon+R_c=0                                    (8)
```

in every monomial of degree at most `N-1`.  It has `binom(N+2,3)` rows.  There
are no omitted higher terms: differentiating a degree-`N` polynomial produces
degree at most `N-1`.

Stack `(8)` under `(6)`.  Exact rational row reduction gives

| `N` | pullback shape | closedness shape | `rank P_N` | closed rank | augmented rank | solution dimension |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | `3 x 3` | `0 x 3` | 2 | 2 | 2 | 1 |
| 1 | `9 x 12` | `1 x 12` | 7 | 8 | 8 | 4 |
| 2 | `18 x 30` | `4 x 30` | 15 | 19 | 19 | 11 |
| 3 | `30 x 60` | `10 x 60` | 26 | 36 | 36 | 24 |
| 4 | `45 x 105` | `20 x 105` | 40 | 60 | 60 | 45 |
| 5 | `63 x 168` | `35 x 168` | 57 | 92 | 92 | 76 |

Every closedness row is independent, and at each degree its rowspace is
transverse to the pullback rowspace:

```text
rank[ P_N ; d ]=rank P_N+binom(N+2,3).                  (9)
```

Equality of the last two rank columns proves consistency over `Q`.  In
particular there is a closed polynomial two-form `Omega_5` of target degree at
most five whose three pullbacks equal

```text
12 d xi wedge dt + terms of source degree at least 6.   (10)
```

The same `Omega_5`, after truncating its pullback to lower source order,
supplies a solution for every earlier `N` in the table.  The `76` in the final
column is an exact affine dimension, not a numerical nullity estimate.

## 3. These are genuine local target-pair jets

At source degree zero, the three right-hand-side entries in `(6)` are

```text
(12,12,12).                                             (11)
```

Consequently every solution realizing `(10)` has `Omega_5(0)!=0` at the
common target point.  A standard local presymplectic Darboux argument now
applies.  For completeness, if `X` spans the one-dimensional kernel of a
nonvanishing closed two-form, then

```text
L_X Omega=d(i_X Omega)+i_X(d Omega)=0.                  (12)
```

Choose a transverse surface, put its nonvanishing area form into local
coordinates `dF wedge dG`, and extend `F,G` constantly along the flow of `X`.
Equation `(12)` gives

```text
Omega_5=dF wedge dG                                     (13)
```

as analytic, hence formal, germs at the common target.  Thus `(10)` is
not merely an arbitrary-form artifact: one local target-pair germ works on all
three source branches through degree five.

This argument does not make `F,G` global or polynomial.  It is a local formal
survival theorem only.

## 4. Degree six is the exact fixed-hostile boundary

THM-3627 proves for the same `Q_h` that the larger arbitrary-two-form system at
degree six has

```text
rank P_6=77,                     rank[P_6|tau_6]=78.     (14)
```

It supplies a sparse exact left-cokernel certificate annihilating all `252`
target columns but not `tau_6`.  Therefore no target two-form at all can
satisfy `(6)` at `N=6`; a fortiori no closed or decomposable one can.

Combining the positive table with `(14)` proves that, for this fixed hostile,

```text
closed/formal target-pair jets survive through N=5,
every arbitrary target two-form fails at N=6.           (15)
```

The word *first* in this boundary refers to total source degree for the fixed
polynomial `(1)`.  It is not a uniform degree bound for all non-even folds;
THM-3630 reserves precisely that separate all-family question.

## 5. Scope and reproduction

Reproduce with

```bash
python3 04-computation/jc2_russell_cylinder_noneven_closed_form_order_five_thm3631.py
python3 -O 04-computation/jc2_russell_cylinder_noneven_closed_form_order_five_thm3631.py
```

The companion pins the promoted THM-3624 and THM-3627 theorem, script, and
output packages; reconstructs every pullback matrix through degree five using
the pinned THM-3624 exact constructor; and imposes all coefficients of `(8)`.
It contains no Python assertion statement.  The independent hostile audit
supplied a separate reconstruction through the THM-3627 sparse-series engine.

This candidate proves no global regular pair on the Russell cylinder, no
polynomial Darboux decomposition, no Keller map, and no counterexample to the
planar Jacobian conjecture.  Its sharp contribution is the local boundary:
closedness leaks no earlier than the degree-six arbitrary-form obstruction for
this fixed non-even hostile.

The independent hostile audit closes the promotion gate.  **QED.**
