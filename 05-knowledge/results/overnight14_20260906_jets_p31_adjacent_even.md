# Adjacent even ideals close an all-depth pair of p31 Smith factors

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED; FINITE-EXACT controls.** No theorem ID is claimed. This note is separate
from the frozen thirteenth result and does not change its statements or hashes.

## 1. Inheritance and the target actually closed

The [thirteenth theorem](overnight13_20260906_jets_p31_intermediate.md)
and its [independent audit](overnight13_20260906_jets_p31_intermediate_audit.md)
prove, for the complete 16-jet Hasse bank on nodes
`0,31^e,31^e a`, with `a,a-1` units in Z_31 and e>=1,

```
E29=631e+1+kappa,
kappa=[a mod31 in {3,11,15,17,21,29}].
```

Here E_r is the sum of the first r residual Smith exponents after the
16 node-zero unit factors. It equals the valuation of the full ideal D_(16+r).
The thirteenth result transported smaller Hermite cofactors into two minimum
odd minors and supplied the next-band divisibility needed at depth one.
Its audit identified the adjacent even ideals as a precise next target.

The canonical hostile is a=3 versus a=4: their metric, ordinary Deuring
bit, AP bit and largest factor agree, while intermediate factors and finite
kernels differ. The corrected near miss is the metric-only extrapolation
recorded in [MISTAKE-547](../../01-canon/MISTAKES.md); that correction also
warns against extrapolating a shallow minimum while omitting another weight
slope. A targeted search of the current theorem/result surface found the
present pair formulas only in the thirteenth audit's explicitly unproved
next-step paragraph. No external priority claim is made.

The live concept board is: complete Hasse banks; even full determinants;
degree/derivative shifts; integral content; next-weight support; and the
kernel precision window. The least-used sidecar is the unique minimum
even minor. It has no nontrivial residue packet, so its content and the next
weight gap can settle an adjacent ideal without another large polynomial
bank.

## 2. Statement and equality boundary

For the same complete observer and every e>=1, every admissible p-adic lift a,

```
E28=v31(D44)=588e+2,
E30=v31(D46)=675e+1.                                (1)
```

Writing lambda_r for the r-th residual Smith exponent, the inherited E29
then gives

```
lambda29=43e-1+kappa,
lambda30=44e-kappa.                                 (2)
```

These are positions 45 and 46 in the full ordered 48-factor list. The two
factors are equal exactly when e=1 and kappa=1. Their sum is always 87e-1;
their gap is e+1-2*kappa. At e=0 the observer is unimodular and both factors
are zero. Formulas (1)--(2) are not applied at that boundary.

This is an all-depth and all-lift two-factor law, not a complete p31 Smith
partition. It asserts nothing about whether the other factors depend on
additional residue data.

## 3. General minimum even minor and the scalar it retains

After clearing the m rows at node zero, use normalized residual entries

```
binom(c,j) x^(c-j),
x in {1,a}, 0<=j<m, m<=c<3m.
```

For rank 2s, where 1<=s<=m, the unique least column sum uses
`c=m,...,m+2s-1`; the unique largest derivative-order sum uses both copies
of `j=m-s,...,m-1`. Put h=m-s. The resulting determinant is exactly

```
C_(m,s) [a(a-1)]^(s^2),
C_(m,s)=product_(c=s)^(3s-1) [(c+h)!/c!]
                / product_(j=0)^(s-1) [((j+h)!/j!)^2].       (3)
```

To prove this, first take m=s. The selected residual matrix is the full
residual observer, whose determinant is the full confluent determinant
`a^(s^2)(a-1)^(s^2)` after clearing the identity block at node zero. For
general m, lower each selected degree and derivative order by h, using

```
binom(c,j) = [(c)_h/(j)_h] binom(c-h,j-h).
```

The power c-j is unchanged. Factoring the rational scalars by columns and
rows yields (3). This is a determinant identity; it does not assert that
nonunit derivative scalings preserve the Smith form. The monomial in (3)
is primitive, and the determinant is an integer polynomial, so Gauss's
lemma also proves C_(m,s) is integral.

On the physical scaled nodes, the minor has the additional factor
`31^(3s^2 e)`. Indeed its degree sum minus derivative-order sum is

```
sum_(c=m)^(m+2s-1)c - 2 sum_(j=m-s)^(m-1)j = 3s^2.
```

These factors are extracted formally inside each determinant. The lower
bounds below retain both this weight and the integral coefficient content.

## 4. Rank 28: two content factors and the complete next band

Take m=16,s=14,h=2. The unique minimum rows are j=2,...,15 at both nodes;
the unique minimum columns are 16,...,43. Their sums are 238 and 826,
so the weight is 588. Formula (3) gives

```
C28=43!42!/[(15!)^3(14!)^3]
   =57292925200157771047992775131648000,
v31(C28)=2.                                         (4)
```

Each numerator factorial contains one factor 31, and neither denominator
factorial contains one. The unit monomial is [a(a-1)]^196. Hence this minor
has valuation exactly 588e+2 for every admissible lift.

Every other minor has nonnegative integral column excess plus row deficit.
There are exactly five minors of total excess one. Four replace one of
the two retained order-two rows by one of the two omitted order-one rows.
The fifth replaces the last retained column 43 by 44. There are no other
choices: both row and column losses are nonnegative integers.

Every one of these five row sets has only j>=1, and every selected column
set contains column 31. For 1<=j<=15,
`binom(31,j)=0 mod31`; therefore column 31 vanishes identically over
F31[a]. Each normalized next-band determinant has all coefficients
divisible by 31. Its valuation is at least 589e+1, which is at least
588e+2 for every e>=1. Higher bands have valuation at least 590e, also
at least 588e+2. The unique minimum attains that value. This proves the
first equality in (1), including the shallow e=1 boundary.

This content step is necessary for the proof: a weight-only estimate of
a next-band minor would give 589 at e=1, below the attained minimum 590.
The actual column support supplies the missing factor 31.

## 5. Rank 30: the weight gap already suffices

Now take m=16,s=15,h=1. The unique minimum rows are j=1,...,15 at both
nodes; the unique minimum columns are 16,...,45. Their sums are 240 and
915, yielding weight 675. The exact scalar and content are

```
C30=45!/(15!)^3=53494979785374631680,
v31(C30)=1.                                         (5)
```

The normalized determinant is C30[a(a-1)]^225, so the minimum minor has
valuation exactly 675e+1. Every other minor has weight at least 676 and
integral coefficients. Thus its valuation is at least 676e>=675e+1 for
every e>=1. This proves the second equality in (1) without assuming any
extra divisibility in the next band.

The complete next band again contains five minors: four replace an
order-one row by an order-zero row and one raises the last column by one.
The rank-28 zero-column proof does not automatically transfer, since some
of these rows have j=0. The present proof needs only the ordinary weight
gap, so it does not make that false structural extrapolation.

## 6. Exact two-factor kernel window

Subtracting adjacent ideal valuations establishes (2). For integer
precision b>=0, write the contribution of this pair to the logarithmic
kernel size as

```
K_pair(b,kappa)=min(b,43e-1+kappa)+min(b,44e-kappa).
```

Directly from the two breakpoints,

```
K_pair(b,1)-K_pair(b,0)
   = 1  if 43e<=b<=44e-1,
   = 0  otherwise.                                  (6)
```

At e=1 this window consists of b=43, recovering the finite 43,43 versus
42,44 control. At larger depth the window has exactly e integer precisions.
The pair sum, and therefore its determinant contribution, is unchanged.
The unit class transfers one valuation unit from the larger factor to the
smaller one, and capped kernel sums detect this transfer precisely in (6).

Equation (6) concerns the two identified factors. It does not imply that
the entire kernel orders for arbitrary inputs in the two classes differ
by exactly a factor 31: other factors have not been classified here. When
a full Smith comparison establishes that all other factors coincide, it
does give that full-kernel consequence. The thirteenth finite a3/a4
depth-one matrices are such a control.

The connection contract is explicit: adjacent determinantal ideals map to
ordered factors by subtraction; ordered factors map to kernel exponents
by the function lambda -> min(b,lambda). The first map requires the full
ideals, not selected minors. The second retains each factor individually;
replacing the pair by its sum destroys exactly the precision window above.

## 7. Exact certificate, scope, and reproduction

The standalone companion has 126,363 always-active gates and imports no
previous producer. It verifies 44 signed literal even determinants: 36 at
s=1,...,6 and m=s,s+1,s+3, and eight large normalized rank-28/rank-30
determinants at a=2,3,4,-29. The exact scalar products agree with the closed
factorial formulas, their integer valuations, and an independent Legendre
factorial-valuation computation.

It enumerates all row and column index sets separately: 35,960 at rank 28
and 496 at rank 30. This verifies the unique minimum and complete next
band for both ranks. The corresponding index-pair universes contain
1,293,121,600 and 246,016 minors; their determinant polynomials are not
enumerated. The proof reduces coefficient work to the unique minimum and,
at rank 28 only, five next-band support obstructions.

Ten literal original 48-by-48 Hasse matrices verify E28,E29,E30 and the
actual factor pair. They use a=3 and a=4 at e=0,1,2,5, plus a=879 at e=3
and the negative lift a=-20 at e=4. The deep individual-packet cancellation
at 879 remains a hostile to replacing the common odd ideal by one packet.
Depth zero is checked separately. The finite matrix elimination uses
fixed-precision least-valued DVR pivots and checks the complete determinant
valuation 768e. The two-factor kernel window is also checked exactly for
e=1,...,20 and all relevant integer precisions; formula (6) itself is proved
by its two breakpoints, not extrapolated from that bank.

```
python -B 04-computation/overnight14_20260906_jets_p31_adjacent_even.py
python -B -O 04-computation/overnight14_20260906_jets_p31_adjacent_even.py
```

Both modes emit identical LF bytes directly, with stdout explicitly configured
for LF; no capture-time newline replacement is needed. Source SHA256:
`e0eb5d1b403187706ccb9902747a58060e157885d5e8cfbefea12d757eb67bfb`.
Output SHA256:
`9c7cd02d21bfea5ffdb80dac65cf572fd9904909b18e35a2d7a1d6ce575e1288`.

The generalization suggested by the mechanism is to pair the universal
odd packet with its two neighboring even determinants and audit only the
weight levels that can beat their coefficient contents. That is a route
to further factor blocks. It is not asserted to succeed uniformly at every
prime, and the full p31 partition remains outside this note.

The [independent audit](overnight14_20260906_jets_p31_adjacent_even_audit.md)
accepts the full proof with132,821 gates, including120 fresh residual-layer
matrices, all minimal/next bands, and signed rational determinant controls.

**Filing:** root read these proofs and audits and integrated the fourteenth
checkpoint. Reproduction commands are relative to the repository root.
