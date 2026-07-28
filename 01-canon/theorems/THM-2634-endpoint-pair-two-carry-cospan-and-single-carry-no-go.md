---
id: THM-2634
title: "Endpoint-pair two-carry cospan and single-carry chronology no-go"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  THM-2625's allocated endpoint determinant is a tensor product of two
  separately integrated signed endpoint sums.  Splitting each sum before
  its endpoint DFT by c=floor(169{x}) mod 13 therefore gives a lawful
  two-carry cospan (c_L,c_R,q,Delta), but no canonical single carry.  The
  carry-difference transform factors into thirteen endpoint-character
  contractions; k=0 is exactly the old THM-2625 determinant bank, while
  difference zero is an endpoint-matched Gram contraction.  Separate
  endpoint totals and singleton norms do not determine that contraction.
  A transverse carry-equivariant physical endpoint section would make
  Delta=beta+alpha c and trivialize the carry torsor exactly when alpha is
  nonzero, but canon supplies neither that section nor the endpoint-matched
  pair twist / same-shell reference needed to form the diagonal.  No
  chronology, row exclusion, or LRC(14) conclusion follows.
source: deep-energy-audit-2026-07-28-endpoint-two-carry
depends_on:
  - THM-2620-endpoint-pair-parabolic-transvection-and-translation-gauge-boundary
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
related:
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2420-affine-shell-cross-reference-composition-and-complete-zero-reference-hostile
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2630-old-wall-affine-clutching-and-successor-sector-no-go
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
script: 04-computation/lrc14_endpoint_pair_two_carry_cospan_thm2634.py
output: 05-knowledge/results/lrc14_endpoint_pair_two_carry_cospan_thm2634.out
script_sha256: 282dded26f7e1f8084f095043eb41df16bf06f87e94ef8a288b3f3af5401ed39
output_sha256: 5bf8bf63cd3b7487216138b0d9076d84ace77be4b30eee2f90a9455c2e8f1287
hash_basis: LF-normalized bytes
---

# THM-2634 -- the endpoint determinant has two carries, not one

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2630 isolates the missing predecessor/successor digit on a physical
root ancestry.  THM-2625 supplies a rich coefficient-side determinant
coordinate.  The tempting next step is to correlate that determinant with
one carry `c=floor(Rx) mod 13`.  That step is ill-typed: the two allocated
endpoint factors have already been integrated separately and need not come
from the same physical `x`.

The repair is a cospan rather than a forced identification.  Retain one carry
on each endpoint factor, prove the exact carry-character factorization, and
state the additional physical section that would make the two sides meet.

## 1. Carry-resolved endpoint factors

Use THM-2625's canonical control, notation, and left-deep allocation:

```text
R=13^2=169,                    N=R*T_DEN,
(X,m,Y)=(13,1,742599),
q=l-r,                         Delta=det(l,r).             (1)
```

For an endpoint point represented by `n/N`, define the right-continuous
endpoint digit

```text
c(n)=floor(R*{n/N}) mod 13
    =floor((n mod N)/T_DEN) mod 13.                         (2)
```

The convention in (2) matters only for the finite signed boundary measure;
it does not assert a positive interior root sheet.  Under a quotient-gauge
translation by `s/13`, the numerator moves by `13s*T_DEN`, so (2) is
unchanged.  Also `R/13=13`, both `X,Y` are divisible by thirteen, and the
extra `m ell_b` phase is unchanged because `c3=0 mod 13`.  In the canonical
clock both endpoint factors therefore split representative-independently.

Write THM-2625's signed endpoint sums as

```text
P_ell=sum_(c in F_13) P_(ell,c),
Q_ell=sum_(c in F_13) Q_(ell,c),                            (3)
```

where `P_(ell,c)` retains endpoints of
`E^ell intersect T^(-2)Q_(1,{a})` with digit `c`, while `Q_(ell,c)` retains
endpoints of `E^ell` with digit `c`.  These are different endpoint lists.
Their endpoint DFTs are

```text
L_c(l)=sum_ell P_(ell,c) zeta_13^(-<tau(ell),l>),
R_c(r)=sum_ell Q_(ell,c) zeta_13^(+<tau(ell),r>).            (4)
```

Consequently

```text
sum_c L_c(l)=Lstar(l),              sum_c R_c(r)=Rstar(r).   (5)
```

The minimal lawful refinement of the allocated endpoint current is

```text
J(l,r;c_L,c_R)
 =C0/169^2 L_(c_L)(l)R_(c_R)(r),                            (6)

S(q,Delta;c_L,c_R)
 =sum_(r:det(q,r)=Delta) L_(c_L)(r+q)R_(c_R)(r).             (7)
```

Summing both carry coordinates in (6)--(7) recovers exactly THM-2625's
`J(l,r)` and `Sstar(q,Delta)`.  This gives the honest cospan

```text
signed left endpoint  <-  (c_L,c_R,q,Delta)  ->  signed right endpoint.
                                                                    (8)
```

It does not give a map from the middle object to one physical carry.

## 2. Exact single-carry no-go

The multiplication in THM-2625 occurs after the two endpoint integrations:

```text
(sum_c P_c)(sum_d Q_d)=sum_(c,d)P_cQ_d.                    (9)
```

Replacing (9) by

```text
sum_c P_cQ_c                                                   (10)
```

is a new diagonal contraction.  It is not an algebraic regrouping of (9).
Even the two separate totals and singleton norms do not determine it.  On
`F_13`, compare

```text
P=delta_0, Q=delta_0,             and
P'=delta_0, Q'=delta_1.                                      (11)
```

In both cases the left and right totals equal one and both singleton squared
norms equal one.  The diagonal contraction is respectively one and zero.
Thus neither full support of the separate endpoint transforms nor survival
of every determinant sector supplies `c_L=c_R`.

The hostile lifts to the entire endpoint plane.  Put

```text
L_0(l)=1 for every l,              R_1(r)=1 for every r,
all other carry bins zero.                                 (11a)
```

Then the collapsed joint array is one on all `28,561` endpoint pairs.  Its
determinant pushforward is nonzero on all `2,185` admissible sectors: a
nonzero-`q` fibre has thirteen terms, while `(q,Delta)=(0,0)` has `169`.
Nevertheless its matched-carry contraction is identically zero.  Moving the
right array from carry one to carry zero leaves every collapsed endpoint
coefficient, determinant sector, total, and singleton norm unchanged, but
makes every admissible matched sector positive.  Thus even THM-2625's exact
full-support conclusion is logically insufficient for the diagonal.

This is the first failed implication:

```text
two endpoint sums in one coefficient formula
   does not imply
two endpoint events on one physical ancestry.                (12)
```

The missing datum is precisely an endpoint-matched polarized Gram observable,
not another statistic of the already collapsed determinant bank.

## 3. Carry-difference Fourier factorization

The strongest unconditional contraction keeps the carry difference.  Put

```text
C_d(q,Delta)
 =1/13 sum_(r:det(q,r)=Delta) sum_c
       L_(c+d)(r+q)R_c(r),                                  (13)

Lhat_k(l)=1/13 sum_c L_c(l) zeta_13^(-kc),
Rhat_k(r)=1/13 sum_c R_c(r) zeta_13^(+kc).                   (14)
```

Taking the normalized transform of (13) in `d` gives the exact identity

```text
Chat_k(q,Delta)
 =sum_(r:det(q,r)=Delta) Lhat_k(r+q)Rhat_k(r).               (15)
```

For `k=0`, equation (15) is `Sstar(q,Delta)/169`.  Thus the old determinant
bank is the neutral carry character.  The matched-carry contraction is

```text
C_0(q,Delta)=sum_k Chat_k(q,Delta),                          (16)
```

with the normalization of (13)--(15).  Equation (16) is mathematically
defined coefficient-side, but it is not automatically a physical current.

The THM-2380 pair-twist mechanism states exactly what would realize it.  If
one can place the two endpoint arrays in one common physical gauge and form

```text
E_(q,Delta,d)(t)
 =1/13 sum_(r:det(q,r)=Delta) sum_c
   |L_(c+d)(r+q)+zeta_13^t conj(R_c(r))|^2.                 (17)
```

Here THM-2625's `R_c` is already the conjugated endpoint factor, so the
physical second amplitude in (17) is `conj(R_c)`.  After subtracting the
singleton norms, the `t=-1` Fourier mode recovers `C_d(q,Delta)`.  One nontrivial
quadrature plus the norms suffices.  Separate
intensities do not.  THM-2420 gives the parallel same-shell statement: a
residue-zero reference can preserve the correlation, but such a reference
is not free and must share the endpoint gauge.

There is a sharp conditional positivity criterion before complex
cancellation.  At one genuinely matched endpoint pair, if the two carry
arrays are nonnegative, then

```text
C_0>0  iff  supp(L) intersect supp(conj(R)) is nonempty.    (17a)
```

In particular `|supp(L)|+|supp(R)|>13` forces `C_0>0`; the threshold is
sharp by two disjoint supports whose sizes sum to thirteen.  THM-2625's
endpoint amplitudes are signed and become complex after the DFT, and its two
lists are separately integrated, so (17a) is a sidecar criterion rather than
an available conclusion.

## 4. The transverse-section criterion

There is a clean conditional closure.  Fix `q!=0`.  A carry-equivariant
physical endpoint section has the affine form

```text
r(c)=r_0+c v,                   l(c)=r(c)+q.                (18)
```

Then

```text
Delta(c)=det(q,r(c))
        =beta+alpha c,

beta=det(q,r_0),               alpha=det(q,v).              (19)
```

Hence `Delta` trivializes the carry torsor if and only if the section is
transverse:

```text
alpha!=0,                      c=alpha^(-1)(Delta-beta).     (20)
```

If `v` is parallel to `q`, then `alpha=0` and `Delta` sees no carry.  A
transverse direction always exists algebraically: for `q=(q_1,q_2)`, take
`v=e_2` when `q_1!=0`, or `v=e_1` when `q_2!=0`.  What is missing is not the
linear algebra; it is a physical rule proving that the endpoint events of
THM-2625 lie on that same section.

More generally, a lawful two-endpoint clutching

```text
Delta=alpha(c_L-c_R)+beta(q,z),       alpha!=0              (21)
```

would identify the matched sector with `Delta=beta`.  Under THM-2620's
common endpoint translation, `beta` changes by `det(q,t)`, so it is section
data rather than a target-difference invariant.  Full nonvanishing of all
THM-2625 determinant sectors does not choose it.

The holotopy square is therefore

```text
physical endpoint pair  ----------------->  (c_L,c_R)
         |                                      |
         | allocated section                    | difference
         v                                      v
      (l,r)  --------------------------->  Delta.           (22)
```

Only the bottom and left maps are presently canonical.  A pair twist supplies
the top matching test; a transverse physical section supplies the right-hand
clutching.  Neither follows from determinant support.

There is a useful but strictly conditional derangement analogy.  Every
nonzero translation of a `C_13` torsor is fixed-point-free.  If each charged
carry component were already known to carry one positive common-root branch
and the transition preserved that same branch, nonzero holonomy would be
impossible; THM-2622 would then give a section.  Full endpoint/current support
is not this branch-survival hypothesis.  The matched pair twist in (17) is a
candidate sidecar for it, not a proof.  This is the `C_13` analogue of
THM-2633's proved derangement-character gate, but the analogy is conditional:
THM-2633 supplies no endpoint pairing, positivity, or branch-survival premise
for this LRC packet.

## 5. Exact canonical control

The companion works over the exact field `F_53`, where `16` has order
thirteen.  It performs `6,591` exhaustive endpoint-digit/gauge checks.  On a
dense deterministic two-carry endpoint plane it computes (13) directly,
computes (15) independently, and checks exact DFT and inverse DFT equality.
It also checks that the neutral character is exactly the old determinant
bank divided by `169`.  The exact support counts are

```text
old determinant bank:       2,143,
carry differences:         27,888,
carry characters:          28,363,
matched difference zero:    2,144.                         (23)
```

These dense-control counts are not claims about THM-2625's canonical
amplitudes; they are hostile arithmetic for the universal identities.

The script separately realizes (11a).  Its collapsed endpoint support is
`28,561`, its collapsed determinant support is all `2,185` admissible
sectors, and its shifted-carry character support is the maximum admissible

```text
2,185*13=28,405,                                       (24)
```

while its matched support is zero.  The carry-zero version has matched
support `2,185`.  Finally, the companion checks `53` exact pair-twist
quadratures and all

```text
168*169*13=369,096                                      (25)
```

affine-section/carry instances.  For every `alpha in F_13`, exactly
`168*13=2,184` pairs `(q,v)` have `det(q,v)=alpha`, confirming the parallel
versus transverse boundary.

Run

```bash
python 04-computation/lrc14_endpoint_pair_two_carry_cospan_thm2634.py
python -O 04-computation/lrc14_endpoint_pair_two_carry_cospan_thm2634.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_endpoint_pair_two_carry_cospan_thm2634.out.
```

## 6. Exact boundary

This theorem proves a typed information-loss boundary and the minimal lawful
two-carry repair on one canonical control.  It does not prove

```text
a common physical endpoint ancestry;
an endpoint-matched pair twist or same-shell reference;
a transverse physical section;
positivity of the signed endpoint boundary measures;
compatibility with THM-2630's successor digit;
one canonical semantic owner/current;
a row exclusion or LRC(14).
```

The next decisive test is no longer another determinant census.  It is to
construct or obstruct one endpoint-matched phase ratio together with a
transverse section on the same physical packet.

QED (candidate; independent hostile audit pending).
