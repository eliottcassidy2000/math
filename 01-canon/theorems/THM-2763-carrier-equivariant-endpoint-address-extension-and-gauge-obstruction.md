---
id: THM-2763
title: "Carrier-equivariant endpoint-address extension and gauge obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The exact
  carrier-address lattice and its one- and two-sided quotient gauges are
  proved algebraically.  The fixed-triangle descent obstruction and one
  gamma-fibre support statements are finite-exact certificates in the first
  cyclotomic field of THM-2625, replayed byte-identically under ordinary and
  optimized Python.  No endpoint current, row exclusion, or LRC(14) follows.
source: root/lrc-endpoint-carrier-gauge-2026-07-28
audit: >
  extended-address-and-gamma-slice-referee-2026-07-28
  (INDEPENDENTLY HOSTILE-AUDITED: reconstructed W.r=l-k and both quotient
  duals/signs; verified 13^3/13^4 dimensions and the e_0 section; replayed
  hardened gauge/gamma scripts under normal/-O with byte-identical stored
  outputs; checked the fixed-carrier full P/Q/H hostile, transported positive
  control, raw-dual {0,2,10} versus full primal DFT support, and one-/two-sided
  scope: PASS)
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
related:
  - MISTAKE-313
script: 04-computation/lrc14_fixed_triangle_carrier_gauge_obstruction.py
output: 05-knowledge/results/lrc14_fixed_triangle_carrier_gauge_obstruction.out
script_sha256: c4c90f10b67734ba6f011390a2b40d431476b722f743da474515a5322cd31835
output_sha256: 64d19d06f265486e779f62b8f526b946a56013fdfaf20dd8009dbe5c0619bf12
secondary_script: 04-computation/lrc14_extended_carrier_gamma_slice.py
secondary_output: 05-knowledge/results/lrc14_extended_carrier_gamma_slice.out
secondary_script_sha256: e5356e2ad0d77c73629efabd1131b09aac708ad5f1f528722d9c2a0dc406ce0e
secondary_output_sha256: c5fb1873901de8e1809ea4067c0d4867437b07c61a0ff0e46e263c4eebbd10bb
hash_basis: LF-normalized bytes
---

# THM-2763 -- carrier-equivariant endpoint-address extension and gauge obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2749's common carrier cannot be inserted into THM-2334 by multiplying
the old `169` representative values and then forgetting where the carrier
harmonics came from.  The failure occurs before target difference, endpoint
determinant, or Radon aggregation: the resulting fixed-carrier coefficient
does not descend under the old `<W>` representative gauge.

There is a canonical repair.  If `k` and `l` are the total carrier harmonics
at the source and target endpoints, then the exact relation address is not
an element of `Lambda=ker W`; it lies in

```text
Lambda_tilde={(r,k,l) in Z^9 x Z^2 : r.W+k-l=0}.       (1)
```

Modulo thirteen, the faithful two-sided target quotient has `13^4=28,561`
elements.  Its quotient retaining only the imbalance `delta=k-l` has
`13^3=2,197` elements; that coarser object is faithful at the collapsed
total-harmonic level for a one-sided carrier insertion, but it forgets the
common harmonic in a two-sided insertion.  On one exact diagonal fibre, only
three dual-twist evaluations are nonzero, while Fourier inversion gives
nonzero coefficients in all
thirteen primal imbalance residues.  None of these statements completes a
THM-2334 current or excludes an LRC row.

## 1. Fixed data and inherited relation packet

Use the proved THM-2749 rail-eight, clock-one, `(s,t)=(0,4)` common carrier
and the deep-owner speed row

```text
p=13,
R=13^6=4826809,
W=(1,14,27,40,53,66,13,13^3,2*13^5),
(X,m,Y)=(12*13^4,1,38*13^4)
       =(342732,1,1085318),
Y-X=W_8=2*13^5=742586.                                (2)
```

Let `L_13` be THM-2334/2309's six-dimensional deep-owner relation packet in

```text
K_13=ker(W:F_13^9 -> F_13).
```

Thus `dim K_13=8`, `dim L_13=6`, and the old target quotient
`K_13/L_13` has `13^2=169` elements.  Its dual is

```text
L_13^perp/<W>.                                        (3)
```

Since `v_13(X)=v_13(Y)=4`, every endpoint phase at this triangle is
representable in the ambient cyclotomic order

```text
R*T/13^4=13^2*T,                                     (4)
```

which is exactly the order already certified in THM-2625.  Minimality of this
conductor is neither needed nor asserted; no new root-order or primality claim
is used.

## 2. The collapsed 169 bank does not descend

For an old character representative `ell`, replacing `ell` by `ell+W`
translates every old present endpoint factor by

```text
-T/13=-22910530602960                                (5)
```

on the common `T`-grid.  The bare target coefficient is unchanged because
`13|Y`.  If the THM-2749 carrier is held fixed, however, exact reduction in
the first THM-2625 field

```text
F_q,  q=352341050142921841,
```

gives the following unequal representative values:

```text
target address  endpoint   value at ell          value at ell+W
(0,0)           source     287659712270709994     0
(0,0)           target     248235870634784933     0
(1,0)           source     287659712270709994     0
(1,0)           target     248235870634784933     0.             (6)
```

The nonzero difference modulo a certified prime proves inequality in
characteristic zero.  These four rows are a hostile certificate; the theorem
does not infer characteristic-zero vanishing from the displayed zero images.

At `(0,0)` the stronger separately allocated `P/Q/H` hostile is

```text
carrier  fixed ell: (P,Q,H,mass)                  ell+W, carrier fixed
source   (189041250036777056,287659712270709994,
          65280867241115379,6320326320)            (0,0,0,0)
target   (218344733173586894,248235870634784933,
          209108233808250489,6320326320)           (0,0,0,0).    (6a)
```

Transporting the carrier by the same `-1/13` physical shift restores the
coefficient identity at all eight tested source/target rows.  More generally,
that covariance is uniform: change variables by `-1/13`; the endpoint phases
are trivial because `13|X,Y`, the delayed word is unchanged because `13|R`,
and carrier translation supplies exactly the phase omitted by the fixed
carrier.  Therefore the fixed collapsed carrier gives a section-dependent
array of 169 numbers, not a function on (3).

## 3. The faithful two-sided exact-address lattice

Expand the source and target carriers into total Fourier harmonics `k` and
`l`.  An atomic two-sided endpoint term satisfies

```text
(u+R beta).W+k=X,
v.W+l=Y.                                             (7)
```

Keep THM-2334's old raw expression

```text
r=u+R beta+m e_8-v.                                  (8)
```

Using `Y=X+mW_8`, equations (7)--(8) give

```text
r.W=(X-k)+(Y-X)-(Y-l)=l-k,
r.W+k-l=0.                                           (9)
```

This proves (1).  The old current embeds as the `k=l=0` slice.  The larger
diagonal locus is

```text
{(r,k,k): r in Lambda, k in Z};                      (10)
```

it must not be conflated with the old lattice.

Reduce (1) modulo thirteen:

```text
K_full=ker(W,1,-1) subset F_13^(9+2).                (11)
```

The normal has a unit coordinate, so reduction from `Lambda_tilde` onto
`K_full` is surjective.  Since (11) has dimension ten and `(L_13,0,0)` has
dimension six,

```text
G_full=K_full/(L_13,0,0),
dim G_full=4,
|G_full|=13^4=28561.                                 (12)
```

Writing `ell in L_13^perp`, its character group has the faithful description

```text
G_full^ =
 {(ell,a,b): ell in L_13^perp, a,b in F_13}
 / <(W,1,-1)>,                                       (13)

(ell,a,b) ~ (ell+sW,a+s,b-s).                        (14)
```

The `a` and `b` coordinates retain the separate source and target carrier
harmonics.  Equation (14) is the missing representative-gauge law.  Setting
`a=b=0` after choosing representatives destroys it, exactly as (6) detects.

A Bezout vector `z.W=1` gives

```text
r_sharp=r-(l-k)z in Lambda.                          (15)
```

For this canonical row `W_0=1`, so `z=e_0` is a distinguished typed section.
Retaining provenance gives the explicit coordinate isomorphism

```text
G_full -> (K_13/L_13) direct_sum F_13^2,
[(r,k,l)] |-> ([r-(l-k)e_0],k,l).                  (15a)
```

Its inverse sends `([r_sharp],k,l)` to
`[(r_sharp+(l-k)e_0,k,l)]`, so the signs in (15a) are forced by (9).
But (15) folds carrier imbalance into the guard coordinate if the `(k,l)`
provenance is discarded.  It forgets
whether that coordinate came from the old endpoint word or the carrier and
therefore does not preserve the atomic-factor or all-`91`-unit masks.  Even
the `28,561` quotient retains only total endpoint harmonics `k,l`; a labelled
carrier-factor allocation remains a necessary sidecar for a factorwise
THM-2625 transplant.

The equality `|G_full|=28,561` happens to match the size of THM-2625's joint
endpoint plane.  After choosing a target basis, (15a) is abstractly a copy of
`G direct_sum G`; it is not a canonical identification with THM-2625's
physical `(left,right)` endpoint addresses.  No transplant map is proved by
the cardinality coincidence.

## 4. The 2,197-element imbalance quotient

Put

```text
delta=k-l.
```

Forgetting the common harmonic maps (11) onto

```text
K_delta=ker(W,1) subset F_13^(9+1),
G_delta=K_delta/(L_13,0).                            (16)
```

Its kernel in `G_full` is the common-harmonic line
`{(0,t,t):t in F_13}`.  Consequently

```text
dim G_delta=3,
|G_delta|=13^3=2197,                                 (17)
```

and

```text
G_delta^ =
 {(ell,gamma): ell in L_13^perp, gamma in F_13}
 / <(W,1)>,

(ell,gamma) ~ (ell+sW,gamma+s).                      (18)
```

Under the dual injection `G_delta^ -> G_full^`,

```text
(ell,gamma) |-> (ell,gamma,-gamma).                  (19)
```

Thus (16) is the diagonal two-sided character slice, or equivalently the
quotient that forgets common carrier harmonic.  If a carrier is inserted
only at the left endpoint (`l=0`), then `delta=k` and (16) is the full
faithful collapsed-total-harmonic residue extension.  It is not the full
allocation object when both endpoint carriers are inserted separately, and
even one-sided use needs labelled atomic refinement for factorwise or
all-`91`-unit masks.

## 5. One diagonal gamma fibre is fully supported after inversion

Fix `ell=(0,0)` in the old target-coordinate section and, for
`gamma=0,...,12`, translate both THM-2749 endpoint carriers by `-gamma/13`.
Let `P_gamma` be the source-carrier/full-terminal coefficient and `Q_gamma`
the target-carrier coefficient with no terminal factor.  This realizes the
diagonal character `(a,b)=(gamma,-gamma)` in (19).

The physical intersection masses are exactly

```text
source: gamma 0 -> 6320326320,
        gamma 2 -> 4072511520,
        gamma 10 -> 6320326320;

target: gamma 0 -> 6320326320,
        gamma 2 -> 4046066640,
        gamma 10 -> 6320326320,                      (20)
```

and zero for every other gamma.  The three overlap numerators are

```text
374977060157700,
241617017842200,
374977060157700.                                    (21)
```

At all three supported twists, `P_gamma`, `Q_gamma`, and
`H_gamma=P_gamma Q_gamma` have nonzero images in `F_q`.  Hence their raw
dual-twist support is exactly

```text
supp P=supp Q=supp H={0,2,10}.                       (22)
```

Now take the thirteen-point DFT in `gamma`.  Every one of the thirteen DFT
entries is nonzero for each of `P`, `Q`, and `H`.  Fourier inversion therefore
shows that every primal imbalance residue occurs nontrivially on this fibre.
This does **not** say that all thirteen raw dual twists are nonzero; (22) is
the exact opposite statement.

For a genuinely one-sided insertion, the uncarried right endpoint at
`ell=(0,0)` has certified image

```text
Q_bare=272457584061297438 !=0.                       (23)
```

Multiplying the full-support DFT of `P_gamma` by this constant shows that all
thirteen primal imbalance residues also occur in the one-sided current on
this fibre.

## 6. Exact boundary

This theorem proves precisely:

```text
fixed collapsed carrier + old 169 representatives
    -/-> a quotient-invariant THM-2334 bank;

separate total harmonics (k,l)
    -> a lawful 28,561-element two-sided residue quotient;

imbalance delta=k-l
    -> a lawful 2,197-element coarse quotient,
       faithful for one-sided insertion;

one diagonal ell-fibre
    -> raw twist support {0,2,10}
       and full thirteen-residue inverse-transform support.           (24)
```

It does not yet construct the factor-labelled 27-plus-carrier atomic word,
identify one old exact relation address, prove an all-`91`-unit coefficient,
form a nonzero endpoint determinant sector, survive target difference or
Radon aggregation, complete a THM-2334 current, exclude a scalar row, or prove
LRC(14).

MISTAKE-313 is a correction boundary, not a dependency: no provisional
THM-2751 wing value is used here.  The actual source-clocked two-sided
THM-2749 carrier is used throughout.

## 7. Reproduction and evidence addresses

From the repository root:

```bash
python3 04-computation/lrc14_fixed_triangle_carrier_gauge_obstruction.py
python3 -O 04-computation/lrc14_fixed_triangle_carrier_gauge_obstruction.py

python3 04-computation/lrc14_extended_carrier_gamma_slice.py
python3 -O 04-computation/lrc14_extended_carrier_gamma_slice.py
```

Each ordinary/optimized pair byte-matches its corresponding transcript in
`05-knowledge/results`.  The scripts use explicit `require` gates, contain no
Python `assert` nodes, and pin the shared exact helper
`04-computation/lrc14_extended_carrier_endpoint_lib.py` at LF-normalized
SHA-256
`4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57`.
The script and output addresses are pinned in the front matter.
