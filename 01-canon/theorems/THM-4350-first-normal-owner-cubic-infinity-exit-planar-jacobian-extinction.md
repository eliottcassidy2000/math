---
id: THM-4350
title: "First-normal-owner cubic infinity-exit planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327/4344 + VERIFIED-EXACT + TWICE
  INDEPENDENTLY AUDITED. In the inherited reduced (2,3), exact-weight-twelve
  seam, Z=beta_11=zeta_3=W=xi_10=0, U*K!=0, and
  (alpha_11,Theta)!=(0,0) are extinct. The 104 exact support patterns have
  exactly three fans. Each has graph genus eleven and one automatically
  smooth genus-three replacement face of positive Keller-form order. A
  512-configuration Boolean atlas has only 336 distinct supports and is not
  a realizability census. The alpha_11=Theta=0 corner, U=0, K=0, seam entry,
  JC(2), and DC(2) are not claimed here.
source: root + cubic-root/fan/global hostile-referee agents / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4344-clean-cubic-infinity-exit-planar-jacobian-extinction
related:
  - THM-4339-clean-interior-cubic-edge-planar-jacobian-extinction
  - THM-4341-odd-self-similar-cusp-reciprocal-tail-duality
  - THM-4342-clean-cubic-zero-exit-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_thm4350.py
primary_output: 05-knowledge/results/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_thm4350.out
primary_script_sha256: b9257f907386398a7b52b2136de3c7501901d41ed76202a0e133d3ea8af5f640
primary_output_sha256: 9da80a441b041cf61fc06b3a2f7132d441978838f64ce887eefa22d9d5e45b33
local_referee_script: 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_local_referee_thm4350.py
local_referee_output: 05-knowledge/results/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_local_referee_thm4350.out
local_referee_script_sha256: 14c8a1e02565eccdaf35000ec2308cea31d1d2b5576275ca6eb0d6e29e12d72c
local_referee_output_sha256: c6260cefe6bacedc8a6b6ce049a91ec58b53dec0af99c4ec5282643e8dd0f750
global_referee_script: 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_global_referee_thm4350.py
global_referee_output: 05-knowledge/results/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_global_referee_thm4350.out
global_referee_script_sha256: 6b740591a5558caad53d8f2a7417128e3775a864be25fc32b160fe3419011aa2
global_referee_output_sha256: 0753f4db583c9310e75dff74cfbc06bae8444d20a8e04c7d9bdff0186071c48f
hash_basis: raw LF bytes
audit: >
  PASS AFTER THE SECOND MASK-QUOTIENT REPAIR. The 894-check primary rebuilds
  the complete source, all four owner fans (including the excluded corner),
  exact faces, carrier normalizations, edges, packets, graphs, form orders,
  and the next even-A5 sidecar. The 213-check local referee independently
  reconstructs all three claimed primitive charts. The import-free global
  referee distinguishes 512 labelled Boolean configurations from 336
  distinct supports, embeds all 128 exact supports, and checks every face,
  edge, Pick/Riemann--Hurwitz ledger, and owner transition. Normal and
  optimized streams byte-match all three frozen outputs.
---

# THM-4350 -- First-normal-owner cubic infinity-exit planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327/4344 + VERIFIED-EXACT + TWICE
INDEPENDENTLY AUDITED. THE DISPLAYED FIRST-NORMAL-OWNER UNION IS EXTINCT.
THE DOUBLE-OWNER-ZERO CORNER, `U=0`, `K=0`, SEAM ENTRY, `JC(2)`, AND `DC(2)`
REMAIN OUTSIDE THIS THEOREM.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)` seam, retaining the
complete sixteen-term source of THM-4230 and

```text
e=-1376/135,                    K=2848/45-(7/6)Delta.       (1)
```

Write `u=upsilon_5` and `alpha=alpha_11`. Impose

```text
Z=beta_11=zeta_3=W=xi_10=0,
U*K!=0,                         (alpha,Theta)!=(0,0).       (2)
```

Then no nonautomorphic planar Keller pair lies on `(1)--(2)`. Every lower
coefficient not displayed in `(2)` is arbitrary, subject only to the
inherited seam. The conclusion is relative to the inherited toroidal and
proper-flat target interfaces; it proves neither entry into this seam nor
`JC(2)` or `DC(2)`.

The inheritance pass was:

- closest proved mechanism: THM-4344's replacement-face principle and
  positive-order proper-flat extinction;
- canonical hostile: deleting a leading owner can expose a face of *larger*
  genus, so genus is not monotone under support deletion;
- corrected near miss: independently toggled row and cancellation labels
  form a configuration over-atlas, not a census of realizable or even
  distinct supports;
- least-used sidecar: THM-4341's attachment warning. It is dormant on `(2)`
  because every carrier is squarefree, but it detects the exact next wall.

The Anchor / Niche / Wildcard board was

```text
first normal owner | support quotient | face merger | carrier genus
primitive base | edge word | graph genus | lawful natural index.             (3)
```

## 2. Exact support quotient

On `(2)`, the four decisive support coefficients are uniquely owned:

```text
(2,6,1)=-U,       (4,2,1)=-K,
(3,5,1)=-alpha,   (4,3,1)=-Theta.                       (4)
```

The relevant multiply visible coefficients are

```text
(2,3,1)=K+1376/135,
(2,4,1)=Theta-Delta,
(2,5,1)=-u.                                            (5)
```

The last item is a coupled source coefficient, not an independently
deletable aggregate. The seam gives

```text
K=0 iff Delta=5696/105,
K+1376/135=0 iff Delta=3968/63.                         (6)
```

Away from `K=0`, the status tuple

```text
(Delta, K+1376/135, Theta, Theta-Delta)
```

has eight realizable classes. Combining it with the independent status bits
of `u,Phi,eta,alpha` gives exactly 128 distinct realizable supports. Their
owner census is

| `(alpha,Theta)` status | exact supports | lower fan |
|---|---:|---|
| `(1,1)` | 40 | `M,D6,E11,T` |
| `(1,0)` | 24 | `M,D6,E10` |
| `(0,1)` | 40 | `M,E01,T` |
| `(0,0)` | 24 | `M,E00` |

Here `1` means nonzero and `0` means zero. Thus the theorem union `(2)`
contains exactly

```text
40+24+40=104                                             (7)
```

realizable supports, and all 104 have one of three stated fans.

For hostile coverage, independently toggling six optional source rows and
three deletion labels produces `2^9=512` *keyed configurations*, 128 per
owner fan. Quotienting equal active point sets leaves only 336 distinct
supports. Of these, 128 are exact and 208 are additional synthetic supports.
The canonical coefficient-coupling map embeds every exact support into an
identical Boolean configuration. The other 384 keyed configurations are
coupling-illegal; 88 nevertheless duplicate an exact active support, while
the remaining 296 collapse to the 208 synthetic supports. These three
numbers answer different questions and are not interchangeable.

## 3. The three replacement fans and their merger law

For a plane `(a,b,c)`, the lifted height is `a*r+b*l+c`. Up to invertible
torus monomials, the complete claimed face list is:

| face | plane | polygon vertices | initial equation | `(2A,B,I)` |
|---|---|---|---|---:|
| `M` | `(1/12,1/6,-1/6)` | `(0,1),(2,0),(2,6),(0,7)` | `(S^2-P)(1-UP^6)` | `(24,14,6)` |
| `D6` | `(1/6,1/6,-1/3)` | `(2,0),(2,6),(3,5)` | `S^2(1-UP^6-alpha SP^5)` | `(6,8,0)` |
| `E11` | `(2/7,1/7,-4/7)` | `(2,0),(3,5),(4,3)` | `S^2(1-alpha SP^5-Theta S^2P^3)` | `(7,3,3)` |
| `T` | `(1/2,0,-1)` | `(2,0),(4,2),(4,3)` | `S^2(1-S^2P^2(K+Theta P))` | `(2,4,0)` |
| `E10` | `(3/8,1/8,-3/4)` | `(2,0),(3,5),(4,2)` | `S^2(1-alpha SP^5-KS^2P^2)` | `(8,4,3)` |
| `E01` | `(1/4,1/6,-1/2)` | `(2,0),(2,6),(4,3)` | `S^2(1-UP^6-Theta S^2P^3)` | `(12,8,3)` |

The local referee's scale aliases are `D7=E11`, `D8=E10`, and `D12=E01`.
They name the same faces, not three additional components.

The equations are reconstructed from the literal source, not inferred from
the polygon alone. Exhaustion of every exact support and all 512 hostile
configurations gives the same fan within each owner class.

The wall crossings have a useful exact merger description:

```text
alpha -> 0, Theta!=0:       D6 + E11  -> E01,       T survives;
Theta -> 0, alpha!=0:       E11 + T   -> E10,       D6 survives.   (8)
```

Thus deletion of an owner contracts two adjacent faces into a larger
replacement face. It does not erase their arithmetic contribution.

The global Pick and component ledgers are

| owner case | global `(2A,B,I)` | `(V,E,b1)` | carrier genus | total genus |
|---|---:|---:|---:|---:|
| `11` | `(39,13,14)` | `(10,20,11)` | 3 | 14 |
| `10` | `(38,12,14)` | `(9,19,11)` | 3 | 14 |
| `01` | `(38,12,14)` | `(9,19,11)` | 3 | 14 |

The equality of total genera across `(8)` is therefore a proved graph-plus-
normalization identity, not a visual inference from the fan.

## 4. Exact normalizations and positive form orders

The `D6` equation is rational by solving linearly for `S`; when present, `T`
is rational after writing `z^2=K+Theta P`. The sole positive-genus face has
one of the exact hyperelliptic normalizations

```text
E11: y=2Theta*S*P^2+alpha*P^4,
     y^2=P(4Theta+alpha^2P^7);                          (9)

E10: y=2K*S*P+alpha*P^4,
     y^2=4K+alpha^2P^8;                                (10)

E01: y=Theta*S*P^2,
     y^2=Theta*P(1-UP^6).                              (11)
```

The branch degrees are respectively `8,8,7`. On their owner gates the first
polynomial has the simple root `P=0` and seven distinct nonzero roots; the
second has eight distinct roots because `K*alpha!=0`; the third has the
simple root `P=0` and six distinct nonzero roots because `U*Theta!=0`.
Consequently all three smooth projective normalizations have genus three.
There is no residual discriminant inside `(2)` and hence no critical-value
series or collision tail to resolve.

For completeness, the literal primitive source charts used by the
certificates are

```text
D6:  Q=sigma^6,  s=sigma^-1 S, p=sigma^-1 P, G=sigma^2 F_Q;
E11: Q=sigma^7,  s=sigma^-2 S, p=sigma^-1 P, G=sigma^4 F_Q;
E10: Q=sigma^8,  s=sigma^-3 S, p=sigma^-1 P, G=sigma^6 F_Q;
E01: Q=sigma^12, s=sigma^-3 S, p=sigma^-2 P, G=sigma^6 F_Q.       (12)
```

For example, the least familiar `E01` chart expands without ellipsis as

```text
G=(S^2-sigma^4P)
  [1-UP^6-Theta*S^2P^3-sigma*eta*SP^4
     -sigma^2(uP^5+K*S^2P^2)-sigma^3Phi*SP^3
     -sigma^4Delta*P^4-e*sigma^6P^3
     -(8/3)sigma^8P^2+3sigma^10P]
  -sigma^12S^2/2.                                      (13)
```

Equation `(13)` proves both the displayed face and the absence of a hidden
lower-order owner at the `alpha=0` transition.

For a supporting plane `(a,b,c)`, the inherited Keller-form density is

```text
5/6-(a+b+c).                                            (14)
```

Clearing both source and target denominators gives the positive orders

| face | least common base | form order |
|---|---:|---:|
| `M` | 12 | 9 |
| `D6` | 6 | 5 |
| `E11` | 42 | 41 |
| `T` | 6 | 8 |
| `E10` | 24 | 26 |
| `E01` | 12 | 11 |

In particular, the order `26` for `E10` is on the common base 24; its density
is `13/12`. This prevents the common error of reporting denominator 12 as an
actual primitive base for the full source-target chart.

## 5. Edges, packets, and component completeness

Up to units and coordinate reversal, the outer boundary schemes are

```text
11: X-1; 1-KX^2; K+Theta X; Theta+alpha X;
    alpha+UX; X-1; U-X^6.

10: X-1; 1-KX^2; K+alpha X; alpha+UX; X-1; U-X^6.

01: X-1; 1-KX^2; K+Theta X; Theta+UX; X-1; U-X^6.       (15)
```

The internal schemes are

```text
11: 1-UX^6; 1-alpha X; 1-Theta X;
10: 1-UX^6; 1-alpha X;
01: 1-UX^6; 1-Theta X.                                 (16)
```

Every scheme in `(15)--(16)` is reduced in characteristic zero under its
displayed owner conditions and `U*K!=0`. Hence the toric intersections are
transverse and the face normalizations listed above are complete.

The source-infinity packets are

```text
11=(11,8,6,3,2,2,1),
10=(11,10,6,2,2,1),
01=(13,11,3,2,2,1).                                    (17)
```

Each satisfies `sum(e_i-1)=26=2*14-2`. This is an independent
Riemann--Hurwitz check, not the proof of completeness.

The `M` face consists of seven rational components with twelve internal
nodes. Case `11` adds six `M--D6` attachments and one attachment at each of
`D6--E11` and `E11--T`, giving `V=10,E=20,b1=11`. Case `10` has
`V=9,E=19,b1=11`; case `01` has the same ledger after the merger `(8)`.
Adding the unique genus-three carrier gives total genus fourteen in all
three cases, as required by the global polygons.

## 6. Proper-flat extinction

After a finite common base change and toric regularization, every special
component on `(2)` is either

1. a rational component from `M`, `D6`, `T`, or a toric chain; or
2. the smooth genus-three component `E11`, `E10`, or `E01`, whose pulled-back
   Keller differential has strictly positive order by `(14)`.

A rational curve cannot map nontrivially to the good elliptic target, and
the positive-order differential makes the map from the genus-three carrier
constant. Retaining the actual positive component multiplicities on a
common dominating base, the inherited proper-flat identity gives, for a
positive-degree target line bundle `L`,

```text
deg(phi_generic^*L)=sum_Gamma m_Gamma deg(phi_Gamma^*L)=0.       (18)
```

This contradicts the positive generic response degree of a nonautomorphic
Keller pair. Therefore the gate `(1)--(2)` is extinct. **QED.**

## 7. Unexpected connections and lawful natural indexing

The exact owner diagram is the nonzero part of a Boolean square. If

```text
a=1_[alpha!=0],                 t=1_[Theta!=0],
n=1+2(1-a)+(1-t),                                      (19)
```

then cases `11,10,01,00` receive the natural labels `1,2,3,4`. Equivalently,
one may call `(2n-1)^2` the `n`th odd-square address and store only `n` as
the ordinal. This is an exact indexing of owner statuses on the fixed
`xi_10=W=0`, `U*K!=0` block.

It is not a sufficient geometric invariant. The integer `n` forgets the
supporting plane, which pair of faces merged in `(8)`, the primitive base,
the edge word, and the carrier equation. The lawful compressed object is

```text
(n; owner divisor, merger arrow, plane, edge word, primitive base).       (20)
```

This makes precise the same lesson seen in the tournament and rational-tree
work: an orientation or ordinal can index a discrete choice, but a tie must
retain its resolving polynomial and incidence address. Here valuation
comparison is a preorder, not a tournament. Away from equality the lower
valuation wins; at equality, the initial polynomial `(9)--(11)` carries
genus. Arbitrarily orienting that tie would destroy the invariant one is
trying to encode.

The surprising reusable result is therefore a nonmonotonicity principle:
support deletion can raise the genus of the unique carrier while preserving
total genus through a face contraction. THM-4344 has a genus-two `D6`
carrier when `xi_10!=0`; setting `xi_10=0` exposes a genus-three carrier in
every first-normal-owner direction. Any future owner-descent algorithm must
recompute the lower hull and graph rather than inherit the preceding genus.

## 8. Reproduction and sharp next boundary

Run from the repository root:

```bash
python3 -B 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_thm4350.py
python3 -B -O 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_thm4350.py
python3 -B 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_local_referee_thm4350.py
python3 -B -O 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_local_referee_thm4350.py
python3 -B 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_global_referee_thm4350.py
python3 -B -O 04-computation/jc2_m12_first_normal_owner_cubic_infinity_exit_extinction_global_referee_thm4350.py
```

All three normal/optimized pairs byte-match their frozen outputs. The primary
contains 894 exact checks, the local referee 213, and the global referee an
independent import-free reconstruction.

The sharp omitted corner is

```text
alpha=Theta=0.                                           (21)
```

Its fan is `M,E00`; generically `E00` has genus two, while its sole candidate
collision is `eta^2=4KU`. The scripts verify an exact even-`A5` sidecar there,
but `(21)` is deliberately not imported into this theorem. It requires its
own promoted chart and proper-flat proof. Beyond it remain `U=0`, `K=0`, seam
entry, and the global `JC(2)`/`DC(2)` interfaces.
