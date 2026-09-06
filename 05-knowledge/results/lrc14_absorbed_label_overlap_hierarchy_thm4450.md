# Dyadic five-ray bank: overlap hierarchy and exact component decoder

**Status: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED; LRC(14) OPEN.**
The equal-radius overlap formula is **RECOVERED** from
THM-739/LEM-042. Its four sharp tenth-order constants, the resulting `q=2`
absorbed-body mass floors, the mixed-radius `q=4` hybrid floor, and the
component-address decoder below are **PROVED**.
The four displayed structured bodies and their twenty ray/scale checks are
**FINITE-EXACT** method hostiles. The restriction to the five primitive rays
and the inequalities `t L_H<beta_(p,q)` are **INHERITED**, specifically the
odd-3-unit specialization of THM-4153 (with THM-2060/2061 supplying the
dyadic seam and arithmetic sieves); they are not a new closure here. No
entire ray is closed, no universal body census is claimed, and LRC(14)
remains **OPEN**.

This is the detailed result note for THM-4450.

## 1. Inheritance pass and live objects

Write

```text
D_s={y in R/Z: ||sy||<1/14},
G_A=(R/Z) minus union_(a in A) D_a.
```

The two structured absorbed bodies are

```text
q=4: H=2C union {r},
q=2 with one even original tail: H=C union {r}.
```

The closest proved mechanism is THM-2061's equivalence between failure of
`2H union {a,b}` and containment of the complete closed set `G_H` in the
open two-tail cross-comb. THM-2060 adds primitivity, while THM-2061 adds
divisor-completeness through 14, owner-parity, and the determinant/width
box. The corrected near miss is that a covered phase or one covered
component is not a row counterexample. The least-used sidecar is the
address of every closed safe component, including isolated equality points.

The live board was: Haar mass; reduced ratio; common scale; component
address; endpoint openness; absorbed-label residue class modulo 6; and
structured fibre provenance.

## 2. A sharp four-tier overlap order statistic

Let `c,r` be distinct positive integers and reduce

```text
c/r=p/q,  gcd(p,q)=1.
```

With `B_2(u)=u^2-u+1/6`, normalized Haar measure satisfies

```text
mu(D_c intersect D_r)
 =1/49+[B_2({(q-p)/14})-B_2({(q+p)/14})]/(pq).          (1)
```

### Recovered proof of (1)

For `A(y)=1_(||y||<1/14)`, the Fourier coefficients are

```text
Ahat(0)=1/7,
Ahat(k)=sin(pi k/7)/(pi k),  k nonzero.
```

Orthogonality in `integral A(py)A(qy)dy` leaves frequencies `(qk,-pk)`.
Product-to-sum and

```text
sum_(k nonzero) cos(2 pi k u)/k^2=2 pi^2 B_2({u})
```

give (1). Common dilation disappears because multiplication on the circle
preserves Haar measure.

Since the oscillation of `B_2` is `1/4`,

```text
mu(D_c intersect D_r)>=1/49-1/(4pq).                    (2)
```

The orientation matters even though the measure is symmetric: because
`q|r`, the reduced denominator can contain only those factors of 2 and 3
already present in `r`. Exact enumeration below the product cutoff in (2)
gives the complete class-uniform table. `Below` and `equal` list the oriented
ratios `c/r`; no height cutoff remains. Here `lambda_g` is the best constant
uniform over all `r` with `gcd(r,6)=g`; a particular fixed `r` may force a
larger constant because it has fewer available denominators.

| `gcd(r,6)` | allowed reduced `q` | `lambda_g` | cutoff `pq` | ratios below `lambda_g` | ratios equal to `lambda_g` |
|---:|:---|---:|---:|:---|:---|
| `1` | `gcd(q,6)=1` | `1/63` | `55` | `1/13,1/11,2/11,3/11,10,11,12,13` | `9/5,9,27` |
| `2` | `3` does not divide `q` | `1/70` | `40` | `1/13,1/11,2/11,3/11,11/2,11,12,13` | `1/10,3/10,10` |
| `3` | `q` odd | `1/70` | `40` | `1/13,1/11,2/11,3/11,11/3,11,12,13` | `10/3,10` |
| `6` | unrestricted | `1/77` | `33` | `1/13,1/12,12,13` | `1/11,2/11,3/11,11/3,11/2,11` |

For example, above the first cutoff,
`1/49-1/(4*56)>1/63`; analogous first excluded products are `41` and
`34`. Thus the finite lists are exhaustive, not bounded evidence.

If `C` consists of ten distinct labels and `r` is not in `C`, the ten ratios
`c/r` are distinct. Each row of the table has fewer than ten subcritical
ratios. Therefore

> **Sharp tenth-overlap theorem.** For every ten-set `C` disjoint from
> `{r}`,
>
> ```text
> max_(c in C) mu(D_c intersect D_r) >= lambda_r,        (3)
> ```
>
> where the class-uniform `lambda_g` is determined by `g=gcd(r,6)` as in
> the table.

The constants in (3) are sharp already with primitive `C`. The following
controls have `|C|=10`, `gcd(C)=1`, `r notin C`, and maximum overlap exactly
`lambda_g`.

| `gcd(r,6)` | `r` | sharp `C` |
|---:|---:|:---|
| `1` | `715` | `55,9295,8580,65,130,195,7865,7150,1287,6435` |
| `2` | `1430` | `110,130,260,390,7865,15730,17160,18590,143,429` |
| `3` | `429` | `33,39,78,117,1573,4719,5148,5577,1430,4290` |
| `6` | `1716` | `132,143,20592,22308,156,312,468,6292,9438,18876` |

Sharpness here is deliberately typed: it is sharp for the maximum of the
ten pair correlations. It does not assert sharpness of the union overlap or
of the structured safe-mass floor derived next.

## 3. New `q=2` absorbed-body floors and entry thresholds

For `H=C union {r}` one has the exact identity

```text
mu(G_H)
 =mu(G_C)-1/7+mu(D_r intersect union_(c in C)D_c).       (4)
```

Combining (3) and (4) gives

> **Absorbed-label floor.** If `|C|=10` and `r notin C`, then
>
> ```text
> gcd(r,6)=1:  mu(G_H)>=mu(G_C)-8/63;
> gcd(r,6)=2 or 3: mu(G_H)>=mu(G_C)-9/70;
> gcd(r,6)=6:  mu(G_H)>=mu(G_C)-10/77.                  (5)
> ```

This improves the previous uniform subtraction of `1/7` in every residue
class. It covers all possible absorbed labels, not only odd 3-units.

THM-2061 gives the general odd-pair obstruction cap `4/63`. For a pair of
odd 3-unit tails, the inherited sharper cap is `4/77`. Compact containment
in a proper open cross-comb excludes equality. Hence (5) yields these
sufficient entry levels on the original ten-body:

| `gcd(r,6)` | loss in (5) | general pair (`4/63`) | odd-3-unit pair (`4/77`) | five-ray localization (`4/91`) |
|---:|---:|---:|---:|---:|
| `1` | `8/63` | `4/21` | `124/693` | `20/117` |
| `2` or `3` | `9/70` | `121/630` | `139/770` | `157/910` |
| `6` | `10/77` | `134/693` | `2/11` | `174/1001` |

For the actual odd-3-unit absorbed label, the threshold `124/693` improves
the earlier `15/77`, and five-ray localization begins at `20/117` rather
than `17/91`. Once localized, the ratio-specific `gcd(r,6)=1` pair-entry
levels are

```text
(1,11): 124/693;       (1,23): 256/1449;
(5,11): 86/495;        (1,37): 404/2331;
(1,25): 272/1575.                                      (6)
```

These are sufficient mass gates. No universal assertion that an arbitrary
ten-body enters one of them is made.

## 4. A `q=4` mixed-radius refinement, but no residual closure

For `H=2C union {r}`, disintegrate over `y=2x`. For `y in G_C`, both lifts
are safe for `2C`; since `r` is odd, at most one is `r`-dangerous. Put

```text
E_r={y:||ry||<1/7}.
```

The exact fibre identity is

```text
mu(G_(2C union {r}))
 =mu(G_C)-(1/2)mu(G_C intersect E_r).                   (7)
```

It immediately recovers the inherited pointwise floor `mu(G_H)>=mu(G_C)/2`.
There is also a mixed-radius analogue of (1): if `c/r=p/q` is reduced,

```text
mu(D_c intersect E_r)
 =2/49+[B_2({(q-2p)/14})-B_2({(q+2p)/14})]/(pq).        (8)
```

For odd 3-unit `r`, the global tenth order statistic of (8) is `1/28`.
There are only seven oriented ratios below it,

```text
1/25,1/13,1/11,2/11,5,6,13,
```

and equality occurs at `4/5,4,12,20`. It is sharp at

```text
r=3575,
C=(21450,325,17875,650,275,46475,143,2860,14300,42900).
```

Equations (7)--(8) therefore give the proved hybrid floor

```text
mu(G_(2C union {r}))
 >=max(mu(G_C)/2, mu(G_C)-1/8).                         (9)
```

In the live residual range `mu(G_C)<1/4`, the first term already dominates,
so (9) does not improve the `q=4` entry level `8/91` or close a ray. It is a
strictly stronger global structured-body statement, not an LRC closure.

The inherited ratio-sensitive `q=4` pair-entry levels remain

```text
(1,11):8/77; (1,23):16/161; (5,11):36/385;
(1,37):24/259; (1,25):16/175.                           (10)
```

## 5. Inherited five-ray localization, explicitly separated

This section records an inheritance, not a new theorem. All cells and widths
in this section use quotient-`y` coordinates. Specialize
THM-4153 to odd 3-unit tails and suppose

```text
mu(G_H)>=4/91.
```

If `2H union {a,b}` failed, write `a=tp,b=tq` with primitive `p<q`. Then

```text
(p,q) in {(1,11),(1,23),(5,11),(1,37),(1,25)},          (11)
t is an odd 3-unit, and t L_H<beta_(p,q),                (12)
```

where `L_H` is the longest positive-length component of `G_H`. The exact
primitive open components, useful for the new address decoder, are:

```text
1:11  (6/77,8/77), (69/77,71/77)
1:23  (6/161,8/161), (20/161,22/161),
      (139/161,141/161), (153/161,155/161)
5:11  (6/35,15/77), (62/77,29/35)
1:37  (6/259,8/259), (20/259,22/259), (34/259,36/259),
      (223/259,225/259), (237/259,239/259),
      (251/259,253/259)
1:25  (6/175,8/175), (4/35,22/175),
      (153/175,31/35), (167/175,169/175).               (13)
```

Their largest quotient-`y` widths are, in the same order,

```text
2/77, 2/161, 9/385, 2/259, 2/175.                       (14)
```

Physical-`x` pullback cells are twice as numerous and half as wide.
Equations (11)--(14), including the five ratios and the beta restriction,
are the inherited THM-4153 specialization confirmed by the prior THM-4449
candidate audit. The new work begins with the exact address vector below.

## 6. Exact finite component-address reduction

Let the complete closed components of `G_H`, including zero-length isolated
safe points, be `J_1,...,J_m`. Let `(A_k,B_k)` range over the primitive open
cells in (13). At a retained scale `t`, define the pulled cells

```text
U_(n,k,t)=((n+A_k)/t,(n+B_k)/t),  0<=n<t.               (15)
```

> **Exact decoder.** Strict failure at ratio `(p,q)` and scale `t` is
> equivalent to: for every `i`, there exists one address `(n,k)` such that
>
> ```text
> J_i strictly contained in U_(n,k,t).                  (16)
> ```

Indeed, the pulled cells are precisely the components of
`m_t^(-1)(Sigma_(p,q))`, and a connected closed component contained in their
open union lies in exactly one cell. This also handles isolated equality
points correctly. Omitting them can turn a boundary escape into a false
containment.

A cheap necessary prefilter is that every endpoint of every pulled cell in
(15) must be `H`-unsafe. A safe endpoint lies outside the open cross-comb and
immediately witnesses noncontainment. The finite hostiles in Section 8 show
that the converse is false, and that inspecting one arbitrary component is
also false. The whole address vector in (16) is the smallest faithful state
exhibited here.

The structured body itself need not be recomputed from scratch:

```text
q=2: G_(C union {r})=G_C intersect G_{r};
q=4: take both half-scale preimages of every component of G_C,
     then intersect with G_{r}.                          (17)
```

These interval identities are exact for closed weak-safe sets.

## 7. Explicit finite scale caps

LRC for eleven speeds gives a phase with body clearance at least `1/12`.
If `R_H=max H`, the minimum-clearance function is `R_H`-Lipschitz, so
`G_H` contains an interval of radius `1/(84R_H)` and hence

```text
L_H>=1/(42R_H).                                         (18)
```

Combining (12), (14), and (18) gives the coarse all-body caps

| ratio | necessary scale inequality |
|---:|---:|
| `1:11` | `t<12R_H/11` |
| `1:23` | `t<12R_H/23` |
| `5:11` | `t<54R_H/55` |
| `1:37` | `t<12R_H/37` |
| `1:25` | `t<12R_H/25` |

The `5:11` line is sharper than the generic THM-2061 tail-height box.
Together with `t` odd and `3` not dividing `t`, this is an explicit finite
scale bank for every fixed body. It is not a uniform height-free body
enumeration.

## 8. Finite-exact structured hostiles

Every body below has eleven distinct speeds, gcd one, is divisor-complete
through 14, has exact structured provenance, and lies in the residual band

```text
4/91 <= mu(G_H) < 4/77.
```

For every displayed body and every one of the five rays, the bank from
(12) is the singleton `{1}`. In the two right columns, entries follow the
ray order `1:11,1:23,5:11,1:37,1:25`.

| body | exact provenance | `mu(G_H)` | components; `L_H` | trapped components | safe pulled endpoints |
|:---|:---|---:|:---|:---|:---|
| `q4_A` | `H=(2,12,14,16,18,20,22,25,26,34,38)=2(1,6,7,8,9,10,11,13,17,19) union {25}` | `1151237/25865840` | `44;17/3192` | `2,2,4,0,2` | `0,0,0,0,0` |
| `q4_B` | `H=(2,6,8,14,20,22,23,26,32,34,36)=2(1,3,4,7,10,11,13,16,17,18) union {23}` | `174683549/3945221280` | `46;1/180` | `4,0,4,2,2` | `2,0,0,0,0` |
| `q2_A` | `H=(2,6,10,14,16,17,18,22,23,26,60)=C union {17}` | `12707741/246576330` | `36;23/3808` | `2,2,0,2,2` | `0,0,0,0,0` |
| `q2_B` | `H=(2,3,10,12,13,14,16,17,18,19,22)=C union {13}` | `10204829/203693490` | `20;1/154` | `0,4,2,0,2` | `0,2,0,2,0` |

Across the two bodies of each provenance class, every ray has a trapped
component in a case where all pulled endpoints are unsafe. Nevertheless,
every one of the twenty checks also has escaping components and an exact
physical full-row witness with clearance at least `1/14`; the audit records
the witness phase and lift. Thus these are hostile to endpoint-only and
single-component arguments, not counterexamples to LRC(14).

## 9. Precise finite next obligation

No complete ratio ray closed in this session. After the inherited
mass/ratio/scale filters, the exact remaining obligation for each proposed
body is:

1. Enforce THM-2060/2061 primitivity, divisor-completeness, owner parity,
   and determinant/width constraints.
2. Construct every closed component of `G_H`, retaining isolated points.
3. For each of the five ratios, enumerate the odd 3-unit scales satisfying
   `tL_H<beta_(p,q)` (or the coarse caps in Section 7).
4. Reject a scale immediately if any pulled endpoint is body-safe.
5. Otherwise compute the address of every component by (15)--(16). One
   escaping component closes that row; only a full address vector remains
   hostile.

This is finite for each body and uses no arbitrary tail-height scan. A
uniform bound on the bodies or a theorem forcing an escape on one whole ray
is still missing.

The information map is:

```text
source: complete closed structured safe set G_H
target: five primitive ratio rays with finite common-scale banks
map: Haar mass -> ratio filter -> width filter -> component address vector
preserved: strict/open type, ratio, scale, all component endpoints
destroyed by scalar stages: circular addresses, isolated points, owners
sidecars restored: endpoint safety, full addresses, THM-2061 parity/determinant
cheapest decisive test: address every component at every retained scale. (19)
```

## 10. Exact replay

The self-contained audit uses `Fraction` arithmetic and explicit `need`
checks, so optimized Python cannot erase its gates. It verifies (1) and (8)
against literal wall integration on `3,309` coprime pairs, reconstructs the
four sharp order-statistic controls, proves the finite product atlases,
reconstructs all four structured safe sets, and checks all twenty finite
ray banks.

```text
python -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450.py
python -O -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450.py
```

Normal and optimized outputs byte-match after newline normalization and end
with

```text
status=PASS;literal_overlap_checks=3309;finite_bank_checks=20;LRC14=OPEN
```

```text
04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450.py
  28e95f61464adcdf127d5b1d3e3d56949f871e644835e677c56b6326ab9e58fa
05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450.out
  a0555be8fe4db9122cc2267f4a14b2de9191eee00ccdd4fa49b32f29b3e995db
```

The clean-room decoder referee and orthogonal arithmetic referee replay as

```text
python -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.py
python -O -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.py
python -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent.py
python -O -B 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent.py
```

```text
04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.py
  832d3012823d374046964affd3028275bd125e259e4f8183ed769255856079e3
05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.out
  40bd7e00e62e0e2a87c591116e217d95a8e2087373202494f1bc6a4fc1743f98
04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent.py
  174c5f2d0db6552595949283e9817c92aea9dd04a2c5517c3611e54f7839b940
05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent.out
  ff175e08f3f0c2b4c230bb202c0fab4d53f6e5300bb5fb9afcfabd92c1df0b36
```

Canonical dependencies:

- `01-canon/theorems/THM-739-pairwise-coprime-bad-overlap-exact-bernoulli-closed-form.md`
- `01-canon/theorems/LEM-042-pair-overlap-law.md`
- `01-canon/theorems/THM-2060-crt-tail-coset-saturation.md`
- `01-canon/theorems/THM-2061-lrc14-dyadic-two-tail-folded-seam.md`
- `01-canon/theorems/THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer.md`
- `01-canon/theorems/THM-4449-lrc14-dyadic-seventh-rounding-energy-and-residual-haar-entry.md`
