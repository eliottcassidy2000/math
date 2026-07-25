---
id: THM-2182
title: "Endpoint-grid product and the joint tail-overlap sidecar"
status: >
  PROVED + VERIFIED-EXACT. If L is divisible by 14e for every speed in
  a finite core E, then every scaled tail factors exactly:
  measure(G_(E union LC))=measure(G_E)measure(G_C). Hence every aligned
  seven-core/six-tail row is strictly LRC(14)-safe, quantitatively by
  1/(49e_2). The tails {3,12} and {4,6} prove that phase-zero labels,
  zero one-comb endpoint currents, equal one-comb masses, and equal
  reciprocal sums still do not determine a two-comb continuation; the
  missing coordinate is their joint normalized overlap law.
source: codex-2026-07-24-knot-relations
depends_on: []
related:
  - THM-1091
  - THM-2094
  - THM-2162
  - THM-2174
script: 04-computation/lrc14_endpoint_grid_product_thm2182.py
output: 05-knowledge/results/lrc14_endpoint_grid_product_thm2182.out
script_sha256: df197a01ca4f882591348a9ec8011c696982cbfe3ee11569d562ae5cea274c41
output_sha256: b1edb2391936b6d8a55393a7e744e242a987c23d6dffd796dcb65e8facf999c3
hash_basis: working-tree bytes (LF)
---

# THM-2182 -- endpoint-grid products retain the whole tail

At radius `1/14`, put

```text
D_v={t in R/Z:||vt||<1/14},
G_E=intersection_(e in E) ((R/Z)\D_e).               (1)
```

All statements about measure are unchanged if either inequality is made
weak at the finitely many boundary points. For an integer `L` and a set
`C`, write

```text
LC={Lc:c in C}.                                      (2)
```

## 1. Exact endpoint-grid product

> **Endpoint-grid product theorem.** Let `E` and `C` be finite sets of
> positive integers. If
>
> ```text
> 14e divides L             for every e in E,         (3)
> ```
>
> then
>
> ```text
> measure(G_(E union LC))
>   =measure(G_E) measure(G_C).                       (4)
> ```

**Proof.** Every boundary point of `D_e` has the form

```text
(14q+/-1)/(14e).                                     (5)
```

Condition (3) puts all these points in `(1/L)Z/Z`. Therefore the indicator
of `G_E` is constant almost everywhere on every open grid cell

```text
I_j=(j/L,(j+1)/L),              j=0,...,L-1.         (6)
```

The complete tail indicator is

```text
B_C(t)=1_(G_C)(Lt)
      =product_(c in C) 1_((R/Z)\D_(Lc))(t).          (7)
```

It is `1/L`-periodic. On every cell, the substitution `x=Lt-j` gives

```text
integral_(I_j) B_C(t)dt=measure(G_C)/L.               (8)
```

If `N` grid cells belong to `G_E`, then

```text
measure(G_E)=N/L,                                    (9)
```

and summing (8) over exactly those cells proves

```text
integral 1_(G_E)(t)B_C(t)dt
  =N measure(G_C)/L
  =measure(G_E)measure(G_C).
```

The left side is (4). QED.

This is an exact whole-tail Haar packet, not an approximation and not a
termwise independence assertion. The multiplication map `t |-> Lt` forgets
the occupied cell but preserves the complete joint law of the normalized
tail `C`.

## 2. An infinite strict LRC(14) family

Let

```text
E={e_1<e_2<...<e_7},        C={c_1,...,c_6},         (10)
```

where both sets have distinct positive entries, and take any `L` satisfying
(3). Because `L>max E`, the row `E union LC` has thirteen distinct positive
relative speeds.

Each danger set has measure `1/7`. The circular interval

```text
(-1/(14e_2),1/(14e_2)) mod 1                        (11)
```

lies in `D_(e_1) intersection D_(e_2)`, so

```text
measure(D_(e_1) intersection D_(e_2))>=1/(7e_2).
                                                               (12)
```

First join those two danger sets, then apply the union bound to the other
five. This gives

```text
measure(G_E)
 =1-measure(union_(e in E)D_e)
 >=1/(7e_2).                                         (13)
```

The ordinary union bound for six danger sets gives

```text
measure(G_C)>=1-6/7=1/7.                             (14)
```

Equation (4) now yields the uniform strict certificate

```text
measure(G_(E union LC))>=1/(49e_2)>0.                (15)
```

Thus **every aligned `7+6` endpoint-grid split satisfies strict LRC(14)**.
This closes an infinite structured family, but does not assert that an
arbitrary thirteen-speed row admits such a split.

## 3. Why marginal endpoint data cannot replace the tail law

Let `K` be any finite positive-speed core with `measure(G_K)>0`, and let
`L` satisfy

```text
14k divides L             for every k in K.          (16)
```

Compare the two normalized tails

```text
C_1={3,12},                 C_2={4,6}.               (17)
```

Every actual tail speed is zero modulo the endpoint grid `L`. Applying
(4) to a singleton tail shows, for all four speeds,

```text
measure(G_(K union {Lc}))=(6/7)measure(G_K).          (18)
```

Equivalently, all four one-comb endpoint currents in the notation of
THM-2162 vanish. The two tails therefore have the same:

```text
endpoint phase word:                  (0,0);
two one-comb continuation masses:     (6mu/7,6mu/7);
two one-comb endpoint currents:       (0,0);
reciprocal-speed sum:                 5/(12L),        (19)
```

where `mu=measure(G_K)`.

Nevertheless their joint laws differ. Multiplication by an integer preserves
Haar measure, so

```text
measure(G_(C_1))=measure(G_{1,4}),
measure(G_(C_2))=measure(G_{2,3}).                    (20)
```

For `{1,4}`, the two danger sets overlap only in the interval about zero of
length `1/28`. For `{2,3}`, their only overlap has length `1/21`. Hence

```text
measure(G_(C_1))=1-2/7+1/28=3/4,
measure(G_(C_2))=1-2/7+1/21=16/21.                   (21)
```

By (4),

```text
measure(G_(K union LC_1))=(3/4)mu,
measure(G_(K union LC_2))=(16/21)mu,                 (22)
```

whose difference is

```text
mu/84>0.                                             (23)
```

Thus the tuple in (19), despite retaining phase, scale through the reciprocal
sum, and every one-tail continuation, is not continuation-congruent for a
two-tail insertion. The first missing coordinate is the ratio/gcd-sensitive
**joint normalized-tail overlap law**.

## 4. Exact thirteen-speed witness

Take

```text
K={1,2,3,4,5,6,8,9,10,11,12},
L=55440=lcm(14k:k in K).                             (24)
```

Two independent exact interval evaluators give

```text
measure(G_K)=883/6930.                               (25)
```

Consequently the two thirteen-speed rows in (22) have measures

```text
883/9240,
7064/72765,                                          (26)
```

respectively, with exact gap

```text
883/582120.                                          (27)
```

The companion script computes the core and normalized-tail measures by both
danger-interval union and complete boundary-cell enumeration, then checks
every fraction in (21), (25)--(27).

## 5. Transfer boundary

The typed transfer is:

```text
source:       THM-2162 endpoint denominator;
map:          t=(j+x)/L cell disintegration;
preserved:    the entire tail product;
destroyed by marginalization:
              tail ratio/gcd overlap;
needed sidecar:
              normalized tail C, or its full safe law.          (28)
```

THM-2094 also conditions on a complete residue law, and THM-1091 gives a
Fourier sheet identity. Neither result states (4), the strict aligned `7+6`
closure (15), or the marginal obstruction (17)--(23).

QED.
