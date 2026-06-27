---
id: HYP-2934
title: LRC14 summand/multiplicand Farey bridge
status: PROOF-INTERFACE / branch separator; not a proof of LRC14
source: codex-2026-06-23-S133
related:
  - HYP-2932
  - HYP-2933
  - HYP-2931
  - HYP-2930
  - HYP-2908
  - HYP-2161
  - HYP-2083
  - HYP-1822
  - HYP-1821
---

# HYP-2934: `2/27` is the C=27 summand-unit branch, `3/41` is the first K33 multiplicand blow-up

S133 reconnects the prompt's Farey-product graph observation with the older
summand/product graph work.  For an LRC14 exact optimum `M(S)=p/q`, there are
three adjacent graph layers that must not be collapsed:

```text
Farey additive pair:     (p,q) is a summand pair at node p+q.
LRC shell denominator:   q may itself be the summand-shell node C=2n-1.
Multiplicand product:    (p,q) is a factor fiber at product node p*q,
                         and K_{p,q} is the incidence blow-up of that fiber.
```

The guardrail is that the LRC summand graph of HYP-2083/HYP-2161 is the
denominator shell `q=27`, not the formal Farey node `p+q`.

## Computation

The script
`04-computation/lrc14_summand_multiplicand_bridge_codex_s133.py` stores output
at
`05-knowledge/results/lrc14_summand_multiplicand_bridge_codex_s133.out`.

It verifies the LRC14 unit-excess chain:

```text
1/13: p+q=14, p*q=13,  K_{1,13}, q-threshold parent.
2/27: p+q=29, p*q=54,  K_{2,27}, C=27 second-gap summand-unit shell.
3/41: p+q=44, p*q=123, K_{3,41}, first K_{3,3} incidence wall.
```

For `n=14`, `C=2n-1=27`, and the denominator-shell summand pairs split as:

```text
unit-visible:    9 shells
nonunit gcd=3:   3 shells
nonunit gcd=9:   1 shell
```

Thus `2/27` is not merely "planar" in the S131 product ledger.  It is exactly
the second-gap value `2/(2n-1)`, living in the `C=27` summand graph acted on by
units.  The two S131 unit-excess petal rows

```text
swap 10->20
swap 13->26
```

both have this bridge tag.

By contrast, `3/41` leaves the `C=27` shell and is the first unit-excess
complete-bipartite incidence blow-up containing a `K_{3,3}` minor.  It is the
known near-miss row `12->36`.

## Product Graph Versus `K_{p,q}`

S133 also separates the ordinary multiplicand graph from the `K_{p,q}`
certificate:

```text
2/27: product node 54 = 2 * 3^3
      factor fibers: (2,27), (3,18), (6,9)
      K_{2,27}: two-block, no K33.

3/41: product node 123 = 3 * 41
      factor fibers: (3,41)
      K_{3,41}: K33-wall.
```

So the product graph and the complete-bipartite certificate are not identical
objects.  The product graph records all factor fibers of `p*q`; `K_{p,q}` is the
incidence expansion of the particular factor fiber `(p,q)`.  This is why
`2/27` can have rich `3`-adic multiplicand fibers but no `K_{3,3}` wall, while
`3/41` has a simple product node and still becomes the first three-owner
incidence obstruction.

## Proof Consequence

Refine HYP-2932's branch split as:

```text
e=0:        AP/GW floor candidates.
e=1,p=1:    coarse q-threshold parent 1/13.
e=1,p=2:    C=27 second-gap summand-unit/petal branch.
e=1,p>=3:   K_{3,3} multiplicand-incidence branch.
```

The `p=2` branch should be attacked with the HYP-2083/HYP-2161 toolkit:
petal/two-block rigidity, unit/nonunit `C=27` shell strata, lift denominators,
CRT conservativity, and owner labels.

The `p>=3` branch should be attacked as a finite three-owner packet: the
`K_{3,3}` incidence wall is not a contradiction by itself, but it marks the
first place a HYP-2908 tournament-state lift or forbidden-H endpoint could have
enough independent owners to be constructed.

## Tournament Analysis

S133 challenges the default vertex choice.  It considers runners, gaps,
denominator shells, Farey pairs, factor fibers, `K_{p,q}` incidence blow-ups,
`K_{3,3}` packets, and proof obligations.

The chosen tournament vertices are proof carriers:

```text
q-binding,
C27-shell,
Farey p+q,
product fiber,
Kpq blowup,
K33 packet,
raw apex iso.
```

The pair observable is a role score:

```text
theorem safety,
additive shell retention,
multiplicand fiber retention,
K33 signal,
LRC specificity,
false-positive resistance.
```

The resulting carrier tournament is transitive:

```text
q-binding > C27-shell > Farey p+q > product fiber
          > Kpq blowup > K33 packet > raw apex iso.
```

A product-only gauge would flip `9` carrier pairs, which is the warning: product
data is useful only after the denominator shell and owner/incidence labels have
been retained.

## Proof Target

The next useful theorem is:

```text
Every remaining q=14 non-AP/GW atom either reduces to the C=27 p=2
summand-unit/petal branch, or crosses the p>=3 K33 packet wall.
```

This does not finish LRC14, but it narrows the handoff.  The `p=2` side should
be discharged by the existing second-gap shell machinery.  The `p>=3` side
should feed a finite three-owner obstruction packet, rather than being treated
as another scalar Farey denominator.
