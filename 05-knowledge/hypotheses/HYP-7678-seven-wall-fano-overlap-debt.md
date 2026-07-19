---
id: HYP-7678
title: The seventh-deletion Kakeya needle carries a Fano/chi7 overlap debt
status: OPEN sharpened target; THM-1153 proves the first-order coefficient is exactly zero and identifies the missing observable, but supplies no positive overlap constant
source: codex-2026-07-18-S75
depends_on: [THM-1153, THM-1025, THM-1149]
---

# HYP-7678 -- the seven-wall overlap debt

Let an actual LRC(14) counterexample be split into a six-speed core `W` and
seven deleted speeds `S`.  Lower-case fattening gives a protected interval

```text
I subset Safe_(1/14)(W),   |I|>=1/(7 max W).
```

The seven danger combs of `S` must cover `I`.  THM-1153 proves that the
first-order fragmentation coefficient is exactly zero here: seven bulk
densities of `1/7` can pay for all of `I`, leaving only endpoint excess
`(1/7) sum_(s in S)1/s`.

The sharpened target is an arithmetic overlap-discrepancy theorem on this
specific needle.  One useful sufficient form is:

> There are absolute `eta>0` and `C<infinity` such that, for every such
> protected interval and every seven distinct integer combs, some spanning
> tree `T` on the seven combs satisfies
>
> ```text
> sum_({a,b} in T) measure(I intersect D_a intersect D_b)
>   >= eta |I| - C sum_(s in S)1/s.                  (F7)
> ```

Hunter's tree inequality and coverage give the reverse bound

```text
sum_({a,b} in T) measure(I intersect D_a intersect D_b)
  <= (1/7) sum_(s in S)1/s.
```

Thus `(F7)` would restore a positive harmonic crown at `r=7`, bound
`min(S)/max(W)`, and extend THM-1153's projective compression below `v7`.

The Fano/`chi_7` language is not decorative.  Bulk duty one forces an almost
partition, while THM-1149's identity says private mass is paid by
multiplicity at least three.  Pairwise runner tournaments can record a tree
but cannot tell whether its overlaps coincide in the same triple chambers.
The likely faithful object has:

```text
vertices = tooth or wall-crossing events on I,
hyperedges = simultaneous owner sets,
weights = overlap length plus endpoint excess,
sidecar = comb modulus and tooth address.
```

Challenge recorded: vertices need not be runners, and the quotient must
preserve coverage of the protected needle.  A transitive runner tournament
destroys exactly the Fano incidence the target needs.

No claim is made that `(F7)` is true with arbitrary `eta,C` as stated.  The
next computation should measure the best scale-normalized pair/tree/triple
deficits on the already-certified hard seven-comb packets and look for a
phase-labelled counterexample before promoting constants.
