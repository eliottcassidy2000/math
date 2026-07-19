---
id: HYP-7678
title: The seventh-deletion Kakeya needle carries a Fano/chi7 overlap debt
status: OPEN uniform target, materially advanced by THM-1156, THM-1166, and THM-1176.  THM-1156 proves the chi7 seam bipartition/private-or-triple alternative.  THM-1166 proves global uncovered mass >=1/12, closes the common-dilate seventh crown G/m<=77/12, and gives exact adaptive forest and Fano triple-gcd necessary inequalities.  THM-1176 proves the complementary slow-gap pressure dichotomy and a finite-depth toothpick ladder.  No arbitrary-packet positive F7 constant or LRC(14) closure is known
source: codex-2026-07-18-S75
depends_on: [THM-1153, THM-1025, THM-1149, THM-1156, THM-1166, THM-1176]
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

## Landed seventh-rung arithmetic

THM-1156 proves that exact two-tooth seams are bipartitioned by the quadratic
character of the 7-primitive speed.  A zero seam inside the closed covered
needle therefore needs third support.  This is qualitative incidence, not a
uniform overlap amount.

THM-1166 supplies the first quantitative consequences beyond the zero
coefficient.  If `rho_ij` is the global pair mass, then every three-speed
packet has pair sum at least `1/24`, hence

```text
R=sum_(i<j)rho_ij>=7/24.
```

The optimal quadratic chamber certificate

```text
Q(C)=C-(2/7)binom(C,2)=C(8-C)/7
```

gives global uncovered mass at least `2R/7>=1/12`.  Consequently, if all
seven deleted speeds share a divisor `G`, periodicity and the protected-needle
length force

```text
G/max(W)<=77/12.                                        (A7)
```

For every labelled Fano plane, with `G_ell` the gcd of the three speeds on a
line, the exact linewise discrepancy decomposition gives

```text
sum_(ell) max(W)/G_ell>=32/231,                         (MF7)
min_(ell) G_ell/max(W)<=1617/32.
```

And every forest `F` obeys the adaptive necessary inequality

```text
sum_({i,j} in F)
 [|I|rho_ij-rho_ij(1-rho_ij)/gcd(si,sj)]
 <=(6/49)sum_i 1/si.                                    (TF7)
```

These close the common-dilate branch and give the Fano probe a numerical
consumer.  The still-open step is to force a contradiction from `(MF7)` or
`(TF7)` for arbitrary mixed-period packets.

THM-1176 attacks the same wall from the slowest-comb gaps rather than global
pair mass.  If the slowest deleted speed is `a`, the other six are
`b_1<...<b_6`, and `m=max(W)`, every counterexample partition obeys

```text
a<13m  or  a sum_i 1/b_i>1,                            (SG7)
```

with `b_1<=6a-3` in the second branch.  The phase-saturation equality is
impossible by an exact parity argument.  At every cardinality, `r` faster
combs covering a full `c`-slow gap force

```text
c sum_i 1/d_i>7-r.                                     (SG-r)
```

Consequently three or fewer faster combs cannot cover such a gap, and
recursive gap nesting gives

```text
b_1/a<13/6  or  b_2/b_1<13/6  or  b_3/b_2<4/3.        (TP7)
```

This is the toothpick self-similarity that the density ledger erased.  It
localizes the unresolved Fano/chi7 work to the harmonic-crowded, mixed-phase
clock region of `(SG7)`; the Fano lines should be placed over slow-gap endpoint
events, not over a fixed runner tournament.

The sharpened target is an arithmetic overlap-discrepancy theorem on this
specific needle.  One useful sufficient form is:

> There are absolute `eta>0` and `C>=0` such that, for every such
> protected interval and every seven distinct integer combs whose union
> covers `I`, some spanning
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
Indeed, writing `H=sum_(s in S)1/s`, the two inequalities give

```text
eta |I| <= (C+1/7)H,
H >= eta/((7C+1) max(W)),
min(S)/max(W) < 7(7C+1)/eta.
```

The last inequality uses distinctness of the seven deleted speeds.

The Fano/`chi_7` language is not decorative. When endpoint excess is small,
bulk duty one forces an almost partition, while THM-1149's identity says private mass is paid by
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
