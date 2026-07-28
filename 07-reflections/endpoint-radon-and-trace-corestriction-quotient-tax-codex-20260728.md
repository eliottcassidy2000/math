# Endpoint Radon and trace corestriction pay the same quotient tax

**Status:** CURRENT SYNTHESIS / RESEARCH GUIDE, not a theorem.  Truth sources
are corrected THM-2621, audited THM-2620/2625, and MISTAKE-301.

Two apparently unrelated frontier moves now have the same exact shape.

```text
LRC current:
  ordered endpoints (L,R)
       -> target difference q=L-R
       loses common endpoint origin and Delta=det(L,R).

Planar inverse map:
  branchwise polynomial action potentials H_i
       -> field trace Phi=sum_i H_i
       loses every trace-zero / anti-invariant branch channel.
```

The quotients are lawful.  THM-2334's Radon marginal is exactly recovered by
summing THM-2625's determinant sectors, and THM-2621's trace differential is
exactly `d Tr(H)`.  The problem is not an incorrect quotient.  It is asking a
consumer-sensitive question after applying a quotient whose kernel the
consumer does not kill.

## The live concept board

| Object | Representation | Preserved predicate | Lost coordinate | Cheapest next test |
|---|---|---|---|---|
| THM-2334 current | `J(L,R)` on `F_13^2 x F_13^2` | exact target marginal | endpoint origin, determinant sector | preserve one nonzero sector through a positive adjacent-clock kernel |
| THM-2620 grammar | parabolic graph at fixed `(q,Delta)` | ordered projective transition | coefficient sign, chronology, ancestry | attach one physical root gate without Radon collapse |
| Degree-four inverse | marked quartic `(f,b)` | Keller PDE and sheet-defect poles | polynomial source primitive | impose the power-sum residue gate on each `D_4` boundary chart |
| `D_4` action | `G_tau=H-tau(H)` | deck anti-invariant channel | trace and matching-quadratic ownership | descend `-G_tau^2` and compare its divisor with the two Kummer parity vectors |
| Positive chronology | common-ancestry cylinders | literal nonempty intersections | principal action / equal masses | test sector coefficients rather than support after the later root gate |

## What the two repairs actually prove

On the LRC side, THM-2625 proves far more support than expected:

```text
all 28,561 endpoint cells survive;
all 2,185 admissible (q,Delta) sums survive;
all 2,016 nondegenerate parabolic sectors survive;
every such thirteen-edge fibre has all thirteen coefficients nonzero.
```

This removes “perhaps every transvection sector cancels” from the live
residual.  It does not produce positivity or time.

On the JC side, MISTAKE-301 moves in the opposite logical direction.  The
attractive traced residue obstruction is impossible for an actual polynomial
Keller map:

```text
theta=x dy-kappa^(-1)P dQ=dH,
omega=Tr(theta)=d Tr(H),
Res_D(omega)=0.
```

The useful data is therefore finer, not larger: each branch separately obeys

```text
Res(x dy)=sum_j j x_(-j)y_j=0,
```

and a `D_4` involution retains `H-tau(H)` even though trace kills it.

## The common mechanism

The successful move has five steps.

1. Name the exact quotient and prove it is lawful.
2. Compute its kernel as a coordinate, not merely as a dimension.
3. Split the object before quotienting.
4. Test the actual coefficient/branch array, including hostile signs.
5. Reapply the quotient and demand exact recovery of the canonical marginal.

THM-2625 passes step 5 through

```text
sum_Delta Sstar(q,Delta)=169 Anum(q).
```

Corrected THM-2621 passes it through

```text
Tr(dH)=d Tr(H).
```

This is stronger than the slogan “keep more information”: it identifies the
minimal split whose recombination is already a proved canonical object.

## Tournament perspective, with its boundary

THM-2620 has an intrinsic directed relation: an ordered endpoint pair at
fixed nondegenerate `(q,Delta)` is the graph of a parabolic transvection.
Here a tournament/projective-transition lens is content-bearing because row
order, endpoint reversal, determinant squareclass, and the lost common origin
are explicit.

The `D_4` branch ledger is not automatically a tournament.  A deck
involution pairs branches, but pole orders can tie and residue zero is forced;
orienting those ties would be cosmetic.  The faithful object is first a
signed/valued deck graph with vertex potentials `H_i` and edge coboundaries
`H_i-H_j`.  Only an intrinsic asymmetric boundary observable would justify an
orientation.

## Cross-frontier next experiment

The strongest honest common experiment is a **quotient-before/after square**.

```text
fine array --physical or polynomial constraint--> fine constrained array
   |                                               |
 quotient                                        quotient
   |                                               |
coarse object -------- known theorem ----------> coarse consequence
```

For LRC, the upper arrow must be positivity plus chronological ancestry; for
JC, it is polynomial exactness plus branchwise valuation.  A proposed bridge
is useful only if the square commutes and the upper arrow preserves a
specific fine coordinate.  This rejects two recurring false moves at once:
using a nonzero marginal as an endpoint edge, and using a zero trace residue
as branchwise exactness.

The new practical rule is: **when an exact aggregate looks too weak, do not
replace it; factor it through the smallest audited fine object and make the
next consumer act there.**
