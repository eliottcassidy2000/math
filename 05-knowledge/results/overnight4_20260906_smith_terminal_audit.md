# Independent audit: terminal clusters determine the largest two-jet factor

**Verdict: PASS.** The all-node analytical argument and the exact one-digit
complete-residue exception in incoming THM-4439 are sound. This audit does
not claim a new theorem or restore metric determination of intermediate
Smith factors. The bounded exact controls also pass independently.

Sources read via `git show origin/main` during the parent-reported
`72a7db79b` incoming update:

- `01-canon/theorems/THM-4439-all-node-twojet-metric-precision-by-terminal-clusters.md`;
- its full proof `05-knowledge/results/hermite-terminal-cluster-precision-overnight-hexagon-sep05.md`;
- THM-4435 and its underlying `four-node-metric-and-hermite-precision-overnight-hexagon-sep05.md`;
- our `overnight3_20260906_smith_triple_single.md` ternary metric hostile.

## Exact target and inverse audit

The observer includes values and first Hasse derivatives at every node
of a finite set X of distinct integers, on the fixed polynomial module
of degree below2|X|. This completeness and degree box are essential.
For a prime p, put S_x=sum_(y!=x)v_p(x-y), f_x=max_(y!=x)v_p(x-y), and
q_x=sum_(y!=x)1/(x-y). The audited conclusion is

    L_p=max_(terminal C) [2S_C+max(0,f_C-[|C|=p])],

with L_p=0 for one node. The complete Hermite inverse supplies

    L_p=max_x {2S_x,2S_x-v_p(2q_x)}.

I checked the underlying denominator statement rather than accepting it
as a headline: the derivative cardinal polynomial has leading coefficient
1/F'(x)^2, and the value cardinal polynomial has leading coefficient
-F''(x)/F'(x)^3. All other coefficients are integral combinations of these
two types. Both are attained, with the F''=0 case correctly retained.
The least denominator of the square inverse is the largest Smith factor.
First Hasse derivative equals ordinary first derivative, while F'' is the
ordinary second derivative; losing its factor2 would break p=2.

## Why the maximum descends to terminal clusters

At a nonterminal nearest ball, the current leaf is a singleton child but
another child has at least two points. Moving to a point in that deeper
child increases S by at least1 and f by at least1. Distances to outside
points are unchanged. Thus T=2S+f increases by at least3. Repeating ends
at a terminal cluster. Its best reciprocal contribution is at least T-1,
so no nonterminal contribution can win the global maximum. The separate
2S term also increases under this move. This proves the descent for both
terms, including shallow f=0 and uneven branches.

## Why the complete-residue exception loses exactly one digit

In a terminal cluster C at depth f, normalized node residues u_x are
distinct modulo p. Put G(U)=product_(x in C)(U-u_x) and Q_x=p^f q_x.
The internal reciprocal part is G''(u_x)/(2G'(u_x)); every G'(u_x) is a
unit. Outside terms belong to p Z_p, and their difference at two u values
is divisible by p^2: the numerator contains p^(2f), whereas the two
outside denominators together have valuation at most2f-2. Therefore
the outside contribution divided by p is constant modulo p.

For odd p and m=|C|<p, G'' has degree m-2 and unit leading coefficient
m(m-1). It cannot vanish on all m distinct residues, so min v_p(Q_x)=0.
For m=p, G mod p is U^p-U. Hence G'' is coefficientwise divisible by p,
and G''/p has degree p-2 with unit leading coefficient p-1. As G'=-1
modulo p, Q_x/p reduces to a nonconstant polynomial of degree p-2 plus
the outside constant. It cannot vanish on all p residues. Consequently
min v_p(Q_x)=1, exactly, regardless of higher unit digits or nearby
outside clusters. This is a simultaneous statement, not a pointwise one.

At p=2 every terminal cluster has two points. Its single internal
reciprocal is a unit and outside terms are even, so min v_2(Q_x)=0.
Here the one-digit subtraction comes from the explicit factor2 in2q_x.
The distinct mechanisms agree in the displayed terminal formula.

Two useful hostiles were reconstructed literally:

- At p=3, X=(0,3,6), the normalized reciprocal valuations are
  (1,infinity,1). Local cancellation can be infinite, although the
  cluster loses exactly one digit. Its Smith list is(0,0,1,3,4,4), so
  omitting the complete-residue exception would wrongly predict5.
- With the nearby outsider1, X=(0,1,3,6), the same terminal cluster has
  normalized reciprocal valuations(2,1,1), yet still largest exponent4.
  The outside constant can deepen cancellation at a selected point but
  cannot eliminate every unit after first division.

The partial p=3 cluster X=(0,3) has largest exponent3, so subtracting a
digit for every terminal cluster would also be false. At depth0 the
baseline2S and max(0,...) correctly keep all exponents nonnegative.

## Exact connection to ternary intermediate-factor blindness

For both X=(0,9,27,81) and X'=(0,9,54,81), the only terminal cluster is
the closest pair{0,81}, with f=4 and S=9. Its size2 is less than p=3,
so THM-4439 gives L_3=2*9+4=22 immediately. Literal integer Smith forms
independently give

    X : (0,0,2,6,7,12,15,22),
    X': (0,0,2,6,7,13,14,22).

The common tree therefore determines the largest factor and the common
determinant valuation64, while the normalized unit residue changes an
intermediate determinantal ideal. Modulo3^13 their kernel orders are
3^53 and3^54. A largest-factor theorem cannot erase that observer
distinction.

Source: the full integral observer. Target: its sharp uniform coefficient
precision loss. Map: form exact inverse columns, maximize their denominator
loss over a complete terminal cluster, then retain the valuation tree.
Preserved predicate: the largest Smith exponent and its uniform precision
threshold. Destroyed information: intermediate Smith factors and finite
precision kernel sizes below that threshold. Needed sidecar for those
stronger targets: the unit-sensitive intermediate-minor residue chi from
our ternary family. The cheapest decisive test is the displayed pair of
literal8-by-8 observers, rather than another largest-factor comparison.

## Independent finite controls

The companion checker imports no repository producer. It builds literal
integer matrices and full Smith forms, computes their rational inverses,
evaluates reciprocal sums, and enumerates terminal balls by all pair
depths. Those four paths agree on24 named observers at p=2,3,5,7 with
one through seven nodes. The bank includes full/partial residue clusters,
depth0/1/2/3, higher unit lifts, outside points at depth f-1, unequal
branches, signed translations, vanishing q, and both dyadic and ternary
intermediate-factor twins. It also checks the confluent determinant
valuation4 sum_(i<j)v_p(x_i-x_j).

Reproduction artifacts are the
[independent source](../../04-computation/overnight4_20260906_smith_terminal_audit.py),
[normal output](overnight4_20260906_smith_terminal_audit.out), and
[optimized output](overnight4_20260906_smith_terminal_audit_optimized.out).
Reproduce from the repository root:

    python 04-computation/overnight4_20260906_smith_terminal_audit.py
    python -O 04-computation/overnight4_20260906_smith_terminal_audit.py

Both pass221 explicit gates. All arithmetic is integer/rational; no
floating-point approximation or sampled inference enters the proof audit.
General higher Hasse jets, incomplete observations, and moving modules
remain outside this theorem's scope.

Source SHA-256:

    b0fc41164d1ccf0b21a7fcbf0f55f029d4516b453ed86e41a7a77fc63699d00d

Normal and optimized outputs are byte-identical after explicit LF
normalization, with SHA-256:

    2f47860b7d6664ae1a86d242b858ab706ec360d7649963757b48c1ea724e4e12
