---
id: HYP-2901
status: PROOF GUARDRAIL / structural synthesis
source: codex-2026-06-22-S114
tags: [lrc14, tournament-recursion, legendre, eisenstein, mobius, exact-period, equidistribution, denominator-wall]
depends_on:
  - HYP-2899
  - HYP-2900
  - THM-523
  - THM-527
  - THM-549
related:
  - HYP-+2876
  - HYP-2886
  - HYP-2890
  - HYP-2898
  - OPEN-Q-108
results:
  - 04-computation/lrc_lcm_committed_denominator_wall_codex_s114.py
  - 05-knowledge/results/lrc_lcm_committed_denominator_wall_codex_s114.out
  - 04-computation/lrc_venn_bonferroni_wide_kps.py
  - 05-knowledge/results/lrc_venn_bonferroni_wide_kps.out
  - 04-computation/lrc_rfar_convergent_tail_kps.py
  - 05-knowledge/results/lrc_rfar_convergent_tail_kps.out
---

# HYP-2901: three-mode composition splits the finite Venn node from the analytic denominator wall

The corrected recursion atlas is:

```text
Mobius mode:     A+B+C-D-E-F+G              all sizes
Eisenstein mode: A+B-C                      even sizes only
Legendre mode:   A+B+D-C-E-F+G              odd sizes only, as a 3-set Venn
```

The letters do **not** name the same subtournament sizes in the three modes.
That is the main correction to preserve.

## Correct Legendre geometry

For odd size `N`, the Legendre half-tiling is a full 3-set Venn over corners:

```text
corners:
  A : N-1
  D : N-2
  B : N-1

pairwise intersections:
  A cap B = C : N-2
  A cap D = E : N-3
  B cap D = F : N-3

triple:
  A cap B cap D = G : N-4
```

Thus the three edge unions are:

```text
A+B-C,  A+D-E,  B+D-F
```

and the center/full union is:

```text
A+B+D-C-E-F+G.
```

The old net-coefficient readout `(2,0,-2,1)` is still numerically correct, but
it hides geometry: `D` and `C` both have size `N-2`, and cancel only after
scalarization.  For LRC work they must remain distinct labels: one is a corner,
the other is an overlap.

## Composition at 14

For `14=2*7`, the even Eisenstein fold is the top operator:

```text
14 even -> A+B-C with A,B of size 13 and C of size 12.
```

This fold exposes the apex `7`.  The odd apex then belongs to the Legendre Venn,
where the corrected corner/edge/center labels carry the finite three-gap or
AP-hull structure.  Mobius is the always-on inclusion-exclusion skeleton and
also the divisor/coprime exact-period packet skeleton (`phi=mu*id`).

So the useful composition is not a single scalar recurrence:

```text
exact-period packets
  -> Mobius divisor/IE labels
  -> Eisenstein even fold 2q -> q
  -> Legendre odd Venn at q
  -> scalar cap/floor only after labels survive.
```

## LCM committed-speed wall

Let

```text
S_X = {1,2,...,11,13,lcm(2,...,X)}.
```

Then for every denominator `D<=X` and every numerator `a`, the committed speed
has phase

```text
lcm(2,...,X) * a / D == 0 mod 1.
```

Hence no fraction `a/D` can be a `1/14` lonely witness for `S_X`.  Therefore
the minimal witness denominator for this family tends to infinity with `X`.

This is a proof-level guardrail:

```text
No fixed finite denominator basis can prove LRC14 for all covering 13-sets.
```

The script confirms the wall and then probes the first denominator above it.
The stronger heuristic `firstD = nextprime(X)` is false:

```text
X=60  firstD=67,  nextprime=61
X=110 firstD=121=11^2, nextprime=113
X=120 firstD=121=11^2, nextprime=127
```

Across `X=14..140`, `72/127` tested rows disagree with `nextprime(X)`.
The right object is a prime-power / divisor-lattice opening plus residue
compatibility with the base row, not primality alone.

## Incoming S45 integration: radical filter before equidistribution

Mac-mini S45 independently refuted the finite-certificate reading of
HYP-+2876 with the same lcm family and added the useful small-prime filter.
If a prime `p<=13` divides no speed in `S`, then `t=1/p` is immediately safe:
every runner has nonzero residue mod `p`, hence distance at least `1/p`, and

```text
1/p >= 1/13 > 1/14.
```

Thus any hard covering row must be prime-covering for

```text
2,3,5,7,11,13
```

and must also kill the composite denominator `14`.  In one-committed-speed
families this says the radical, not the absolute size, is the first filter;
the radical-saturated case `30030 | v` is where the small-prime witness no
longer closes the row and equidistribution becomes necessary.

S114's correction to that filter is that the **minimal** witness denominator
above the wall is not governed by primes alone: `D=121=11^2` appears first for
`X=110,120`.  The effective Node-3 theorem should therefore be stated over
exact-period prime-power packets, with the radical filter as the first easy
branch and torus equidistribution as the saturated branch.

## Incoming S31t integration: the Venn gives a Bonferroni-3 wide target

KPS S31t applies the same corrected Venn geometry to the wide p0/cap branch.
For a bounded base `B` and far runners, the cover expansion is the Newton /
Mobius expansion in the far runners:

```text
p0(B union far) = T_1 + T_2 + T_3 + T_4 + ...
```

where the corrected Venn identifies:

```text
T_1  corners       one-far packet
T_2  edges         two-far / doublet packet
T_3  center        three-far packet
```

The verified pattern is:

```text
T_1 = 0 or small in genuine-wide rows,
T_2 is the binding doublet edge,
T_3 is the first subdominant correction,
T_4,T_5,... alternate with decreasing size in sampled rows.
```

If the uniform sign/decrease statement is proved, then the wide branch has the
Bonferroni upper bound

```text
p0(B union far) <= T_1 + T_2 + T_3,
```

so the corrected Venn reduces the unbounded multi-far cap problem to a
third-order target: doublet plus triple, with the `r>=4` tail nonpositive.  This
fits HYP-2901's split: the Venn supplies the finite packet names, while the
uniform decay/sign theorem is still analytic Dedekind/equidistribution content.

## Proof implication

This separates the two remaining LRC14 nodes cleanly.

**Finite node (Node 2):** prove the cap/extremality statement in the corrected
Legendre Venn geometry.  The live finite route is AP/three-gap rigidity:
consecutive/dilated rows are the only configurations whose sector-gap profile
stays low-complexity under all phases, and the Venn labels keep corner and
overlap contributions distinct until the cap comparison.

**Analytic node (Node 3 / finite Part A):** no finite denominator basis can
replace effective equidistribution.  The lcm family forces witnesses beyond
every bounded denominator wall.  The next theorem should be an effective
prime-power/residue equidistribution lemma: once the committed denominator wall
is passed, enough unit residue packets survive to intersect the `GOOD cap G_P`
floor, with the finite `#arcs/Vmax` loss from THM-527.

## Assumption challenge

The relevant tournament vertices here are not runners.  Useful vertices are:

```text
recursion modes, Venn regions, exact-period denominator packets,
prime-power openings, cap/floor proof obligations, and scalarization guards.
```

This quotient preserves which inclusion-exclusion packet can feed the LRC
floor/cap proof.  It destroys the actual circle-time geometry and must hand
that geometry back to the analytic equidistribution/Part-A node.
