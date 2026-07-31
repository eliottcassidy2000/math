# K3 primitive-unit phase -> carrier spanning-tree bridge

Status: **OPEN RESEARCH NOTE.**  This is a proposed connection and test plan,
not a theorem or dependency.

## Live sources

- THM-2984 retains the primitive high-ray unit `u` and the attained/unattained
  sign bit before taking a denominator maximum.
- The candidate reflected-ray branch
  `origin/codex/session-lrc-ray-continue3` retains carrier-addressed pair
  overlaps and uses the pointwise spanning-tree inequality

  ```text
  sum_i 1_(D_i) - sum_(ij in T) 1_(D_i intersect D_j) <= 1
  ```

  because the edges induced by any active vertex set form a forest.
- The analytic inequality itself is already sorry-free in
  `TournamentH7.LRCTreeHunter.tree_hunter_add_le`; a future package should
  reuse that Lean dependency and formalize only the new carrier/ray inputs.
- The sharp cyclic Cayley lemma in the pending THM-2979/2981 packages says
  cardinality-only torsion is already optimal.  A failure of that gate must
  use residue shape or location.
- THM-2984 now has a second sharp cardinality theorem for the *absolute*
  unit gate.  The strict bad band modulo `d` has exactly
  `beta(d)=2 floor((d-1)/14)+1` residues.  Thus `|S|>beta(d)` clears every
  primitive unit, optimally for a count-only argument.  Below that boundary,
  the relevant object is the multiplicative orbit of the actual shape `S`,
  not merely its size.

## Precise proposed connection

Source object:

```text
one projected k=3 terminal case
=(fixed first label, two finite low labels, high denominator d,
  primitive unit u, fixed body carrier C_E).
```

Target object:

```text
the complete graph on the four literal labels, weighted by their exact
pair-overlap masses on C_E (or on the fixed-low-safe subcarrier).
```

Map:

```text
(case,u) -> high ray z_m=z_0(u)+mL
         -> exact carrier pair weights omega_ij(u,m)
         -> maximum spanning-tree weight.
```

Preserved predicate: pointwise membership of each carrier point in every
strict-open danger comb, hence the forest inequality and an honest upper bound
on union mass.

Destroyed by the old denominator quotient: primitive direction, absolute
cell address, pair co-occurrence, and whether the high-ray zero is attained.

Needed sidecars: THM-2984 reachability trichotomy, exact residue-ray recurrence,
carrier component word, and exact pair-overlap integrals.  The torsion-density
certificate is an optional earlier exit, not a dependency of the tree gate.

## Why this is not a tournament

The intrinsic binary observable is symmetric overlap mass, not an oriented
comparison.  The correct object is a weighted complete graph and its graphic
matroid.  Orienting edges would add a gauge and lose the forest inequality's
meaning.

## Cheapest decisive test

1. Wait for the first terminal `(case,u)` on which THM-2984 finds no absolute
   fixed-safe cell.
2. Reconstruct the exact fixed-low-safe carrier as closed rational intervals.
3. For the four literal combs, compute all six exact pair-overlap masses on
   that carrier at the first scalar-reachable high point of the ray.
4. Enumerate the 16 labelled spanning trees of `K4`; test whether

   ```text
   sum singleton masses - max_T sum_(ij in T) omega_ij < carrier mass.
   ```

5. Repeat along the exact `A/z` ray.  Determine whether every edge weight is
   affine/rational in `1/z`, piecewise rational with finitely many endpoint
   events, or genuinely changes combinatorial type.

A positive result closes the ray without a safe point shared by the coarse
cell model.  A negative result must record the minimizing unit, high height,
active endpoint word, best tree, exact deficit, and first combinatorial event;
that witness then becomes the smallest target for a three-edge triangle or
matroid-rank refinement.

## Multiplicative transporter compression

Before enumerating carrier trees, package the absolute-cell failures as the
transporter

```text
T_d(S,B_d)={u in (Z/dZ)^*:uS subset B_d},
B_d={r:14 min(r,d-r)<d}.
```

For a terminal case the exact obstruction is not "some unit looks bad" but

```text
T_d(S,B_d) intersect U_active is nonempty,
```

where `U_active` is supplied by THM-2984's signed-ray attainment law.  This
description separates geometry (`T_d`) from the scalar wall (`U_active`) and
permits several exact early exits:

1. `|S|>beta(d)` makes the transporter empty sharply.
2. Units preserve every gcd stratum
   `O_g={r:gcd(r,d)=g}`.  If
   `|S intersect O_g|>|B_d intersect O_g|` for any `g|d`, the transporter is
   empty even when total cardinality is inconclusive.
3. If `S` contains a unit `s`, the image `us` determines `u`; only the unit
   residues of `B_d` need be tried.  More generally one can choose the
   residue whose stabilizer in `(Z/dZ)^*` is smallest and intersect the
   resulting candidate cosets recursively.
4. Multiplicative stabilizers compress duplicate work: when `|S|=|B_d|`, a
   successful inclusion is equality, and all transporters (if any) form a
   coset of `Stab(B_d)` after one witness is fixed.

The last stabilizer is now explicit.  If `d>=15`, then
`Stab(B_d)={+1,-1}`.  A stabilizer sends `1` into `B_d`, so its least absolute
representative is `a<=b=floor((d-1)/14)`.  For `a>=2`, the bad residue
`k=floor(b/a)+1` has `b<ak<=2b<d-b`, contradiction.  Hence an equality-sized
nonempty transporter has exactly the two candidates `+-u_0`.  For
`2<=d<=14`, `B_d={0}` and the stabilizer is the whole unit group.

There is also a direct additive-to-multiplicative bridge.  If `S` contains
two residues whose difference has additive order at most seven, then for
every unit their images have circular separation at least `d/7`; they cannot
both lie in the strict radius-`d/14` band.  Thus every located short-order
pair is already a universal unit-cell certificate.  The converse fails at
`d=8,S={1}`.  In particular, a nonempty transporter forces `S` to be an
independent set in THM-2979's Cayley graph `G(d,7)`.  This explains the exact
hierarchy rather than treating the two gates as unrelated:

```text
short-order difference -> transporter empty -> all active unit rays close,
but transporter empty need not supply any pair.
```

The cheap next computation should therefore report, for the first failed
count gate, the gcd-stratum vectors of `S` and `B_d`, transporter size,
stabilizer size, and its intersection with the positive/zero/negative active
unit classes.  This is the correct group-action sidecar; it avoids inventing
a tournament where the intrinsic object is a bipartite incidence relation.

## Recursive extension

The hierarchy is

```text
cardinality-forced low-order pair
  -> beta(d)-forced absolute unit escape
  -> gcd-stratum / transporter test T_d(S,B_d) intersect U_active
  -> primitive-unit carrier spanning tree
  -> unit-indexed graphic-matroid rank / selected higher intersections.
```

Each step retains exactly one coordinate discarded by the previous one and
has a finite hostile probe.  Do not jump to higher intersections before an
actual tree deficit is exhibited.
