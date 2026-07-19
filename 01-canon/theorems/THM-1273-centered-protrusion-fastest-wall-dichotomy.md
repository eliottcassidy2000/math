---
id: THM-1273
title: Centered-protrusion fastest-wall dichotomy
status: PROVED.  In the THM-1267 small-tail branch, four-prefix mass forces a lower-safe/fastest-active point in K while the protruding five-prefix survivor supplies a fastest-safe point in the endpoint tail.  The first intervening fastest wall is either a bare j=4 flood-tail wall outside K or is crossed by a lower tooth and carries an exact 1/(14*lcm) seam quantum.  Paper topology, optimization-safe exact referee, and sorry-free Lean arithmetic/propositional consumer are supplied.  This is a positioned bridge, not a new additive invoice or LRC(14) closure
source: codex-2026-07-19 protrusion/flood bridge audit
depends_on: [THM-1198, THM-1267]
related: [THM-1196, THM-1253, THM-1266, THM-1272, THM-1275]
script: 04-computation/lrc14_centered_protrusion_fastest_wall_thm1273.py
output: 05-knowledge/results/lrc14_centered_protrusion_fastest_wall_thm1273.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCenteredProtrusionFastestWall.lean
script_sha256: fa3618a4b06cbfe673d4461a738094b1442d7d8674be755702f34ce2f9efb50e
output_sha256: 926f86e969f2826c3e4d857b4efa4d52cb422c4cb67798ff9739bd2fa604af3f
formalization_sha256: 8d9718fa443f187596c43e610c8f499937ca74f5995ccf11b3925f6cb11c3380
---

# THM-1273 -- centered-protrusion fastest-wall dichotomy

## 1. Statement

Use the notation and hypotheses of THM-1267.  Thus

```text
G=[(14k+1)/(14c),(14k+13)/(14c)],
c<d1<d2<d3<d4<d5<d6,                                 (1)
```

the six strict danger combs cover `G`, `S` is the complete closed `d1`-safe
component selected by the centered spoke, and

```text
K=S intersect G,                  E=S minus G.        (2)
```

The set `E` is empty or one endpoint interval.  Normalize `S` to `[0,1]`,
write `ell=|E|/|S|`, and let `f` be THM-1198's six-bin probability density

```text
f=(3/4,13/12,7/6,7/6,13/12,3/4).                    (3)
```

Put

```text
h=d6,                   B={d2,d3,d4,d5},             (4)
V=S minus union_(j in B) Dbar_j,
U=S minus union_(j in B union {h}) Dbar_j.            (5)
```

Reflection handles the two tail orientations, so orient `S` with `K` first
and `E` last.  Then one of the following holds.

1. **Large protrusion:** `ell>=1/6`.
2. **Fastest-wall bridge:** `ell<1/6`, and there are regular points

   ```text
   x in int(K) intersect V intersect D_h,
   y in int(E) intersect U,                            (6)
   ```

   with `x<y`.  Let `z` be the first wall of the `h`-tooth containing `x`
   in the direction of `y`.  Then exactly the following useful dichotomy
   holds:

   * **bare `j=4` flood wall:** every `j in B` is non-dangerous at `z`.
     In that case `z` is outside `K`, hence `z in E`; or
   * **lower-crossed wall:** some `j in B` is strictly dangerous at `z`.
     Its tooth crosses the `h` wall and the adjacent positive intersection
     has length

     ```text
     omega_z >=gcd(h,j)/(14hj)=1/[14 lcm(h,j)].       (7)
     ```

At a bare wall, `d1` is safe because `z in S`, the four lower owners are
safe by definition, and `h` is exactly at distance `1/14`.  Thus all six
fast combs are non-dangerous there.  It is a literal terminal `j=4` flood
obligation, not merely a scalar survivor count.

## 2. The two mass obligations

Every normalized load belonging to a speed strictly above `d1` is less than
`7/36`.  Deleting only the four lower owners in (4) therefore gives

```text
integral_V f >1-4(7/36)=2/9.                          (8)
```

In the small-tail branch, `E` lies wholly in one endpoint sixth, where
`f=3/4` almost everywhere.  Hence

```text
integral_E f=(3/4)ell<1/8.                            (9)
```

Since `S=K union E` up to one null interface point, (8)--(9) imply the new
positioned interior obligation

```text
integral_(V intersect K) f>2/9-1/8=7/72.             (10)
```

Choose `x` from this positive-measure set away from the finite collection of
tooth and interval boundaries.  The owner `d1` is non-dangerous throughout
`S`, and all four owners in `B` are closed-safe at `x`.  Because `x in G`
and the six strict combs cover `G`, the fastest owner must be strictly
dangerous at `x`.  This proves the first half of (6).

THM-1267's separated-load refinement independently gives

```text
integral_U f>11/360,                  U subset E.     (11)
```

Choose `y` from `U` away from the same finite boundary set.  It is strictly
safe for all five owners `d2,...,d6`, proving the second half of (6).

This is the key change of carrier.  THM-1267 supplied one endpoint interval;
(10)--(11) put two oppositely typed proof obligations on that oriented
needle: a fastest-active four-prefix survivor on the `K` side and a
fastest-safe five-prefix survivor on the `E` side.

## 3. The first fastest wall

The point `x` lies inside one open `h`-tooth, while `y` is outside the closed
`h`-comb.  The oriented endpoint of the tooth containing `x` therefore lies
strictly between them.  Call it `z`.  It has the exact form

```text
z=(14n+1)/(14h)       or       z=(14n-1)/(14h),       (12)
```

and `||hz||=1/14`, so the strict `h`-comb does not cover `z`.

If a lower owner `j in B` is dangerous at `z`, openness puts an interval on
both sides of `z` inside one `j`-tooth.  On the `h`-active side this interval
meets the adjacent `h`-tooth positively.  The intersection either uses the
whole `h`-tooth or has one `h` endpoint and one `j` endpoint.  In the latter
case, clearing `14hj` gives a positive integer numerator divisible by
`gcd(h,j)`; in the former case its length is `1/(7h)`, which is larger still.
This proves (7).

Suppose instead no lower owner is dangerous at `z`.  If `z` belonged to
`K`, then `z in G`, while `d1`, all four lower owners, and `h` would all be
non-dangerous.  This contradicts the strict cover of `G`.  Therefore
`z notin K`.  Since the complete segment from `x` to `y` lies in `S` and
`S=K union E`, one has `z in E`.  This proves the bare branch.

Closed combs in (5) are used only for the mass survivor.  The regular points
`x,y` are selected off their finitely many boundaries.  At `z` the cover
test is strict: equality for `h`, or for any lower owner, is safe.  Thus no
null-set convention is smuggled into the wall dichotomy.

## 4. What is and is not paid

The crossed branch produces a genuine located overlap quantum.  It is not
automatically a *new* quantum.  If `z in K`, the same overlap may already be
represented among THM-1253's selected chronological seams or THM-1275's
flood/turn invoice.  Adding it again without a selection/disjointness lemma
would double count.  The bare branch is a wall event of measure zero and
also carries no scalar credit by itself.

Accordingly THM-1273 closes a missing **transport statement**, not the final
budget.  It proves that a centered protrusion cannot remain invisible to the
fastest chronology: it lands either on a literal four-prefix flood exit or
on an actual lower-fast overlap.  The remaining quantitative lemma must do
one of the following.

* show that the crossed wall is a skipped-tooth or multi-low occurrence not
  already selected by THM-1275;
* transport a bare wall to a complementary-coverer/Fano obligation while
  retaining its phase and wall address; or
* select several such centered needles with disjoint wall-side intervals.

This precise no-double-count boundary is why another unlocated pair sum or
runner-colour Fano count does not yet close the branch.

## 5. The exact `c=140` positive control

For THM-1266's sharp cell

```text
c=140, d1=254, (d2,d3,d4,d5,h)=(255,256,257,292,1805),
G=[1121/1960,1133/1960],                              (13)
```

the centered `d1` component is

```text
S=[2045/3556,2057/3556],              ell=33/280<1/6. (14)
```

The four-prefix survivor in `S` has exactly five components; adding the
full `1805` comb leaves exactly three.  Take

```text
x=7476011/12938240 in int(K),
y=7425603/12837160 in int(E).                         (15)
```

The point `x` is closed-safe for `255,256,257,292` and strictly dangerous
for `1805`; `y` is closed-safe for all five.  The three intervening fastest
walls are

```text
14603/25270,        2923/5054,        14617/25270.   (16)
```

The first two are crossed by the tooth `256@148`.  Their actual seam lengths
are

```text
213/6469120,                 65/1293824,              (17)
```

respectively `213` and `325` times the primitive
`1/[14(1805)(256)]` quantum.  The last wall has no lower crosser and lies
strictly in `E`.  At that point `1805` has depth exactly `1/14`, while all of
`254,255,256,257,292` have depth strictly larger than `1/14`.

Thus the sharp star realizes **both** outputs of the bridge on one oriented
needle: two crossed walls followed by a bare terminal wall.  It remains a
positive control, not a six-cover.  In particular the bridge does not falsely
reject the row merely because its local five-rung ladder is sharp.

## 6. Kakeya, `j=4`, Fano, and tournament carriers

The normalized component `S` is the Kakeya needle.  The correct vertices for
this theorem are not runners but the ordered pair of survivor obligations
and the intervening fastest wall events.  A wall vertex carries the active
lower-owner subset

```text
A(z)={j in B: ||jz||<1/14}.                           (18)
```

The bare branch is `A(z)=empty`; the crossed branch is `A(z) nonempty` plus
the exact endpoint numerator.  This preserves the predicate proved here.

For the sharp row, ordering its three wall vertices gives the transitive
tournament with score histogram `(0,1,2)`, no directed cycles, three
singleton SCCs, and one Hamiltonian path.  Its active-set labels are

```text
({256},{256},empty).                                  (19)
```

The ordinary tournament on `d2,...,d6` is also transitive, with score
histogram `(0,1,2,3,4)`, no cycles, five singleton SCCs, and one Hamiltonian
path.  It loses (18)--(19), all three wall positions, and which seam numerator
is present.

We challenged runners, selected teeth, fastest-safe gaps, section boundaries,
wall-crossing events, residues, Fano lines, and proof obligations as vertices.
The faithful carrier is

```text
(oriented S; V-components; inside h-active obligation; outside h-safe
 obligation; ordered h walls; active-low subset and endpoint digit per wall). (20)
```

For the `j=4` thread, a bare vertex is exactly the terminal endpoint of the
fastest flood of a four-prefix survivor.  For the Fano/`chi_7` thread, a bare
wall should be coloured by the *complementary obligation which covers that
located phase*, not by a preassigned runner colour.  This agrees with
THM-1260's negative local-surjectivity audit: incidence among several located
events, not one runner-colour fork, is the missing hypergraphic state.

## 7. Verification and scope

The dependency-free exact referee has zero Python `assert` nodes.  It checks
the full rational mass ledger, `198` strict boundary-approach rows, and
`40,504` signed-address lower-crossed wall configurations.  Every cleared
overlap numerator is integral and gcd-divisible, every seam pays (7), and the
minimum quantum ratio is exactly one.  It independently reconstructs the
four- and five-prefix survivor components, oriented points, three wall events,
active-low subsets, bare wall, and two seam lengths in (13)--(19).  Normal and
optimized Python outputs are byte-identical.

The sorry-free Lean module checks the mass constants and strict implications,
the bare-versus-crosser propositional cover logic, the inherited natural
gcd/lcm quantum, and every displayed rational comparison in the sharp row.
The analytic load bounds, positive-measure point selection, interval
orientation, continuity/first-wall selection, and endpoint ownership are the
explicit paper providers.  No theorem-specific axiom or `native_decide` is
used.

THM-1273 proves a phase-located bridge from centered protrusion to fastest
wall geometry.  It does not prove that the crossed seam is new, that the bare
wall is covered by a prescribed complementary owner, uniform six-comb
noncoverage, the empty sporadic branch, or LRC(14).  ∎

Frozen artifact hashes are

```text
source         fa3618a4b06cbfe673d4461a738094b1442d7d8674be755702f34ce2f9efb50e
output         926f86e969f2826c3e4d857b4efa4d52cb422c4cb67798ff9739bd2fa604af3f
formalization  8d9718fa443f187596c43e610c8f499937ca74f5995ccf11b3925f6cb11c3380
```
