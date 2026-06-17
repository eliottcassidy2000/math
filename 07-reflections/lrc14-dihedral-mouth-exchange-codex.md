# LRC14 dihedral mouths: the drop-6 core is a reflected hexagon shadow

**Source:** codex-2026-06-17-S2.  Prompt: spend a long session seeing the deep connection
between dihedral groups and LRC and use it to make proof progress.

This note targets OPEN-Q-108, the uniform fattening lemma for 12-speed cores

`G_C={t in [0,1): ||v t|| > 1/14 for every v in C}`.

## Dihedral endpoint lemma

For one speed `v`, the danger set

`D_v={t: ||v t|| <= 1/14}`

is the pullback of the base moat under `t -> v*t`.  Its endpoints are

`(14*k - 1)/(14*v)` and `(14*k + 1)/(14*v)` modulo `1`.

These endpoints are acted on by the dihedral clock of the regular `v`-gon:
rotations send `k -> k+r`, and reflection sends `k -> -k` while swapping left
and right endpoints.  For a full core `C`, the shared rotation usually vanishes,
but the global reflection `t -> 1-t` always survives.  Thus safe components of
`G_C` come in reflected orbits.

If a safe component starts at the right endpoint of speed `a`, center `k/a`, and
ends at the left endpoint of speed `b`, center `l/b`, then its exact length is

`(a*(14*l - 1) - b*(14*k + 1))/(14*a*b)`.

I am calling the numerator the **dihedral mouth determinant**.  It is the integer
gap between two adjacent endpoint events after they are placed on their common
`14ab` clock.

## Drop-6 decomposition

The previous gauntlet found the current worst 12-core

`C6={1,2,3,4,5,7,8,9,10,11,12,13}`

with

`meas(G_C6)=7/858`.

The new endpoint atlas shows what that number is made of.  The safe set has only
four components, paired by reflection into two dihedral mouth orbits:

| orbit | endpoint mouth | component length | reflected total |
|---|---:|---:|---:|
| outer | `13R2 -> 12L2` and mirror | `1/728` | `1/364` |
| inner | `12R2 -> 11L2` and mirror | `5/1848` | `5/924` |

So

`meas(G_C6)=2*(1/728)+2*(5/1848)=7/858`.

Even more importantly, both reflected pairs lie inside the omitted speed-6 danger
moat, specifically the two opposite teeth centered at `1/6` and `5/6`.  Thus the
worst known core is not merely "AP13 minus 6"; it is a **D6 missing-clock
shadow**, with the other AP speeds covering all but two reflected pieces of that
hexagon moat.

## Perturbation trade

The script `04-computation/lrc14_dihedral_endpoint_atlas_codex.py` scans the same
two-AP-deletion plus one-replacement family through `w<=180`, but now compares
each row to the old drop-6 mouths.

Among the `2004` rows whose missing pair includes `6`:

- `1104` rows damage at least part of the old drop-6 mouths.
- The best total row is still delete `(6,10)`, add `20`, with
  `meas(G_C)=3859/420420 = 7/858 + 1/980`.
- That best row damages none of the old hexagon mouths; it just adds a reflected
  decagon mouth pair inside the missing speed-10 moat.
- The best row that actually damages the old drop-6 mouths is delete `(6,10)`,
  add `40`, with `meas(G_C)=2669/194040`, already much larger.
- The replacement `69` can kill most of the old hexagon mouths, but the new
  mouth mass it forces elsewhere is huge; the largest-loss row in the scan has
  measure `5/77`.

This suggests the proof target:

> **Dihedral mouth-exchange inequality.**  Starting from the AP drop-6 shadow,
> any coordinated replacement that covers old hexagon mouth mass must create at
> least as much reflected mouth mass in another missing-clock shadow, with a
> positive surplus.  In the tested family the smallest surplus is `1/980`, and it
> occurs without covering the old mouths.

If this exchange inequality can be proved in a scale-invariant form, it would
upgrade HYP-2568's computational extremal statement into a real attack on
OPEN-Q-108.

## Tournament Analysis

The script uses reflected safe-component orbits, not runners, as the main
Tournament Analysis vertices.

- Pairwise observable: total measure contributed by a reflected mouth orbit.
- Switch/gauge: orient the larger mouth contribution toward the smaller one,
  with left-endpoint order as the Hamiltonian tie path.
- Fingerprints: drop-6 has two mouth-orbit vertices, one Hamiltonian path, zero
  directed cycles, and zero edge flips between measure ranking and raw determinant
  ranking.
- The AP deletion-compression tournament over omitted speeds is transitive, with
  order `6,12,10,4,2,13,5,3,8,9,11,1,7`; delete `6` is the unique champion.

This is deliberately a different quotient from the speed-load tournament.  It
preserves exact safe measure and endpoint adjacency, but it destroys most runner
history away from the two boundary owners.  That is the right trade for this
proof step: OPEN-Q-108 needs a positive safe interval, and safe intervals are
born at adjacent endpoint events.

## Honest status

This is proof progress, not a proof.  It rigorously explains the current
extremal example and gives a new inequality to aim at, but it does not control
arbitrary coordinated growth of three or more large speeds.  The next useful
steps are:

1. Prove the endpoint-mouth formula as a reusable lemma in the canon.
2. Prove the AP missing-clock shadow lemma: `G_{AP13\{m}}` is contained in the
   omitted `D_m` moat and decomposes into reflected `D_m` teeth.
3. Upgrade the two-delete/one-replacement trade into an infinite inequality:
   replacement speeds that cover old mouths force compensating new mouths.
4. Test whether non-AP primitive 12-cores can be normalized into a mouth-exchange
   sequence without decreasing `meas(G_C)`.

Cross-links: OPEN-Q-108, HYP-2568, HYP-2569, THM-523, THM-522, T836, T837.
