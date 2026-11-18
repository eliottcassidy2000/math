# Dihedral height selection, with an actual rotation group

A targeted change from [sevenfold norm readings](NORM_SEVEN.md): this time the upper witness has a literal seven-rotation and a reflection, and the requested integer chooses the width. A second construction gives exact paired readings on the same arbitrary input. Neither is a universal prime selector.

## A marked dihedral construction

Use the five-vertex quotient with outgoing masks (20,13,8,17,6). The involution (1 2)(3 4), fixing 0, reverses every arc. Put the carousel C on Z/7 with steps {1,2,3} at vertex 0. Put U and its converse at 1 and 2, retaining singleton vertices 3,4. The rotation r acts by x↦x+1 on C and fixes the outside. The signed reflection s acts by x↦−x on C, exchanges corresponding vertices of U and its converse, and swaps the singletons. It reverses all arcs. Thus on Hamilton paths it must be followed by reversal of the path order. The two maps satisfy srs=r^(-1) and give a genuine dihedral action on paths.

Write h=H(U) and v=f_U(2)−2h, the number of its ordered nonempty two-path covers. Exact compressed-word enumeration gives
```
H(upper)/7 = 4573 h² + 7392 h v + 2919 v².
```
The marked quotient on 0,1,2 is transitive and there are just two ordinary singleton positions: cutting a word at them bounds the visits to each mark by three. The enumerator certifies that only one or two runs at either U-mark are actually possible. Weight the carousel runs by its ordered nonempty cover counts divided by 7 (free rotation on the decorated cover). The resulting coefficients for the two U-marks at (1,1),(1,2),(2,1),(2,2) are (4573,3696,3696,2919). Using U and its converse gives the same h,v by path reversal.

For transitive U of width k≥1 we have h=1,v=2^k−2. Given an absolute target p, form Δ=7392²+4·2919(p−4573). Test whether it is a nonnegative square, whether (sqrt(Δ)−7392)/(2·2919) is a nonnegative integer v, and whether v+2 is a power of two. This selects k directly, without a forward scan or a residue window. The polynomial is increasing on the physical cone at fixed h, so the width, when it exists, is unique. For a general block the norm inverse can instead try h≤sqrt(p/4573) and recover an exact v by the same quadratic inversion with the corresponding h-dependent coefficients. Physical realization of (h,v) is then still needed.

| Prescribed quotient p | Width | Literal upper count | Upper order |
|---:|---:|---:|---:|
| 31,033 | 2 | 217,231 | 13 |
| 47,278,009 | 7 | 330,946,063 | 23 |

Both quotients are prime by exact trial division. For the first, a separate lower witness comes from the signed wrapper with D=(80,41,67,5,14,29,58), contact mask 29, and a width-five transitive module. Its readings give 978·32−263=31,033. The literal subset counter verifies both graphs. The ordinary prime is not assumed to inherit the rotation: in fact the full upper automorphism group is C7 in these examples. It contains C7, its order divides 7p by the free action on paths, and an automorphism of order p cannot act on fewer than p vertices. The signed group, which also contains s, is therefore exactly D7 (order 14).

## Where the reflection sits in the absolute count

First divide the path set by r. Reflection fixes a rotation orbit iff that orbit has a representative fixed by s. This representative is unique: if s(P)=r^a(P) under the signed path action, replacing P by r^j(P) reduces the equation to 2j=a on Z/7. Consequently
```
p = F + 2M,
```
where F is the number of paths literally fixed by s and M the number of reflected pairs of seven-orbits. This uses the inherited signed-action mechanism in [THM-4747](../../../01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md), but the following exact size-dependent fold is specific to the new witness.

A fixed path is determined by its first half: choose one from each mirrored vertex pair, make a directed path through those choices, and end with an arc into the center 0. In this example there are four fixed-size choices of pairs (three from the carousel, one from the ordinary singletons), plus k matched U-pairs. Compress visits to the two U copies in the first half. Over all choices from the four fixed-size pairs the certified numbers of templates with 0,1,2 total U-runs are 9,111,93, with none higher. To fill u runs assign the k distinct matched pairs to u labeled nonempty runs; any such assignment dictates the unique increasing or decreasing order in that copy. Hence for k≥1,
```
F = 111 + 93(2^k−2) = 93·2^k−66.
```
For p=31,033 this is F=297 and M=15,368. At p=47,278,009 it is F=11,829 and M=23,633,090. The reflected half is linear in the binary width; the unpaired full paths supply a genuinely quadratic scale. Knowing just the reflection-fixed paths would not place the absolute target.

## A same-input paired affine reader using a genuine seven-rotation

A different five-vertex quotient (16,1,3,7,14), with carousel at 0 and a single arbitrary block U at 1, has divided run counts 1863 for one U-run and 1036 for two. Thus H/7=1036 B−209 h. The lower signed wrapper D=(190,124,120,240,96,64,129,54), contact mask 19, has (ordinary, first-ended, filtered) readings (155,1036,672), hence the very same affine read on arbitrary U. This matches the known slope-259 carousel read with two levels of transitive scaling, but supplies the physical lower partner without dividing the carousel's nonrealizable profile. The rotation is literal in this upper branch; the reflection of the previous section is not asserted here.

If a base input profile (h,b) is chosen, the requested p determines the padding from (p+209h)/(1036b), which must be an exact positive power of two. For p=265,007, use a singleton and seven padding singletons; for p=13,258,501, use the strong count-11 core (16,17,3,7,12) with b=100 and seven padding singletons. These yield 1,855,049 and 92,809,507 respectively on the upper axis. The prime and sevenfold constructions work on the same source.

**What transfers from the group theory:** retain both the reflection-fixed and paired terms, then make their common physical width the output of an exact inverse. The remaining Chan–Pak step is to choose a suitable quotient and actual input profile for each target not on these routes. A dihedral decomposition alone is not that supply theorem.

Replay: `python3 dihedral_replay.py` alongside the earlier packet scripts (stdlib).
