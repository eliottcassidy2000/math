# Literal strong-suffix returns from 170807

Status: **FINITE-EXACT** named realizations, with independently checked native
operations and integer path counts. A reusable all-depth ternary grammar is
**OPEN**. No theorem IDs are claimed by this packet.

## Inheritance and exact interface

The source is the actual order-17 graph in
`../../triangle_return_170807_20261103/interface.json`, field `target`, with
physical names
`[1,2,4,6,7,8,14,21,22,31,32,33,39,40,50,51,54]`.
The inherited zero-ready certificate has root 39, sole inport 31, transitive
prefix `50 -> 51`, and a strong order-14 suffix B on
`[1,2,4,6,7,8,14,21,22,31,32,33,40,54]`.
Its actual endpoint data are

\[
h=H(B)=2195,\qquad F=E_{31}(B,2)=42153,\qquad h+4F=170807.
\]

The closest proved mechanism is the inherited zero-ready identity
`H(T_s)=h+2^s F`, together with the isolated-edge compiler from the inherited
`realize.py`. The hostile control is the unchanged marked suffix: changing a
scalar clock does not give a return rule. The corrected near miss is equating
one numerical hit with a uniformly repeatable affine recurrence. The least-used
sidecar here is the *choice of actual inport inside the strong component*.
The concept board is: actual named suffix; movable inport; transitive-prefix
clock; native path-complement guard; full port relation; repeated grammar entry.

The source/target map flips specified physical pairs inside B, transfers root
39's sole inport, and appends prefix vertices when required. It preserves every
old named vertex and re-enters the same zero-ready structural class. It does
not preserve the old root star or its phase gauge. Those data are explicitly
recomputed in the independent `../history_phase/` audit.

## Both missing next branches occur

Every pair in the following table means reversing its current direction.
Changing inport 31 to w additionally reverses `{39,31}` and `{39,w}`.

| Record in `realized_depth3.json` | B pairs reversed | New inport | Prefix size | h | F | Actual H |
|---|---|---:|---:|---:|---:|---:|
| 0 | `{4,31}`, `{6,22}` | 4 | 4 | 5655 | 31673 | 512423 |
| 1 | `{1,7}`, `{6,31}`, `{6,54}` | 1 | 2 | 39331 | 118273 | 512423 |
| 4 | `{1,31}`, `{14,33}`, `{21,22}` | 40 | 3 | 46073 | 58294 | 512425 |

Thus the literal graph with H=170807 has native realizations of both
`3H+2` and `3H+4`. Record 1 keeps order 17. Record 4 adds only permanent prefix
vertex 55, giving order 18. Record 0 uses fewer native complements but adds
prefix vertices 55 and 56, giving order 19.

Each old edge reversal uses two actual Hamiltonian-path complements with one
fresh temporary helper, then deletes that helper alone. Record 0 uses eight
complements; records 1 and 4 use ten. The full named paths, intervening rows,
birth cuts and deletions are recorded in the JSON. Every completed isolated
repair is strong. Common local frames at ports 50 and 51 are maintained in
each complement pair. Root 39, spare 6, old inport 31 and all other old names
survive, but their relations need not remain unchanged.

`verify.py` checks independently that the target root has the advertised sole
inport, that the prefix really dominates the indicated strong suffix, and that
`E_39(T,2)=H(T)`. The return is therefore an entry into the structural grammar,
not merely a count realized by an unrelated graph.

## Exact finite menu and its limitation

The discovery menu changes zero, one, two or three of the 91 B edges, rejects
non-strong B, considers every one of the 14 physical inports, and changes the
transitive prefix clock. The original recorded range is prefix size 2 through
5. At depth at most two, 4187 suffixes were considered, 3815 were strong, and
only record 0 hit either target. At depth at most three, 125672 suffixes were
considered, 113361 were strong, and six hits were found: five +2 and one +4.
The other +2 rows remain in the exact JSON, without quotienting by isomorphism.

Removing the prefix cutoff gives **eight** hits: the same five +2 witnesses and
three +4 witnesses. The additional +4 rows use prefix size 6, with `(h,F)` equal
to `(6377,7907)` at inport 40 or `(2985,7960)` at inport 21. Both have actual
order 21 and are compiled in
`realized_suffix_depth3_moving1_all_prefixes.json`. Its record indices differ
from the earlier bounded file: the selected bounded indices 0, 1, 4 become
full-menu indices 0, 3, 6. Both files retain their original stable order.

For the six hits, allowing the same source prefix s and the recorded relative
prefix shift k gives `D(s)=H'(s+k)-3H(s)` as follows:

| Records | k | D(s) |
|---|---:|---|
| 0, 2, 3 | 2 | `-930 + 233*2^s` |
| 1 | 0 | `32746 - 8186*2^s` |
| 4 | 1 | `39488 - 9871*2^s` |
| 5 | 2 | `21486 - 5371*2^s` |

Each slope is nonzero. The indicated offset is attained only at s=2 for this
fixed transformation and shift. The appearances of 233 and of the new inports
are exact arithmetic facts; they do not supply a recurrence by themselves.

## A complete second generation

The two selected targets also each admit both next branches, by the same exact
menu with at most three B-edge reversals, all 14 inports, and every prefix size
at least the actual inherited size. No old prefix vertex is removed.

| Source | Chosen record | Target | New inport | Prefix size | h | F | Order |
|---:|---|---:|---:|---:|---:|---:|---:|
| 512423 | seed1, index3 | 1537271 | 1 | 3 | 80399 | 182109 | 18 |
| 512423 | seed1, index5 | 1537273 | 22 | 2 | 65365 | 367977 | 17 |
| 512425 | seed4, index0 | 1537277 | 22 | 3 | 59237 | 184755 | 18 |
| 512425 | seed4, index1 | 1537279 | 40 | 4 | 80319 | 91060 | 19 |

Here `seed1` and `seed4` refer to the stable bounded first-generation indices,
not the renumbered all-prefix menu. The two second-generation files are
`realized_generation2_seed1_depth3_all_prefixes.json` and
`realized_generation2_seed4_depth3_all_prefixes.json`.
All their hits, not just the four selected rows, have native certificates.

The seed1 scan considers 125672 suffixes, of which 119045 are strong, finding
four +2 and two +4 hits. The seed4 scan considers 125672 suffixes, of which
116479 are strong, finding one hit of each kind. These are **FINITE-EXACT**
depth-two closure witnesses, not an infinite grammar theorem.

The four selected B-edit supports are all a two-edge fan together with a
disjoint edge, namely the undirected graph `P3 disjoint-union K2`. Their fan
centers are respectively 21, 14, 2 and 33. The best first-generation +2 witness
has the same support type, with fan center 6; the best first-generation +4
witness instead has three disjoint edges. This is a precise recurring operation
shape, but it forgets the directions and context needed for path counts.

A stronger nested observation also has an exact failure boundary. Both chosen
first-generation targets, and three of the four chosen grandchildren, retain
the interior collar with root 32, sole inport 33 and prefix vertex 40, because
all their B edits lie inside the remaining C12. The +4-then-+4 branch reverses
`{22,32}`, giving root 32 an additional inport 22 and breaking this collar.
In fact this last B14 has **no** zero-ready vertex: its H is 80319 and every
`E_v(B,2)` is strictly greater than H. The other three selected grandchildren
each have exactly one zero-ready B vertex, namely 32. Thus merely choosing a
different B root cannot restore the missing collar in the final branch.
A recursive rule therefore needs an interior root-defect coordinate or a larger
choice of nested carrier; a fan-support label alone is insufficient. The exact
support, collar and all-vertex endpoint check is reproduced by `motif_check.py`
and recorded in `motifs.json`. No third-generation blind expansion was undertaken.

## Reproduction and adversarial checks

From the repository root, with Python and g++ on PATH:

```text
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/scan.py 3 1 --all-prefixes
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/realize.py 05-knowledge/results/zero_ready_strong_switch_20261103/literal/suffix_depth3_moving1_all_prefixes.json
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/verify.py
python -O 05-knowledge/results/zero_ready_strong_switch_20261103/literal/verify.py
```

The verifier imports no historical checker. It validates every native path in
its actual successive tournament and computes all source/target Hamiltonian
counts plus suffix h/F using `count.cpp`'s deletion recurrence. The discovery
code instead propagates partial paths forward. The independent checker guards
every signed 64-bit addition and multiplication, including the order-21 targets;
no overflow occurs. For the B14 discovery menu the factorial bound alone already
suffices, since every endpoint convolution is at most `14!`.
The transitive triple and directed triangle are positive controls. Duplicating
a path vertex, supplying the wrong inport, or corrupting a target row must all
be rejected. Checks use explicit exceptions and remain active under Python -O.

`scan.py --all-prefixes` removes the upper clock cutoff: since F>0, the count
strictly increases with the prefix size and the scan terminates after its first
overshoot. In these recorded universes every target is below 2^21, so the code's
safety maximum of 60 cannot exclude a solution. Second-generation menus use
`--seed 1` and `--seed 4`, and prohibit shrinking the actual inherited prefix.
The exact second-generation reproduction commands are:

```text
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/scan.py 3 1 --seed 1 --all-prefixes
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/scan.py 3 1 --seed 4 --all-prefixes
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/realize.py 05-knowledge/results/zero_ready_strong_switch_20261103/literal/generation2_seed1_depth3_all_prefixes.json
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/realize.py 05-knowledge/results/zero_ready_strong_switch_20261103/literal/generation2_seed4_depth3_all_prefixes.json
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/verify.py 05-knowledge/results/zero_ready_strong_switch_20261103/literal/realized_generation2_seed1_depth3_all_prefixes.json
python -O 05-knowledge/results/zero_ready_strong_switch_20261103/literal/verify.py 05-knowledge/results/zero_ready_strong_switch_20261103/literal/realized_generation2_seed1_depth3_all_prefixes.json
python 05-knowledge/results/zero_ready_strong_switch_20261103/literal/verify.py 05-knowledge/results/zero_ready_strong_switch_20261103/literal/realized_generation2_seed4_depth3_all_prefixes.json
python -O 05-knowledge/results/zero_ready_strong_switch_20261103/literal/verify.py 05-knowledge/results/zero_ready_strong_switch_20261103/literal/realized_generation2_seed4_depth3_all_prefixes.json
```
