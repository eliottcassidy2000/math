# Full complement words force a sharp small sheet-gcd edge

**PROVED reduction to an exact finite graph certificate; FINITE-EXACT exhaustive
certificate and actual sharpness control; INDEPENDENTLY AUDITED.** Every
actual balanced `6+7` decoder entry surviving the full inherited gcd profiles
has an edge in its connected seven-component atlas with sheet gcd at most six.
The constant six is attained by a complete actual entry with all inherited
filters retained. This is an edge-selection theorem, not a claim of strict
failure or a solution of LRC(14).

## 1. Inheritance and the missing coordinate

The primary supplier is the completed inherited-profile computation in
[overnight12, gcd-semigroup decoder descent](overnight12_20260906_lrc_gcd_semigroup.md),
including its full complement words, frozen in
[the profile JSON](../../04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json).
The strict actual atlas, the two-component kernel interpretation, and the
full two orientations of mixed three-support tests are inherited from the
same decoder route and
[the dual-pair consumer](continuing1_20260906_lrc_dual_pair.md).
The immediate new consumer is [the translated grid argument](third_20260906_grid.md):
a selected pair supplies overlap credit, and its sheet gcd multiplies the
component-count error. That consumer retains the actual pair ratio.

The six-object board is: scalar gcd caps; full complement words; sheet gcd
states; actual inert atlas edges; translated grid overlaps; and the mixed
support kernel. The scalar caps permit a connected seven-state hostile with
all edge gcds at least ten. Restoring the full words removes that hostile
and, after exact graph compilation, forces an edge gcd at most six.
The actual sharpness control preserves the kernel and all inherited profiles,
so the successful graph statistic has an actual decoder interpretation.

The source-to-target map sends the physical row `tV union gU` to the seven
integers `d_i=gcd(t,u_i)`. It preserves every sheet edge gcd and every profile
using all six smaller-component labels. It destroys primitive edge ratios,
physical heights, exact phase, and the other inherited profile words. The
complete actual atlas and all physical words are therefore checked separately
on the sharpness control. No converse realization theorem is inferred from a
valid seven-state row. Targeted searches of the inherited profile, graph and
small-edge statements found no existing version of this edge-selection result;
no external priority claim is made.

## 2. Exact reduction from actual entries

Let `Q=91^6`, and suppose an actual primitive thirteen-label decoder entry
has precisely two components

```text
tV union gU,  |V|=6, |U|=7,
gcd(V)=gcd(U)=gcd(t,g)=1.
```

Assume it survives every inherited gcd/complement-word profile. This is a
necessary condition for a hypothetical strict failure. Write

```text
d_i=gcd(t,u_i),                    i=1,...,7.
```

For any nonempty set `I` of `k<=6` tail indices, the subset consisting of
all six physical `tV` labels and the `k` physical `gu_i` labels has gcd

```text
c_I=gcd(d_i:i in I).
```

Indeed `gcd(tV)=t` and `gcd(t,g)=1`. Its remaining complement word is
exactly

```text
sort(gcd(c_I,d_j):j not in I).
```

Thus `(c_I,word)` belongs to the frozen profile level `7-k`. In particular,
each `d_i` is one of the 42 allowed core gcds in level six:

```text
1,2,3,4,5,6,8,9,10,11,12,15,16,17,18,20,22,23,24,25,
27,29,30,32,33,34,36,40,44,45,46,48,50,51,54,58,60,64,
66,72,88,90.
```

Also `gcd(d_1,...,d_7)=1`. The scalar shadow of these profiles gives
bounds `(90,30,9,4,2,1,1)` for gcds of subsets of sizes one through seven,
but the proof below retains the complete words.

If an actual atlas edge joins `u_i,u_j`, its sheet gcd is

```text
e_ij=gcd(t,gcd(u_i,u_j))=gcd(d_i,d_j).
```

The actual seven-component atlas is connected. If all its edges had
`e_ij>=7`, the graph on the seven states with every available edge
`gcd(d_i,d_j)>=7` would be connected too. The exact finite certificate in
the next section excludes precisely this possibility. Therefore some
actual atlas edge has `e_ij<=6`. This conclusion is about at least one
edge; it does not bound all edges or prescribe the selected edge's ratio.

## 3. Exhaustive graph certificate, including relabeling and repeated states

For threshold `s`, define a graph on a seven-entry multiset by putting an
edge between distinct positions when their state gcd is at least `s`.
Repeated state values are allowed: distinct physical speeds can have the
same sheet gcd. Sorting the multiset only quotients permutation symmetry;
it neither identifies positions inside profile words nor deletes repetitions.

For a partial multiset `S` of length `n<=7`, consider each subset `I` of
`k<=min(n,6)` positions. Its currently visible complement word must be a
submultiset of some complete allowed word of level `7-k` at core gcd
`gcd(S_i:i in I)`. All such possible subwords are generated directly from
the full frozen profile table. A partial multiset failing this test cannot
extend to a valid full multiset. At `n=7` the word lengths agree, so the
test is exact full profile membership; the global gcd-one condition is
also imposed.

Starting with every allowed one-state value at least `s`, the compiler
extends every retained sorted multiset by every allowed state connected to
at least one existing state. It retains the extension exactly when the
partial profile test passes. This growth is exhaustive: every connected
finite graph admits deletion of a leaf of a spanning tree while retaining
a connected induced subgraph, and every valid full profile restricts to a
permitted subword on every partial multiset. Induction therefore retains
a growth ordering of every possible full connected valid state multiset.

The exact layer sizes are:

| Minimum allowed edge gcd | size 1 | size 2 | size 3 | size 4 | size 5 | size 6 | size 7 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 7 | 36 | 71 | 153 | 249 | 256 | 111 | **0** |
| 6 | 37 | 111 | 381 | 825 | 1,052 | 575 | **26** |

Every one of the 26 final survivors is checked again by literal complete
profile membership, and the full sorted list is frozen in the output.
There is no state at threshold seven. This proves the asserted finite
certificate, and Section 2 transfers it to actual entries at all heights.
The finite arithmetic is a load-bearing part of this theorem; no analytic
classification of all profile graphs is claimed.

## 4. Why scalar caps and a loose atlas probe are insufficient

The seven states

```text
(14,48,50,63,75,80,84)
```

satisfy all 127 nonempty scalar subset gcd caps. Their graph of edges with
gcd at least ten is the connected tree

```text
14--84(14), 48--80(16), 48--84(12),
50--75(25), 50--80(10), 63--84(21).
```

The parenthesized numbers are edge gcds. At the state 14, however, the
complete complement word is `(1,2,2,2,7,14)`, and this word is absent from
the inherited level-six table. The first failed implication is replacing
full word membership by separate scalar subset ceilings. This example
already fails the exact allowed core-gcd projection, since 14 is absent.

The independent audit found a stronger hostile:

```text
(8,8,9,9,10,60,72).
```

Every subset gcd belongs to the **exact** allowed core-gcd set for its size,
and the threshold-seven graph is connected, with smallest edge gcd eight.
At state 60, however, the full word `(3,3,4,4,10,12)` fails level six.
This proves that even the complete scalar projections lose information needed
by the edge theorem. Both controls live in arithmetic quotients; neither is
claimed to be an actual surviving decoder entry.

A separate exploratory lift initially used membership in the full inert
multiplicative semigroup where the strict decoder requires each prime
exponent at most two. The smallest relevant witness was the proposed sum
`125=5^3`. Its inert prime factors do not make it a valid strict atlas sum.
This unpublished probe was corrected before using an actual control; the
frozen source rebuilds the strict atlas and actively rejects 125 while
accepting 25. The missing coordinate was the prime exponent, and no
mathematical truth surface was promoted from that probe.

## 5. A sharp actual entry with every inherited filter

Take

```text
V=(1,2,3,4,5,6),
U=(12584,14872,117,9999,98890,132990,10296),
g=1,
t=360*(1000Q+1)=204432930734760360.
```

Both shapes are primitive; their physical union has global gcd one and
thirteen distinct positive labels. The scales are coprime because `g=1`.
The six-component scale `t` is a component gcd, and is not claimed to be a
primitive integer. Its coprime multiplier `1000Q+1` does not alter any
sheet state. Their exact list is

```text
(8,8,9,9,10,30,72).
```

The full physical sum is `4293091545430247308<Q^2`. Rebuilding the entire
strict 5,855-pair atlas gives precisely two physical graph components of
sizes six and seven, with exact rational edge rank eleven. The actual
seven-component edges are:

| Tail-index pair (indices start at 0) | Reduced ratio in index order | Sheet gcd |
|---|---|---:|
| 0,6 | 11:9 | 8 |
| 1,6 | 13:9 | 8 |
| 2,6 | 1:88 | 9 |
| 3,6 | 101:104 | 9 |
| 4,5 | 29:39 | 10 |
| 5,6 | 155:12 | **6** |

Their reduced sums are `20,22,89,205,68,167`, all in the strict atlas.
The graph is a tree; its smallest sheet edge gcd is exactly six.

All 231 mixed supports fail the exact bounded signed-box test. The simple
reason is `t>Q` times the sum of the two largest `U` labels. For two `V`
labels and one tail, the smallest possible nonzero distinguished coefficient
is at least `t/max(U)>Q`. In the reverse orientation, every nonzero `tV`
term has magnitude at least `t`, exceeding the sum of the two bounded tail
terms. Internal support relations lie in the two component kernels;
connected decoder edges span those kernels. Together these facts establish
the actual equality `W_(Q,3)=V_dec` of rank eleven.

Every one of the 4,095 full physical inherited profiles passes, including
words not used by the seven-state reduction. Maximum subset gcds for sizes
seven through twelve are `(72,10,9,3,2,1)`. This is therefore an actual
sharpness control for the edge-selection theorem, not merely an abstract
graph realization. Its literal phase `x=1/7` has clearance exactly `1/7`.
It is safe, and no claim is made that other existing safety methods miss it.

## 6. Consumer and reproduction

For the selected edge, the [translated overlap inequality](third_20260906_grid.md)
can replace its general sheet-error multiplier 30 by six. Its reduced
primitive ratio remains unspecified, so any uniform ratio envelope must
still retain the full actual atlas. This is a direct consumer of the full
word information; sharpness prevents a universal improvement to five on
the same actual entry domain.

The [standalone producer](../../04-computation/third_20260906_decoder_profile_graph.py)
imports no repository implementation. Its only data dependency is the frozen
inherited profile table, pinned by raw LF SHA256
`935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.
The [frozen output](third_20260906_decoder_profile_graph.out) records both
complete layer counts, all 26 threshold-six survivors, the scalar hostile,
and the actual sharpness control.

```bash
python3 -B 04-computation/third_20260906_decoder_profile_graph.py
python3 -B -O 04-computation/third_20260906_decoder_profile_graph.py
```

No theorem ID, external priority, universal decoder closure, or LRC(14)
solution is asserted. The exact remaining phase and ratio obligations stay
with their respective consumers.

The [independent audit](third_20260906_profile_edge_audit.md) uses a different
search: retain only exact scalar gcd projections during connected growth,
then test full words on all 104 threshold-seven and 493 threshold-six final
rows. It independently gets zero and 26 survivors, reconstructs the actual
sharpness control, and accepts the analytic reduction. Its stronger scalar
projection hostile is retained in Section 4 and in the primary checker.

Normal and optimized primary outputs agree byte for byte, with 8,118
always-active gates. Raw LF SHA256 values:

```text
source 05bb4e0da7c5a527efa4791819ae19a118a9a8e6f800dd9fa93c1ee645440058
output 4dc33b2ed7d3f1b4e4cd5c4549d3a6dbf07fddf5be1f61f1d2d8ca99563dc24d
semantic 140e1ad43bbccc9e7d0c2732d2fc72ebff9e9c3fb6343b7e87f6b221cf045f45
```
