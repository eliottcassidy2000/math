# Zero first, then repair the count and return a nonzero guard

**PROVED actual construction / FINITE-EXACT programs — [independent audit](DIRECT_AUDIT.md) passed.**
This is a different operation from H-neutral preparation followed by C2.
It uses exactly two permanent births to go from the actual eleven-name
H1429 source to a strong thirteen-name target carrying a positive joint
C2 guard. All old names and the physical root survive.

## 1. The common actual prefix

Start from the literal PREPARED_FOUR output in PROOF §1, on names0,...,10,
with root8 and `(H,F_8)=(1429,1563)`. Its actual auxiliary9 is zero-ready.
Apply the ordinary two-birth C2 there. The retained-name graph on0,...,12 is

```
(7934,7812,7816,7922,7814,8086,7862,7424,7775,4224,6656,512,2048).
```

Both births are strong. Exact direct path counts give
`H=4287=3·1429`, `F_8=4689=3·1563`. The two new permanent names are11,12.

We now change old arcs of this thirteen-name graph. These changes repair
H by+2 or+4 **after** the zero multiplication and simultaneously prepare a
positive auxiliary. They are not H-neutral at the preparation stage and
are not covered by ORDER_OBSTRUCTION's H-neutral-preparation+C2 hypothesis.

## 2. Three literal repairs with the same physical root

Reverse the following pairs, in the displayed orders. Every pair is
compiled by the guarded two-Hamilton-complement helper recipe; no old-pair
reversal is assumed to be a native primitive.

```
A: (2,3),(2,9),(3,8),(6,8),(8,12),(11,12).
B: (1,8),(2,3),(2,9),(3,8),(7,9),(8,12),(11,12).
C: (1,8),(2,3),(2,7),(2,12),(3,8),(4,8),(6,8),(8,11),(9,11).
```

The literal completed targets are:

```
A:
(7934,7812,7296,8182,7814,8086,8118,7424,3607,4228,6656,4608,256)

B:
(7934,8068,7296,8182,7814,8086,7862,7936,3669,4100,6656,4608,256)

C:
(7934,8068,3584,8182,8070,8086,8118,7428,5637,6272,6656,256,2052).
```

Their exact observers are:

|Repair|H|F_8=G|Auxiliary|S|L|K|Start gap|Root defect|
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|A|4289|7653|12|4291|1253|6400|2|3364|
|B|4291|6403|12|4292|1027|5376|1|2112|
|C|4291|8371|11|4293|1715|6656|2|4080|

In each target the auxiliary has old root8 as its sole out-neighbour.
The exact joint identity `G=L+K=F_8` from PROOF §2 therefore applies, and
the deletion-graph counters verify L and K separately. The old physical
root's response is **not** claimed to stay4689 during repair.

The targets also have independent ready-language constructors:

|Repair|Core|Entries|Tail|
|---|---|---|---|
|A|`(6,4,24,19,3)`|{0}|110001|
|B|`(2,4,25,19,3)`|{0}|010101|
|C|`(2,12,9,17,7)`|{0}|000001|

The source reconstructs these via funnel, indegree-one birth, converse,
and the explicit physical permutation in the receipt. That permutation
identifies a target on the retained names; the old-arc program implements
the actual map from the common zero child. No unnamed replacement is used.

## 3. Every native stage and the immediate next nonzero step

The three repairs have6,7,9 old-pair steps. Their complete receipts contain
44 guarded Hamilton complements,22 temporary helpers and helper-only
deletions, and66 helper/native states. Every such state is strong, and
every completed adjacency-stage graph is strong. The compiler is the
audited one in the frozen first-round source, while a separate named-set
interpreter checks each helper cut, both directed guards, both complements,
and deletion of only the helper. All intermediate `(H,F_8,F_8−H)` are
recorded with root8, never the current final array entry.

From the original H1429 graph to each ready target there are only two
permanent births. Temporary helper order is14, followed by deletion back
to13. Applying ordinary C2 immediately at the positive auxiliary now gives

|Repair|Two actual digits from1429|Next H|Next F_8|Next root defect|Order|
|---|---|---:|---:|---:|---:|
|A|+2, then+4|12871|22959|10088|15|
|B|+4, then+2|12875|19209|6334|15|
|C|+4, then+4|12877|25113|12236|15|

Both subsequent births are strong. Thus A supplies a strong version of
the numerical(+2,+4) branch that was weak in CHILDREN, by changing the
operation and retaining only thirteen names at the first child. It does
not refute ORDER_OBSTRUCTION: that theorem requires a strong H-neutral
preparation before each nonzero C2 multiplication. Here the first C2 uses
gap0 and the following repair changes H.

These are finite actual returns to specified positive C2 input interfaces.
The next C2 outputs are not asserted to belong to one class closed under
both future digits. No unbounded grammar is proved by these programs.

## 4. The obvious clock lift fails, with an exact error

There is a cheap natural unbounded test. Replace old physical0, the original
prefix singleton, by a transitive module TT_k, preserving every outside
incidence and keeping the same outside pair operations. Its physical root
is still old8, although its array index changes. All three final targets
remain strong and have the same positive gap and joint identity.

Their finite suffix signatures `(h,p,f,g,j,e,s)` are

```
A: (11,6400,3821,5651,3364,11,11),
B: (9,5376,3197,3587,2112,9,9),
C: (13,6656,4179,6513,4080,13,13).
```

The inherited exact ordered-prefix identity `H=h+2^k f−j`, `F=h+2^k f`
therefore proves, for every k>=1,

```
A: H=3821·2^k−3353, F=3821·2^k+11;
B: H=3197·2^k−2103, F=3197·2^k+9;
C: H=4179·2^k−4067, F=4179·2^k+13.
```

The scalar targets requested by a uniform +2/+4 map from the original
prepared clock are `2268·2^k−247` or `2268·2^k−245`. The errors are

```
A:1553(2^k−2), B:929(2^k−2), C:1911(2^k−2).
```

They vanish only at k=1. The first literal inflated clock k=2 has H11931,
10685,12649 instead of8825,8827,8827. The source directly checks those
three k2 matrices. The all-k formulas come from the exact suffix response,
not from interpolation. The auxiliary still has sole exit r; after
deleting it, a Hamilton path ending at r must put the transitive prefix
first, so its relevant ending count and hence start gap are independent
of k. Strongness follows by expanding the original strong one-vertex
module. This preserves the guard but fails the desired count map.

The failed implication is therefore precise: a small fixed-clock native
repair with a positive output guard does not automatically lift to the
old clock family. The missing coordinate is the full dyadic response,
alongside root defect and retained order.

## 5. Reproduction and status

```
python3 04-computation/nonzero_preparation_20261103/surgery_direct.py --output /tmp/direct.json
python3 -O 04-computation/nonzero_preparation_20261103/surgery_direct.py --output /tmp/direct-O.json
```

The source pins both imported first-round modules and reconstructs every
graph and program. It does not consume old numerical receipts. The full
new receipt and source/output pins are in `DIRECT_MANIFEST.json`. The independent analytic/source/replay audit has passed.
