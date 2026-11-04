# A complete actual second-round fork, with one explicit weak branch

**PROVED — [independent audit](CHILDREN_AUDIT.md) passed.** This supplements the frozen
[first-round proof](PROOF.md) without changing its source or receipts. It
provides a finite second round, not an unbounded returning class. Three
second-round preparations and children stay strong throughout. The fourth
ends with a transitive singleton prefix and is deliberately classified as
weak.

The two actual first children have old physical root8:

```
H=4291, F_8=4479, order13;
H=4289, F_8=5769, order14.
```

All names of each particular child are kept in its continuation. No
preparatory permanent vertex is born in this second round. Each completed
ready target then receives the two ordinary C2 births, with every previous
name still present.

## 1. The four literal targets and their access data

The target constructors use the exact `ready_target(core,ports,word)`
function in the frozen first-round source. It builds the funnel, adds the
indegree-one newborn, and takes converse. Its auxiliary has the designated
root as its sole out-neighbour. Each construction below is followed by an
explicit permutation, recorded in `children.json`, sending its designated
root to physical8. This permutation is an identification of the **new
target graph** on the same name set, not a free relabelling of the old input:
the subsequent literal arc changes implement the difference.

|Input H|gap|Core|Ports|Tail word|Target order|Strong?|
|---:|---:|---|---|---|---:|---|
|4291|1|`(2,4,1)`|{0}|00001111|13|yes|
|4291|2|`(2,12,9,17,7)`|{0}|000001|13|yes|
|4289|1|`(0,)`|{0}|10010111001|14|yes|
|4289|2|`(2,4,9,3)`|{1,2}|0101101, then prepend TT1|14|**no**|

The last ready core has order13. Its extra prefix is not a new permanent
birth or permission to delete an old name: the full fourteen-name target
identifies an existing old name with that prefix. The retained name3
becomes the global source at the final old-pair flip.

The actual target rows, in retained physical label order, are:

```
4291, gap1:
(7934,8068,8064,7926,7846,8070,7862,256,6233,4480,7040,640,2176)

4291, gap2:
(7934,8068,5248,8182,8102,8070,8118,256,7173,6532,4736,1156,2176)

4289, gap1:
(16126,16004,15488,16374,16006,16022,16310,15616,4151,15748,
 14592,8448,10240,256)

4289, gap2, weak:
(16118,16260,16000,16375,16006,16022,16310,3328,12341,15744,
 14592,256,2176,6272)
```

The auxiliaries are respectively physical7,7,13,11. In each target their
outgoing row is exactly the bit for physical8, so the general identity
`G=F_8` from PROOF §2 applies, including the weak target.

## 2. Native programs, observer transport, and the precise weak boundary

The old-pair sequences are:

```
4291, gap1: (1,9),(2,8),(7,12),(8,10),(9,11).
4291, gap2: (1,9),(2,9),(2,11),(3,8),(4,8),(6,8),(7,12),(10,11).
4289, gap1: (2,9),(6,8),(8,13),(12,13).
4289, gap2: (1,8),(6,8),(7,12),(7,13),(11,13),(0,3).
```

Each pair is implemented by the inherited helper/two-Hamilton-complement
compiler, with independent named-set replay. The complete receipt contains
all guard words, helper cuts, intermediate rows, and `(H,F_8,F_8−H)` at
every native and completed adjacency stage. There are23 old-pair flips,
46 guarded complements,23 temporary helpers/deletions, and69 helper/native
states. **All69 are strong.** The three unpadded branches also have every
completed adjacency stage strong, as do both C2 births and final children.

On the padded branch only the final helper deletion loses strongness:
the last flip `(0,3)` produces the displayed weak target with global source3.
Before deleting the temporary helper, that native state is still strong.
After deletion, the other thirteen vertices induce the strong ready core.
The two C2 births leave source3 dominating everything and act on that
strong core. Thus both C2 stages and the final child remain weak, with
exactly the declared source-prefix structure. This is not an all-strong
four-leaf tree.

The exact joint and final observers are:

|Input H|gap|Target F=G|S|L|K|Target defect|Child H|Child F_8|Child defect|
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|4291|1|4563|4292|1491|3072|272|12875|13689|814|
|4291|2|8371|4293|1715|6656|4080|12877|25113|12236|
|4289|1|6983|4290|2887|4096|2694|12869|20949|8080|
|4289|2|12094|4291|3902|8192|7805|12871|36282|23411|

Here `L=F_8(T−a)` and `K=P_{T−{a,8}}(2)` are counted independently on
the deletion graphs. Each `L+K=F_8` is an exact joint test, rather than a
claim inferred merely from the child H. Forward DP independently verifies
the backward-table H/F observations at the source, target, and both C2
stages of every branch. All of these are actual matrices, not response-only
realizability guesses.

## 3. Why padding preserves the needed gap

For any tournament T with retained a,r, prepend a transitive chain Z of
length l, all of whose vertices beat T. Hamilton paths must traverse Z
first, so H is unchanged. A path beginning at a cannot contain a prefix
vertex. Hence every prefix vertex must belong to the P slot in `[P A_a]`,
where it is forced into its unique order at the beginning: S is unchanged.

For `[W_r P]`, each prefix vertex may independently be assigned to either
of the two slots; the vertices assigned to a slot precede its T-part in
their unique order. Thus F is multiplied by `2^l`. The same assignment
argument multiplies L and K by `2^l`: the A_a or Q_ar slot cannot contain
prefix vertices, and the other two slots may receive them freely. Therefore
G scales by `2^l` and the zero joint deficit persists.

This is an exact, unbounded padding operation on the full response, but it
produces a weak tournament whenever l>0. In the only padded continuation
here l=1. It addresses name retention and scalar readiness, not the
independent all-strong output obligation. The explicit sole-out-neighbour
identity is an additional direct check of G=F in the literal target.

## 4. Scope, controls, and remaining work

Combining with the first-round programs gives actual values

```
                     1429
                 /          \
              4289          4291
             /    \        /    \
          12869  12871   12875  12877
          strong  weak   strong strong.
```

Every edge has a full named native implementation followed by C2, and no
old physical name is deleted. The fixed physical root remains8. The root
defect changes at every preparation and is fully carried through the tree.
The tree has finite depth two. It is neither a global integer-coverage
theorem nor a modular-universality theorem nor an unbounded closure claim.

The weak leaf is a positive broader-tournament control and an all-strong
near miss. The first failed implication is that a smaller scalar-ready
representative can receive a larger old name set while staying strong.
The missing coordinate is retained order. Finding a strong ready target
for H4289 of order at least14, or another direct strong next operation,
remains OPEN. A bounded enlarged-core scout did not pay an all-class
absence claim and is not used in this proof.

Reproduce with:

```
python3 04-computation/nonzero_preparation_20261103/surgery_children.py --output /tmp/children.json
python3 -O 04-computation/nonzero_preparation_20261103/surgery_children.py --output /tmp/children-O.json
```

The source checks exact hashes of the frozen first-round source modules.
It does not read their response receipts. The634 always-active gates and
full JSON receipt are pinned in `CHILDREN_MANIFEST.json`; independent
audit passed; see [CHILDREN_AUDIT.md](CHILDREN_AUDIT.md).
