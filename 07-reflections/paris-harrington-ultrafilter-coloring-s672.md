# Paris-Harrington Ultrafilter Coloring, S672

The useful correction is small but important:

```text
ultrafilter side is not enough;
we need extension rank.
```

HYP-2245 already taught us that the divisor-210 picture is exact on a Boolean
cube and leaky on the tournament quotient.  HYP-2246 then made the repair
algorithmic: attach the deleted owner and its `L/M/U` side.  Paris-Harrington
adds the recursive version of the same repair.  A coloring is a side-choice on
tuple atoms, but the real proof object is the tree of bad side-choices under
outer extension.

The exact pair miniature is the cleanest little model I found.

For pairs colored in two colors, with a relatively-large homogeneous triple,
the bad counts are:

```text
N=1..6: 1, 2, 6, 18, 12, 0
```

The tree is better than the sequence.  At `N=4`, only the middle edge-count
shell extends:

```text
edge_count 2: dead
edge_count 3: alive for one more step
edge_count 4: dead
```

At `N=5`, the middle shell itself dies.  That is exactly the user’s
upper/lower picture with one extra coordinate attached.  Blue/black says which
side of a complement pair the coloring lives on; the derivative rank says
whether the side-choice still has coherent children.

The ternary scout keeps me honest.  Random colorings stop finding bad examples
around `N=10`, but a local-repair search quickly finds bad colorings through
`N=20`.  That is not a theorem, but it is the right warning: a random pressure
heuristic can see density collapse while a coherent branch keeps escaping into
the tail.  That is the finite intuition behind the Paris-Harrington witness
function outrunning PA-provably recursive bounds.

So I would now read the repo chain as:

```text
HYP-2245:
  upper/lower filter side lives on a cube, leaks on quotient

HYP-2246:
  attach L/M/U deleted-owner address to repair tournament enumeration

HYP-2247:
  attach extension-rank profile to repair recursive bad-branch reasoning
```

For LRC14, this suggests a concrete proof target.  Define a bad node to be an
owner/carry fiber over the `C=27` shell with no lonely witness.  Define outer
extensions as coherent `+27` carry-owner lifts.  After adding the
owner-private deletion bit from HYP-2241, try to prove a rank drop:

```text
bad-child profile strictly decreases under every local carry extension
that preserves the visible Res_27 floor shadow.
```

That would make the owner/carry label more than a classifier.  It would become
a termination certificate.

For unit distance, the same template says: do not ask only whether the edge
count or unit-spine count is high enough.  Ask whether a construction shadow
can keep postponing every mandatory unit-spine witness into bulk/slab tail
coordinates.  The rank should be a spine-owner deletion profile plus remaining
extension count.

This also clarifies the ultrafilter language.  The tournament metagraph can
look like a divisor lattice, and the divisor lattice can be read as four
principal ultrafilters on `210`.  But recursion asks a sharper question:

```text
Can this side-choice keep extending after the initial segment is named?
```

Paris-Harrington is the theorem-shaped answer: standardly no, but the bound is
too large for PA to prove uniformly.  The repo method should imitate the local
part, not the grand unprovability claim: find the side channel that makes bad
outer extensions lose rank.
